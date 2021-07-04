
create_method = function(method,x_select=NULL,args=list(),name=NULL) {
  
  if (!(is.character(method) & length(method) == 1)) stop("Provide single string to define method.")
  if (!(any(method == c("mean","ols","ridge","plasso","forest_grf","lasso", "lasso_inter",
                        "elastic_net", "xgboost", "neural_net", "nn_keras")))) stop("Provide one of these options c(\"mean\",\"ols\",\"ridge\",\"plasso\", \"lasso_inter\",\"forest_grf\",\"xgboost\", \"nn_keras\" ,\"neural_net\")")
  if (!(is.null(args) | is.list(args))) stop("Provide either NULL or list for args.")
  if (!(is.null(x_select) | is.logical(x_select))) stop("Provide either NULL or logical for x_select.")
  if (!((is.character(name) & length(name) == 1) | is.null(name))) stop("Provide single string to name method.")
  
  list(method=method,args=args,x_select=x_select,name=name,weights=weights)
}



#### This function takes a dataset and strings with variable names as input

design_matrix = function(data,int=NULL,int_p=2,poly=NULL,poly_p=2,log=NULL) {
  
  # Part for interactions
  if (!is.null(int)) {
    int = paste0("+(",paste0(int,collapse="+"),")^",toString(int_p))
  }
  
  # Part for polynomials
  if (!is.null(poly)) {
    poly = paste0("+poly(",paste0(poly,collapse=paste0(",",toString(poly_p),",raw=TRUE)+poly(")),",",toString(poly_p),",raw=TRUE)",
                   "-(",paste0(paste0(poly,collapse="+")),")")
  }
  
  # Part for logs
  if (!is.null(log)) {
    # Check whether some variables can't be logged because not positive
    ind_neg = colMins(data[,log])<=0
    if (sum(ind_neg)>0) {
      cat("\n Following variables not modified to be logged because of non-positive values:",paste(colnames(data[,log])[ind_neg]),"\n" )
      log = log[!ind_neg]
    }
    log = paste0("+log(",paste0(log,collapse=")+log("),")")
  }
  
  # Combine the three parts
  fmla = as.formula(paste("~0",int,poly,log))
  
  # Generate matrix
  data = model.matrix(fmla,data=as.data.frame(data))
  # Clean variable names to make sense
  colnames(data) = gsub("poly\\(","",colnames(data))
  colnames(data) = gsub(paste0(", ",toString(poly_p),", raw = TRUE)"),"",colnames(data))
  colnames(data) = gsub("log\\(","ln_",colnames(data))
  colnames(data) = gsub("\\)","",colnames(data))
  
  return(data)
}



#### This function takes a matrix of data and removes
#### 1. Variables without variation
#### 2. Dummy variables where one group is nearly empty (optional in one of both treatment groups)
#### 3. Redundant (highly correlated variables)

data_screen = function(data,treat=NULL,bin_cut=0.01,corr_cut=0.99,print=TRUE) {
  
  ## Kick out variables with no variation
  # Identify the names
  nm_del = colnames(data)[colSds(data) == 0]
  # Describe identified variables
  if (print==TRUE) {
    cat("\n\n Variables with no variation:",nm_del,"\n\n")
    if (identical(nm_del, character(0)) == FALSE) print(describe(data[,nm_del],fast=TRUE))
  }
  # Remove identified variables
  if (identical(nm_del, character(0)) == FALSE) data = data[,!colnames(data) %in% nm_del]
  
  ## Remove dummy variables lower than threshold in one of the two treatment groups
  # Identify dummies
  bin = apply(data,2,function(x) { all(x %in% 0:1) })

  # Calculate means of all variables and check whether they are potentially close to 0 or 1
  if (is.null(treat)) {
    mean = colMeans(data)
    bel_cut = (mean<bin_cut | mean > (1-bin_cut))
  } else {
    mean1 = colMeans(data[d==1,])
    mean0 = colMeans(data[d==0,])
    bel_cut = (mean1<bin_cut | mean1 > (1-bin_cut) | mean0<bin_cut | mean0 > (1-bin_cut))
  }
  
  # Identify names that are binary and close to 0 and 1
  nm_del = colnames(data)[bin & bel_cut]
  if (print==TRUE) {
    cat("\n\n Dummy variables close to 0 or 1:",nm_del,"\n\n")
    if (identical(nm_del, character(0)) == FALSE) print(describe(data[,nm_del],fast=TRUE))
  }
  
  # Remove identified variables
  if (identical(nm_del, character(0)) == FALSE) data = data[,!colnames(data) %in% nm_del]
  
  ## Remove all redundant (nearly perfectly correlated) variables
  # Calculate correlation matrix and consider only upper diagonal
  cor = (abs(cor(data))>corr_cut)
  cor[lower.tri(cor, diag=TRUE)] = FALSE
  
  # Identify names of redundant variables
  nm_del = colnames(cor)[colSums(cor)>0]
  
  if (print==TRUE) {
    cat("\n\n Variables (nearly) perfectly correlated:",nm_del,"\n\n")
    if (identical(nm_del, character(0)) == FALSE) print(describe(data[,nm_del],fast=TRUE))
  }
  if (identical(nm_del, character(0)) == FALSE) data = data[,colSums(cor)==0]
  
  return(data)
}


#' Matrix for cross-fitting indicators
#'
#' Creates matrix of binary fold indicators (n x # cross-folds)
#'
#' @param n Number of of observations
#' @param cf Number of cross-fitting folds
#' @param cl Optional vector of cluster variable if cross-fitting should account for clusters.
#'
#' @importFrom dplyr ntile
#' @import stats
#'
#' @return n times # cross-folds matrix of binary fold indicators
#'
#' @export

prep_cf_mat = function(n,cf,cl=NULL) {
  
  if (cf == 1) cf_mat = matrix(rep(1,n),ncol=1)
  
  if (!is.null(cl)) {
    fold = ntile(runif(length(unique(cl))),cf)
    fold = factor(fold[match(cl,unique(cl))])
    cf_mat = model.matrix(~0+fold)
  }
  else {
    fold = factor(ntile(runif(n),cf))
    cf_mat = model.matrix(~0+fold)
  }
  colnames(cf_mat) = sprintf("CF %d",1:cf)
  
  return(cf_mat)
}

# interacts all inputted variables with each other and forms second order polynomials

interact.all <- function(input){
  output <- do.call("cbind", lapply(1:ncol(input), FUN=function(z){ 
    do.call("cbind", lapply(z:ncol(input), FUN=function(x){
      tmp <- data.frame(input[,z] * input[,x])
      colnames(tmp)[1] <- paste(colnames(input)[z], ":", colnames(input)[x], sep="") # Multiplication symbol changed to : as per R style instead of x to avoid confusion with variables containing x
      tmp
    }))
  }))
  output
}


#' Model construction for Keras
#' creates model structure as defined with given parameters
#'

build_model <- function(data, 
                        spec, 
                        twolayers = FALSE,
                        units1 = 64, 
                        units2 = 32, 
                        act.fct1 = "relu", 
                        act.fct2 = "relu",
                        act.fctfinal = NA, 
                        loss.fct = "mse", 
                        eval.metric = "mean_squared_error", 
                        l1_1, 
                        learning_rate = 0.01,
                        optimizer = "adagrad") {
  
  input <- layer_input_from_dataset(data %>% dplyr::select(-label))
  if (twolayers == FALSE) {
    if(is.na(act.fctfinal)){
      output <- input %>% 
        layer_dense_features(dense_features(spec)) %>% 
        layer_dense(units = units1, activation = act.fct1, kernel_regularizer = regularizer_l1(l1_1))%>%
        layer_dense(units = 1)
    } else {
      output <- input %>% 
        layer_dense_features(dense_features(spec)) %>% 
        layer_dense(units = units1, activation = act.fct1, kernel_regularizer = regularizer_l1(l1_1))%>%
        layer_dense(units = 1, activation = act.fctfinal) 
    }
  } else {
    if(is.na(act.fctfinal)){
      output <- input %>% 
        layer_dense_features(dense_features(spec)) %>% 
        layer_dense(units = units1, activation = act.fct1, kernel_regularizer = regularizer_l1(l1_1))%>%
        layer_dense(units = units2, activation = act.fct2)%>%
        layer_dense(units = 1)
    } else {
      output <- input %>% 
        layer_dense_features(dense_features(spec)) %>% 
        layer_dense(units = units1, activation = act.fct1, kernel_regularizer = regularizer_l1(l1_1))%>%
        layer_dense(units = units2, activation = act.fct2)%>%
        layer_dense(units = 1, activation = act.fctfinal) 
    }
  }
  
  model <- keras_model(input, output)
  
  if(optimizer == "adam") {
    model %>% 
      compile(
        loss = loss.fct,
        optimizer = optimizer_adam(lr = learning_rate),
        metrics = list(eval.metric)
      )
  } else if (optimizer == "adamax") {
    model %>% 
      compile(
        loss = loss.fct,
        optimizer = optimizer_adamax(lr = learning_rate),
        metrics = list(eval.metric)
      )
  } else if (optimizer == "adagrad") {
    model %>% 
      compile(
        loss = loss.fct,
        optimizer = optimizer_adagrad(lr = learning_rate),
        metrics = list(eval.metric)
      )
  } else if (optimizer == "rmsprop") {
    model %>% 
      compile(
        loss = loss.fct,
        optimizer = optimizer_rmsprop(lr = learning_rate),
        metrics = list(eval.metric)
      )
  } else if (optimizer == "sgd") {
    model %>% 
      compile(
        loss = loss.fct,
        optimizer = optimizer_sgd(lr = learning_rate),
        metrics = list(eval.metric)
      )
  } else {
    print("unknown optimizer!")
  }
  
  
  model
}






