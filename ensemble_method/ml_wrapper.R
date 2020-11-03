# Averaging method (comparison purposes) ----------------------------------

mean_fit = function(x,y) {
  mean = mean(y)
  output = list("mean"=mean,"n"=nrow(x))
  output
}

predict.mean_fit = function(mean_fit,x,y,xnew=NULL,weights=FALSE) {
  
  if (is.null(xnew)) fit = rep(mean_fit$mean,nrow(x))
  else fit = rep(mean_fit$mean,nrow(xnew))
  
  if (isTRUE(weights)) w = matrix(1 / length(y),nrow(xnew),nrow(x))
  else w = NULL
  
  list("prediction"=fit,"weights"=w)
}


ols_fit = function(x,y) {
  
  x = cbind(rep(1,nrow(x)),x)
  
  ols_coef = lm.fit(x,y)$coefficients

  ols_coef
}


# Least squares regression ------------------------------------------------

predict.ols_fit = function(ols_fit,x,y,xnew=NULL,weights=FALSE) {
  
  if (is.null(xnew)) xnew = x
  
  x = cbind(rep(1,nrow(x)),x)
  xnew = cbind(rep(1,nrow(xnew)),xnew)
  
  # Remove variables that were dropped due to collinearity
  x = x[,!is.na(ols_fit)]
  xnew = xnew[,!is.na(ols_fit)]
  
  # Calculate hat matrix
  hat_mat = xnew %*% solve(crossprod(x),tol=2.225074e-308) %*% t(x)
  fit = hat_mat %*% y
  
  if (weights==FALSE) hat_mat = NULL
  
  list("prediction"=fit,"weights"=hat_mat)
}


# Cross-validated ridge regression ----------------------------------------

#' This function estimates cross-validated ridge regression based on the \code{\link{glmnet}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param w vector of weights
#' @param ... Pass \code{\link{glmnet}} options
#' @import glmnet
#'
#' @return An object with S3 class "glmnet"
#' @export

ridge_fit = function(x,y,args=list()) {
  
  ridge = do.call(cv.glmnet,c(list(x=x,y=y,alpha=0),args))
  
  ridge
}

predict.ridge_fit = function(ridge_fit,x,y,xnew=NULL,weights=FALSE) {
  
  if (is.null(xnew)) xnew = x
  
  fit = predict(ridge_fit,newx=xnew,type="response")
  
  if (weights==FALSE) hat_mat = NULL 
  else {
    # Get covariate matrices
    n = nrow(x)
    p = ncol(x)
    x = scale(x)
    x = cbind(rep(1,nrow(x)),x)
    xnew = scale(xnew)
    xnew = cbind(rep(1,nrow(xnew)),xnew)
    
    # Calculate hat matrix, see also (https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat)
    hat_mat = xnew %*% solve(crossprod(x) + ridge_fit$lambda.min  * n / sd(y) * diag(x = c(0, rep(1,p)))) %*% t(x)
    fit = hat_mat %*% y
  }
  
  list("prediction"=fit,"weights"=hat_mat)
}


# Cross-validated LASSO regression ----------------------------------------

#' This function estimates cross-validated lasso regression based on the \code{\link{glmnet}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param w vector of weights
#' @param ... Pass \code{\link{glmnet}} options
#' @import glmnet
#'
#' @return An object with S3 class "glmnet"
#' @export

lasso_fit = function(x,y,args=list()) {
  
  lasso = do.call(cv.glmnet,c(list(x=x,y=y),args = args))
  lasso
}

predict.lasso_fit = function(lasso_fit,x,y,xnew=NULL,weights=FALSE) {
  
  if (isTRUE(weights)) stop("No weighted representation of Lasso available.")
  if (is.null(xnew)) xnew = x
  
  fit = predict(lasso_fit,newx=xnew,type="response",s="lambda.min")
  
  list("prediction"=fit,"weights"="No weighted representation of Lasso available.")
}

# Cross-validated LASSO regression with first order interaction terms ----------------------------------------

#' This function estimates cross-validated lasso regression based on the \code{\link{glmnet}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param w vector of weights
#' @param ... Pass \code{\link{glmnet}} options
#' @import glmnet
#'
#' @return An object with S3 class "glmnet"
#' @export

lasso_inter_fit = function(x,y,args=list()) {
  # Final version
  inter_poly = interact.all(x)
  tmp = as.matrix(cbind(x, inter_poly))

  
  lasso_inter = do.call(cv.glmnet,c(list(x=tmp,y=y),args = args))
  lasso_inter
}

predict.lasso_inter_fit = function(lasso_inter_fit,x,y,xnew=NULL,weights=FALSE) {
  
  if (isTRUE(weights)) stop("No weighted representation of Lasso available.")
  if (is.null(xnew)) {
    inter_poly = interact.all(x)
    tmp = as.matrix(cbind(x, inter_poly))
    } else {
      inter_poly = interact.all(xnew)
      tmp = as.matrix(cbind(xnew, inter_poly))
      }
  
  fit = predict(lasso_inter_fit, newx = tmp, type = "response", s = "lambda.1se")
  
  list("prediction"=fit,"weights"="No weighted representation of Lasso available.")
}

# Random Forest -----------------------------------------------------------

forest_grf_fit = function(x,y,args=list()) {

  rf = do.call(regression_forest,c(list(X=x,Y=y),args))

  rf
}


predict.forest_grf_fit = function(forest_grf_fit,x,y,xnew=NULL,weights=FALSE) {
  
  if (is.null(xnew)) xnew = x
  
  fit = predict(forest_grf_fit,newdata=xnew)$prediction
  
  if (weights==TRUE) w = get_sample_weights(forest_grf_fit,newdata=xnew)
  else w = NULL
  
  list("prediction"=fit,"weights"=w)
}


# Post-LASSO regression ---------------------------------------------------

plasso_fit = function(x,y,args=list()) {
  plasso = do.call(plasso,c(list(x=x,y=y),args))
  
  plasso
}

predict.plasso_fit = function(plasso_fit,x,y,xnew=NULL,weights=FALSE) {

  if (is.null(xnew)) xnew = x
  x = add_intercept(x)
  xnew = add_intercept(xnew)
  
  # Fitted values for post lasso
  nm_act = names(coef(plasso_fit$lasso_full)[,plasso_fit$ind_min_pl])[which(coef(plasso_fit$lasso_full)[,plasso_fit$ind_min_pl] != 0)]
  
  xact = x[,nm_act,drop=F]
  xactnew = xnew[,nm_act,drop=F]
  
  # Remove potentially collinear variables
  coef = lm.fit(xact,y)$coefficients
  xact = xact[,!is.na(coef)]
  xactnew = xactnew[,!is.na(coef)]
  
  hat_mat = xactnew %*% solve(crossprod(xact),tol=2.225074e-308) %*% t(xact)
  fit_plasso = hat_mat %*% y
  if (weights==FALSE) hat_mat = NULL
  
  list("prediction"=fit_plasso,"weights"=hat_mat)
}
  


# XGBoost (eXtreme Gradient Boosting) -------------------------------------

xgboost_fit = function(x,y,args=list()) {
  
  dx = xgb.DMatrix(data = x, label = y) 
  
  xgb = do.call(xgb.train, c(list(data = dx),args))
  
  xgb
}


predict.xgboost_fit = function(xgb_fit,x,y,xnew=NULL,weights=FALSE) {
  
  if (is.null(xnew)) xnew = x
  
  fit = predict(xgb_fit,newdata=xnew)

  list("prediction"=fit, "weights" = "No weighted representation of XGBoost available")
}


# Neural Network implementation with the neuralnet package ----------------------------------------------------------

neural_net_fit = function(x,y,args=list()) {
  
  x_nn = as.data.frame(x)
  y_nn = as.data.frame(y)
  
  colnames(y_nn) = c("Y1")
  
  tempdata = cbind(y_nn, x_nn)
  
  nn_formula = as.formula(paste("Y1 ~", paste(colnames(x_nn), collapse = " + ")))
  assign("failed", FALSE, env=globalenv())
  tryCatch({nnet = do.call(neuralnet, c(list(formula = nn_formula, data = tempdata), args))
  state = FALSE
  return(list("model" = nnet, "state" = state))
  }, warning = function(w){w
    assign("failed", TRUE, env=globalenv())
  }, 
  finally = {
    if (get("failed", env=globalenv())){
      state = TRUE
      return(list("model" = NULL, "state" = state))
    } else {
      state = FALSE
      return(list("model" = nnet, "state" = state))
    }
  })
}


predict.neural_net_fit = function(neural_net_fit,x,y,xnew=NULL,weights=FALSE) {
  
  if (is.null(xnew)) xnew = x
  
  xnew_nn = as.data.frame(xnew)
  
  if (neural_net_fit$state){
    preds = rep(mean(y), nrow(xnew_nn))
  } else {
    preds = predict(neural_net_fit$model, newdata = xnew_nn)
  }
  
  list("prediction"= preds, "weights" = "No weighted representation of the neural network available")
}


# Neural Network implemented with Keras and Tensorflow --------------------

nn_keras_fit = function(x,y,args=list()) {
  column_names = rep(NA, ncol(x))
  for(i in 1:ncol(x)){
    column_names[i] = paste0("x",i)
  }
  
  data = x %>% as_tibble(.name_repair = "minimal") %>%
    setNames(column_names) %>%
    mutate(label = y)
  
  spec = feature_spec(data, label ~.) %>%
    step_numeric_column(all_numeric()) %>%
    fit()
  
  layer = layer_dense_features(
    feature_columns = dense_features(spec),
    dtype = tf$float32
  )
  
  model = build_model(data, 
                      spec, 
                      units1 = args$units1, 
                      units2 = args$units2, 
                      act.fct1 = args$act.fct1, 
                      act.fct2 = args$act.fct2,
                      act.fctfinal = args$act.fctfinal, 
                      loss.fct = args$loss.fct, 
                      eval.metric = args$eval.metric)
  
  early_stop <- callback_early_stopping(monitor = "val_loss", patience = 30)
  
  model %>% fit(
    x = data %>% dplyr::select(-label),
    y = data$label,
    epochs = 500,
    validation_split = 0.2,
    verbose = 0,
    callbacks = list(early_stop)
  )
  return(model)
}


predict.nn_keras_fit = function(nn_keras_fit,x,y,xnew=NULL,weights=FALSE) {
  
  if (is.null(xnew)) xnew = x
  
  column_names = rep(NA, ncol(xnew))
  for(i in 1:ncol(xnew)){
    column_names[i] = paste0("x",i)
  }
  
  newdata = xnew %>% as_tibble(.name_repair = "minimal") %>%
    setNames(column_names)
  
  preds = nn_keras_fit %>% predict(newdata)
  
  list("prediction"= preds, "weights" = "No weighted representation of the neural network available")
}
