dml_ens_ml = function(y, d, x, ps_methods, oc_methods, ml_methods, ens_folds = 2) {
  # Double Machine Learning estimator using the efficient score function as introduced by Chernozhukov et al. (2018)
  # Implementation inspired by Knaus DML for multiple treatments https://github.com/MCKnaus/dmlmt 
  # check if the given treatment vector is binary
  binary <- all(d %in% 0:1)
  if (binary){
    # split data set randomly into two halves for cross-fitting
    n = length(y)
    fold = model.matrix(~0 + factor(ntile(runif(n),2)))[,1]
    
    x_main <- x[!fold, ]
    x_aux <- x[as.logical(fold), ]
    
    y_main <- y[!fold]
    y_aux <- y[as.logical(fold)]
    
    d_main <- d[!fold]
    d_aux <- d[as.logical(fold)]
    
    # get dummy for treated (d) and controls (1-d)
    treatment_main = cbind(d_main, 1-d_main)
    treatment_aux = cbind(d_aux, 1-d_aux)
    
    # estimate the propensity score p(x) and 1-p(x) for inverse probability weighting
    # create empty list of matrices to store values
    p_hat_main = lapply(seq_along(ml_methods), function(f){f = matrix(NA, nrow(x_main), 2)})
    p_hat_aux = lapply(seq_along(ml_methods), function(f){f = matrix(NA, nrow(x_aux), 2)})
    names(p_hat_main) = ml_methods
    names(p_hat_aux) = ml_methods
    
    # run the ensemble with main data for the model and aux data for the evaluation
    p_hat_aux_ensemble = ensemble(ps_methods, x = x_main, y = d_main, xnew = x_aux, nfolds = ens_folds) # predict the propensity score with main model, aux data
    # extract predictions for ensemble
    p_hat_aux[[1]][,1] = p_hat_aux_ensemble$ensemble # extract prediction values
    p_hat_aux[[1]][,2] = 1 - p_hat_aux[[1]][,1] # retrieve 1-p(x) for main model and aux data
    
    # run the ensemble with switched data sets
    p_hat_main_ensemble = ensemble(ps_methods, x = x_aux, y = d_aux, xnew = x_main, nfolds = ens_folds) # predict the propensity score with aux model, main data
    p_hat_main[[1]][,1] = p_hat_main_ensemble$ensemble # extract prediction values
    p_hat_main[[1]][,2] = 1 - p_hat_main[[1]][,1] # retrieve 1-p(x) for aux model and main data
    
    # extract predictions for each ml method
    for(i in 1:(length(ml_methods) - 1)){
      i_loc = i + 1
      p_hat_main[[i_loc]][,1] = p_hat_main_ensemble$fit_full$predictions[,i]
      p_hat_main[[i_loc]][,2] = 1 - p_hat_main[[i_loc]][,1]
    }
    for(i in 1:(length(ml_methods) - 1)){
      i_loc = i + 1
      p_hat_aux[[i_loc]][,1] = p_hat_main_ensemble$fit_full$predictions[,i]
      p_hat_aux[[i_loc]][,2] = 1 - p_hat_main[[i_loc]][,1]
    }
    
    # estimate the potential outcome 
    # create empty list of matrices to store values
    y_hat_main = lapply(seq_along(ml_methods), function(f){f = matrix(NA, nrow(x_main), 2)})
    y_hat_aux = lapply(seq_along(ml_methods), function(f){f = matrix(NA, nrow(x_aux), 2)})
    names(y_hat_main) = ml_methods
    names(y_hat_aux) = ml_methods
    
    
    for (i in 1:length(unique(d))){
      # only use either treated or non-treated as groups (i = 1,2)
      # derive the model with the main data set and predict the outcomes with the auxiliary data set
      # the prediction will be made with data of both treated and controls to receive potential outcomes
      y_hat_aux_ensemble = ensemble(oc_methods, x = x_main[treatment_main[,i] == 1, ], y = y_main[treatment_main[,i] == 1], xnew = x_aux, nfolds = ens_folds)
      y_hat_aux[[1]][,i] = y_hat_aux_ensemble$ensemble
      # same procedure as above but with switched data sets for cross-fitting purposes
      y_hat_main_ensemble = ensemble(oc_methods, x = x_aux[treatment_aux[,i] == 1, ], y = y_aux[treatment_aux[,i] == 1], xnew = x_main, nfolds = ens_folds)
      y_hat_main[[1]][,i] = y_hat_main_ensemble$ensemble
      
      for(j in 1:(length(ml_methods) - 1)){
        j_loc = j + 1
        y_hat_main[[j_loc]][,i] = y_hat_main_ensemble$fit_full$predictions[,j]
      }
      
      for(j in 1:(length(ml_methods) - 1)){
        j_loc = j + 1
        y_hat_aux[[j_loc]][,i] = y_hat_aux_ensemble$fit_full$predictions[,j]
      }
    }
    
    # estimations for the average treatment effect ate
    # calculate efficient score E[(y - mu_hat)/p(x) + mu_hat]
    # both scores combined fulfill Neyman orthogonality condition
    # this doubly robust estimator remains consistent with either p(x) or mu(x) not precisely specified
    # see Chernozhukov et al (2018), Hahn (1998), Robins and Rotnitzky (1995)
    ipw_mat_main = lapply(seq_along(ml_methods), function(f){f = matrix(NA, nrow(x_main), 2)})
    ipw_mat_aux = lapply(seq_along(ml_methods), function(f){f = matrix(NA, nrow(x_aux), 2)})
    names(ipw_mat_main) = ml_methods
    names(ipw_mat_aux) = ml_methods
    
    mu_mat_main = lapply(seq_along(ml_methods), function(f){f = matrix(NA, nrow(x_main), 2)})
    mu_mat_aux = lapply(seq_along(ml_methods), function(f){f = matrix(NA, nrow(x_aux), 2)})
    names(mu_mat_main) = ml_methods
    names(mu_mat_aux) = ml_methods
    
    
    # calculate the inverse probability weights d/p(x) and (1-d)/(1-p(x))
    for (i in 1:length(unique(d))){
      # predictions of all applications
      # main
      for(j in 1:length(ml_methods)){
        ipw_mat_main[[j]][,i] = treatment_main[,i] / p_hat_main[[j]][,i]
      }
      # auxiliary
      for(j in 1:length(ml_methods)){
        ipw_mat_aux[[j]][,i] = treatment_aux[,i] / p_hat_aux[[j]][,i]
      }
    }
    
    # finalize the score calculation
    for (i in 1:length(unique(d))){
      # main
      for(j in 1:length(ml_methods)){
        mu_mat_main[[j]][,i] = (y_main - y_hat_main[[j]][,i])*ipw_mat_main[[j]][,i] + y_hat_main[[j]][,i]
      }
      # auxiliary
      for(j in 1:length(ml_methods)){
        mu_mat_aux[[j]][,i] = (y_aux - y_hat_aux[[j]][,i])*ipw_mat_aux[[j]][,i] + y_hat_aux[[j]][,i]
      }
    }
    
    # calculate standard errors SE
    se_mu_main = matrix(NA, nrow = length(ml_methods), ncol = length(unique(d)))
    se_mu_aux = matrix(NA, nrow = length(ml_methods), ncol = length(unique(d)))
    for (i in 1:length(unique(d))){
      for (j in 1:length(ml_methods)) {
        se_mu_main[j,i] = sqrt(sum((mu_mat_main[[j]][,i] - mean(y_main))^2)/nrow(x_main))
      }
      for (j in 1:length(ml_methods)) {
        se_mu_aux[j,i] = sqrt(sum((mu_mat_aux[[j]][,i] - mean(y_aux))^2)/nrow(x_main))
      }
    }
    
    se_mu = (se_mu_main + se_mu_aux)/2
    colnames(se_mu) = c("mu_1", "mu_0")
    rownames(se_mu) = c("ensemble", "lasso", "xgb", "nn")
    
    
    # calculate te(D, theta, eta) = (g(1,X) - g(0,X)) + D(Y - g(1,X))/p(x) + (1-D)(Y - g(0,X))/(1-p(x)) - theta
    # which with moment and orthogonality condition fulfilled the difference between the two scores equals the treatment effect theta
    te_mat_main = lapply(seq_along(ml_methods), function(f){f = matrix(NA, nrow(x_main), 1)})
    te_mat_aux = lapply(seq_along(ml_methods), function(f){f = matrix(NA, nrow(x_aux), 1)})
    names(te_mat_main) = ml_methods
    names(te_mat_aux) = ml_methods
    # main
    for(j in 1:length(ml_methods)){
      te_mat_main[[j]] = mu_mat_main[[j]][,1] - mu_mat_main[[j]][,2]
    }
    # auxiliary
    for(j in 1:length(ml_methods)){
      te_mat_aux[[j]] = mu_mat_aux[[j]][,1] - mu_mat_aux[[j]][,2]
    }
    
    se_te_main = matrix(NA, nrow = length(ml_methods), 1)
    se_te_aux = matrix(NA, nrow = length(ml_methods), 1)
    se_te_main = do.call(rbind, lapply(te_mat_main, FUN = function(x){sd(x)}))
    se_te_aux = do.call(rbind, lapply(te_mat_aux, FUN = function(x){sd(x)}))

    se_te = (se_te_main + se_te_aux)/2
    colnames(se_te) = c("se_te")
    rownames(se_te) = c("ensemble", "lasso", "xgb", "nn")

    # get the cross-fitted average treatment effect
    ate_main = lapply(te_mat_main, FUN = function(x){mean(x)})
    ate_aux = lapply(te_mat_aux, FUN = function(x){mean(x)})
    
    # merge cross-fitted results
    ate = mapply(function(x,y) mean(c(x,y)), ate_main, ate_aux, SIMPLIFY = FALSE)
    te = lapply(Map(cbind,te_mat_main, te_mat_aux), rowMeans)
    
    # extract weights of ensemble
    w_ens_ps = colMeans(rbind(p_hat_main_ensemble$nnls_weights, p_hat_aux_ensemble$nnls_weights))
    w_ens_oc = colMeans(rbind(y_hat_main_ensemble$nnls_weights, y_hat_aux_ensemble$nnls_weights))
    
    # define single output
    output = list("ate" = ate, "te" = te, "w_ens_ps" = w_ens_ps, "w_ens_oc" = w_ens_oc, "se_te" = se_te, "se_po" = se_mu)
    
    return(output)
    
  } else {
    print("Non-binary treatments, please check your data")
  }
}