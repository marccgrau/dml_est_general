dml_est_cf_ensemble = function(y, d, x, ps_methods, oc_methods) {
  # Double Machine Learning estimator using the efficient score function as in Hahn 1998
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
    p_hat_main = matrix(NA, nrow(x), 2) # create empty matrix to store values
    p_hat_aux = matrix(NA, nrow(x), 2)
    
    p_hat_aux_ensemble = ensemble(ps_methods, x = x_main, y = d_main, xnew = x_aux)
    p_hat_aux[,1] = p_hat_aux_ensemble$fit_full$predictions # predict the propensity score with main model, aux data
    p_hat_aux[,2] = 1 - p_hat_aux[,1] # retrieve 1-p(x) for main model and aux data
    
    p_hat_main_ensemble = ensemble(ps_methods, x = x_aux, y = d_aux, xnew = x_main)
    p_hat_main[,1] = p_hat_main_ensemble$fit_full$predictions # predict the propensity score with aux model, main data
    p_hat_main[,2] = 1 - p_hat_main[,1] # retrieve 1-p(x) for aux model and main data
    
    # estimate the potential outcome 
    y_hat_main = matrix(NA, nrow(x), 2) # create empty matrix to store values
    y_hat_aux = matrix(NA, nrow(x), 2)
    
    for (i in 1:2){
      y_hat_aux_ensemble = ensemble(oc_methods, x = x_main[treatment_main[,i] == 1, ], y = y_main[treatment_main[,i] == 1], xnew = x_aux)
      y_hat_aux[,i] = y_hat_aux_ensemble$fit_full$predictions
      
      y_hat_main_ensemble = ensemble(oc_methods, x = x_aux[treatment_aux[,i] == 1, ], y = y_aux[treatment_aux[,i] == 1], xnew = x_main)
      y_hat_main[,i] = y_hat_main_ensemble$fit_full$predictions
    }
    
    # calculate efficient score E[(y - mu_hat)/p(x) + mu_hat]
    # this doubly robust estimator remains consistent with either p(x) or mu(x) misspecified
    mu_mat_main = matrix(NA, nrow(x), 2)
    ipw_mat_main = matrix(NA, nrow(x), 2)
    
    mu_mat_aux = matrix(NA, nrow(x), 2)
    ipw_mat_aux = matrix(NA, nrow(x), 2)
    
    # calculate the inverse probability weights d/p(x) and (1-d)/(1-p(x))
    for (i in 1:2){
      ipw_mat_main[,i] = treatment_main[,i] / p_hat_main[,i]
      
      ipw_mat_aux[,i] = treatment_aux[,i] / p_hat_aux[,i]
    }
    
    for (i in 1:2){
      mu_mat_main[,i] = (y_main - y_hat_main[,i])*ipw_mat_main[,i] + y_hat_main[,i]
      
      mu_mat_aux[,i] = (y_aux - y_hat_aux[,i])*ipw_mat_aux[,i] + y_hat_aux[,i]
    }
    
    # calculate psi(D, theta, eta) = (g(1,X) - g(0,X)) + D(Y - g(1,X))/p(x) + (1-D)(Y - g(0,X))/(1-p(x)) - theta
    # which with moment and orthogonality condition fulfilled the difference between the two scores equals the treatment effect theta
    psi_main = mu_mat_main[,1] - mu_mat_main[,2]
    psi_aux = mu_mat_aux[,1] - mu_mat_aux[,2]
    
    # average all the treatment effects to get the ATE
    ate_main = mean(psi_main)
    ate_aux = mean(psi_aux)
    
    ate = mean(c(ate_main, ate_aux))
    
    # extract weights of ensemble
    w_ens_ps = colMeans(rbind(p_hat_main_ensemble, p_hat_aux_ensemble))
    w_ens_oc = colMeans(rbind(y_hat_main_ensemble, y_hat_aux_ensemble))
    
    return("ate" = ate, "w_ens_ps" = w_ens_ps, "w_ens_oc" = w_ens_oc)
    
  } else {
    print("Non-binary treatments, please check your data")
  }
}