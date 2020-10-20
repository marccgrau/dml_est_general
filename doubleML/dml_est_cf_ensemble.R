dml_est_cf_ensemble = function(y, d, x, ps_methods, oc_methods) {
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
    # create empty matrices to store values
    p_hat_main_ens = matrix(NA, nrow(x), 2) 
    p_hat_aux_ens = matrix(NA, nrow(x), 2)
    
    p_hat_aux_lasso = matrix(NA, nrow(x), 2)
    p_hat_aux_xgb = matrix(NA, nrow(x), 2)
    p_hat_aux_nn = matrix(NA, nrow(x), 2)
    p_hat_main_lasso = matrix(NA, nrow(x), 2)
    p_hat_main_xgb = matrix(NA, nrow(x), 2)
    p_hat_main_nn = matrix(NA, nrow(x), 2)
    
    # run the ensemble with main data for the model and aux data for the evaluation
    p_hat_aux_ensemble = ensemble(ps_methods, x = x_main, y = d_main, xnew = x_aux) # predict the propensity score with main model, aux data
    # extract predictions for ensemble
    p_hat_aux_ens[,1] = p_hat_aux_ensemble$ensemble # extract prediction values
    p_hat_aux_ens[,2] = 1 - p_hat_aux_ens[,1] # retrieve 1-p(x) for main model and aux data
    
    # extract predictions for each ml method
    p_hat_aux_lasso[,1] = p_hat_aux_ensemble$fit_full$predictions[,1]
    p_hat_aux_lasso[,2] = 1 - p_hat_aux_lasso[,1]
    p_hat_aux_xgb[,1] = p_hat_aux_ensemble$fit_full$predictions[,2]
    p_hat_aux_xgb[,2] = 1 - p_hat_aux_xgb[,1]
    p_hat_aux_nn[,1] = p_hat_aux_ensemble$fit_full$predictions[,3]
    p_hat_aux_nn[,2] = 1 - p_hat_aux_nn[,1]

    # run the ensemble with switched data sets
    p_hat_main_ensemble = ensemble(ps_methods, x = x_aux, y = d_aux, xnew = x_main) # predict the propensity score with aux model, main data
    p_hat_main_ens[,1] = p_hat_main_ensemble$ensemble # extract prediction values
    p_hat_main_ens[,2] = 1 - p_hat_main_ens[,1] # retrieve 1-p(x) for aux model and main data
    
    # extract predictions for each ml method
    p_hat_main_lasso[,1] = p_hat_main_ensemble$fit_full$predictions[,1]
    p_hat_main_lasso[,2] = 1 - p_hat_main_lasso[,1]
    p_hat_main_xgb[,1] = p_hat_main_ensemble$fit_full$predictions[,2]
    p_hat_main_xgb[,2] = 1 - p_hat_main_xgb[,1]
    p_hat_main_nn[,1] = p_hat_main_ensemble$fit_full$predictions[,3]
    p_hat_main_nn[,2] = 1 - p_hat_main_nn[,1]
    
    # estimate the potential outcome 
    y_hat_main_ens = matrix(NA, nrow(x), 2) # create empty matrix to store values
    y_hat_aux_ens = matrix(NA, nrow(x), 2)
    
    y_hat_main_lasso = matrix(NA, nrow(x), 2)
    y_hat_aux_lasso = matrix(NA, nrow(x), 2)
    y_hat_main_xgb = matrix(NA, nrow(x), 2)
    y_hat_aux_xgb = matrix(NA, nrow(x), 2)
    y_hat_main_nn = matrix(NA, nrow(x), 2)
    y_hat_aux_nn = matrix(NA, nrow(x), 2)
    
    for (i in 1:2){
      # only use either treated or non-treated as groups (i = 1,2)
      # derive the model with the main data set and predict the outcomes with the auxiliary data set
      # the prediction will be made with data of both treated and controls to receive potential outcomes
      y_hat_aux_ensemble = ensemble(oc_methods, x = x_main[treatment_main[,i] == 1, ], y = y_main[treatment_main[,i] == 1], xnew = x_aux)
      y_hat_aux_ens[,i] = y_hat_aux_ensemble$ensemble
      # extract predictions for each ml method
      y_hat_aux_lasso[,i] = y_hat_aux_ensemble$fit_full$predictions[,1]
      y_hat_aux_xgb[,i] = y_hat_aux_ensemble$fit_full$predictions[,2]
      y_hat_aux_nn[,i] = y_hat_aux_ensemble$fit_full$predictions[,3]
      
      # same procedure as above but with switched data sets for cross-fitting purposes
      y_hat_main_ensemble = ensemble(oc_methods, x = x_aux[treatment_aux[,i] == 1, ], y = y_aux[treatment_aux[,i] == 1], xnew = x_main)
      y_hat_main_ens[,i] = y_hat_main_ensemble$ensemble
      # extract predictions for each ml method
      y_hat_main_lasso[,i] = y_hat_main_ensemble$fit_full$predictions[,1]
      y_hat_main_xgb[,i] = y_hat_main_ensemble$fit_full$predictions[,2]
      y_hat_main_nn[,i] = y_hat_main_ensemble$fit_full$predictions[,3]
      
    }
    
    # calculate efficient score E[(y - mu_hat)/p(x) + mu_hat]
    # both scores combined fulfill Neyman orthogonality condition
    # this doubly robust estimator remains consistent with either p(x) or mu(x) not precisely specified
    # see Chernozhukov et al (2018), Hahn (1998), Robins and Rotnitzky (1995)
    mu_mat_main_ens = matrix(NA, nrow(x), 2)
    mu_mat_aux_ens = matrix(NA, nrow(x), 2)
    ipw_mat_main_ens = matrix(NA, nrow(x), 2)
    ipw_mat_aux_ens = matrix(NA, nrow(x), 2)
    
    mu_mat_main_lasso = matrix(NA, nrow(x), 2)
    mu_mat_main_xgb = matrix(NA, nrow(x), 2)
    mu_mat_main_nn = matrix(NA, nrow(x), 2)
    mu_mat_aux_lasso = matrix(NA, nrow(x), 2)
    mu_mat_aux_xgb = matrix(NA, nrow(x), 2)
    mu_mat_aux_nn = matrix(NA, nrow(x), 2)
    ipw_mat_main_lasso = matrix(NA, nrow(x), 2)
    ipw_mat_main_xgb = matrix(NA, nrow(x), 2)
    ipw_mat_main_nn = matrix(NA, nrow(x), 2)
    ipw_mat_aux_lasso = matrix(NA, nrow(x), 2)
    ipw_mat_aux_xgb = matrix(NA, nrow(x), 2)
    ipw_mat_aux_nn = matrix(NA, nrow(x), 2)
    
    # calculate the inverse probability weights d/p(x) and (1-d)/(1-p(x))
    for (i in 1:2){
      # ensemble predictions
      ipw_mat_main_ens[,i] = treatment_main[,i] / p_hat_main_ens[,i]
      # predictions of each ml method
      ipw_mat_main_lasso[,i] = treatment_main[,i] / p_hat_main_lasso[,i]
      ipw_mat_main_xgb[,i] = treatment_main[,i] / p_hat_main_xgb[,i]
      ipw_mat_main_nn[,i] = treatment_main[,i] / p_hat_main_nn[,i]
      # ensemble prediction switched data sets
      ipw_mat_aux_ens[,i] = treatment_aux_ens[,i] / p_hat_aux_ens[,i]
      # predictions of each ml method
      ipw_mat_aux_lasso[,i] = treatment_aux[,i] / p_hat_aux_lasso[,i]
      ipw_mat_aux_xgb[,i] = treatment_aux[,i] / p_hat_aux_xgb[,i]
      ipw_mat_aux_nn[,i] = treatment_aux[,i] / p_hat_aux_nn[,i]
    }
    
    # finalize the score calculation
    for (i in 1:2){
      # ensemble prediction
      mu_mat_main_ens[,i] = (y_main - y_hat_main_ens[,i])*ipw_mat_main_ens[,i] + y_hat_main_ens[,i]
      # predictions for each ml method
      mu_mat_main_lasso[,i] = (y_main - y_hat_main_lasso[,i])*ipw_mat_main_lasso[,i] + y_hat_main_lasso[,i]
      mu_mat_main_xgb[,i] = (y_main - y_hat_main_xgb[,i])*ipw_mat_main_xgb[,i] + y_hat_main_xgb[,i]
      mu_mat_main_nn[,i] = (y_main - y_hat_main_nn[,i])*ipw_mat_main_nn[,i] + y_hat_main_nn[,i]
      # predictions of ensemble switched data sets
      mu_mat_aux_ens[,i] = (y_aux - y_hat_aux_ens[,i])*ipw_mat_aux_ens[,i] + y_hat_aux_ens[,i]
      # predictions of each ml method
      mu_mat_aux_lasso[,i] = (y_aux - y_hat_aux_lasso[,i])*ipw_mat_aux_lasso[,i] + y_hat_aux_lasso[,i]
      mu_mat_aux_xgb[,i] = (y_aux - y_hat_aux_xgb[,i])*ipw_mat_aux_xgb[,i] + y_hat_aux_xgb[,i]
      mu_mat_aux_nn[,i] = (y_aux - y_hat_aux_nn[,i])*ipw_mat_aux_nn[,i] + y_hat_aux_nn[,i]
    }
    
    # calculate psi(D, theta, eta) = (g(1,X) - g(0,X)) + D(Y - g(1,X))/p(x) + (1-D)(Y - g(0,X))/(1-p(x)) - theta
    # which with moment and orthogonality condition fulfilled the difference between the two scores equals the treatment effect theta
    psi_main_ens = mu_mat_main_ens[,1] - mu_mat_main_ens[,2]
    psi_main_lasso = mu_mat_main_lasso[,1] - mu_mat_main_lasso[,2]
    psi_main_xgb = mu_mat_main_xgb[,1] - mu_mat_main_xgb[,2]
    psi_main_nn = mu_mat_main_nn[,1] - mu_mat_main_nn[,2]
    
    psi_aux_ens = mu_mat_aux_ens[,1] - mu_mat_aux_ens[,2]
    psi_aux_lasso = mu_mat_aux_lasso[,1] - mu_mat_aux_lasso[,2]
    psi_aux_xgb = mu_mat_aux_xgb[,1] - mu_mat_aux_xgb[,2]
    psi_aux_nn = mu_mat_aux_nn[,1] - mu_mat_aux_nn[,2]
    
    # average all the treatment effects to get the ATE
    ate_main_ens = mean(psi_main_ens)
    ate_main_lasso = mean(psi_main_lasso)
    ate_main_xgb = mean(psi_main_xgb)
    ate_main_nn = mean(psi_main_nn)
    
    ate_aux_ens = mean(psi_aux_ens)
    ate_aux_lasso = mean(psi_aux_lasso)
    ate_aux_xgb = mean(psi_aux_xgb)
    ate_aux_nn = mean(psi_aux_nn)
    
    # get the cross-fitted average treatment effect
    ate_ens = mean(c(ate_main_ens, ate_aux_ens))
    
    ate_lasso = mean(c(ate_main_lasso, ate_aux_lasso))
    ate_xgb = mean(c(ate_main_xgb, ate_aux_xgb))
    ate_nn = mean(c(ate_main_nn, ate_aux_nn))
    
    # extract weights of ensemble
    w_ens_ps = colMeans(rbind(p_hat_main_ensemble$nnls_weights, p_hat_aux_ensemble$nnls_weights))
    w_ens_oc = colMeans(rbind(y_hat_main_ensemble$nnls_weights, y_hat_aux_ensemble$nnls_weights))
    
    # define single output
    output = list("ate_ens" = ate_ens, "ate_lasso" = ate_lasso, "ate_xgb" = ate_xgb, "ate_nn" = ate_nn, "w_ens_ps" = w_ens_ps, "w_ens_oc" = w_ens_oc)
    
    return(output)
    
  } else {
    print("Non-binary treatments, please check your data")
  }
}