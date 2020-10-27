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
    
    # estimations for the average treatment effect ate
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
      ipw_mat_aux_ens[,i] = treatment_aux[,i] / p_hat_aux_ens[,i]
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
    
    # calculate standard errors SE
    se_mu = matrix(NA, nrow = (length(oc_methods) + 1), ncol = ncol(mu_mat_main_ens))

    for (i in 1:ncol(mu_mat_main_ens)){
      se_mu[1,i] = sqrt(sum((rowMeans(cbind(mu_mat_main_ens[,i], mu_mat_aux_ens[,i])) - mean(rowMeans(cbind(mu_mat_main_ens[,i], mu_mat_aux_ens[,i]))))^2) / length(mu_mat_main_ens[,i]))
      se_mu[2,i] = sqrt(sum((rowMeans(cbind(mu_mat_main_lasso[,i], mu_mat_aux_lasso[,i])) - mean(rowMeans(cbind(mu_mat_main_lasso[,i], mu_mat_aux_lasso[,i]))))^2) / length(mu_mat_main_lasso[,i]))
      se_mu[3,i] = sqrt(sum((rowMeans(cbind(mu_mat_main_xgb[,i], mu_mat_aux_xgb[,i])) - mean(rowMeans(cbind(mu_mat_main_xgb[,i], mu_mat_aux_xgb[,i]))))^2) / length(mu_mat_main_lasso[,i]))
      se_mu[4,i] = sqrt(sum((rowMeans(cbind(mu_mat_main_nn[,i], mu_mat_aux_nn[,i])) - mean(rowMeans(cbind(mu_mat_main_nn[,i], mu_mat_aux_nn[,i]))))^2) / length(mu_mat_main_lasso[,i]))
    }
    
    colnames(se_mu) = c("mu_1", "mu_0")
    rownames(se_mu) = c("ensemble", "lasso", "xgb", "nn")
    
    
    # calculate te(D, theta, eta) = (g(1,X) - g(0,X)) + D(Y - g(1,X))/p(x) + (1-D)(Y - g(0,X))/(1-p(x)) - theta
    # which with moment and orthogonality condition fulfilled the difference between the two scores equals the treatment effect theta
    te_main_ens = mu_mat_main_ens[,1] - mu_mat_main_ens[,2]
    te_main_lasso = mu_mat_main_lasso[,1] - mu_mat_main_lasso[,2]
    te_main_xgb = mu_mat_main_xgb[,1] - mu_mat_main_xgb[,2]
    te_main_nn = mu_mat_main_nn[,1] - mu_mat_main_nn[,2]
    
    te_aux_ens = mu_mat_aux_ens[,1] - mu_mat_aux_ens[,2]
    te_aux_lasso = mu_mat_aux_lasso[,1] - mu_mat_aux_lasso[,2]
    te_aux_xgb = mu_mat_aux_xgb[,1] - mu_mat_aux_xgb[,2]
    te_aux_nn = mu_mat_aux_nn[,1] - mu_mat_aux_nn[,2]
    
    te_ens = rowMeans(cbind(te_main_ens, te_aux_ens))
    te_lasso = rowMeans(cbind(te_main_lasso, te_aux_lasso))
    te_xgb = rowMeans(cbind(te_main_xgb, te_aux_xgb))
    te_nn = rowMeans(cbind(te_main_nn, te_aux_nn))
    
    se_te = matrix(NA, nrow = (length(oc_methods) + 1), ncol = 1)
    
    se_te[1,] = sqrt(sum((te_ens - mean(te_ens))^2) / length(te_ens))
    se_te[2,] = sqrt(sum((te_lasso - mean(te_lasso))^2) / length(te_lasso))
    se_te[3,] = sqrt(sum((te_xgb - mean(te_xgb))^2) / length(te_xgb))
    se_te[4,] = sqrt(sum((te_nn - mean(te_nn))^2) / length(te_nn))
    
    colnames(se_te) = c("se_te")
    rownames(se_te) = c("ensemble", "lasso", "xgb", "nn")
    
    # get the cross-fitted average treatment effect
    ate_ens = mean(te_ens)
    ate_lasso = mean(te_lasso)
    ate_xgb = mean(te_xgb)
    ate_nn = mean(te_nn)
    
    # extract weights of ensemble
    w_ens_ps = colMeans(rbind(p_hat_main_ensemble$nnls_weights, p_hat_aux_ensemble$nnls_weights))
    w_ens_oc = colMeans(rbind(y_hat_main_ensemble$nnls_weights, y_hat_aux_ensemble$nnls_weights))
    
    # define single output
    output = list("ate_ens" = ate_ens, "ate_lasso" = ate_lasso, "ate_xgb" = ate_xgb, "ate_nn" = ate_nn,
                  "te_ens" = te_ens, "te_lasso" = te_lasso, "te_xgb" = te_xgb, "te_nn" = te_nn, "w_ens_ps" = w_ens_ps, "w_ens_oc" = w_ens_oc, 
                  "se_te" = se_te, "se_po" = se_mu)
    
    return(output)
    
  } else {
    print("Non-binary treatments, please check your data")
  }
}