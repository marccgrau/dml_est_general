hyperparam_tuning_1 = function(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma, cv_folds, w){
  ## Data simulation for cross-validation of ml methods to select hyperparameters
  # Using the general DGP of Powers et al. (2018)
  data_cv = generalDGP(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma, w)
  Y_cv = data_cv[[1]]
  D_cv = data_cv[[2]]
  X_cv = data_cv[[3]]
  
  ## Lasso hyperparameters
  ### Lasso hyperparameters are computationally less expensive to estimate
  ### Nontheless, the approximate region of the best lambda minimzing the deviance is determined before the Monte Carlo simulation
  ### Potential Outcome: Cross-validated lambda
  params_lasso_oc = hyperparam_lasso_oc_1(Y_cv, X_cv)
  
  ### Propensity score: Cross-validated lambda
  params_lasso_ps = hyperparam_lasso_ps_1(D_cv, X_cv)
  
  ## XGBoost hyperparameters
  ### following the idea of: https://towardsdatascience.com/getting-to-an-hyperparameter-tuned-xgboost-model-in-no-time-a9560f8eb54b
  ### Potential outcome: Random Search Algorithm
  params_xgb_oc = hyperparam_xgboost_oc_1(Y_cv, X_cv, cv_folds)
  
  ### Propensity Score: Random Search Algorithm
  params_xgb_ps = hyperparam_xgboost_ps_1(D_cv, X_cv, cv_folds)
  
  
  ## Neural Network Hyperparameters
  ### Potential outcome: Grid search algorithm
  params_nn_oc = hyperparam_nnet_oc_1(Y_cv, X_cv)
  
  ### Propensity score: Grid search algorithm
  params_nn_ps = hyperparam_nnet_ps_1(D_cv, X_cv)
  
  
  # Setup the ml methods used in the ensemble for the estimation of the nuisance parameters
  # ML methods used for propensity score estimation
  lasso_ps = create_method("lasso", name = "Lasso_ps_1", args = params_lasso_ps)
  xgb_ps = create_method("xgboost", name = "XGBoost_ps_1", args = params_xgb_ps)
  nnet_ps = create_method("neural_net", name = "NeuralNet_ps_1", args = params_nn_ps)
  
  # ML methods used for potential outcome estimation
  lasso_oc = create_method("lasso", name = "Lasso_oc_1", args = params_lasso_oc)
  xgb_oc = create_method("xgboost", name = "XGBoost_oc_1", args = params_xgb_oc)
  nnet_oc = create_method("neural_net", name = "NeuralNet_oc_1", args = params_nn_oc)
  
  # list the respective methods for each ensemble
  ps_methods = list(lasso_ps, xgb_ps, nnet_ps)
  oc_methods = list(lasso_oc, xgb_oc, nnet_oc)
  
  output = list("ps_methods" = ps_methods, "oc_methods" = oc_methods)
  
  return(output)
}