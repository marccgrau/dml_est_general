# load necessary source code for individual hyperparam tunings
source("hyperparam_tuning/hyperparam_lasso.R")
source("hyperparam_tuning/hyperparam_xgboost.R")
source("hyperparam_tuning/randomsearch_xgb.R")
source("hyperparam_tuning/hyperparam_nnet.R")
source("hyperparam_tuning/gridframes_nn.R")


hyperparam_tuning_1 = function(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma = 1, cv_folds = 2, w = 0){
  ## Data simulation for cross-validation of ml methods to select hyperparameters
  # Using the general DGP of Powers et al. (2018)
  data = generalDGP(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma, w)
  Y = data[[1]]
  D = data[[2]]
  X = data[[3]]
  
  ## Lasso hyperparameters
  ### Lasso hyperparameters are computationally less expensive to estimate
  ### Nontheless, the approximate region of the best lambda minimzing the deviance is determined before the Monte Carlo simulation
  ### Potential Outcome: Cross-validated lambda
  params_lasso_oc = list(family = "gaussian")
  
  ### Propensity score: Cross-validated lambda
  params_lasso_ps = list(family = "binomial")
  
  ## XGBoost hyperparameters
  ### Potential outcome: Random Search Algorithm
  # get random search grid for DGP 1
  grid_xgb_oc = randomsearch_xgb_oc_1()
  params_xgb_oc = hyperparam_xgboost(Y, X, cv_folds, grid_xgb_oc, ps = FALSE)
  
  ### Propensity Score: Random Search Algorithm
  grid_xgb_ps = randomsearch_xgb_ps_1()
  params_xgb_ps = hyperparam_xgboost(D, X, cv_folds, grid_xgb_ps, ps = TRUE)
  
  
  ## Neural Network Hyperparameters
  ### Potential outcome: Grid search algorithm
  grid_nn_oc = gridframe_oc_1()
  params_nn_oc = hyperparam_nnet(Y, X, grid_nn_oc)
  
  ### Propensity score: Grid search algorithm
  grid_nn_ps = gridframe_ps_1()
  params_nn_ps = hyperparam_nnet(D, X, grid_nn_ps)
  
  # Setup the ml methods used in the ensemble for the estimation of the nuisance parameters
  # ML methods used for propensity score estimation
  lasso_ps = create_method("lasso_inter", name = "Lasso_ps_1", args = params_lasso_ps)
  xgb_ps = create_method("xgboost", name = "XGBoost_ps_1", args = params_xgb_ps)
  nnet_ps = create_method("neural_net", name = "NeuralNet_ps_1", args = params_nn_ps)
  
  # ML methods used for potential outcome estimation
  lasso_oc = create_method("lasso_inter", name = "Lasso_oc_1", args = params_lasso_oc)
  xgb_oc = create_method("xgboost", name = "XGBoost_oc_1", args = params_xgb_oc)
  nnet_oc = create_method("neural_net", name = "NeuralNet_oc_1", args = params_nn_oc)
  
  # list the respective methods for each ensemble
  ps_methods = list(lasso_ps, xgb_ps, nnet_ps)
  oc_methods = list(lasso_oc, xgb_oc, nnet_oc)
  
  output = list("ps_methods" = ps_methods, "oc_methods" = oc_methods)
  
  return(output)
}

hyperparam_tuning_2 = function(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma = 1, cv_folds = 2, w = 0){
  ## Data simulation for cross-validation of ml methods to select hyperparameters
  # Using the general DGP of Powers et al. (2018)
  data = generalDGP(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma, w)
  Y = data[[1]]
  D = data[[2]]
  X = data[[3]]
  
  ## Lasso hyperparameters
  ### Lasso hyperparameters are computationally less expensive to estimate
  ### Nontheless, the approximate region of the best lambda minimzing the deviance is determined before the Monte Carlo simulation
  ### Potential Outcome: Cross-validated lambda
  params_lasso_oc = hyperparam_lasso(Y, X, family = "gaussian")
  
  ### Propensity score: Cross-validated lambda
  params_lasso_ps = hyperparam_lasso(D, X, family = "binomial")
  
  ## XGBoost hyperparameters
  ### Potential outcome: Random Search Algorithm
  # get random search grid for DGP 1
  grid_xgb_oc = randomsearch_xgb_oc_2()
  params_xgb_oc = hyperparam_xgboost(Y, X, cv_folds, grid_xgb_oc, ps = FALSE)
  
  ### Propensity Score: Random Search Algorithm
  grid_xgb_ps = randomsearch_xgb_ps_2()
  params_xgb_ps = hyperparam_xgboost(D, X, cv_folds, grid_xgb_ps, ps = TRUE)
  
  
  ## Neural Network Hyperparameters
  ### Potential outcome: Grid search algorithm
  grid_nn_oc = gridframe_oc_2()
  params_nn_oc = hyperparam_nnet(Y, X, grid_nn_oc)
  
  ### Propensity score: Grid search algorithm
  grid_nn_ps = gridframe_ps_2()
  params_nn_ps = hyperparam_nnet(D, X, grid_nn_ps)
  
  # Setup the ml methods used in the ensemble for the estimation of the nuisance parameters
  # ML methods used for propensity score estimation
  lasso_ps = create_method("lasso_inter", name = "Lasso_ps_1", args = params_lasso_ps)
  xgb_ps = create_method("xgboost", name = "XGBoost_ps_1", args = params_xgb_ps)
  nnet_ps = create_method("neural_net", name = "NeuralNet_ps_1", args = params_nn_ps)
  
  # ML methods used for potential outcome estimation
  lasso_oc = create_method("lasso_inter", name = "Lasso_oc_1", args = params_lasso_oc)
  xgb_oc = create_method("xgboost", name = "XGBoost_oc_1", args = params_xgb_oc)
  nnet_oc = create_method("neural_net", name = "NeuralNet_oc_1", args = params_nn_oc)
  
  # list the respective methods for each ensemble
  ps_methods = list(lasso_ps, xgb_ps, nnet_ps)
  oc_methods = list(lasso_oc, xgb_oc, nnet_oc)
  
  output = list("ps_methods" = ps_methods, "oc_methods" = oc_methods)
  
  return(output)
}

hyperparam_tuning_3 = function(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma = 1, cv_folds = 2, w = 0){
  ## Data simulation for cross-validation of ml methods to select hyperparameters
  # Using the general DGP of Powers et al. (2018)
  data = generalDGP(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma, w)
  Y = data[[1]]
  D = data[[2]]
  X = data[[3]]
  
  ## Lasso hyperparameters
  ### Lasso hyperparameters are computationally less expensive to estimate
  ### Nontheless, the approximate region of the best lambda minimzing the deviance is determined before the Monte Carlo simulation
  ### Potential Outcome: Cross-validated lambda
  params_lasso_oc = hyperparam_lasso(Y, X, family = "gaussian")
  
  ### Propensity score: Cross-validated lambda
  params_lasso_ps = hyperparam_lasso(D, X, family = "binomial")
  
  ## XGBoost hyperparameters
  ### Potential outcome: Random Search Algorithm
  # get random search grid for DGP 1
  grid_xgb_oc = randomsearch_xgb_oc_3()
  params_xgb_oc = hyperparam_xgboost(Y, X, cv_folds, grid_xgb_oc, ps = FALSE)
  
  ### Propensity Score: Random Search Algorithm
  grid_xgb_ps = randomsearch_xgb_ps_3()
  params_xgb_ps = hyperparam_xgboost(D, X, cv_folds, grid_xgb_ps, ps = TRUE)
  
  
  ## Neural Network Hyperparameters
  ### Potential outcome: Grid search algorithm
  grid_nn_oc = gridframe_oc_3()
  params_nn_oc = hyperparam_nnet(Y, X, grid_nn_oc)
  
  ### Propensity score: Grid search algorithm
  grid_nn_ps = gridframe_ps_3()
  params_nn_ps = hyperparam_nnet(D, X, grid_nn_ps)
  
  # Setup the ml methods used in the ensemble for the estimation of the nuisance parameters
  # ML methods used for propensity score estimation
  lasso_ps = create_method("lasso_inter", name = "Lasso_ps_1", args = params_lasso_ps)
  xgb_ps = create_method("xgboost", name = "XGBoost_ps_1", args = params_xgb_ps)
  nnet_ps = create_method("neural_net", name = "NeuralNet_ps_1", args = params_nn_ps)
  
  # ML methods used for potential outcome estimation
  lasso_oc = create_method("lasso_inter", name = "Lasso_oc_1", args = params_lasso_oc)
  xgb_oc = create_method("xgboost", name = "XGBoost_oc_1", args = params_xgb_oc)
  nnet_oc = create_method("neural_net", name = "NeuralNet_oc_1", args = params_nn_oc)
  
  # list the respective methods for each ensemble
  ps_methods = list(lasso_ps, xgb_ps, nnet_ps)
  oc_methods = list(lasso_oc, xgb_oc, nnet_oc)
  
  output = list("ps_methods" = ps_methods, "oc_methods" = oc_methods)
  
  return(output)
}