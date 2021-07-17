# load necessary source code for individual hyperparam tunings
source("hyperparam_tuning/hyperparam_lasso.R")
source("hyperparam_tuning/hyperparam_xgboost.R")
source("hyperparam_tuning/randomsearch_xgb.R")
source("hyperparam_tuning/hyperparam_nnet_keras.R")
source("hyperparam_tuning/gridframes_nn_keras.R")


hyperparam_tuning_1_50 = function(n_covariates, n_observations, mufunc, taufunc, psfunc, 
                               sigma = 1, cv_folds = 2, smalltreat = FALSE){
  ## Data simulation for cross-validation of ml methods to select hyperparameters
  # Using the general DGP of Powers et al. (2018)
  data = generalDGP(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma, smalltreat)
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
  grid_xgb_oc = randomsearch_xgb_oc()
  params_xgb_oc = hyperparam_xgboost(Y, X, cv_folds, grid_xgb_oc, ps = FALSE)
  
  ### Propensity Score: Random Search Algorithm
  grid_xgb_ps = randomsearch_xgb_ps()
  params_xgb_ps = hyperparam_xgboost(D, X, cv_folds, grid_xgb_ps, ps = TRUE)
  
  
  ## Neural Network Hyperparameters
  ### Potential outcome: Grid search algorithm
  grid_nn_oc = grid_keras_oc_1_50()
  params_nn_oc = hyperparam_nnet_keras(Y, X, grid_nn_oc, twolayers = TRUE, optimizer = "adamax", learning_rate = 0.01)
  
  ### Propensity score: Grid search algorithm
  grid_nn_ps = grid_keras_ps_1_50()
  params_nn_ps = hyperparam_nnet_keras(D, X, grid_nn_ps, twolayers = FALSE, optimizer = "adam", learning_rate = 0.01)
  
  # Setup the ml methods used in the ensemble for the estimation of the nuisance parameters
  # ML methods used for propensity score estimation
  lasso_ps = create_method("lasso_inter", name = "Lasso_ps", args = params_lasso_ps)
  xgb_ps = create_method("xgboost", name = "XGBoost_ps", args = params_xgb_ps)
  nnet_ps = create_method("nn_keras", name = "nn_ps", args = params_nn_ps)
  
  # ML methods used for potential outcome estimation
  lasso_oc = create_method("lasso_inter", name = "Lasso_oc", args = params_lasso_oc)
  xgb_oc = create_method("xgboost", name = "XGBoost_oc", args = params_xgb_oc)
  nnet_oc = create_method("nn_keras", name = "nn_oc", args = params_nn_oc)
  
  # list the respective methods for each ensemble
  ps_methods = list(lasso_ps, xgb_ps, nnet_ps)
  oc_methods = list(lasso_oc, xgb_oc, nnet_oc)
  
  output = list("ps_methods" = ps_methods, "oc_methods" = oc_methods)
  
  return(output)
}

hyperparam_tuning_1_10 = function(n_covariates, n_observations, mufunc, taufunc, psfunc, 
                               sigma = 1, cv_folds = 2, smalltreat = FALSE){
  ## Data simulation for cross-validation of ml methods to select hyperparameters
  # Using the general DGP of Powers et al. (2018)
  data = generalDGP(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma, smalltreat)
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
  # get random search grid for DGP 2
  grid_xgb_oc = randomsearch_xgb_oc()
  params_xgb_oc = hyperparam_xgboost(Y, X, cv_folds, grid_xgb_oc, ps = FALSE)
  
  ### Propensity Score: Random Search Algorithm
  grid_xgb_ps = randomsearch_xgb_ps()
  params_xgb_ps = hyperparam_xgboost(D, X, cv_folds, grid_xgb_ps, ps = TRUE)
  
  
  ## Neural Network Hyperparameters
  ### Potential outcome: Grid search algorithm
  grid_nn_oc = grid_keras_oc_1_10()
  params_nn_oc = hyperparam_nnet_keras(Y, X, grid_nn_oc, twolayers = FALSE, optimizer = "sgd", learning_rate = 0.01)
  
  ### Propensity score: Grid search algorithm
  grid_nn_ps = grid_keras_ps_1_10()
  params_nn_ps = hyperparam_nnet_keras(D, X, grid_nn_ps, twolayers = FALSE, optimizer = "adam", learning_rate = 0.01)
  
  # Setup the ml methods used in the ensemble for the estimation of the nuisance parameters
  # ML methods used for propensity score estimation
  lasso_ps = create_method("lasso_inter", name = "Lasso_ps", args = params_lasso_ps)
  xgb_ps = create_method("xgboost", name = "XGBoost_ps", args = params_xgb_ps)
  nnet_ps = create_method("nn_keras", name = "NeuralNet_ps", args = params_nn_ps)
  
  # ML methods used for potential outcome estimation
  lasso_oc = create_method("lasso_inter", name = "Lasso_oc", args = params_lasso_oc)
  xgb_oc = create_method("xgboost", name = "XGBoost_oc", args = params_xgb_oc)
  nnet_oc = create_method("nn_keras", name = "NeuralNet_oc", args = params_nn_oc)
  
  # list the respective methods for each ensemble
  ps_methods = list(lasso_ps, xgb_ps, nnet_ps)
  oc_methods = list(lasso_oc, xgb_oc, nnet_oc)
  
  output = list("ps_methods" = ps_methods, "oc_methods" = oc_methods)
  
  return(output)
}

hyperparam_tuning_2_50 = function(n_covariates, n_observations, mufunc, taufunc, psfunc, 
                                  sigma = 1, cv_folds = 2, smalltreat = FALSE){
  ## Data simulation for cross-validation of ml methods to select hyperparameters
  # Using the general DGP of Powers et al. (2018)
  data = generalDGP(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma, smalltreat)
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
  grid_xgb_oc = randomsearch_xgb_oc()
  params_xgb_oc = hyperparam_xgboost(Y, X, cv_folds, grid_xgb_oc, ps = FALSE)
  
  ### Propensity Score: Random Search Algorithm
  grid_xgb_ps = randomsearch_xgb_ps()
  params_xgb_ps = hyperparam_xgboost(D, X, cv_folds, grid_xgb_ps, ps = TRUE)
  
  
  ## Neural Network Hyperparameters
  ### Potential outcome: Grid search algorithm
  grid_nn_oc = grid_keras_oc_2_50()
  params_nn_oc = hyperparam_nnet_keras(Y, X, grid_nn_oc, twolayers = TRUE, optimizer = "adamax", learning_rate = 0.01)
  
  ### Propensity score: Grid search algorithm
  grid_nn_ps = grid_keras_ps_2_50()
  params_nn_ps = hyperparam_nnet_keras(D, X, grid_nn_ps, twolayers = FALSE, optimizer = "adamax", learning_rate = 0.01)
  
  # Setup the ml methods used in the ensemble for the estimation of the nuisance parameters
  # ML methods used for propensity score estimation
  lasso_ps = create_method("lasso_inter", name = "Lasso_ps", args = params_lasso_ps)
  xgb_ps = create_method("xgboost", name = "XGBoost_ps", args = params_xgb_ps)
  nnet_ps = create_method("nn_keras", name = "nn_ps", args = params_nn_ps)
  
  # ML methods used for potential outcome estimation
  lasso_oc = create_method("lasso_inter", name = "Lasso_oc", args = params_lasso_oc)
  xgb_oc = create_method("xgboost", name = "XGBoost_oc", args = params_xgb_oc)
  nnet_oc = create_method("nn_keras", name = "nn_oc", args = params_nn_oc)
  
  # list the respective methods for each ensemble
  ps_methods = list(lasso_ps, xgb_ps, nnet_ps)
  oc_methods = list(lasso_oc, xgb_oc, nnet_oc)
  
  output = list("ps_methods" = ps_methods, "oc_methods" = oc_methods)
  
  return(output)
}

hyperparam_tuning_2_10 = function(n_covariates, n_observations, mufunc, taufunc, psfunc, 
                                  sigma = 1, cv_folds = 2, smalltreat = FALSE){
  ## Data simulation for cross-validation of ml methods to select hyperparameters
  # Using the general DGP of Powers et al. (2018)
  data = generalDGP(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma, smalltreat)
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
  grid_xgb_oc = randomsearch_xgb_oc()
  params_xgb_oc = hyperparam_xgboost(Y, X, cv_folds, grid_xgb_oc, ps = FALSE)
  
  ### Propensity Score: Random Search Algorithm
  grid_xgb_ps = randomsearch_xgb_ps()
  params_xgb_ps = hyperparam_xgboost(D, X, cv_folds, grid_xgb_ps, ps = TRUE)
  
  
  ## Neural Network Hyperparameters
  ### Potential outcome: Grid search algorithm
  grid_nn_oc = grid_keras_oc_2_10()
  params_nn_oc = hyperparam_nnet_keras(Y, X, grid_nn_oc, twolayers = FALSE, optimizer = "adamax", learning_rate = 0.01)
  
  ### Propensity score: Grid search algorithm
  grid_nn_ps = grid_keras_ps_2_10()
  params_nn_ps = hyperparam_nnet_keras(D, X, grid_nn_ps, twolayers = FALSE, optimizer = "adamax", learning_rate = 0.01)
  
  # Setup the ml methods used in the ensemble for the estimation of the nuisance parameters
  # ML methods used for propensity score estimation
  lasso_ps = create_method("lasso_inter", name = "Lasso_ps", args = params_lasso_ps)
  xgb_ps = create_method("xgboost", name = "XGBoost_ps", args = params_xgb_ps)
  nnet_ps = create_method("nn_keras", name = "nn_ps", args = params_nn_ps)
  
  # ML methods used for potential outcome estimation
  lasso_oc = create_method("lasso_inter", name = "Lasso_oc", args = params_lasso_oc)
  xgb_oc = create_method("xgboost", name = "XGBoost_oc", args = params_xgb_oc)
  nnet_oc = create_method("nn_keras", name = "nn_oc", args = params_nn_oc)
  
  # list the respective methods for each ensemble
  ps_methods = list(lasso_ps, xgb_ps, nnet_ps)
  oc_methods = list(lasso_oc, xgb_oc, nnet_oc)
  
  output = list("ps_methods" = ps_methods, "oc_methods" = oc_methods)
  
  return(output)
}

hyperparam_tuning_3_50 = function(n_covariates, n_observations, mufunc, taufunc, psfunc, 
                                  sigma = 1, cv_folds = 2, smalltreat = FALSE){
  ## Data simulation for cross-validation of ml methods to select hyperparameters
  # Using the general DGP of Powers et al. (2018)
  data = generalDGP(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma, smalltreat)
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
  grid_xgb_oc = randomsearch_xgb_oc()
  params_xgb_oc = hyperparam_xgboost(Y, X, cv_folds, grid_xgb_oc, ps = FALSE)
  
  ### Propensity Score: Random Search Algorithm
  grid_xgb_ps = randomsearch_xgb_ps()
  params_xgb_ps = hyperparam_xgboost(D, X, cv_folds, grid_xgb_ps, ps = TRUE)
  
  
  ## Neural Network Hyperparameters
  ### Potential outcome: Grid search algorithm
  grid_nn_oc = grid_keras_oc_3_50()
  params_nn_oc = hyperparam_nnet_keras(Y, X, grid_nn_oc, twolayers = FALSE, optimizer = "adagrad", learning_rate = 0.1)
  
  ### Propensity score: Grid search algorithm
  grid_nn_ps = grid_keras_ps_3_50()
  params_nn_ps = hyperparam_nnet_keras(D, X, grid_nn_ps, twolayers = FALSE, optimizer = "adagrad", learning_rate = 0.1)
  
  # Setup the ml methods used in the ensemble for the estimation of the nuisance parameters
  # ML methods used for propensity score estimation
  lasso_ps = create_method("lasso_inter", name = "Lasso_ps", args = params_lasso_ps)
  xgb_ps = create_method("xgboost", name = "XGBoost_ps", args = params_xgb_ps)
  nnet_ps = create_method("nn_keras", name = "nn_ps", args = params_nn_ps)
  
  # ML methods used for potential outcome estimation
  lasso_oc = create_method("lasso_inter", name = "Lasso_oc", args = params_lasso_oc)
  xgb_oc = create_method("xgboost", name = "XGBoost_oc", args = params_xgb_oc)
  nnet_oc = create_method("nn_keras", name = "nn_oc", args = params_nn_oc)
  
  # list the respective methods for each ensemble
  ps_methods = list(lasso_ps, xgb_ps, nnet_ps)
  oc_methods = list(lasso_oc, xgb_oc, nnet_oc)
  
  output = list("ps_methods" = ps_methods, "oc_methods" = oc_methods)
  
  return(output)
}

hyperparam_tuning_3_10 = function(n_covariates, n_observations, mufunc, taufunc, psfunc, 
                                  sigma = 1, cv_folds = 2, smalltreat = FALSE){
  ## Data simulation for cross-validation of ml methods to select hyperparameters
  # Using the general DGP of Powers et al. (2018)
  data = generalDGP(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma, smalltreat)
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
  grid_xgb_oc = randomsearch_xgb_oc()
  params_xgb_oc = hyperparam_xgboost(Y, X, cv_folds, grid_xgb_oc, ps = FALSE)
  
  ### Propensity Score: Random Search Algorithm
  grid_xgb_ps = randomsearch_xgb_ps()
  params_xgb_ps = hyperparam_xgboost(D, X, cv_folds, grid_xgb_ps, ps = TRUE)
  
  
  ## Neural Network Hyperparameters
  ### Potential outcome: Grid search algorithm
  grid_nn_oc = grid_keras_oc_3_10()
  params_nn_oc = hyperparam_nnet_keras(Y, X, grid_nn_oc, twolayers = FALSE, optimizer = "sgd", learning_rate = 0.01)
  
  ### Propensity score: Grid search algorithm
  grid_nn_ps = grid_keras_ps_3_10()
  params_nn_ps = hyperparam_nnet_keras(D, X, grid_nn_ps, twolayers = TRUE, optimizer = "rmsprop", learning_rate = 0.01)
  
  # Setup the ml methods used in the ensemble for the estimation of the nuisance parameters
  # ML methods used for propensity score estimation
  lasso_ps = create_method("lasso_inter", name = "Lasso_ps", args = params_lasso_ps)
  xgb_ps = create_method("xgboost", name = "XGBoost_ps", args = params_xgb_ps)
  nnet_ps = create_method("nn_keras", name = "nn_ps", args = params_nn_ps)
  
  # ML methods used for potential outcome estimation
  lasso_oc = create_method("lasso_inter", name = "Lasso_oc", args = params_lasso_oc)
  xgb_oc = create_method("xgboost", name = "XGBoost_oc", args = params_xgb_oc)
  nnet_oc = create_method("nn_keras", name = "nn_oc", args = params_nn_oc)
  
  # list the respective methods for each ensemble
  ps_methods = list(lasso_ps, xgb_ps, nnet_ps)
  oc_methods = list(lasso_oc, xgb_oc, nnet_oc)
  
  output = list("ps_methods" = ps_methods, "oc_methods" = oc_methods)
  
  return(output)
}
