#' Starts the whole simulation and creates all graphical outputs
#' 
#' parameters necessary for Monte Carlo simulation
#' @param n_simulations number of simulations for the Monte Carlo study
#' 
#' parameters necessary for data generating processes
#' @param n_covariates number of covariates
#'
#' parameters necessary for DML estimator
#' @param k_folds folds used for cross-fitting


# Preliminaries -----------------------------------------------------------

## Load necessary packages, set working directory and seed, remove previously stored variables

toload <- c("grf", "tidyverse", "hdm", "glmnet", "nnls", "Matrix", "matrixStats", "xgboost", "neuralnet")
toinstall <- toload[which(toload %in% installed.packages()[,1] == F)]
lapply(toinstall, install.packages, character.only = TRUE)
lapply(toload, require, character.only = TRUE)

directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

rm(list = ls())
set.seed(123)

## Load source functions
# General Functions
source("general_functions/general_utils.R")

# DGPs
source("DGP/DGP1.R")
source("DGP/DGP2.R")

# Hyperparameter Tuning
source("hyperparam_tuning/hyperparam_lasso.R")
source("hyperparam_tuning/hyperparam_xgboost.R")
source("hyperparam_tuning/hyperparam_nnet.R")

# DML estimator
source("doubleML/dml_est.R")
source("doubleML/dml_est_cf.R")

# Ensemble learner
source("ensemble_method/ensemble.R")
source("ensemble_method/ml_wrapper.R")
source("ensemble_method/utils_ensemble.R")


# Parameters --------------------------------------------------------------

### Define necessary parameters
## Monte Carlo Simulation
n_simulations = 10                  # Number of simulation rounds for Monte Carlo Study

## Data
n_covariates = 10                    # Number of confounders
n_observations = 100000              # Number of observations in simulated dataset
effect = 10                          # True value for effect
beta = seq(1, n_covariates, 1)/10   # Coefficients for confounders in DGP

## Ensemble method
cv_folds = 2                        # Number of folds for cross-validation of used ML methods in the ensemble method


# Simulation 1: Linear Case -----------------------------------------------

# Hyperparameter Tuning for DGP 1
## Data simulation for cross-validation of ml methods to select hyperparameters
data_cv = DGP1(n_simulations = n_simulations,n_covariates = n_covariates, n_observations = n_observations, beta = beta, effect = effect)
Y_cv = data_cv[[1]]
D_cv = data_cv[[2]]
X_cv = data_cv[[3]]

## Lasso hyperparameters
### Lasso hyperparameters are computationally less expensive to estimate
### Nontheless, the approximate region of the best lambda minimzing the deviance is determined before the Monte Carlo simulation
### Potential Outcome: Cross-validated lambda
lambdas_oc = hyperparam_lasso(Y_cv, X_cv)

### Propensity score: Cross-validated lambda
lambdas_ps = hyperparam_lasso(D_cv, X_cv)

## XGBoost hyperparameters
### following the idea of: https://towardsdatascience.com/getting-to-an-hyperparameter-tuned-xgboost-model-in-no-time-a9560f8eb54b
### Potential outcome: Random Search Algorithm
params_xgb_oc = hyperparam_xgboost(Y_cv, X_cv, cvfold)

### Propensity Score: Random Search Algorithm
params_xgb_ps = hyperparam_xgboost(D_cv, X_cv, cvfold)


## Neural Network Hyperparameters
### Potential outcome: Grid search algorithm
params_nn_oc = hyperparam_nnet(Y_cv_train, X_cv_train, Y_cv_test, X_cv_test)

### Propensity score: Grid search algorithm
params_nn_ps = hyperparam_nnet(D_cv_train, X_cv_train, D_cv_test, X_cv_test)


# Setup the ml methods used in the ensemble for the estimation of the nuisance parameters
# ML methods used for propensity score estimation
lasso_bin_ps_1 = create_method("lasso", name = "Lasso ps 1", args = list(family = "binomial", lambda = lambdas_ps))
xgb_ps_1 = create_method("xgboost", name = "XGBoost ps", args = params_xgb_ps)
nnet_ps_1 = create_method("neural_net", name = "NeuralNet oc", args = list(hidden = c(5), linear.output = FALSE, stepmax = 20000, threshold = 0.4))

# ML methods used for potential outcome estimation
lasso_oc_1 = create_method("lasso", name = "Lasso oc 1", args = list(family = "gaussian", lambdas_oc))
xgb_oc_1 = create_method("xgboost", name = "XGBoost oc", args = params_xgb_oc)
nnet_oc_1 = create_method("neural_net", name = "NeuralNet oc", args = list(hidden = c(5), linear.output = FALSE, stepmax = 20000, threshold = 0.4))

# list the respective methods for each ensemble
ps_methods_1 = list(lasso_bin_ps_1, xgb_ps_1, nnet_ps_1)
oc_methods_1 = list(lasso_oc_1, xgb_oc_1, nnet_oc_1)

# create folds for cross-fitting

#theta_cf = rep(NA, k_folds)
theta = rep(NA, n_simulations)

oc_ensemble = matrix(NA, n_simulations, length(oc_methods_1))
ps_ensemble = matrix(NA, n_simulations, length(ps_methods_1)) 

for (j in 1:n_simulations) {
  
  # simulate data
  data = DGP1(n_simulations = n_simulations,n_covariates = n_covariates, n_observations = n_observations, beta = beta, effect = effect)
  Y = data[[1]]
  D = data[[2]]
  X = data[[3]]
  n_obs = seq(1, nrow(X), 1)
  
  # run the DML estimator, cross-fitting is done within the algorithm, ate and average weights of ensemble are extracted
  dml_estimator  = dml_est_cf_ensemble(Y, D, X, ps_methods_1, oc_methods_1)
  theta_est = dml_estimator$ate                   # extract the average treatment effect
  weights_ensemble_ps = dml_estimator$w_ens_ps    # extract the ensemble weights for each ml method for the propensity score
  weights_ensemble_oc = dml_estimator$w_ens_oc    # extract the ensemble weights for each ml method for the potential outcome
  
  # update list of estimates for current simulation round
  theta[j] = theta_est                            # estimated effect theta in current simulation round
  ps_ensemble[j,] = weights_ensemble_ps           # store weights of current simulation round
  oc_ensemble[j,] = weights_ensemble_oc           # store weights of current simulation round
  
}

# Averaging over all simulations
# Average treatment effect
avg_effect = mean(theta)                            # average effect over all simulation rounds

# Ensemble weights of E[Y|X]
oc_ensemble_weights = as.data.frame(t(colMeans(oc_ensemble)))
for (i in 1:length(oc_methods)) {
  if (!is.null(oc_methods[[i]]$name)) colnames(oc_ensemble_weights)[i] = oc_methods[[i]]$name
  oc_ensemble = as.data.frame(oc_ensemble)
  colnames(oc_ensemble)[i] = oc_methods[[i]]$name
}


# Ensemble weights of E[D|X]
ps_ensemble_weights = as.data.frame(t(colMeans(ps_ensemble)))
for (i in 1:length(ps_methods)) {
  if (!is.null(ps_methods[[i]]$name)) colnames(ps_ensemble_weights)[i] = ps_methods[[i]]$name
  ps_ensemble = as.data.frame(ps_ensemble)
  colnames(ps_ensemble)[i] = ps_methods[[i]]$name
}

# Print the results
paste("Average treatment effect:", round(avg_effect, 3))
paste(sprintf("Ensemble weight E[Y|X] %s:",colnames(oc_ensemble_weights)), round(oc_ensemble_weights, 3))
paste(sprintf("Ensemble weight E[Y|X] %s:",colnames(ps_ensemble_weights)), round(ps_ensemble_weights, 3))


