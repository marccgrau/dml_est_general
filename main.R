#' Main file of master thesis
#' sources all additional code
#' starts the whole simulation and creates all graphical outputs
#' 
#' parameters necessary for Monte Carlo simulation
#' @param n_simulations number of simulations for the Monte Carlo study
#' 
#' parameters necessary for data generating processes
#' @param n_covariates number of covariates
#' @param n_observations number of rows in each simulation round
#' @param Y outcome variable
#' @param D treatment variable
#' @param X confounders
#' @param beta coefficients of controls
#'
#' input parameters necessary for ensemble and DML
#' @param cv_folds folds used for cross-fitting in ensemble
#' @param oc_methods ml methods applied for outcome nuisance parameter
#' @param ps_methods ml methods applied for propensity score parameter
#' 
#' hyperparameter tuning
#' @param params_lasso definition of hyperparameters to use within simulation for lasso
#' @param params_xgb definition of hyperparameters to use within simulation for xgboost
#' @param params_nn definition of hyperparameters to use within simulation for neural network
#' 
#' Results
#' @param theta vector of estimated average treatment effects of each simulation round
#' @param avg_effect average of all thetas
#' @param oc_ensemble_weights weights of each ml method used in ensemble for outcome nuisance parameter
#' @param ps_ensemble_weights weights of each ml method used in ensemble for propensity score nuisance parameter
#' 
#' The main file calls every procedure necessary for the Monte Carlo Study. 
#' It is separated into the following parts
#' 1. Clean workspace, setwd, load all packages and functions
#' 2. Definition of necessary parameters
#' 3. Case 1
#'     a. Evaluation of hyperparameters
#'     b. Monte Carlo simulation
#'     c. Storage of results
#' 4. Case 2
#'     a. Evaluation of hyperparameters
#'     b. Monte Carlo simulation
#'     c. Storage of results
#' 5. Case 3
#'     a. Evaluation of hyperparameters
#'     b. Monte Carlo simulation
#'     c. Storage of results
#' 6. Results
#'     a. Preparation/Processing
#'     b. Graphical Illustration


# Preliminaries -----------------------------------------------------------

## Load necessary packages, set working directory and seed, remove previously stored variables

toload <- c("grf", "tidyverse", "hdm", "glmnet", "nnls", "Matrix", "matrixStats", "xgboost", "neuralnet")
toinstall <- toload[which(toload %in% installed.packages()[,1] == F)]
lapply(toinstall, install.packages, character.only = TRUE)
lapply(toload, require, character.only = TRUE)

directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

rm(list = ls())
set.seed(12345)

## Load source functions
# General Functions
source("general_functions/general_utils.R")

# DGPs
source("DGP/Powers2018_nuisancefunctions.R")
source("DGP/Powers2018_DGPfunctions.R")
source("DGP/Powers2018_generalDGP.R")
source("DGP/DGP1.R")
source("DGP/DGP1_bounded.R")
source("DGP/DGP2.R")
source("DGP/DGP2_bounded.R")
source("DGP/DGP3.R")
source("DGP/DGP3_bounded.R")

# Hyperparameter Tuning
source("hyperparam_tuning/hyperparam_lasso_ps.R")
source("hyperparam_tuning/hyperparam_lasso_oc.R")
source("hyperparam_tuning/hyperparam_xgboost_ps.R")
source("hyperparam_tuning/hyperparam_xgboost_oc.R")
source("hyperparam_tuning/hyperparam_nnet_ps.R")
source("hyperparam_tuning/hyperparam_nnet_oc.R")

# DML estimator
source("doubleML/dml_est_cf_ensemble.R")

# Ensemble learner
source("ensemble_method/ensemble.R")
source("ensemble_method/ml_wrapper.R")
source("ensemble_method/utils_ensemble.R")


# Parameters --------------------------------------------------------------

### Define necessary parameters
## Monte Carlo Simulation
n_simulations = 3                  # Number of simulation rounds for Monte Carlo Study

## Data
n_covariates = 10                    # Number of confounders
n_observations = 1000               # Number of observations in simulated dataset
effect = 2                          # True value for effect
beta = seq(1, n_covariates, 1)/10   # Coefficients for confounders in DGP

## Ensemble method
cv_folds = 2                        # Number of folds for cross-validation of used ML methods in the ensemble method


# Simulation 1: Linear Case -----------------------------------------------

# Hyperparameter Tuning for DGP 1
## Data simulation for cross-validation of ml methods to select hyperparameters
# Using the general DGP of Powers et al. (2018)
data_cv = generalDGP(n_covariates, n_observations, f3, f5, 1)
Y_cv = data_cv[[1]]
D_cv = data_cv[[2]]
X_cv = data_cv[[3]]

## Lasso hyperparameters
### Lasso hyperparameters are computationally less expensive to estimate
### Nontheless, the approximate region of the best lambda minimzing the deviance is determined before the Monte Carlo simulation
### Potential Outcome: Cross-validated lambda
params_lasso_oc = hyperparam_lasso_oc(Y_cv, X_cv)

### Propensity score: Cross-validated lambda
params_lasso_ps = hyperparam_lasso_ps(D_cv, X_cv)

## XGBoost hyperparameters
### following the idea of: https://towardsdatascience.com/getting-to-an-hyperparameter-tuned-xgboost-model-in-no-time-a9560f8eb54b
### Potential outcome: Random Search Algorithm
params_xgb_oc = hyperparam_xgboost_oc(Y_cv, X_cv, cv_folds)

### Propensity Score: Random Search Algorithm
params_xgb_ps = hyperparam_xgboost_ps(D_cv, X_cv, cv_folds)


## Neural Network Hyperparameters
### Potential outcome: Grid search algorithm
params_nn_oc = hyperparam_nnet_oc(Y_cv, X_cv)

### Propensity score: Grid search algorithm
params_nn_ps = hyperparam_nnet_ps(D_cv, X_cv)


# Setup the ml methods used in the ensemble for the estimation of the nuisance parameters
# ML methods used for propensity score estimation
lasso_ps_1 = create_method("lasso", name = "Lasso_ps_1", args = params_lasso_ps)
xgb_ps_1 = create_method("xgboost", name = "XGBoost_ps_1", args = params_xgb_ps)
nnet_ps_1 = create_method("neural_net", name = "NeuralNet_ps_1", args = params_nn_ps)

# ML methods used for potential outcome estimation
lasso_oc_1 = create_method("lasso", name = "Lasso_oc_1", args = params_lasso_oc)
xgb_oc_1 = create_method("xgboost", name = "XGBoost_oc_1", args = params_xgb_oc)
nnet_oc_1 = create_method("neural_net", name = "NeuralNet_oc_1", args = params_nn_oc)

# list the respective methods for each ensemble
ps_methods_1 = list(lasso_ps_1, xgb_ps_1, nnet_ps_1)
oc_methods_1 = list(lasso_oc_1, xgb_oc_1, nnet_oc_1)

# create empty matrices to fill throughout the simulation
ate_ens = rep(NA, n_simulations)
ate_lasso = rep(NA, n_simulations)
ate_xgb = rep(NA, n_simulations)
ate_nn = rep(NA, n_simulations)

te_ens = matrix(NA, n_observations, n_simulations)
te_lasso = matrix(NA, n_observations, n_simulations)
te_xgb = matrix(NA, n_observations, n_simulations)
te_nn = matrix(NA, n_observations, n_simulations)

se_po_ens = matrix(NA, n_simulations, 2)
se_po_lasso = matrix(NA, n_simulations, 2)
se_po_xgb = matrix(NA, n_simulations, 2)
se_po_nn = matrix(NA, n_simulations, 2)

se_te_ens = rep(NA, n_simulations)
se_te_lasso = rep(NA, n_simulations)
se_te_xgb = rep(NA, n_simulations)
se_te_nn = rep(NA, n_simulations)

oc_ensemble = matrix(NA, n_simulations, length(oc_methods_1))
ps_ensemble = matrix(NA, n_simulations, length(ps_methods_1)) 



for (j in 1:n_simulations) {
  
  # simulate data
  data = generalDGP(n_covariates, n_observations, f3, f5, 1)
  Y = data[[1]]
  D = data[[2]]
  X = data[[3]]
  n_obs = seq(1, nrow(X), 1)
  
  # run the DML estimator, cross-fitting is done within the algorithm, ate and average weights of ensemble are extracted
  dml_estimator  = dml_est_cf_ensemble(Y, D, X, ps_methods_1, oc_methods_1)
  
  # update list of estimates for current simulation round
  ate_ens[j] = dml_estimator$ate_ens                      # estimated effect theta in current simulation round
  ate_lasso[j] = dml_estimator$ate_lasso
  ate_xgb[j] = dml_estimator$ate_xgb
  ate_nn[j] = dml_estimator$ate_nn
  te_ens[,j] = dml_estimator$te_ens
  te_lasso[,j] = dml_estimator$te_lasso
  te_xgb[,j] = dml_estimator$te_xgb
  te_nn[,j] = dml_estimator$te_nn
  se_po_ens[j,] = dml_estimator$se_mu[1,]
  se_po_lasso[j,] = dml_estimator$se_mu[2,]
  se_po_xgb[j,] = dml_estimator$se_mu[3,]
  se_po_nn[j,] = dml_estimator$se_mu[4,]
  se_te_ens[j] = dml_estimator$se_te[1]
  se_te_lasso[j] = dml_estimator$se_te[2]
  se_te_xgb[j] = dml_estimator$se_te[3]
  se_te_nn[j] = dml_estimator$se_te[4]
  ps_ensemble[j,] = dml_estimator$w_ens_ps           # store weights of current simulation round
  oc_ensemble[j,] = dml_estimator$w_ens_oc          # store weights of current simulation round
  
  print(paste("Simulation round", j))
}

# Averaging over all simulations
# Average treatment effect
avg_effect_ens = mean(ate_ens)                          # average effect over all simulation rounds
avg_effect_lasso = mean(ate_lasso)
avg_effect_xgb = mean(ate_xgb)
avg_effect_nn = mean(ate_nn)

# Ensemble weights of E[Y|X]
oc_ensemble_weights = as.data.frame(t(colMeans(oc_ensemble)))
for (i in 1:length(oc_methods_1)) {
  if (!is.null(oc_methods_1[[i]]$name)) colnames(oc_ensemble_weights)[i] = oc_methods_1[[i]]$name
  oc_ensemble = as.data.frame(oc_ensemble)
  colnames(oc_ensemble)[i] = oc_methods_1[[i]]$name
}


# Ensemble weights of E[D|X]
ps_ensemble_weights = as.data.frame(t(colMeans(ps_ensemble)))
for (i in 1:length(ps_methods_1)) {
  if (!is.null(ps_methods_1[[i]]$name)) colnames(ps_ensemble_weights)[i] = ps_methods_1[[i]]$name
  ps_ensemble = as.data.frame(ps_ensemble)
  colnames(ps_ensemble)[i] = ps_methods_1[[i]]$name
}

# Print the results
paste("Average treatment effect:", round(avg_effect_ens, 3))
paste(sprintf("Ensemble weight E[Y|X] %s:",colnames(oc_ensemble_weights)), round(oc_ensemble_weights, 3))
paste(sprintf("Ensemble weight E[D|X] %s:",colnames(ps_ensemble_weights)), round(ps_ensemble_weights, 3))


