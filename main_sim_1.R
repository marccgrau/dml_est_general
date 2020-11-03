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

toload <- c("grf", "tidyverse", "hdm", "glmnet", "nnls", "Matrix", "matrixStats", "xgboost", "neuralnet", "MASS", "MLmetrics", "keras", "tfdatasets")
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
source("DGP/GeneralDGP.R")
source("DGP/nuisance_functions.R")
source("DGP/DGPfunctions.R")

# Hyperparameter Tuning
source("hyperparam_tuning/hyperparam_tuning.R")


# DML estimator
source("doubleML/dml_ens_ml.R")

# Ensemble learner
source("ensemble_method/ensemble.R")
source("ensemble_method/ml_wrapper.R")
source("ensemble_method/utils_ensemble.R")


# Parameters --------------------------------------------------------------

### Define necessary parameters
## Monte Carlo Simulation
n_simulations = 5                  # Number of simulation rounds for Monte Carlo Study

## Data
n_covariates = 15                    # Number of confounders
n_observations = 1000               # Number of observations in simulated dataset

## Ensemble method
cv_folds = 2                        # Number of folds for cross-validation of used ML methods in the ensemble method

# declare which ml methods are going to be used by providing short versions of names
# will be further used to create necessary variables
ml_methods = c("ens", "lasso", "xgb", "nn")

# Simulation 1; 50/50 -----------------------------------------------

# Hyperparameter Tuning for DGP 1
hyperparams_1 = hyperparam_tuning_1(n_covariates, n_observations, mu1, tau1, pi1, sigma = 1, cv_folds, w = 0)
ps_methods_1 = hyperparams_1$ps_methods
oc_methods_1 = hyperparams_1$oc_methods

## create empty matrices to fill throughout the simulation
# average treatment effect
for(i in 1:length(ml_methods)) { 
  nam <- paste("ate_", ml_methods[i], "_1", sep = "")
  assign(nam, rep(NA, n_simulations))
}
# treatment effects
for(i in 1:length(ml_methods)) { 
  nam <- paste("te_", ml_methods[i], "_1", sep = "")
  assign(nam, matrix(NA, n_observations, n_simulations))
}
# standard error potential outcomes
for(i in 1:length(ml_methods)) { 
  nam <- paste("se_po_", ml_methods[i], "_1", sep = "")
  assign(nam, matrix(NA, n_simulations, 2))
}

for(i in 1:length(ml_methods)) { 
  nam <- paste("se_te_", ml_methods[i], "_1", sep = "")
  assign(nam, rep(NA, n_simulations))
}


oc_ensemble_1 = matrix(NA, n_simulations, length(oc_methods_1))
ps_ensemble_1 = matrix(NA, n_simulations, length(ps_methods_1)) 

# empty matrices to store effective values of treatment
te_t = matrix(NA, n_observations, n_simulations)
ate_t = rep(NA, n_simulations)
se_te_t = rep(NA, n_simulations)

ate_ens_1 = matrix(NA, n_simulations,length(ml_methods))
colnames(ate_ens_1) = ml_methods

for (j in 1:n_simulations) {
  ###########
  ## 50/50 ##
  ###########
  # simulate data
  data = generalDGP(n_covariates, n_observations, mu1, tau1, pi1, sigma = 1, w = 0)
  Y = data[[1]]
  D = data[[2]]
  X = data[[3]]
  true_te = data[[4]]
  te_t[,j] = true_te
  ate_t[j] = mean(true_te)
  se_te_t[j] = sd(true_te)
  
  n_obs = seq(1, nrow(X), 1)
  
  # run the DML estimator, cross-fitting is done within the algorithm, ate and average weights of ensemble are extracted
  dml_estimator  = dml_ens_ml(Y, D, X, ps_methods_1, oc_methods_1, ml_methods)
   
  # update list of estimates for current simulation round
  ate_ens_1[j,] = do.call(cbind, dml_estimator$ate)                   # estimated effect theta in current simulation round
  ate_lasso_1[j] = dml_estimator$ate$lasso
  ate_xgb_1[j] = dml_estimator$ate$xgb
  ate_nn_1[j] = dml_estimator$ate$nn
  te_ens_1[,j] = dml_estimator$te$ens
  te_lasso_1[,j] = dml_estimator$te$lasso
  te_xgb_1[,j] = dml_estimator$te$xgb
  te_nn_1[,j] = dml_estimator$te$nn
  se_po_ens_1[j,] = dml_estimator$se_po[1,]
  se_po_lasso_1[j,] = dml_estimator$se_po[2,]
  se_po_xgb_1[j,] = dml_estimator$se_po[3,]
  se_po_nn_1[j,] = dml_estimator$se_po[4,]
  se_te_ens_1[j] = dml_estimator$se_te[1]
  se_te_lasso_1[j] = dml_estimator$se_te[2]
  se_te_xgb_1[j] = dml_estimator$se_te[3]
  se_te_nn_1[j] = dml_estimator$se_te[4]
  ps_ensemble_1[j,] = dml_estimator$w_ens_ps           # store weights of current simulation round
  oc_ensemble_1[j,] = dml_estimator$w_ens_oc          # store weights of current simulation round
  
  print(paste("Simulation 1, 50/50, round: ", j))
}

#### Results
# Average treatment effect 50/50
avg_effect_ens_1 = mean(ate_ens_1)                          
avg_effect_lasso_1 = mean(ate_lasso_1)
avg_effect_xgb_1 = mean(ate_xgb_1)
avg_effect_nn_1 = mean(ate_nn_1)

# Ensemble weights of E[Y|X]
oc_ensemble_weights_1 = as.data.frame(t(colMeans(oc_ensemble_1)))
for (i in 1:length(oc_methods_1)) {
  if (!is.null(oc_methods_1[[i]]$name)) colnames(oc_ensemble_weights_1)[i] = oc_methods_1[[i]]$name
  oc_ensemble_1 = as.data.frame(oc_ensemble_1)
  colnames(oc_ensemble_1)[i] = oc_methods_1[[i]]$name
}

# Ensemble weights of E[D|X]
ps_ensemble_weights_1 = as.data.frame(t(colMeans(ps_ensemble_1)))
for (i in 1:length(ps_methods_1)) {
  if (!is.null(ps_methods_1[[i]]$name)) colnames(ps_ensemble_weights_1)[i] = ps_methods_1[[i]]$name
  ps_ensemble_1 = as.data.frame(ps_ensemble_1)
  colnames(ps_ensemble_1)[i] = ps_methods_1[[i]]$name
}

# Print the results
paste("Average treatment effect:", round(avg_effect_ens_1, 3))
paste(sprintf("Ensemble weight E[Y|X] %s:",colnames(oc_ensemble_weights_1)), round(oc_ensemble_weights_1, 3))
paste(sprintf("Ensemble weight E[D|X] %s:",colnames(ps_ensemble_weights_1)), round(ps_ensemble_weights_1, 3))

#### Evaluation



# Simulation 1; 10/90 -----------------------------------------------

# Hyperparameter Tuning for DGP 1
hyperparams_small_1 = hyperparam_tuning_1(n_covariates, n_observations, mu1, tau1, pi1, sigma = 1, cv_folds, w = -3.2)
ps_methods_small_1 = hyperparams_small_1$ps_methods
oc_methods_small_1 = hyperparams_small_1$oc_methods

for (j in 1:n_simulations) {
  ###########
  ## 10/90 ##
  ###########
  # simulate data
  data_small = generalDGP(n_covariates, n_observations, mu1, tau1, pi1, sigma = 1, w = -3.2)
  Y_small = data[[1]]
  D_small = data[[2]]
  X_small = data[[3]]
  n_obs = seq(1, nrow(X), 1)
  
  # run the DML estimator, cross-fitting is done within the algorithm, ate and average weights of ensemble are extracted
  dml_estimator_small  = dml_ens_ml(Y, D, X, ps_methods_1, oc_methods_1, ml_methods)
  
  # update list of estimates for current simulation round
  ate_ens_small_1[j] = dml_estimator_small$ate_ens                      # estimated effect theta in current simulation round
  ate_lasso_small_1[j] = dml_estimator_small$ate_lasso
  ate_xgb_small_1[j] = dml_estimator_small$ate_xgb
  ate_nn_small_1[j] = dml_estimator_small$ate_nn
  te_ens_small_1[,j] = dml_estimator_small$te_ens
  te_lasso_small_1[,j] = dml_estimator_small$te_lasso
  te_xgb_small_1[,j] = dml_estimator_small$te_xgb
  te_nn_small_1[,j] = dml_estimator_small$te_nn
  se_po_ens_small_1[j,] = dml_estimator_small$se_po[1,]
  se_po_lasso_small_1[j,] = dml_estimator_small$se_po[2,]
  se_po_xgb_small_1[j,] = dml_estimator_small$se_po[3,]
  se_po_nn_small_1[j,] = dml_estimator_small$se_po[4,]
  se_te_ens_small_1[j] = dml_estimator_small$se_te[1]
  se_te_lasso_small_1[j] = dml_estimator_small$se_te[2]
  se_te_xgb_small_1[j] = dml_estimator_small$se_te[3]
  se_te_nn_small_1[j] = dml_estimator_small$se_te[4]
  ps_ensemble_small_1[j,] = dml_estimator_small$w_ens_ps           # store weights of current simulation round
  oc_ensemble_small_1[j,] = dml_estimator_small$w_ens_oc          # store weights of current simulation round
  
  print(paste("Simulation 1, 10/90, round: ", j))
}

# Average treatment effect 10/90
avg_effect_ens_small_1 = mean(ate_ens)                          
avg_effect_lasso_small_1 = mean(ate_lasso)
avg_effect_xgb_small_1 = mean(ate_xgb)
avg_effect_nn_small_1 = mean(ate_nn)

# Ensemble weights of E[Y|X]
oc_ensemble_weights_small_1 = as.data.frame(t(colMeans(oc_ensemble_small_1)))
for (i in 1:length(oc_methods_small_1)) {
  if (!is.null(oc_methods_small_1[[i]]$name)) colnames(oc_ensemble_weights_small_1)[i] = oc_methods_small_1[[i]]$name
  oc_ensemble_small_1 = as.data.frame(oc_ensemble_small_1)
  colnames(oc_ensemble_small_1)[i] = oc_methods_small_1[[i]]$name
}


# Ensemble weights of E[D|X]
ps_ensemble_weights_small_1 = as.data.frame(t(colMeans(ps_ensemble_small_1)))
for (i in 1:length(ps_methods_small_1)) {
  if (!is.null(ps_methods_small_1[[i]]$name)) colnames(ps_ensemble_weights_small_1)[i] = ps_methods_small_1[[i]]$name
  ps_ensemble_small_1 = as.data.frame(ps_ensemble_small_1)
  colnames(ps_ensemble_small_1)[i] = ps_methods_small_1[[i]]$name
}

# Print the results
paste("Average treatment effect:", round(avg_effect_ens_small_1, 3))
paste(sprintf("Ensemble weight E[Y|X] %s:",colnames(oc_ensemble_weights_small_1)), round(oc_ensemble_weights_small_1, 3))
paste(sprintf("Ensemble weight E[D|X] %s:",colnames(ps_ensemble_weights_small_1)), round(ps_ensemble_weights_small_1, 3))

