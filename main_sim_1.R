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

toload <- c("grf", "tidyverse", "hdm", "glmnet", "nnls", "Matrix", 
            "matrixStats", "xgboost", "neuralnet", "MASS", "MLmetrics", 
            "keras", "tfdatasets", "data.table", "lessR", "ggthemes", "tictoc")
toinstall <- toload[which(toload %in% installed.packages()[,1] == F)]
lapply(toinstall, install.packages, character.only = TRUE)
lapply(toload, require, character.only = TRUE)

rm(list = ls())
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

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
source("doubleML/dml_ens_trim.R")

# Ensemble learner
source("ensemble_method/ensemble.R")
source("ensemble_method/ml_wrapper.R")
source("ensemble_method/utils_ensemble.R")


# Parameters --------------------------------------------------------------

### Define necessary parameters
## Monte Carlo Simulation
n_simulations = 20                  # Number of simulation rounds for Monte Carlo Study

## Data
n_covariates = 15                    # Number of confounders
n_observations = 2000               # Number of observations in simulated dataset

## DML estimator
cv_folds = 2                        # Number of folds for cross-validation of used ML methods in the ensemble method
trimmer = 0.05                      # quantile which gets deducted from each side to avoid extreme values in the score estimation

# declare which ml methods are going to be used by providing short versions of names
# will be further used to create necessary variables
ml_methods = c("ens", "lasso", "xgb", "nn")

# Simulation 1; 50/50 -----------------------------------------------

# Hyperparameter Tuning for DGP 1
hyperparams = hyperparam_tuning_1(n_covariates, n_observations, mu1, tau1, pi1, sigma = 1, cv_folds, w = 0)
ps_methods = hyperparams$ps_methods
oc_methods = hyperparams$oc_methods

## create empty matrices to fill throughout the simulation
# average treatment effect
for(i in 1:length(ml_methods)) { 
  nam <- paste("ate_", ml_methods[i], sep = "")
  assign(nam, rep(NA, n_simulations))
}
# treatment effects
for(i in 1:length(ml_methods)) { 
  nam <- paste("te_", ml_methods[i], sep = "")
  assign(nam, matrix(NA, (n_observations*(1-(2*trimmer))/2), n_simulations))
}
# standard error potential outcomes
for(i in 1:length(ml_methods)) { 
  nam <- paste("se_po_", ml_methods[i], sep = "")
  assign(nam, matrix(NA, n_simulations, 2))
}

for(i in 1:length(ml_methods)) { 
  nam <- paste("se_te_", ml_methods[i], sep = "")
  assign(nam, rep(NA, n_simulations))
}


oc_ensemble = matrix(NA, n_simulations, length(oc_methods))
ps_ensemble = matrix(NA, n_simulations, length(ps_methods)) 

# empty matrices to store effective values of treatment
te_t = matrix(NA, n_observations, n_simulations)
p_t = matrix(NA, n_observations, n_simulations)
ate_t = rep(NA, n_simulations)
se_te_t = rep(NA, n_simulations)

ate = matrix(NA, n_simulations,length(ml_methods))
colnames(ate) = ml_methods
tic()
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
  true_p = data[[5]]
  te_t[,j] = true_te
  p_t[,j] = true_p
  ate_t[j] = mean(true_te)
  se_te_t[j] = sd(true_te)
  
  n_obs = seq(1, nrow(X), 1)
  rm(data)
  
  # run the DML estimator, cross-fitting is done within the algorithm, ate and average weights of ensemble are extracted
  dml_estimator  = dml_ens_trim(Y, D, X, ps_methods, oc_methods, ml_methods, ens_folds = 2, trim = trimmer)
   
  # update list of estimates for current simulation round
  ate[j,] = do.call(cbind, dml_estimator$ate)                   # estimated effect theta in current simulation round
  ate_ens[j] = dml_estimator$ate$ens
  ate_lasso[j] = dml_estimator$ate$lasso
  ate_xgb[j] = dml_estimator$ate$xgb
  ate_nn[j] = dml_estimator$ate$nn
  te_ens[,j] = dml_estimator$te$ens
  te_lasso[,j] = dml_estimator$te$lasso
  te_xgb[,j] = dml_estimator$te$xgb
  te_nn[,j] = dml_estimator$te$nn
  se_po_ens[j,] = dml_estimator$se_po[1,]
  se_po_lasso[j,] = dml_estimator$se_po[2,]
  se_po_xgb[j,] = dml_estimator$se_po[3,]
  se_po_nn[j,] = dml_estimator$se_po[4,]
  se_te_ens[j] = dml_estimator$se_te[1]
  
  se_te_lasso[j] = dml_estimator$se_te[2]
  se_te_xgb[j] = dml_estimator$se_te[3]
  se_te_nn[j] = dml_estimator$se_te[4]
  ps_ensemble[j,] = dml_estimator$w_ens_ps           # store weights of current simulation round
  oc_ensemble[j,] = dml_estimator$w_ens_oc          # store weights of current simulation round
  
  print(paste("Simulation 1, 50/50, round: ", j))
  rm(dml_estimator)
  gc()
}
toc()

#### Results
# Average treatment effect 50/50
avg_effect_t = mean(ate_t)
avg_effect_ens = mean(ate_ens)                          
avg_effect_lasso = mean(ate_lasso)
avg_effect_xgb = mean(ate_xgb)
avg_effect_nn = mean(ate_nn)

# Ensemble weights of E[Y|X]
oc_ensemble_weights = as.data.frame(t(colMeans(oc_ensemble)))
colnames(oc_ensemble) = ml_methods[-1]

# Ensemble weights of E[D|X]
ps_ensemble_weights = as.data.frame(t(colMeans(ps_ensemble)))
colnames(ps_ensemble) = ml_methods[-1]

# Print the results
paste("Average treatment effect:", round(avg_effect_ens, 3))
paste(sprintf("Ensemble weight E[Y|X] %s:",colnames(oc_ensemble_weights)), round(oc_ensemble_weights, 3))
paste(sprintf("Ensemble weight E[D|X] %s:",colnames(ps_ensemble_weights)), round(ps_ensemble_weights, 3))

#### Store the results as csv-files
folder = "output/sim_50_1" # set folder to store values
# store ate values
output_ate = cbind(ate_t, ate)
fwrite(as.data.table(output_ate), file = file.path(directory_path, folder, "ate.csv"))

# store te values, btw assign each simulation the respective column name
colnames_sim = to("sim_", until = n_simulations, from = 1)
for (i in 1:length(ml_methods)){
  temp = as.data.table(get(paste0("te_", ml_methods[i])))
  colnames(temp) = colnames_sim
  fwrite(temp, file = file.path(directory_path, folder, paste0("te_", ml_methods[i], ".csv")))
}

# standard errors of score
output_se_te = cbind(se_te_t, se_te_ens, se_te_lasso, se_te_xgb, se_te_nn)
fwrite(as.data.table(output_se_te), file = file.path(directory_path, folder, "se_te.csv"))

# Each ensembles weights
fwrite(as.data.table(ps_ensemble), file = file.path(directory_path, folder, "ps_ensemble.csv"))
fwrite(as.data.table(oc_ensemble), file = file.path(directory_path, folder, "oc_ensemble.csv"))

# Simulation 1; 10/90 -----------------------------------------------

# Hyperparameter Tuning for DGP 1
hyperparams_small = hyperparam_tuning_1(n_covariates, n_observations, mu1, tau1, pi1, sigma = 1, cv_folds, w = -3.2)
ps_methods_small = hyperparams_small$ps_methods
oc_methods_small = hyperparams_small$oc_methods

## create empty matrices to fill throughout the simulation
# average treatment effect
for(i in 1:length(ml_methods)) { 
  nam <- paste("ate_", ml_methods[i], "_small", sep = "")
  assign(nam, rep(NA, n_simulations))
}
# treatment effects
for(i in 1:length(ml_methods)) { 
  nam <- paste("te_", ml_methods[i], "_small", sep = "")
  assign(nam, matrix(NA, (n_observations*(1-(2*trimmer))/2), n_simulations))
}
# standard error potential outcomes
for(i in 1:length(ml_methods)) { 
  nam <- paste("se_po_", ml_methods[i], "_small", sep = "")
  assign(nam, matrix(NA, n_simulations, 2))
}

for(i in 1:length(ml_methods)) { 
  nam <- paste("se_te_", ml_methods[i], "_small", sep = "")
  assign(nam, rep(NA, n_simulations))
}

oc_ensemble_small = matrix(NA, n_simulations, length(oc_methods_small))
ps_ensemble_small = matrix(NA, n_simulations, length(ps_methods_small)) 

# empty matrices to store effective values of treatment
te_t_small = matrix(NA, n_observations, n_simulations)
p_t_small = matrix(NA, n_observations, n_simulations)
ate_t_small = rep(NA, n_simulations)
se_te_t_small = rep(NA, n_simulations)

ate_small = matrix(NA, n_simulations,length(ml_methods))
colnames(ate_small) = ml_methods

tic()
for (j in 1:n_simulations) {
  ###########
  ## 50/50 ##
  ###########
  # simulate data
  data = generalDGP(n_covariates, n_observations, mu1, tau1, pi1, sigma = 1, w = -3.2)
  Y = data[[1]]
  D = data[[2]]
  X = data[[3]]
  true_te_small = data[[4]]
  true_p_small = data[[5]]
  te_t_small[,j] = true_te_small
  p_t_small[,j] = true_p_small
  ate_t_small[j] = mean(true_te_small)
  se_te_t_small[j] = sd(true_te_small)
  
  n_obs = seq(1, nrow(X), 1)
  rm(data)
  
  # run the DML estimator, cross-fitting is done within the algorithm, ate and average weights of ensemble are extracted
  dml_estimator_small  = dml_ens_trim(Y, D, X, ps_methods_small, oc_methods_small, ml_methods, ens_folds = 2, trim = trimmer)
  
  # update list of estimates for current simulation round
  ate_small[j,] = do.call(cbind, dml_estimator_small$ate)                   # estimated effect theta in current simulation round
  ate_ens_small[j] = dml_estimator_small$ate$ens
  ate_lasso_small[j] = dml_estimator_small$ate$lasso
  ate_xgb_small[j] = dml_estimator_small$ate$xgb
  ate_nn_small[j] = dml_estimator_small$ate$nn
  te_ens_small[,j] = dml_estimator_small$te$ens
  te_lasso_small[,j] = dml_estimator_small$te$lasso
  te_xgb_small[,j] = dml_estimator_small$te$xgb
  te_nn_small[,j] = dml_estimator_small$te$nn
  se_po_ens_small[j,] = dml_estimator_small$se_po[1,]
  se_po_lasso_small[j,] = dml_estimator_small$se_po[2,]
  se_po_xgb_small[j,] = dml_estimator_small$se_po[3,]
  se_po_nn_small[j,] = dml_estimator_small$se_po[4,]
  se_te_ens_small[j] = dml_estimator_small$se_te[1]
  
  se_te_lasso_small[j] = dml_estimator_small$se_te[2]
  se_te_xgb_small[j] = dml_estimator_small$se_te[3]
  se_te_nn_small[j] = dml_estimator_small$se_te[4]
  ps_ensemble_small[j,] = dml_estimator_small$w_ens_ps           # store weights of current simulation round
  oc_ensemble_small[j,] = dml_estimator_small$w_ens_oc          # store weights of current simulation round
  
  print(paste("Simulation 1, 50/50, round: ", j))
  rm(dml_estimator_small)
  gc()
}
toc()

#### Results
# Average treatment effect 50/50
avg_effect_t_small = mean(ate_t_small)
avg_effect_ens_small = mean(ate_ens_small)                          
avg_effect_lasso_small = mean(ate_lasso_small)
avg_effect_xgb_small = mean(ate_xgb_small)
avg_effect_nn_small = mean(ate_nn_small)

# Ensemble weights of E[Y|X]
oc_ensemble_weights_small = as.data.frame(t(colMeans(oc_ensemble_small)))
colnames(oc_ensemble_small) = ml_methods[-1]

# Ensemble weights of E[D|X]
ps_ensemble_weights_small = as.data.frame(t(colMeans(ps_ensemble_small)))
colnames(ps_ensemble_small) = ml_methods[-1]

# Print the results
paste("Average treatment effect:", round(avg_effect_ens_small, 3))
paste(sprintf("Ensemble weight E[Y|X] %s:",colnames(oc_ensemble_weights_small)), round(oc_ensemble_weights_small, 3))
paste(sprintf("Ensemble weight E[D|X] %s:",colnames(ps_ensemble_weights_small)), round(ps_ensemble_weights_small, 3))

#### Store the results as csv-files
folder_small = "output/sim_10_1" # set folder to store values
# store ate values
output_ate_small = cbind(ate_t_small, ate_small)
fwrite(as.data.table(output_ate_small), file = file.path(directory_path, folder_small, "ate.csv"))

# store te values, btw assign each simulation the respective column name
colnames_sim = to("sim_", until = n_simulations, from = 1)
for (i in 1:length(ml_methods)){
  temp = as.data.table(get(paste0("te_", ml_methods[i], "_small")))
  colnames(temp) = colnames_sim
  fwrite(temp, file = file.path(directory_path, folder_small, paste0("te_", ml_methods[i], ".csv")))
}

# standard errors of score
output_se_te_small = cbind(se_te_t_small, se_te_ens_small, se_te_lasso_small, se_te_xgb_small, se_te_nn_small)
fwrite(as.data.table(output_se_te_small), file = file.path(directory_path, folder_small, "se_te.csv"))

# Each ensembles weights
fwrite(as.data.table(ps_ensemble_small), file = file.path(directory_path, folder_small, "ps_ensemble.csv"))
fwrite(as.data.table(oc_ensemble_small), file = file.path(directory_path, folder_small, "oc_ensemble.csv"))

