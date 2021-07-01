# Preliminaries -----------------------------------------------------------

## Load necessary packages, set working directory and seed, remove previously stored variables

toload <- c("grf", "tidyverse", "hdm", "glmnet", "nnls", "Matrix", "matrixStats", "xgboost", "neuralnet", "profvis", "MASS")
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
source("hyperparam_tuning/hyperparam_lasso.R")
source("hyperparam_tuning/hyperparam_xgboost.R")
source("hyperparam_tuning/randomsearch_xgb.R")
source("hyperparam_tuning/hyperparam_nnet.R")
source("hyperparam_tuning/gridframes_nn.R")

# DML estimator
source("doubleML/dml_est_cf_ensemble.R")

# Ensemble learner
source("ensemble_method/ensemble.R")
source("ensemble_method/ml_wrapper.R")
source("ensemble_method/utils_ensemble.R")


### Data Simulation
# parameters
n_covariates = 15                    # Number of confounders
n_observations = 1000               # Number of observations in simulated dataset

# simulation
test = profvis(
{data = generalDGP(n_covariates, n_observations, mu1, tau1, pi1, sigma = 1, w = 0)}
)

Y = data[[1]]
D = data[[2]]
X = data[[3]]
n_obs = seq(1, nrow(X), 1)

profvis()
