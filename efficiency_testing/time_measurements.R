## Intro ----
local({r <- getOption("repos")
r["CRAN"] <- "http://cloud.r-project.org" 
options(repos=r)
})


toload <- c("grf", "tidyverse", "hdm", "glmnet", "nnls", "Matrix", 
            "matrixStats", "xgboost", "neuralnet", "MASS", "MLmetrics", 
            "keras", "tfdatasets", "data.table", "lessR", "ggthemes", "tictoc",
            "RSQLite", "RColorBrewer")
toinstall <- toload[which(toload %in% installed.packages()[,1] == F)]
lapply(toinstall, install.packages, character.only = TRUE)
lapply(toload, require, character.only = TRUE)
rm(list = ls())

thisPath <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  if (length(grep("^-f$", cmdArgs)) > 0) {
    # R console option
    normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
  } else if (length(grep("^--file=", cmdArgs)) > 0) {
    # Rscript/R console option
    scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
  } else if (Sys.getenv("RSTUDIO") == "1") {
    # RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else if (is.null(attr(stub, "srcref")) == FALSE) {
    # 'source'd via R console
    dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
  } else {
    stop("Cannot find file path")
  }
}

directory_path <- dirname(thisPath())
setwd(directory_path)

set.seed(12345)

## Load source functions
# General Functions
source(file.path(directory_path, "general_functions/general_utils.R"))

# DGPs
source(file.path(directory_path, "DGP/GeneralDGP.R"))
source(file.path(directory_path, "DGP/nuisance_functions.R"))
source(file.path(directory_path, "DGP/DGPfunctions.R"))

# Hyperparameter Tuning
source(file.path(directory_path, "hyperparam_tuning/hyperparam_tuning.R"))


# DML estimator
source(file.path(directory_path, "doubleML/dml_ens.R"))

# Ensemble learner
source(file.path(directory_path, "ensemble_method/ensemble.R"))
source(file.path(directory_path, "ensemble_method/ml_wrapper.R"))
source(file.path(directory_path, "ensemble_method/utils_ensemble.R"))

# Parameters --------------------------------------------------------------

### Define necessary parameters
folders = list.dirs("output", recursive = F)
folders = folders[-7]

## Monte Carlo Simulation
n_simulations = 10                  # Number of simulation rounds for Monte Carlo Study

## Data
n_covariates = 15                    # Number of confounders
n_observations = 2000               # Number of observations in simulated dataset

## DML estimator
cv_folds = 2                        # Number of folds for cross-validation of used ML methods in the ensemble method

ml_methods = c("ens", "lasso", "xgb", "nn")

# DGP 1, 10/90 ----
data = generalDGP(n_covariates, n_observations, mu1, tau1, pi1, sigma = 1, smalltreat = TRUE)
Y = data[[1]]
D = data[[2]]
X = data[[3]]
true_te = data[[4]]
true_p = data[[5]]
ate_t = mean(true_te)
se_te_t = sd(true_te)

n_obs = seq(1, nrow(X), 1)
rm(data)

time_vec_ens = rep(NA, n_simulations)
time_vec_lasso = rep(NA, n_simulations)
time_vec_xgb = rep(NA, n_simulations)
time_vec_nn = rep(NA, n_simulations)

ps_methods = readRDS(file = file.path(folders[1], "hyperparams", "ps_methods.rds"))
oc_methods = readRDS(file = file.path(folders[1], "hyperparams", "oc_methods.rds"))

## Ensemble
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods, oc_methods, ml_methods, ens_folds = 2)
  time = toc()
  time_vec_ens[i] = round(time$toc - time$tic, 3)
}
## Lasso
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[1], oc_methods[1], ml_methods[c(1,2)], ens_folds = 2)
  time = toc()
  time_vec_lasso[i] = round(time$toc - time$tic, 3)
}

## XGB
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[2], oc_methods[2], ml_methods[c(1,3)], ens_folds = 2)
  time = toc()
  time_vec_xgb[i] = round(time$toc - time$tic, 3)
}
## Neural Network
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[3], oc_methods[3], ml_methods[c(1,4)], ens_folds = 2)
  time = toc()
  time_vec_nn[i] = round(time$toc - time$tic, 3)
}

times_1_10 = c(mean(time_vec_ens), mean(time_vec_lasso), mean(time_vec_xgb), mean(time_vec_nn))

# DGP 1, 50/50 ----
data = generalDGP(n_covariates, n_observations, mu1, tau1, pi1, sigma = 1, smalltreat = FALSE)
Y = data[[1]]
D = data[[2]]
X = data[[3]]
true_te = data[[4]]
true_p = data[[5]]
ate_t = mean(true_te)
se_te_t = sd(true_te)

n_obs = seq(1, nrow(X), 1)
rm(data)

time_vec_ens = rep(NA, n_simulations)
time_vec_lasso = rep(NA, n_simulations)
time_vec_xgb = rep(NA, n_simulations)
time_vec_nn = rep(NA, n_simulations)

ps_methods = readRDS(file = file.path(folders[2], "hyperparams", "ps_methods.rds"))
oc_methods = readRDS(file = file.path(folder[2], "hyperparams", "oc_methods.rds"))

## Ensemble
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods, oc_methods, ml_methods, ens_folds = 2)
  time = toc()
  time_vec_ens[i] = round(time$toc - time$tic, 3)
}
## Lasso
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[1], oc_methods[1], ml_methods[c(1,2)], ens_folds = 2)
  time = toc()
  time_vec_lasso[i] = round(time$toc - time$tic, 3)
}

## XGB
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[2], oc_methods[2], ml_methods[c(1,3)], ens_folds = 2)
  time = toc()
  time_vec_xgb[i] = round(time$toc - time$tic, 3)
}
## Neural Network
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[3], oc_methods[3], ml_methods[c(1,4)], ens_folds = 2)
  time = toc()
  time_vec_nn[i] = round(time$toc - time$tic, 3)
}

times_1_50 = c(mean(time_vec_ens), mean(time_vec_lasso), mean(time_vec_xgb), mean(time_vec_nn))

# DGP 2, 10/90 ----
data = generalDGP(n_covariates, n_observations, mu2, tau2, pi2, sigma = 1, smalltreat = TRUE)
Y = data[[1]]
D = data[[2]]
X = data[[3]]
true_te = data[[4]]
true_p = data[[5]]
ate_t = mean(true_te)
se_te_t = sd(true_te)

n_obs = seq(1, nrow(X), 1)
rm(data)

time_vec_ens = rep(NA, n_simulations)
time_vec_lasso = rep(NA, n_simulations)
time_vec_xgb = rep(NA, n_simulations)
time_vec_nn = rep(NA, n_simulations)

ps_methods = readRDS(file = file.path(folders[3], "hyperparams", "ps_methods.rds"))
oc_methods = readRDS(file = file.path(folders[3], "hyperparams", "oc_methods.rds"))

## Ensemble
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods, oc_methods, ml_methods, ens_folds = 2)
  time = toc()
  time_vec_ens[i] = round(time$toc - time$tic, 3)
}
## Lasso
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[1], oc_methods[1], ml_methods[c(1,2)], ens_folds = 2)
  time = toc()
  time_vec_lasso[i] = round(time$toc - time$tic, 3)
}

## XGB
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[2], oc_methods[2], ml_methods[c(1,3)], ens_folds = 2)
  time = toc()
  time_vec_xgb[i] = round(time$toc - time$tic, 3)
}
## Neural Network
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[3], oc_methods[3], ml_methods[c(1,4)], ens_folds = 2)
  time = toc()
  time_vec_nn[i] = round(time$toc - time$tic, 3)
}

times_2_10 = c(mean(time_vec_ens), mean(time_vec_lasso), mean(time_vec_xgb), mean(time_vec_nn))

# DGP 2, 50/50 ----
data = generalDGP(n_covariates, n_observations, mu2, tau2, pi2, sigma = 1, smalltreat = FALSE)
Y = data[[1]]
D = data[[2]]
X = data[[3]]
true_te = data[[4]]
true_p = data[[5]]
ate_t = mean(true_te)
se_te_t = sd(true_te)

n_obs = seq(1, nrow(X), 1)
rm(data)

time_vec_ens = rep(NA, n_simulations)
time_vec_lasso = rep(NA, n_simulations)
time_vec_xgb = rep(NA, n_simulations)
time_vec_nn = rep(NA, n_simulations)

ps_methods = readRDS(file = file.path(folders[4], "hyperparams", "ps_methods.rds"))
oc_methods = readRDS(file = file.path(folder[4], "hyperparams", "oc_methods.rds"))

## Ensemble
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods, oc_methods, ml_methods, ens_folds = 2)
  time = toc()
  time_vec_ens[i] = round(time$toc - time$tic, 3)
}
## Lasso
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[1], oc_methods[1], ml_methods[c(1,2)], ens_folds = 2)
  time = toc()
  time_vec_lasso[i] = round(time$toc - time$tic, 3)
}

## XGB
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[2], oc_methods[2], ml_methods[c(1,3)], ens_folds = 2)
  time = toc()
  time_vec_xgb[i] = round(time$toc - time$tic, 3)
}
## Neural Network
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[3], oc_methods[3], ml_methods[c(1,4)], ens_folds = 2)
  time = toc()
  time_vec_nn[i] = round(time$toc - time$tic, 3)
}

times_2_50 = c(mean(time_vec_ens), mean(time_vec_lasso), mean(time_vec_xgb), mean(time_vec_nn))

# DGP 3, 10/90 ----
data = generalDGP(n_covariates, n_observations, mu3, tau3, pi3, sigma = 1, smalltreat = TRUE)
Y = data[[1]]
D = data[[2]]
X = data[[3]]
true_te = data[[4]]
true_p = data[[5]]
ate_t = mean(true_te)
se_te_t = sd(true_te)

n_obs = seq(1, nrow(X), 1)
rm(data)

time_vec_ens = rep(NA, n_simulations)
time_vec_lasso = rep(NA, n_simulations)
time_vec_xgb = rep(NA, n_simulations)
time_vec_nn = rep(NA, n_simulations)

ps_methods = readRDS(file = file.path(folders[5], "hyperparams", "ps_methods.rds"))
oc_methods = readRDS(file = file.path(folders[5], "hyperparams", "oc_methods.rds"))

## Ensemble
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods, oc_methods, ml_methods, ens_folds = 2)
  time = toc()
  time_vec_ens[i] = round(time$toc - time$tic, 3)
}
## Lasso
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[1], oc_methods[1], ml_methods[c(1,2)], ens_folds = 2)
  time = toc()
  time_vec_lasso[i] = round(time$toc - time$tic, 3)
}

## XGB
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[2], oc_methods[2], ml_methods[c(1,3)], ens_folds = 2)
  time = toc()
  time_vec_xgb[i] = round(time$toc - time$tic, 3)
}
## Neural Network
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[3], oc_methods[3], ml_methods[c(1,4)], ens_folds = 2)
  time = toc()
  time_vec_nn[i] = round(time$toc - time$tic, 3)
}

times_3_10 = c(mean(time_vec_ens), mean(time_vec_lasso), mean(time_vec_xgb), mean(time_vec_nn))


# DGP 3, 50/50 ----
data = generalDGP(n_covariates, n_observations, mu3, tau3, pi3, sigma = 1, smalltreat = FALSE)
Y = data[[1]]
D = data[[2]]
X = data[[3]]
true_te = data[[4]]
true_p = data[[5]]
ate_t = mean(true_te)
se_te_t = sd(true_te)

n_obs = seq(1, nrow(X), 1)
rm(data)

time_vec_ens = rep(NA, n_simulations)
time_vec_lasso = rep(NA, n_simulations)
time_vec_xgb = rep(NA, n_simulations)
time_vec_nn = rep(NA, n_simulations)

ps_methods = readRDS(file = file.path(folders[6], "hyperparams", "ps_methods.rds"))
oc_methods = readRDS(file = file.path(folder[6], "hyperparams", "oc_methods.rds"))

## Ensemble
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods, oc_methods, ml_methods, ens_folds = 2)
  time = toc()
  time_vec_ens[i] = round(time$toc - time$tic, 3)
}
## Lasso
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[1], oc_methods[1], ml_methods[c(1,2)], ens_folds = 2)
  time = toc()
  time_vec_lasso[i] = round(time$toc - time$tic, 3)
}

## XGB
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[2], oc_methods[2], ml_methods[c(1,3)], ens_folds = 2)
  time = toc()
  time_vec_xgb[i] = round(time$toc - time$tic, 3)
}
## Neural Network
for (i in 1:n_simulations){
  tic()
  dml_estimator  = dml_ens(Y, D, X, true_p, true_te, ps_methods[3], oc_methods[3], ml_methods[c(1,4)], ens_folds = 2)
  time = toc()
  time_vec_nn[i] = round(time$toc - time$tic, 3)
}

times_3_50 = c(mean(time_vec_ens), mean(time_vec_lasso), mean(time_vec_xgb), mean(time_vec_nn))



# Output ----
all_times = rbind.data.frame(times_1_10, times_1_50, times_2_10, times_2_50, times_3_10, times_3_50)
output = cbind.data.frame(folders, all_times)
colnames(output) = c("simulation", ml_methods)

