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

directory_path <- dirname(dirname(thisPath()))
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
source(file.path(directory_path, "doubleML/dml_ens_ml.R"))
source(file.path(directory_path, "doubleML/dml_ens_trim.R"))

# Ensemble learner
source(file.path(directory_path, "ensemble_method/ensemble.R"))
source(file.path(directory_path, "ensemble_method/ml_wrapper.R"))
source(file.path(directory_path, "ensemble_method/utils_ensemble.R"))

# Parameters --------------------------------------------------------------

### Define necessary parameters
# set folder again
folder = "output/sim_50_1"
## Monte Carlo Simulation
n_simulations = 50                  # Number of simulation rounds for Monte Carlo Study

## Data
n_covariates = 15                    # Number of confounders
n_observations = 2000               # Number of observations in simulated dataset

## DML estimator
cv_folds = 2                        # Number of folds for cross-validation of used ML methods in the ensemble method
trimmer = 0.05                      # quantile which gets deducted from each side to avoid extreme values in the score estimation

# declare which ml methods are going to be used by providing short versions of names
# will be further used to create necessary variables
ml_methods = c("ens", "lasso", "xgb", "nn")

ps_methods = readRDS(file = file.path(folder, "hyperparams", "ps_methods.rds"))
oc_methods = readRDS(file = file.path(folder, "hyperparams", "oc_methods.rds"))

# Simulation --------------------------------------------------------------

for (j in 1:n_simulations) {
  
  db50_1 = dbConnect(SQLite(), dbname = file.path(folder, "db50_1"))
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
  ate_t = mean(true_te)
  se_te_t = sd(true_te)
  
  n_obs = seq(1, nrow(X), 1)
  rm(data)
  
  # run the DML estimator, cross-fitting is done within the algorithm, ate and average weights of ensemble are extracted
  dml_estimator  = dml_ens_trim(Y, D, X, true_p, true_te, ps_methods, oc_methods, ml_methods, ens_folds = 2, trim = trimmer)
  
  # update list of estimates for current simulation round
  ate = do.call(cbind, dml_estimator$ate)                   # estimated effect theta in current simulation round
  ate_ens = dml_estimator$ate$ens
  ate_lasso = dml_estimator$ate$lasso
  ate_xgb = dml_estimator$ate$xgb
  ate_nn = dml_estimator$ate$nn
  te_ens = dml_estimator$te$ens
  te_lasso = dml_estimator$te$lasso
  te_xgb = dml_estimator$te$xgb
  te_nn = dml_estimator$te$nn
  se_po_ens = dml_estimator$se_po[1,]
  se_po_lasso = dml_estimator$se_po[2,]
  se_po_xgb = dml_estimator$se_po[3,]
  se_po_nn = dml_estimator$se_po[4,]
  se_te_ens = dml_estimator$se_te[1]
  
  se_te_lasso = dml_estimator$se_te[2]
  se_te_xgb = dml_estimator$se_te[3]
  se_te_nn = dml_estimator$se_te[4]
  ps_ensemble = dml_estimator$w_ens_ps           # store weights of current simulation round
  oc_ensemble = dml_estimator$w_ens_oc          # store weights of current simulation round
  
  y_ens = dml_estimator$y_ens
  p_ens = dml_estimator$p_ens
  
  # get trimmed true values for comparison
  p_t_trim = dml_estimator$p_t_trim
  y_t_trim = dml_estimator$y_t_trim
  te_t_trim = dml_estimator$te_t_trim
  
  ## save results in csv-files
  # prepare consolidate output
  output_ate = cbind(ate_t, ate)
  output_se_te = cbind(se_te_t, se_te_ens, se_te_lasso, se_te_xgb, se_te_nn)
  output_se_po = cbind(se_po_ens, se_po_lasso, se_po_xgb, se_po_nn)
  
  ## write to database
  # average treatment effects
  dbWriteTable(conn = db50_1, name = "average_te", value = as.data.table(output_ate), row.names = FALSE, header = TRUE,
               overwrite = FALSE, append = TRUE)
  
  # treatment effects
  for (i in 1:length(ml_methods)){
    temp = as.data.table(t(c(get(paste0("te_", ml_methods[i])))))
    dbWriteTable(conn = db50_1, name = paste0("treatmenteffect_", ml_methods[i]), value = temp, row.names = FALSE, header = FALSE,
                 overwrite = FALSE, append = TRUE)
    rm(temp)
  }
  
  # standard errors of score
  dbWriteTable(conn = db50_1, name = "standerror_te", value = as.data.table(output_se_te), row.names = FALSE, header = TRUE,
               overwrite = FALSE, append = TRUE)
  dbWriteTable(conn = db50_1, name = "standerror_po", value = as.data.table(output_se_po), row.names = FALSE, header = TRUE,
               overwrite = FALSE, append = TRUE)
  
  # Each ensembles weights
  dbWriteTable(conn = db50_1, name = "propensity_ensemble", value = as.data.table(t(ps_ensemble)), row.names = FALSE, header = TRUE,
               overwrite = FALSE, append = TRUE)
  dbWriteTable(conn = db50_1, name = "outcome_ensemble", value = as.data.table(t(oc_ensemble)), row.names = FALSE, header = TRUE,
               overwrite = FALSE, append = TRUE)
  
  # estimated values and true values for treatment effects, propensity score and outcome
  dbWriteTable(conn = db50_1, name = "y_ensemble", value = as.data.table(t(c(y_ens))), row.names = FALSE, header = FALSE,
               overwrite = FALSE, append = TRUE)
  dbWriteTable(conn = db50_1, name = "p_ensemble", value = as.data.table(t(c(p_ens))), row.names = FALSE, header = FALSE,
               overwrite = FALSE, append = TRUE)
  dbWriteTable(conn = db50_1, name = "p_true_trim", value = as.data.table(t(c(p_t_trim))), row.names = FALSE, header = FALSE,
               overwrite = FALSE, append = TRUE)
  dbWriteTable(conn = db50_1, name = "te_true_trim", value = as.data.table(t(c(te_t_trim))), row.names = FALSE, header = FALSE,
               overwrite = FALSE, append = TRUE)
  
  print(paste("Simulation 1, 50/50, round: ", j))
  rm(dml_estimator)
  dbDisconnect(db50_1)
  rm(db50_1)
  gc()
}
