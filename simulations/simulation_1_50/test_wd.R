## Default repo
local({r <- getOption("repos")
r["CRAN"] <- "http://cloud.r-project.org" 
options(repos=r)
})

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
print(directory_path)

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


### Define necessary parameters
folder = "output/sim_1_50" # set folder to store values
## Monte Carlo Simulation
n_simulations = 2                  # Number of simulation rounds for Monte Carlo Study

## Data
n_covariates = 15                    # Number of confounders
n_observations = 2000               # Number of observations in simulated dataset

## DML estimator
cv_folds = 2                        # Number of folds for cross-validation of used ML methods in the ensemble method

# declare which ml methods are going to be used by providing short versions of names
# will be further used to create necessary variables
ml_methods = c("ens", "lasso", "xgb", "nn")

# Simulation 1; 50/50 -----------------------------------------------

# create first simulation entry to later append all further simulations
# necessary due to computational restrictions
# store te values, btw assign each simulation the respective column name
colnames_sim = to("sim_", until = n_simulations, from = 1)

# start db for 50/50 Simulation
db1_50 = dbConnect(SQLite(), dbname = file.path(folder, "db1_50"))