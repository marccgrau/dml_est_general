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