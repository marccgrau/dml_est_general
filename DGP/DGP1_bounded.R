#' Data Generating Process I
#' Linear transformations M(.) and G(.)
#' Simplest case
#' 
#' Use the generalised logistic function in order to get more flexible upper and lower bounds of the distribution
#' 
#' parameters necessary for data generating processes
#' @param n_covariates number of covariates
#' @param n_observations number of rows in each simulation round
#' @param beta coefficients of controls
#' 
#' Output
#' @param Y outcome variable
#' @param D treatment variable
#' @param X confounders

DGP1 = function(n_covariates, n_observations, beta, effect) {
  
  # construct the correlation matrix
  init_matrix = matrix(rnorm(n_covariates*n_covariates,mean=0,sd=1), n_covariates, n_covariates)
  cov_matrix = init_matrix %*% t(init_matrix)
  corr_matrix = t(chol(cov_matrix))
  
  # simulate data
  random.normal = matrix(rnorm(n_covariates*n_observations,0,1), nrow=n_covariates, ncol=n_observations)
  init_X = corr_matrix %*% random.normal
  X = t(init_X)
  
  # noise
  epsilon = rnorm(n_observations, mean = 0, sd = 1)
  eta = rnorm(n_observations, mean = 0, sd = 1)
  
  # construct nuisance functions
  M = X %*% beta
  G = X %*% beta^(-1)
  
  # generalised logistic function
  p <- 0.2/(1 + exp(-(M + epsilon))) # K = 0.2 as upper asymptote; A = 0 as lower asymptote
  
  # construct treatment vector with binary treatments
  D <- rbinom(n_observations, 1, p)
  
  # construct outcomes 
  Y = effect * D + G + eta
  
  return(list(Y, D, X))
}

