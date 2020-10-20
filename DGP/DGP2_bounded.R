#' Data Generating Process II
#' Non-linear form with polynomials
#' medium complex case
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

g_func = function(X){
  for (i in 1:ncol(X)){
    X[,i] = X[,i]^(i)
  }
}

m_func = function(X, nu = 0, gamma = 1) {
  for (i in 1:ncol(X)){
    X[,i] = X[,i]^(sample(1:5,1,TRUE))
  }
}

DGP2 = function(n_covariates, n_observations, beta, effect) {
  
  # construct the correlation matrix
  init_matrix = matrix( rnorm(n_covariates*n_covariates,mean=0,sd=1), n_covariates, n_covariates)
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
  M = m_func(X %*% beta)
  G = g_func(X %*% beta)
  
  # logistic function to retrieve value in probability bounds [0,1]
  p <- 0.2/(1 + exp(-(M + epsilon)))
  
  # construct treatment vector with binary treatments
  D <- rbinom(n_observations, 1, p)
  
  # construct outcomes 
  Y = effect * D + G + eta
  
  return(list(Y, D, X))
}

