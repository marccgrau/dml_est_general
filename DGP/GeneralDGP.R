#' Data Generating Process
#' 
#' All DGPs are following the paper of Powers et al. (2018)
#' DGP 1 as standard version of Powers et al.
#' 
#' parameters necessary for data generating processes
#' @param n_covariates number of covariates
#' @param n_observations number of rows in each simulation round
#' @param beta coefficients of controls
#' 
#' The data generating model is determined by the following variables
#' @param X_i feature vectors distributed according to distribution D, odd numbered features iid from a standard normal, even numbered features from a Bernoulli
#' @param D_i treatment vector D_i ~ind Bernoulli(pi(X_i))
#' @param Y_i outcome vector Y_i ~ind Normal(mu(X_i) + (D_i - 1/2)tau(X_i), sigma_Y^2 )

generalDGP = function(n_covariates, n_observations, mufunc, taufunc, psfunc, sigma = 1, w = 0) {
  
  # construct the correlation matrix
  init_matrix = qr.Q(qr(matrix(rnorm(n_covariates^2), n_covariates)))
  cov_matrix = crossprod(init_matrix, init_matrix*(n_covariates:1))
  corr_matrix = cov2cor(cov_matrix)
  
  copX = mvrnorm(n_observations, mu = rep(0, n_covariates), Sigma = corr_matrix)
  uniforms = pnorm(copX)
  

  X = matrix(NA, nrow = n_observations, ncol = n_covariates)
    for (i in 1:n_covariates){
      if(i %% 2 == 0) {
        X[,i] = qbinom(uniforms[,i], 1, 0.5)
      } else {
        X[,i] = qnorm(uniforms[,i], 0, 1)
      }
    }
    
    ps = propensity_func(X, psfunc, w)
    
    D = rbinom(n_observations, 1, ps)
    
    Y = y_func(n_observations, X, D, mufunc, taufunc, sigma)
    
    true_te = taufunc(X)
    
  return(list(Y, D, X, true_te))
}

