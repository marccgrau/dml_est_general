#' Data Generating Process I
#' Adapted from Powers (2018); Some methods for heterogeneous treatment effect estimation in high dimensions
#' 
#' 
#' @param n observations in the data set
#' @param mu conditional mean function (nuisance function)
#' @param pi propensity function (nuisance function)
#' @param tau treatment effect function
#' 
#' The data generating model is determined by the following variables
#' @param X_i feature vectors distributed according to distribution D
#' @param D_i treatment vector D_i ~ind Bernoulli(pi(X_i))
#' @param Y_i outcome vector Y_i ~ind Normal(mu(X_i) + (D_i - 1/2)tau(X_i), sigma_Y^2 )

# conditional mean function 
# mu_1 = mu(X) + tau(X)/2
# mu_0 = mu(X) - tau(X)/2
condmean_func = function(X, D, mufunc, taufunc){
  mu = D*(mufunc(X) + taufunc(X)/2) + (1 - D)*(mufunc(X) - taufunc(X)/2)
  return(mu)
}

# as DML is destined to be used for observational studies we simulate observational data
# persons with greater mean effect are more likely to receive the treatment
# this represents a usual bias in observational studies as mentioned by Powers (2018)
propensity_func = function(X, D, mufunc, taufunc, treatment_func){
  pi = (exp(condmean_func(X, D, mufunc, taufunc) - treatment_func(X)/2))/(1 + exp(condmean_func(X, D, mufunc, taufunc) - treatment_func(X)/2))
  return(pi)
}

Y_func = function(n_obs, X, D, mufunc, taufunc, sigma = 1) {
  y = rnorm(n_obs, mean = mufunc(X) + (D - (1/2))*taufunc(X), sd = sigma)
  return(y)
}
