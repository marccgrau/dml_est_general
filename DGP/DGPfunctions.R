#' Data Generating Process I
#' Adapted from Powers (2018); Some methods for heterogeneous treatment effect estimation in high dimensions
#' 
#' 
#' @param mu conditional mean function (nuisance function)
#' @param pi propensity function (nuisance function)
#' @param tau treatment effect function



# conditional mean function 
# mu_1 = mu(X) + tau(X)/2
# mu_0 = mu(X) - tau(X)/2
# for determining the true value of the ate 
#condmean_func = function(X, D, mufunc, taufunc){
#  mu = D*(mufunc(X) + taufunc(X)/2) + (1 - D)*(mufunc(X) - taufunc(X)/2)
#  return(mu)
#}

propensity_func = function(X, psfunc, w = 0){
  pi = 1/(1 + exp(-psfunc(X,w)))
  return(pi)
}

y_func = function(n_obs, X, D, mufunc, taufunc, sigma = 1) {
  y = rnorm(n_obs, mean = mufunc(X) + (D - (1/2))*taufunc(X), sd = sigma)
  return(y)
}




