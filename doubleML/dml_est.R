dml_est = function(y, d, x, ps_methods, oc_methods) {
  # Double Machine Learning estimator using the efficient score function as in Hahn 1998
  # Implementation inspired by Knaus DML for multiple treatments https://github.com/MCKnaus/dmlmt 
  # check if the given treatment vector is binary
  binary <- all(d %in% 0:1)
  if (binary){
    
    # get dummy for treated (d) and controls (1-d)
    treatment = cbind(d, 1-d)
    
    # estimate the propensity score p(x) and 1-p(x) for inverse probability weighting
    p_hat = matrix(NA, nrow(x), 2) # create empty matrix to store values
    
    p_fit = regression_forest(x, d) # fit model d ~ x
    p_hat[,1] = predict(p_fit, newdata = x)$predictions # predict the propensity score
    p_hat[,2] = 1 - p_hat[,1] # retrieve 1-p(x)
    
    # estimate the potential outcome 
    y_hat = matrix(NA, nrow(x), 2) # create empty matrix to store values
    
    for (i in 1:2){
      y_fit = regression_forest(x[treatment[,i] == 1, ], y[treatment[,i] == 1]) # fit potential outcome model for both treated and controls
      y_hat[,i] = predict(y_fit, newdata = x)$predictions # predict potential outcomes y_hat for treated and controls
    }
    
    # calculate efficient score E[(y - mu_hat)/p(x) + mu_hat]
    # this doubly robust estimator remains consistent with either p(x) or mu(x) misspecified
    mu_mat = matrix(NA, nrow(x), 2)
    ipw_mat = matrix(NA, nrow(x), 2)
    
    # calculate the inverse probability weights d/p(x) and (1-d)/(1-p(x))
    for (i in 1:2){
      ipw_mat[,i] = treatment[,i] / p_hat[,i]
    }
    
    for (i in 1:2){
      mu_mat[,i] = (y - y_hat[,i])*ipw_mat[,i] + y_hat[,i]
    }
    
    # calculate psi(D, theta, eta) = (g(1,X) - g(0,X)) + D(Y - g(1,X))/p(x) + (1-D)(Y - g(0,X))/(1-p(x)) - theta
    # which with moment and orthogonality condition fulfilled the difference between the two scores equals the treatment effect theta
    psi = mu_mat[,1] - mu_mat[,2]
    
    # average all the treatment effects to get the ATE
    ate = mean(psi)
    
    return("ATE" = ate)
    
  } else {
    print("Non-binary treatments, please check your data")
  }
}