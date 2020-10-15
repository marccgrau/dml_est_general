hyperparam_lasso_ps = function(y, x) {
  ## this function performs cross-validation for lasso regressions
  ## the output is a sequence of lambdas, which gives the approximate region of the optimal lambda for the respective DGP
  
  # define sequence of lambdas for cross-validation
  #seq_lambda_test = 10^seq(1, -2, by = -.01)
  
  # train cross-validated models
  lambda_cv = cv.glmnet(x, y, nfolds = 10, family = "binomial", alpha = 1)$lambda.min
  
  seq_lambda_final = if(lambda_cv - 0.05 < 0) {seq(0, lambda_cv + 0.05, 0.001)} else{seq(lambda_cv - 0.05, lambda_cv + 0.05, 0.001)}
  
  args = list("lambda" = seq_lambda_final, "family" = "binomial", "alpha" = 1)
  
  return(args)
}