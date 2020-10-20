hyperparam_lasso_ps = function(y, x) {
  # propensity score nuisance parameter
  # this function performs cross-validation for lasso regression to get valid values for the hyperparameters
  # the result is a list of arguments which will further be used for the Monte Carlo simulation
  # to not only pass a single best value of lambda, but allow for some variation, a sequence will be constructed

  # train cross-validated models
  # the algorithm automatically selects an optimal region of lambdas to start the evaluation with
  # as the propensity score is a probability the results have to be within [0,1], therefore we use a binomial version
  # set alpha = 1 to only implement lasso regression
  lambda_cv = cv.glmnet(x, y, nfolds = 10, family = "binomial", alpha = 1)$lambda.min
  
  # construct sequence around optimal lambda
  seq_lambda_final = if(lambda_cv - 0.05 < 0) {seq(0, lambda_cv + 0.05, 0.001)} else{seq(lambda_cv - 0.05, lambda_cv + 0.05, 0.001)}
  
  args = list("lambda" = seq_lambda_final, "family" = "binomial", "alpha" = 1, "type.measure" = "mse")
  
  return(args)
}