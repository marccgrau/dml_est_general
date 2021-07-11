randomsearch_xgb_oc = function(){
  # randomly pick parameter values within defined boundaries
  parameters_list = list()
  for (i in 1:500){
    param <- list(booster = "gbtree", # to account for non-linearities a tree based estimator is used
                  objective = "reg:squarederror",
                  max_depth = sample(3:10, 1),
                  eta = runif(1, 0.01, 0.4),
                  subsample = runif(1, 0.7, 1),
                  colsample_bytree = runif(1, 0.6, 1),
                  min_child_weight = sample(0:10, 1),
                  alpha = sample(0:6, 1),
                  nrounds= 300,
                  eval_metric = "rmse",
                  early_stopping_rounds= 30
    )
    parameters <- as.data.frame(param)
    parameters_list[[i]] <- parameters
  }
  
  # create dataframe of all randomly simulated parameter sets
  params_df = do.call(rbind, parameters_list)
  return(params_df)
}

randomsearch_xgb_ps = function(){
  # randomly pick parameter values within defined boundaries
  parameters_list = list()
  for (i in 1:500){
    param <- list(booster = "gbtree", # to account for non-linearities a tree based estimator is used
                  objective = "binary:logistic",
                  max_depth = sample(3:10, 1),
                  eta = runif(1, 0.01, 0.4),
                  subsample = runif(1, 0.7, 1),
                  colsample_bytree = runif(1, 0.6, 1),
                  min_child_weight = sample(0:10, 1),
                  alpha = sample(0:6, 1),
                  nrounds= 300,
                  eval_metric = "error",
                  early_stopping_rounds= 30
    )
    parameters <- as.data.frame(param)
    parameters_list[[i]] <- parameters
  }
  
  # create dataframe of all randomly simulated parameter sets
  params_df = do.call(rbind, parameters_list)
  return(params_df)
}
