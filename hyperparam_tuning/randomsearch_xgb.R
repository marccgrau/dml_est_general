randomsearch_xgb_1 = function(){
  # randomly pick parameter values within defined boundaries
  parameters_list = list()
  for (i in 1:100){
    param <- list(booster = "gbtree", # to account for non-linearities a tree based estimator is used
                  objective = "reg:squarederror",
                  max_depth = sample(3:10, 1),
                  eta = runif(1, .01, .3),
                  subsample = runif(1, .7, 1),
                  colsample_bytree = runif(1, .6, 1),
                  min_child_weight = sample(0:10, 1),
                  lambda = sample(0:5, 1),
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