hyperparam_xgboost = function(y, x, cvfold) {
  # implementation of a random search algorithm to find the best possible combination of hyperparameters for xgboost
  
  parameters_list = list()
  
  for (i in 1:10){
    param <- list(booster = "gbtree",
                  objective = "reg:squarederror",
                  max_depth = sample(3:10, 1),
                  eta = runif(1, .01, .3),
                  subsample = runif(1, .7, 1),
                  colsample_bytree = runif(1, .6, 1),
                  min_child_weight = sample(0:10, 1),
                  lambda = sample(0:5, 1)
    )
    parameters <- as.data.frame(param)
    parameters_list[[i]] <- parameters
  }
  
  ### Create object that contains all randomly created hyperparameters
  params_df = do.call(rbind, parameters_list)
  
  dt_cv = xgb.DMatrix(data = x[as.logical(cvfold), ], label = y[as.logical(cvfold)])
  dval_cv = xgb.DMatrix(data = x[!cvfold, ], label = y[!cvfold])
  lowest_error_list = list()
  
  ### Use randomly created parameters to create 10,000 XGBoost-models
  for (row in 1:nrow(params_df)){
    model_cv <- xgb.train(data=dt_cv,
                        booster = "gbtree",
                        objective = "reg:squarederror",
                        max_depth = params_df$max_depth[row],
                        eta = params_df$eta[row],
                        subsample = params_df$subsample[row],
                        colsample_bytree = params_df$colsample_bytree[row],
                        min_child_weight = params_df$min_child_weight[row],
                        lambda = params_df$lambda,
                        nrounds= 300,
                        eval_metric = "rmse",
                        early_stopping_rounds= 30,
                        print_every_n = 100,
                        watchlist = list(train = dt_cv, val = dval_cv)
    )
    lowest_error_list[[row]] <- as.data.frame(1 - min(model_cv$evaluation_log$train_rmse))
  }
  
  ### Create object that contains all accuracy's
  lowest_error_df = do.call(rbind, lowest_error_list)
  
  ### Bind columns of accuracy values and random hyperparameter values
  randomsearch = cbind(lowest_error_df, params_df)
  
  ### Quickly display highest accuracy
  bestparams = randomsearch[which.max(randomsearch$`1 - min(model_cv$evaluation_log$train_rmse)`), ]
  
  finalparams = list(booster = bestparams$booster, 
                       objective = bestparams$objective,
                       max_depth = bestparams$max_depth,
                       eta = bestparams$eta,
                       subsample = bestparams$subsample,
                       colsample_bytree = bestparams$colsample_bytree,
                       min_child_weight = bestparams$min_child_weight,
                       lambda = bestparams$lambda,
                       nrounds= 300,
                       eval_metric = "rmse")
  
  return(finalparams)
}