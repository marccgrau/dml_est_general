hyperparam_xgboost = function(y, x, cvfold, params_df) {
  # potential outcome
  # implementation of a random search algorithm to find the best possible combination of hyperparameters for xgboost
  # to save computational resources the random search is preferable over the grid search algorithm
  # especially in the case of many potential hyperparameters
  
  # split into training and testing set for out of sample evaluation
  n = nrow(x)
  fold = sample(seq_len(nrow(x)), size = n*0.75)
  
  x_tr = x[fold, ]
  x_te = x[-fold, ]
  
  y_tr = y[fold ]
  y_te = y[-fold]
  
  # bring to xgb-format
  dt_cv = xgb.DMatrix(data = x_tr, label = y_tr)
  dval_cv = xgb.DMatrix(data = x_te, label = y_te)
  lowest_errors = list()
  
  # simulate for each parameter set an xgboost model
  for (row in 1:nrow(params_df)){
    model_cv <- xgb.train(data=dt_cv,
                        booster = params_df$booster[row],
                        objective = params_df$objective[row],
                        max_depth = params_df$max_depth[row],
                        eta = params_df$eta[row],
                        subsample = params_df$subsample[row],
                        colsample_bytree = params_df$colsample_bytree[row],
                        min_child_weight = params_df$min_child_weight[row],
                        lambda = params_df$lambda,
                        nrounds= params_df$nrounds,
                        eval_metric = params_df$eval_metric,
                        early_stopping_rounds= params_df$early_stopping_rounds,
                        print_every_n = 100,
                        watchlist = list(train = dt_cv, val = dval_cv)
    )
    lowest_errors[[row]] <- as.data.frame(min(model_cv$evaluation_log$val_rmse))
  }
  
  # extract the rmse of each parameter set
  lowest_error_df = do.call(rbind, lowest_errors)
  
  # associate each rmse with its parameter set
  randomsearch = cbind(lowest_error_df, params_df)
  
  # extract lowest rmse parameter set
  bestparams = randomsearch[which.min(randomsearch$`min(model_cv$evaluation_log$val_rmse)`), ]
  
  # output the best performing hyperparameters
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