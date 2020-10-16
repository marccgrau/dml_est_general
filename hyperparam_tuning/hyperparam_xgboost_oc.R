hyperparam_xgboost_oc = function(y, x, cvfold) {
  # potential outcome
  # implementation of a random search algorithm to find the best possible combination of hyperparameters for xgboost
  # to save computational resources the random search is preferable over the grid search algorithm
  # especially in the case of many potential hyperparameters
  
  # split into training and testing set for out of sample evaluation
  parameters_list = list()
  n = nrow(x)
  fold = sample(seq_len(nrow(x)), size = n*0.75)
  
  x_tr = x[fold, ]
  x_te = x[-fold, ]
  
  y_tr = y[fold ]
  y_te = y[-fold]
  
  # randomly pick parameter values within defined boundaries
  for (i in 1:100){
    param <- list(booster = "gbtree", # to account for non-linearities a tree based estimator is used
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
  
  # create dataframe of all randomly simulated parameter sets
  params_df = do.call(rbind, parameters_list)
  
  # bring to xgb-format
  dt_cv = xgb.DMatrix(data = x_tr, label = y_tr)
  dval_cv = xgb.DMatrix(data = x_te, label = y_te)
  lowest_error_list = list()
  
  # simulate for each parameter set an xgboost model
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
    lowest_error_list[[row]] <- as.data.frame(min(model_cv$evaluation_log$val_rmse))
  }
  
  # extract the rmse of each parameter set
  lowest_error_df = do.call(rbind, lowest_error_list)
  
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