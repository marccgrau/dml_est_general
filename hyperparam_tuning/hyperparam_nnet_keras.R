hyperparam_nnet_keras = function(y, x, grid_nn) {
  
  column_names = rep(NA, ncol(x))
  for(i in 1:ncol(x)){
    column_names[i] = paste0("x",i)
  }
  
  n = length(y)
  fold = sample(seq_len(nrow(x)), size = n*0.75)
  
  x_tr = x[fold, ]
  x_te = x[-fold, ]
  
  y_tr = y[fold ]
  y_te = y[-fold]
  
  df_tr = x_tr %>% as_tibble(.name_repair = "minimal") %>%
    setNames(column_names) %>%
    mutate(label = y_tr)
  
  df_te = x_te %>% as_tibble(.name_repair = "minimal") %>%
    setNames(column_names) %>%
    mutate(label = y_te)
  
  spec = feature_spec(df_tr, label ~.) %>%
    step_numeric_column(all_numeric()) %>%
    fit()
  
  layer = layer_dense_features(
    feature_columns = dense_features(spec),
    dtype = tf$float32
  )
  
  early_stop <- callback_early_stopping(monitor = "val_loss", patience = 15)
  
  lowest_errors = rep(NA, nrow(grid_nn))
  
  for (i in 1:nrow(grid_nn)){
    model = build_model(df_tr, spec, 
                        units1 = grid_nn$units1[i], 
                        act.fct1 = grid_nn$act.fct1[i], 
                        act.fctfinal = grid_nn$act.fctfinal[i], 
                        loss.fct = grid_nn$loss.fct[i], 
                        eval.metric = grid_nn$eval.metric[i],
                        l1_1 = grid_nn$l1_1[i]
    )
    
    model %>% fit(
      x = df_tr %>% dplyr::select(-label),
      y = df_tr$label,
      epochs = 500,
      validation_split = 0.2,
      batch_size = 4,
      verbose = 0,
      callbacks = list(early_stop)
    )
    
    preds = model %>% predict(df_te %>% dplyr::select(-label))
    if (is.null(grid_nn$act.fctfinal[i])){
      lowest_errors[i] = rmse_calc(y_te, preds)
    } else {
      lowest_errors[i] = LogLoss(preds, y_te)
    }
  }
  
  # combine the errors with the respective parameter set
  gridsearch = cbind(lowest_errors, grid_nn)
  
  # evaluate best performing set
  bestparams = gridsearch[which.min(lowest_errors), ]
  
  # output the best hyperparameter set
  finalparams = list(units1 = bestparams$units1,
                     act.fct1 = bestparams$act.fct1, 
                     act.fctfinal = bestparams$act.fctfinal,  
                     loss.fct = bestparams$loss.fct, 
                     eval.metric = bestparams$eval.metric,
                     l1_1 = bestparams$l1_1
                     )
  
  return(finalparams)
}