hyperparam_nnet = function(y_tr, x_tr, y_te, x_te) {
  
  names_nn = colnames(as.data.frame(x_tr))
  train_nn = as.data.frame(cbind(y_tr, x_tr))
  test_X_nn = as.data.frame(x_te)
  test_Y_nn = as.data.frame(y_te)
  
  colnames(train_nn) = c("Y", names_nn)
  nn_formula = as.formula(paste("Y ~", paste(names_nn, collapse = " + ")))
  
  params_nn = list(
    act.fct = c("tanh", "logistic"),
    neurons = c(5:8),
    threshold = c(0.02, 0.03),
    err.fct = "sse",
    stepmax = 100000,
    linear.output = TRUE,
    rep = c(1:3)
  )
  
  grid_frame_nn = expand.grid(params_nn)
  
  for (row in 1:nrow(grid_frame_nn)) {
    nncv <- neuralnet(formula = nn_formula, 
                        data=train_nn,
                        act.fct = grid_frame_nn$act.fct[row],
                        hidden = grid_frame_nn$neurons[row],
                        stepmax = grid_frame_nn$stepmax[row],
                        linear.output = grid_frame_nn$linear.output[row],
                        err.fct = grid_frame_nn$err.fct[row],
                        threshold = grid_frame_nn$threshold[row],
                        rep = grid_frame_nn$rep[row]
    )
    pred_nn = predict(nncv, newdata = test_X_nn)
    lowest_error_list_Y_nn[[row]] = rmse_calc(test_Y_nn, preds_nn)
  }
  
  ### Create object that contains all accuracy's
  lowest_error_df = do.call(rbind, lowest_error_list)
  
  ### Bind columns of accuracy values and random hyperparameter values
  gridsearch = cbind(lowest_error_df, grid_frame_nn)
  
  ### Quickly display highest accuracy
  bestparams = gridsearch[which.min(lowest_error_df), ]
  
  finalparams = list(act.fct = bestparams$act.fct, 
                     hidden = bestparams$hidden,
                     stepmax = bestparams$stepmax,
                     linear.output = bestparams$linear.output,
                     err.fct = bestparams$err.fct,
                     threshold = bestparams$threshold,
                     rep = bestparams$rep)
  
  return(finalparams)
}