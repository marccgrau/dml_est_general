hyperparam_nnet_oc = function(y, x) {
  
  n = nrow(x)
  fold = sample(seq_len(nrow(x)), size = n*0.75)
  
  x_tr = x[fold, ]
  x_te = x[-fold, ]
  
  y_tr = y[fold ]
  y_te = y[-fold]
  
  names_nn = colnames(as.data.frame(x_tr))
  train_nn = as.data.frame(cbind(y_tr, x_tr))
  test_X_nn = as.data.frame(x_te)
  test_Y_nn = as.data.frame(y_te)
  
  colnames(train_nn) = c("Y", names_nn)
  nn_formula = as.formula(paste("Y ~", paste(names_nn, collapse = " + ")))
  
  params_nn = list(
    act.fct = c("tanh", "logistic"),
    algorithm = c("rprop+"),
    neurons = c(4:6),
    threshold = 0.8,
    learningrate.limit = c(0.01, 0.02, 0.03),
    err.fct = "sse",
    stepmax = 200000,
    linear.output = TRUE,
    rep = c(1)
  )
  
  grid_frame_nn = expand.grid(params_nn)
  lowest_error_list = list()
  
  for (row in 1:nrow(grid_frame_nn)) {
    nncv <- neuralnet(formula = nn_formula, 
                        data=train_nn,
                        act.fct = grid_frame_nn$act.fct[row],
                        hidden = grid_frame_nn$neurons[row],
                        stepmax = grid_frame_nn$stepmax[row],
                        linear.output = grid_frame_nn$linear.output[row],
                        err.fct = grid_frame_nn$err.fct[row],
                        threshold = grid_frame_nn$threshold[row],
                        learningrate = grid_frame_nn$learningrate[row],
                        rep = grid_frame_nn$rep[row],
                        algorithm = grid_frame_nn$algorithm[row],
                        lifesign = "full",
                        lifesign.step = 10000
    )
    preds_nn = predict(nncv, newdata = test_X_nn)
    lowest_error_list[[row]] = rmse_calc(test_Y_nn, preds_nn)
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
                     learningrate = bestparams$learningrate,
                     rep = bestparams$rep,
                     algorithm = bestparams$algorithm)
  
  return(finalparams)
}