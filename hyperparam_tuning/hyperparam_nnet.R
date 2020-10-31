hyperparam_nnet = function(y, x, grid_nn) {
  # potential outcome
  # hyperparameter tuning for neural network
  
  # split data for out of sample evaluation
  n = nrow(x)
  fold = sample(seq_len(nrow(x)), size = n*0.75)

  x_tr = x[fold, ]
  x_te = x[-fold, ]
  
  y_tr = y[fold ]
  y_te = y[-fold]
  
  # define optimal data types for implementation of neuralnet()
  names_nn = colnames(as.data.frame(x_tr))
  train_nn = as.data.frame(cbind(y_tr, x_tr))
  test_X_nn = as.data.frame(x_te)
  test_Y_nn = as.data.frame(y_te)
  
  # define neural network formula
  colnames(train_nn) = c("Y", names_nn)
  nn_formula = as.formula(paste("Y ~", paste(names_nn, collapse = " + ")))
  
  lowest_error_list = rep(NA, nrow(grid_nn))
  
  # calculate a neural network for each set of hyperparameters as defined in the grid
  for (row in 1:nrow(grid_nn)) {
    assign("skip_to_next", FALSE, env=globalenv())
    tryCatch({nncv = neuralnet(formula = nn_formula, 
                        data=train_nn,
                        act.fct = grid_nn$act.fct[row],
                        hidden = grid_nn$neurons[row],
                        stepmax = grid_nn$stepmax[row],
                        linear.output = grid_nn$linear.output[row],
                        err.fct = grid_nn$err.fct[row],
                        threshold = grid_nn$threshold[row],
                        rep = grid_nn$rep[row],
                        algorithm = grid_nn$algorithm[row],
                        lifesign = "full",
                        lifesign.step = 10000
    )}, warning = function(w){w
      assign("skip_to_next", TRUE, env=globalenv())
      }, 
    finally = {
      if (get("skip_to_next", env=globalenv())){
        preds_nn = rep(mean(y_tr), nrow(test_Y_nn))
        lowest_error_list[row] = rmse_calc(test_Y_nn, preds_nn)
        next
      } else {
        preds_nn = predict(nncv, newdata = test_X_nn)
        lowest_error_list[row] = rmse_calc(test_Y_nn, preds_nn)
        next
      }
    })
  }
  
  # combine the errors with the respective parameter set
  gridsearch = cbind(lowest_error_list, grid_nn)
  
  # evaluate best performing set
  bestparams = gridsearch[which.min(lowest_error_list), ]
  
  # output the best hyperparameter set
  finalparams = list(act.fct = bestparams$act.fct, 
                     hidden = bestparams$neurons,
                     stepmax = bestparams$stepmax,
                     linear.output = bestparams$linear.output,
                     err.fct = bestparams$err.fct,
                     threshold = bestparams$threshold,
                     rep = bestparams$rep,
                     algorithm = bestparams$algorithm)
  
  return(finalparams)
}