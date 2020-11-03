grid_keras_oc_1 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    units1 = c(50, 70, 90),
    units2 = c(50, 70, 90),
    act.fct1 = "sigmoid",
    act.fct2 = "sigmoid",
    act.fctfinal = "sigmoid",
    act.output = FALSE,
    loss.fct = "mse",
    eval.metric = "mean_squared_error"
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}

grid_keras_ps_1 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    units1 = c(50, 100),
    units2 = c(50, 100),
    act.fct1 = c("sigmoid"),
    act.fct2 = c("sigmoid"),
    act.fctfinal = "sigmoid",
    act.output = TRUE,
    loss.fct = "binary_crossentropy",
    eval.metric = "binary_crossentropy"
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}