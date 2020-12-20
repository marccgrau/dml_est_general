grid_keras_oc_1 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    units1 = c(32),
    units2 = c(32),
    act.fct1 = c("relu"),
    act.fct2 = c("sigmoid"),
    act.fctfinal = NA,
    act.output = FALSE,
    loss.fct = "mse",
    eval.metric = "mean_absolute_error",
    l1_1 = c(0.0001),
    l1_2 = c(0.0001)
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}

grid_keras_ps_1 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    units1 = c(32),
    units2 = c(32),
    act.fct1 = c("relu"),
    act.fct2 = c("sigmoid"),
    act.fctfinal = c("sigmoid"),
    loss.fct = "binary_crossentropy",
    eval.metric = "binary_crossentropy",
    l1_1 = c(0.0001),
    l1_2 = c(0.0001)
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}