grid_keras_oc_1_50 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    units1 = c(16),
    units2 = c(16),
    act.fct1 = c("sigmoid"),
    act.fct2 = c("relu"),
    act.fctfinal = NA,
    act.output = FALSE,
    loss.fct = "mse",
    eval.metric = "mean_squared_error",
    l1_1 = c(0, 0.0001, 0.001, 0.01, 0.1)
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}

grid_keras_ps_1_50 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    units1 = c(8),
    act.fct1 = c("relu"),
    act.fctfinal = c("sigmoid"),
    loss.fct = "binary_crossentropy",
    eval.metric = "binary_crossentropy",
    l1_1 = c(0, 0.0001, 0.001, 0.01, 0.1)
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}

grid_keras_oc_1_10 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    units1 = c(64),
    act.fct1 = c("relu"),
    act.fctfinal = NA,
    act.output = FALSE,
    loss.fct = "mse",
    eval.metric = "mean_squared_error",
    l1_1 = c(0, 0.0001, 0.001, 0.01, 0.1)
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}

grid_keras_ps_1_10 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    units1 = c(64),
    act.fct1 = c("sigmoid"),
    act.fctfinal = c("sigmoid"),
    loss.fct = "binary_crossentropy",
    eval.metric = "binary_crossentropy",
    l1_1 = c(0, 0.0001, 0.001, 0.01, 0.1)
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}

grid_keras_oc_2_50 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    units1 = c(16),
    units2 = c(8),
    act.fct1 = c("sigmoid"),
    act.fct2 = c("relu"),
    act.fctfinal = NA,
    act.output = FALSE,
    loss.fct = "mse",
    eval.metric = "mean_squared_error",
    l1_1 = c(0, 0.0001, 0.001, 0.01, 0.1)
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}

grid_keras_ps_2_50 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    units1 = c(32),
    act.fct1 = c("relu"),
    act.fctfinal = c("sigmoid"),
    loss.fct = "binary_crossentropy",
    eval.metric = "binary_crossentropy",
    l1_1 = c(0, 0.0001, 0.001, 0.01, 0.1)
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}

grid_keras_oc_2_10 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    units1 = c(16),
    units2 = c(8),
    act.fct1 = c("sigmoid"),
    act.fct2 = c("relu"),
    act.fctfinal = NA,
    act.output = FALSE,
    loss.fct = "mse",
    eval.metric = "mean_squared_error",
    l1_1 = c(0, 0.0001, 0.001, 0.01, 0.1)
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}

grid_keras_ps_2_10 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    units1 = c(32),
    act.fct1 = c("relu"),
    act.fctfinal = c("sigmoid"),
    loss.fct = "binary_crossentropy",
    eval.metric = "binary_crossentropy",
    l1_1 = c(0, 0.0001, 0.001, 0.01, 0.1)
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}

