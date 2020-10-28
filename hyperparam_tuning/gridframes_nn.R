gridframe_oc_1 = function(){
  # define set of possible hyperparameters
  params_nn = list(
    act.fct = c("logistic"),
    algorithm = c("rprop-"),
    neurons = c(6:10),
    threshold = c(0.2, 0.1),
    err.fct = "sse",
    stepmax = 200000,
    linear.output = TRUE,
    rep = c(1)
  )
  
  # expand the grid according to defined hyperparameters
  grid_frame_nn = expand.grid(params_nn)
  return(grid_frame_nn)
}

