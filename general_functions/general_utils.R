# calculate the root mean square error of a model's predictions
rmse_calc = function(true_value, predictions) {
  rmse = sqrt(sum((predictions - true_value)^2)/nrow(true_value))
  return(rmse)
}