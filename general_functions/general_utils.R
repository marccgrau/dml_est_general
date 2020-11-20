# calculate the root mean square error of a model's predictions
rmse_calc = function(true_value, predictions) {
  rmse = sqrt(sum((predictions - true_value)^2)/length(true_value))
  return(rmse)
}

ensemble_plot = function(ensemble_values){
  df_bc = as.data.table(round(colMeans(ensemble_values), 3)*100)
  types = as.data.table(c("Lasso", "XGBoost", "Neural Network"))
  ord = seq(1:3)
  
  # geom_text_repel(aes(x = type, label = paste0(value, "%")), position = position_stack(vjust = 0.5), size = 4)
  
  dat = cbind(df_bc, types) %>% gather(type, value)
  dat = cbind(dat, ord)
  colnames(dat) = c("Methods", "type", "value", "ord")
  
  h_bar = ggplot(dat, aes(fill = Methods, y = value, group = order(ord, decreasing = T), order = ord)) +
    theme_tufte(base_family = "serif") +
    theme(axis.text.y=element_blank(), panel.grid.minor = element_blank(), 
          panel.grid.major = element_line(colour = "black"), axis.line.y = element_line("black")) +
    geom_col(aes(x = type), position = "stack", width = 1.17) +
    geom_text(aes(x = type, label = paste0(value, "%")), position = position_stack(vjust = 0.5), size = 3) +
    theme(legend.position="bottom", legend.title = element_blank()) + 
    labs(y = "Percentage", x = "Ensemble") +
    scale_fill_brewer(palette="Greys", direction = 1, labels = c("Lasso", "XGBoost", "Neural Network")) +
    scale_x_discrete(breaks = NULL) +
    scale_y_continuous(labels = scales::percent_format(scale = 1), breaks = c(0, 25, 50, 75, 100)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 100) +
    coord_flip()
  
  h_bar
}

