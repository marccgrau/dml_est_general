# calculate the root mean square error of a model's predictions
rmse_calc = function(true_value, predictions) {
  rmse = sqrt(sum((predictions - true_value)^2)/length(true_value))
  return(rmse)
}

ensemble_plot = function(ensemble_values){
  df_bc = as.data.frame(t(round(colMeans(ensemble_values), 3)))
  colnames(df_bc) = c("Lasso", "XGBoost", "Neural Network")
  dat = melt(df_bc)
  data = cbind(dat, rep(1, length(ensemble_values)))
  colnames(data) = c("variable", "value", "col_pos")
  data$variable = factor(data$variable, levels = unique(data$variable))

  h_bar = ggplot(data, aes(fill = fct_rev(variable), y = value)) +
    theme_tufte(base_family = "serif") +
    theme(axis.text.y=element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_line(colour = "darkgrey"), 
          axis.line.y = element_line("black"),
          legend.position="bottom",
          legend.title = element_blank(), 
          axis.ticks.x = element_blank()) +
    geom_bar(aes(x = col_pos, y = value), position = "stack", stat = "identity") +
    geom_text(aes(x = col_pos, label = scales::percent(value, accuracy = 0.1)), 
              position = position_stack(vjust = 0.5), size = 3) +
    labs(y = "Percentage", x = "Ensemble") +
    scale_fill_manual(values = rev(c("#D6D6CEFF", "#FFB547FF", "#ADB17DFF")), guide = guide_legend(reverse = TRUE)) +
    scale_x_discrete(breaks = NULL) +
    scale_y_continuous(labels = scales::percent, breaks = c(0, .25, .50, .75, 1),
                       expand = expansion(mult = c(0,0))) +
    geom_hline(yintercept = 0, aes(colour = "black")) +
    geom_hline(yintercept = 1, aes(colour = "black")) +
    coord_flip()
  
  h_bar
}

