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
              position = position_stack(vjust = 0.5), size = 2.5) +
    labs(y = NULL, x = "Ensemble") +
    scale_fill_manual(values = rev(c("#D6D6CEFF", "#FFB547FF", "#ADB17DFF")), guide = guide_legend(reverse = TRUE)) +
    scale_x_discrete(breaks = NULL) +
    scale_y_continuous(labels = scales::percent, breaks = c(0, .25, .50, .75, 1),
                       expand = expansion(mult = c(0,0))) +
    geom_hline(yintercept = 0, aes(colour = "black")) +
    geom_hline(yintercept = 1, aes(colour = "black")) +
    coord_flip()
  
  h_bar
}

require(ggplot2)

label_size = 4.5

theme_MGthesis <- function(base_size = 22.5, base_family = "",
                       base_line_size = 0.5,
                       base_rect_size = 0.5, 
                       percentage = FALSE) {
  
  font = "serif"
  # Starts with theme_grey and then modify some parts
  
  
  theme_grey(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace%
    theme(
      # white background and dark border
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border     = element_rect(fill = NA, colour = "grey20"),
      # make gridlines dark, same contrast with white as in theme_grey
      panel.grid = element_line(colour = "grey92"),
      panel.grid.minor = element_line(size = rel(0.5)),
      # contour strips to match panel contour
      strip.background = element_rect(fill = "grey85", colour = "grey20"),
      # match legend key to background
      legend.key       = element_rect(fill = "white", colour = NA),
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 14,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0.5,                #left align
        vjust = 1,
        margin = margin(0,0,15,0)),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 16),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 16,                 #font size
        hjust = 1),               #right align
      
      axis.line = element_line(color = "darkgrey",
                               size = 0.5),
      
      axis.ticks = element_line(color = "darkgrey"),
      
      axis.title.x = element_text(             #axis titles
        family = font,            #font family
        size = 16,
        margin = margin(15,0,0,0), 
        color = "black"),               #font size
      
      axis.title.y = element_text(             #axis titles
        family = font,            #font family
        size = 16,
        margin = margin(0,15,0,0), 
        color = "black"),
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 16, 
        color = "black"),                #font size
      
      strip.text.x = element_text(              #axis text
        family = font,            #axis famuly
        size = 16, 
        color = "black"),
      
      legend.position = "bottom", 
      legend.text = element_text(size = 16),
      legend.title=element_text(size=16, face = 'bold'),
      
      text = element_text(family = "serif")
      
    )
  
  
}
