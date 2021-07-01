toload <- c("grf", "tidyverse", "hdm", "glmnet", "nnls", "Matrix", 
            "matrixStats", "xgboost", "neuralnet", "MASS", "MLmetrics", 
            "keras", "tfdatasets", "data.table", "lessR", "ggthemes", "tictoc",
            "RSQLite")
toinstall <- toload[which(toload %in% installed.packages()[,1] == F)]
lapply(toinstall, install.packages, character.only = TRUE)
lapply(toload, require, character.only = TRUE)

rm(list = ls())
directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

set.seed(12345)

source("../general_functions/general_utils.R")

x = as.matrix(fread(file.path(directory_path, "x_data.csv")))
d = as.matrix(fread(file.path(directory_path, "d_data.csv")))

y = d

# data
column_names = rep(NA, ncol(x))
for(i in 1:ncol(x)){
  column_names[i] = paste0("x",i)
}

fold = sample(seq_len(nrow(x)), size = nrow(x)*0.75)

x_tr = x[fold, ]
x_te = x[-fold, ]

y_tr = y[fold]
y_te = y[-fold]

data = x_tr %>% as_tibble(.name_repair = "minimal") %>%
  setNames(column_names) %>%
  mutate(label = y_tr)

# grid to test on
single_layer_params = list(
  units1 = c(8, 16, 32, 64),
  act.fct1 = c("relu", "sigmoid"),
  act.fctfinal = c("sigmoid")
)

double_layer_params = list(
  units1 = c(8, 16, 32, 64),
  units2 = c(8, 16, 32, 64),
  act.fct1 = c("relu", "sigmoid"),
  act.fct2 = c("relu", "sigmoid"),
  act.fctfinal = c("sigmoid")
)

# expand the grid according to defined hyperparameters
grid_frame_single_layer = expand.grid(single_layer_params)
grid_frame_double_layer = expand.grid(double_layer_params)

# simulations
reps = 10

# result matrices
single_layer_res = matrix(NA, nrow = nrow(grid_frame_single_layer), ncol = 2)
double_layer_res = matrix(NA, nrow = nrow(grid_frame_double_layer), ncol = 2)

# single layer evaluation
for (j in 1:nrow(grid_frame_single_layer)){
  time_vec = rep(NA, reps)
  LogLoss_vec = rep(NA, reps)
  for (i in 1:reps){
    spec = feature_spec(data, label ~.) %>%
      step_numeric_column(all_numeric()) %>%
      fit()
    
    layer = layer_dense_features(
      feature_columns = dense_features(spec),
      dtype = tf$float32
    )
    
    input <- layer_input_from_dataset(data %>% dplyr::select(-label))
    
    output <- input %>% 
      layer_dense_features(dense_features(spec)) %>% 
      layer_dense(units = grid_frame_single_layer$units1[j], activation = grid_frame_single_layer$act.fct1[j])%>%
      layer_dense(units = 1, activation = grid_frame_single_layer$act.fctfinal[j])
    
    model <- keras_model(input, output)
    
    model %>% 
      compile(
        loss = "binary_crossentropy",
        optimizer = optimizer_sgd(lr = 0.01),
        metrics = list("binary_crossentropy")
      )
    
    early_stop <- callback_early_stopping(monitor = "val_loss", patience = 30)
    
    tic()
    
    model %>% fit(
      x = data %>% dplyr::select(-label),
      y = data$label,
      epochs = 500,
      validation_split = 0.2,
      batch_size = 4,
      verbose = 0,
      callbacks = list(early_stop)
    )
    
    temp = toc()
    time_vec[i] = round(temp$toc - temp$tic, 3)
    
    df_te = x_te %>% as_tibble(.name_repair = "minimal") %>%
      setNames(column_names) %>%
      mutate(label = y_te)
    
    preds = model %>% predict(df_te %>% dplyr::select(-label))
    LogLoss_vec[i] = LogLoss(preds, y_te)
    
    average_time = round(mean(time_vec), 2)
    average_LogLoss = round(mean(LogLoss_vec), 2)
    
    single_layer_res[j,1] = average_time
    single_layer_res[j,2] = average_LogLoss
    
    rm(model)
    gc()
  }
}


