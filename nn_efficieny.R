column_names = rep(NA, ncol(x))
for(i in 1:ncol(x)){
  column_names[i] = paste0("x",i)
}

data = x %>% as_tibble(.name_repair = "minimal") %>%
  setNames(column_names) %>%
  mutate(label = y)

spec = feature_spec(data, label ~.) %>%
  step_numeric_column(all_numeric()) %>%
  fit()

layer = layer_dense_features(
  feature_columns = dense_features(spec),
  dtype = tf$float32
)
tic()
model = build_model(data, 
                    spec, 
                    units1 = 128, 
                    units2 = 128, 
                    act.fct1 = "relu", 
                    act.fct2 = "sigmoid",
                    loss.fct = "mse", 
                    eval.metric = "mae",
                    l1_1 = 0,
                    l1_2 = 0)

early_stop <- callback_early_stopping(monitor = "val_mae", patience = 15)

model %>% fit(
  x = data %>% dplyr::select(-label),
  y = data$label,
  epochs = 500,
  validation_split = 0.2,
  batch_size = 4,
  verbose = 2,
  callbacks = list(early_stop)
)

toc()