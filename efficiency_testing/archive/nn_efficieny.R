column_names = rep(NA, ncol(x))
for(i in 1:ncol(x)){
  column_names[i] = paste0("x",i)
}

fold = sample(seq_len(nrow(x)), size = nrow(x)*0.75)

x_tr = x[fold, ]
x_te = x[-fold, ]

y_tr = y[fold ]
y_te = y[-fold]

data = x_tr %>% as_tibble(.name_repair = "minimal") %>%
  setNames(column_names) %>%
  mutate(label = y_tr)

rm(model)
gc()

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
  layer_dense(units = 32, activation = "relu", kernel_regularizer = regularizer_l1(0.001))%>%
  layer_dense(units = 32, activation = "sigmoid", kernel_regularizer = regularizer_l1(0)) %>%
  layer_dense(units = 1, activation = "sigmoid")

model <- keras_model(input, output)

model %>% 
  compile(
    loss = "binary_crossentropy",
    optimizer = optimizer_adamax(lr = 0.1),
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

toc()

df_te = x_te %>% as_tibble(.name_repair = "minimal") %>%
  setNames(column_names) %>%
  mutate(label = y_te)

preds = model %>% predict(df_te %>% dplyr::select(-label))
LogLoss(preds, y_te)
