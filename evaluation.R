#' Evaluation of all simulated results
#' Additionally, the necessary graphical outputs are constructed
#' 

## Load necessary packages, set working directory and seed, remove previously stored variables

toload <- c("grf", "tidyverse", "hdm", "glmnet", "nnls", "Matrix", 
            "matrixStats", "xgboost", "neuralnet", "MASS", "MLmetrics", 
            "keras", "tfdatasets", "data.table", "ggthemes")
toinstall <- toload[which(toload %in% installed.packages()[,1] == F)]
lapply(toinstall, install.packages, character.only = TRUE)
lapply(toload, require, character.only = TRUE)

directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

set.seed(12345)

ml_names = c("Ensemble", "Lasso", "XGBoost", "Neural Network")

# Results simulation 1: -----------------------------------------------------------
# read in the results
folder = "output/sim_50_1"
ate_1 = fread(file = file.path(directory_path, folder, "ate.csv"))
te_ens_1 = fread(file = file.path(directory_path, folder, "te_ens.csv"))
te_lasso_1 = fread(file = file.path(directory_path, folder, "te_lasso.csv"))
te_xgb_1 = fread(file = file.path(directory_path, folder, "te_xgb.csv"))
te_nn_1 = fread(file = file.path(directory_path, folder, "te_nn.csv"))
se_te_1 = fread(file = file.path(directory_path, folder, "se_te.csv"))
ps_ensemble_1 = fread(file = file.path(directory_path, folder, "ps_ensemble.csv"))
oc_ensemble_1 = fread(file = file.path(directory_path, folder, "oc_ensemble.csv"))

## Ensemble bar chart
figure_path = "output/sim_50_1/figures"
# propensity score
ps_plot_1 = ensemble_plot(ps_ensemble_1)
oc_plot_1 = ensemble_plot(oc_ensemble_1)

ggsave(filename = "ps_ensemble_1.png", device = "png", width = 760, height = 200, units = "mm", scale = 0.25, plot = ps_plot_1, dpi = 300, path = file.path(directory_path, figure_path))
ggsave(filename = "oc_ensemble_1.png", device = "png", width = 760, height = 200, units = "mm", scale = 0.25, plot = oc_plot_1, dpi = 300, path = file.path(directory_path, figure_path))

# rmse
rmse_ens_1 = round(rmse_calc(ate_1$ate_t, ate_1$ens), 2)
rmse_lasso_1 = round(rmse_calc(ate_1$ate_t, ate_1$lasso), 2)
rmse_xgb_1 = round(rmse_calc(ate_1$ate_t, ate_1$xgb), 2)
rmse_nn_1 = round(rmse_calc(ate_1$ate_t, ate_1$nn), 2)


rmse_1 = as.data.frame(cbind(rbind(rmse_ens_1, rmse_lasso_1, rmse_xgb_1, rmse_nn_1), ml_names))
rmse_1$V1 = as.numeric(rmse_1$V1)
colnames(rmse_1) = c("RMSE", "Method")
rownames(rmse_1) = ml_names

rmse_plot_1 = ggplot(rmse_1, aes(x = Method, y = RMSE, fill = Method)) +
  theme_tufte(base_family = "serif") +
  geom_col(aes(y = RMSE), stat = "fill") + 
  labs(x = "Methods") +
  coord_cartesian(ylim = c(0,5)) +
  geom_label(aes(x = Method, label = RMSE), position = position_stack(vjust = 1), size = 3) +
  scale_fill_brewer(palette="Greys") +
  geom_hline(yintercept = 0) +
  guides(fill = FALSE)

rmse_plot_1
# coverage probability
# construct the confidence interval
CI_95 = CI(ate$ate_t, ci = .95)

