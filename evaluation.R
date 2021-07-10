#' Evaluation of all simulated results
#' Additionally, the necessary graphical outputs are constructed
#' 

## Load necessary packages, set working directory and seed, remove previously stored variables

toload <- c("grf", "tidyverse", "hdm", "glmnet", "nnls", "Matrix", 
            "matrixStats", "xgboost", "neuralnet", "MASS", "MLmetrics", 
            "keras", "tfdatasets", "data.table", "lessR", "ggthemes", "tictoc",
            "RSQLite", "RColorBrewer", "ggsci")
toinstall <- toload[which(toload %in% installed.packages()[,1] == F)]
lapply(toinstall, install.packages, character.only = TRUE)
lapply(toload, require, character.only = TRUE)

directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

set.seed(12345)

ml_names = c("Ensemble", "Lasso", "XGBoost", "Neural Network")

## Load source functions
# General Functions
source("general_functions/general_utils.R")

# Results simulation 1: -----------------------------------------------------------
# read in the results
folder = "output/sim_1_50"
# start db for 50/50 Simulation
db1_50 = dbConnect(SQLite(), dbname = file.path(folder, "db1_50"))

#exemplary dataset
x_example = fread(file = file.path(folder, "exemplary", "X_data.csv"))

# Treatment effects and prop score
ate_1 = dbGetQuery(db1_50, "SELECT * FROM average_te")
te_ens_1 = as.matrix(dbGetQuery(db1_50, "SELECT * FROM treatmenteffect_ens"))
te_lasso_1 = as.matrix(dbGetQuery(db1_50, "SELECT * FROM treatmenteffect_lasso"))
te_xgb_1 = as.matrix(dbGetQuery(db1_50, "SELECT * FROM treatmenteffect_xgb"))
te_nn_1 = as.matrix(dbGetQuery(db1_50, "SELECT * FROM treatmenteffect_nn"))
true_te_1 = as.matrix(dbGetQuery(db1_50, "SELECT * FROM te_true"))
true_p_1 = as.matrix(dbGetQuery(db1_50, "SELECT * FROM p_true"))
# standard errors
se_te_1 = dbGetQuery(db1_50, "SELECT * FROM standerror_te")
se_po_1 = dbGetQuery(db1_50, "SELECT * FROM standerror_po")
#ensemble weights
ps_ensemble_1 = dbGetQuery(db1_50, "SELECT * FROM propensity_ensemble")
oc_ensemble_1 = dbGetQuery(db1_50, "SELECT * FROM outcome_ensemble")
# Confidence Intervals
CI_up_ens_1 = dbGetQuery(db1_50, "SELECT * FROM CI_up_ens")
CI_down_ens_1 = dbGetQuery(db1_50, "SELECT * FROM CI_down_ens")
CI_up_lasso_1 = dbGetQuery(db1_50, "SELECT * FROM CI_up_lasso")
CI_down_lasso_1 = dbGetQuery(db1_50, "SELECT * FROM CI_down_lasso")
CI_up_xgb_1 = dbGetQuery(db1_50, "SELECT * FROM CI_up_xgb")
CI_down_xgb_1 = dbGetQuery(db1_50, "SELECT * FROM CI_down_xgb")
CI_up_nn_1 = dbGetQuery(db1_50, "SELECT * FROM CI_up_nn")
CI_down_nn_1 = dbGetQuery(db1_50, "SELECT * FROM CI_down_nn")

## Ensemble bar chart
figure_path = "output/sim_1_50/figures"
# propensity score
ps_plot_1 = ensemble_plot(ps_ensemble_1)
oc_plot_1 = ensemble_plot(oc_ensemble_1)

ggsave(filename = "ps_ensemble_1.png", device = "png", width = 760, height = 200, units = "mm", scale = 0.25, plot = ps_plot_1, dpi = 300, path = file.path(directory_path, figure_path))
ggsave(filename = "oc_ensemble_1.png", device = "png", width = 760, height = 200, units = "mm", scale = 0.25, plot = oc_plot_1, dpi = 300, path = file.path(directory_path, figure_path))

# rmse of ate
rmseATE_ens_1 = round(rmse_calc(ate_1$ate_t, ate_1$ens), 2)
rmseATE_lasso_1 = round(rmse_calc(ate_1$ate_t, ate_1$lasso), 2)
rmseATE_xgb_1 = round(rmse_calc(ate_1$ate_t, ate_1$xgb), 2)
rmseATE_nn_1 = round(rmse_calc(ate_1$ate_t, ate_1$nn), 2)


rmseATE_1 = as.data.frame(cbind(rbind(rmseATE_ens_1, rmseATE_lasso_1, rmseATE_xgb_1, rmseATE_nn_1), ml_names))
rmseATE_1$V1 = as.numeric(rmseATE_1$V1)
colnames(rmseATE_1) = c("RMSE", "Method")
rownames(rmseATE_1) = ml_names

rmseATE_plot_1 = ggplot(rmseATE_1, aes(x = Method, y = RMSE, fill = Method)) +
  theme_tufte(base_family = "serif") +
  geom_col(aes(y = RMSE), stat = "fill") + 
  labs(x = "Methods") +
  coord_cartesian(ylim = c(0,0.3)) +
  geom_label(aes(x = Method, label = RMSE), position = position_stack(vjust = 1), size = 3) +
  scale_fill_uchicago(palette = "light") +
  geom_hline(yintercept = 0) +
  guides(fill = FALSE)

rmseATE_plot_1

# bias of ate
biasATE_ens_1 = round(abs(mean(ate_1$ate_t - ate_1$ens)), 2)
biasATE_lasso_1 = round(abs(mean(ate_1$ate_t - ate_1$lasso)), 2)
biasATE_xgb_1 = round(abs(mean(ate_1$ate_t - ate_1$xgb)), 2)
biasATE_nn_1 = round(abs(mean(ate_1$ate_t - ate_1$nn)), 2)

biasATE_1 = as.data.frame(cbind(rbind(biasATE_ens_1, biasATE_lasso_1, biasATE_xgb_1, biasATE_nn_1), ml_names))
biasATE_1$V1 = as.numeric(biasATE_1$V1)
colnames(biasATE_1) = c("Bias", "Method")
rownames(biasATE_1) = ml_names

biasATE_plot_1 = ggplot(biasATE_1, aes(x = Method, y = Bias, fill = Method)) +
  theme_tufte(base_family = "serif") +
  geom_col(aes(y = Bias), stat = "fill") + 
  labs(x = "Methods") +
  coord_cartesian(ylim = c(0,0.1)) +
  geom_label(aes(x = Method, label = Bias), position = position_stack(vjust = 1), size = 3) +
  scale_fill_uchicago(palette = "light") +
  geom_hline(yintercept = 0) +
  guides(fill = FALSE)

biasATE_plot_1

# average rmse of simulations
rmse_ens_1 = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  rmse_ens_1[i] = sqrt(sum((data.table(te_ens_1[i,]) - data.table(true_te_1[i,]))^2)/nrow(data.table(true_te_1[i,])))
}
rmse_lasso_1 = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  rmse_lasso_1[i] = sqrt(sum((data.table(te_lasso_1[i,]) - data.table(true_te_1[i,]))^2)/nrow(data.table(true_te_1[i,])))
}
rmse_xgb_1 = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  rmse_xgb_1[i] = sqrt(sum((data.table(te_xgb_1[i,]) - data.table(true_te_1[i,]))^2)/nrow(data.table(true_te_1[i,])))
}
rmse_nn_1 = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  rmse_nn_1[i] = sqrt(sum((data.table(te_nn_1[i,]) - data.table(true_te_1[i,]))^2)/nrow(data.table(true_te_1[i,])))
}

avgrmse_ens_1 = round(mean(rmse_ens_1), 2)
avgrmse_lasso_1 = round(mean(rmse_lasso_1), 2)
avgrmse_xgb_1 = round(mean(rmse_xgb_1), 2)
avgrmse_nn_1 = round(mean(rmse_nn_1), 2)


avgrmse_1 = as.data.frame(cbind(rbind(avgrmse_ens_1, avgrmse_lasso_1, avgrmse_xgb_1, avgrmse_nn_1), ml_names))
avgrmse_1$V1 = as.numeric(avgrmse_1$V1)
colnames(avgrmse_1) = c("RMSE", "Method")
rownames(avgrmse_1) = ml_names

avgrmse_plot_1 = ggplot(avgrmse_1, aes(x = Method, y = RMSE, fill = Method)) +
  theme_tufte(base_family = "serif") +
  geom_col(aes(y = RMSE), stat = "fill") + 
  labs(x = "Methods") +
  coord_cartesian(ylim = c(0,10)) +
  geom_label(aes(x = Method, label = RMSE), position = position_stack(vjust = 1), size = 3) +
  scale_fill_uchicago(palette = "light") +
  geom_hline(yintercept = 0) +
  guides(fill = FALSE)

avgrmse_plot_1

# development of ensemble
## propensity score
index = seq(1, nrow(ps_ensemble_1), 1)
data_ps_development = melt(ps_ensemble_1)
data_ps_development = cbind.data.frame(index, data_ps_development)

plot_ps_development = ggplot(data_ps_development, aes(x=index, y=value, fill=variable)) +
  theme_tufte(base_family = "serif") +
  theme(legend.position = "bottom") +
  scale_fill_uchicago(palette = "light") +
  geom_area() +
  labs(x = "Simulations", y = "Share in ensemble")

## conditional outcome
index = seq(1, nrow(oc_ensemble_1), 1)
data_oc_development = melt(oc_ensemble_1)
data_oc_development = cbind.data.frame(index, data_oc_development)

plot_oc_development = ggplot(data_oc_development, aes(x=index, y=value, fill=variable)) +
  theme_tufte(base_family = "serif") +
  theme(legend.position = "bottom") +
  scale_fill_uchicago(palette = "light") +
  geom_area() +
  labs(x = "Simulations", y = "Share in ensemble")

# show correlation between input factors and true te
data_corr_map = cbind.data.frame(data.table(true_te_1[1,]), x_example)
colnames(data_corr_map) = c("TE", "X1", "X2", "X3", "X4", "X5", "X6", "X7",
                            "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15")

mtcars %>%
  gather(-mpg, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = mpg)) +
  facet_wrap(~ var, scales = "free") +
  geom_point() +
  stat_smooth()

data_corr_map %>%
  gather(-TE, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = TE)) +
  facet_wrap(~ var, scales = "free") +
  geom_point() +
  stat_smooth()

r <- cor(data_corr_map, use="complete.obs")
r = round(r,2)

ggcorrplot(r, 
           hc.order = TRUE, 
           type = "lower",
           lab = TRUE)

# coverage probability
# construct the confidence interval
#CR_ens %>% 
#  mutate(in_CI = te_ens_1 >= CI_down_ens & te_ens_1 <= CI_up_ens) 
