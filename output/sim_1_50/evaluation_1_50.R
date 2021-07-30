#' Evaluation of all simulated results
#' Additionally, the necessary graphical outputs are constructed
#' 

## Load necessary packages, set working directory and seed, remove previously stored variables

toload <- c("grf", "tidyverse", "hdm", "glmnet", "nnls", "Matrix", 
            "matrixStats", "xgboost", "neuralnet", "MASS", "MLmetrics", 
            "keras", "tfdatasets", "data.table", "lessR", "ggthemes", "tictoc",
            "RSQLite", "RColorBrewer", "ggsci", "ggcorrplot", "grid", "gridExtra",
            "extrafont", "normtest", "stargazer")
toinstall <- toload[which(toload %in% installed.packages()[,1] == F)]
lapply(toinstall, install.packages, character.only = TRUE)
lapply(toload, require, character.only = TRUE)

directory_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(directory_path)

set.seed(12345)

ml_names = c("Ensemble", "Lasso", "XGBoost", "Neural Network")

## Load source functions
# General Functions
source(file.path(dirname(dirname(directory_path)), "general_functions/general_utils.R"))

# Results simulation 1: -----------------------------------------------------------
# read in the results
# start db for 50/50 Simulation
db1_50 = dbConnect(SQLite(), dbname = "db1_50")

n_simulations = 500

figure_path = "figures"

#exemplary dataset
x_example = fread(file = file.path("exemplary", "X_data.csv"))
y_example = fread(file = file.path("exemplary", "Y_data.csv"))
d_example = fread(file = file.path("exemplary", "D_data.csv"))

# Treatment effects and prop score
ate = dbGetQuery(db1_50, "SELECT * FROM average_te")
te_ens = as.matrix(dbGetQuery(db1_50, "SELECT * FROM treatmenteffect_ens"))
te_lasso = as.matrix(dbGetQuery(db1_50, "SELECT * FROM treatmenteffect_lasso"))
te_xgb = as.matrix(dbGetQuery(db1_50, "SELECT * FROM treatmenteffect_xgb"))
te_nn = as.matrix(dbGetQuery(db1_50, "SELECT * FROM treatmenteffect_nn"))
true_te = as.matrix(dbGetQuery(db1_50, "SELECT * FROM te_true"))
true_p = as.matrix(dbGetQuery(db1_50, "SELECT * FROM p_true"))
# standard errors
se_te = dbGetQuery(db1_50, "SELECT * FROM standerror_te")
se_po = dbGetQuery(db1_50, "SELECT * FROM standerror_po")
#ensemble weights
ps_ensemble = dbGetQuery(db1_50, "SELECT * FROM propensity_ensemble")
oc_ensemble = dbGetQuery(db1_50, "SELECT * FROM outcome_ensemble")
# Confidence Intervals
CI_up_ens = dbGetQuery(db1_50, "SELECT * FROM CI_up_ens")
CI_down_ens = dbGetQuery(db1_50, "SELECT * FROM CI_down_ens")
CI_up_lasso = dbGetQuery(db1_50, "SELECT * FROM CI_up_lasso")
CI_down_lasso = dbGetQuery(db1_50, "SELECT * FROM CI_down_lasso")
CI_up_xgb = dbGetQuery(db1_50, "SELECT * FROM CI_up_xgb")
CI_down_xgb = dbGetQuery(db1_50, "SELECT * FROM CI_down_xgb")
CI_up_nn = dbGetQuery(db1_50, "SELECT * FROM CI_up_nn")
CI_down_nn = dbGetQuery(db1_50, "SELECT * FROM CI_down_nn")

# overlap
data_overlap = cbind.data.frame(d_example, data.table(true_p[1,]))
colnames(data_overlap) = c("group", "propensity")
data_overlap$group = factor(data_overlap$group, levels = unique(data_overlap$group))

plot_overlap = ggplot(data_overlap, aes(x = propensity, fill = group)) +
  theme_tufte(base_family = "serif") +
  theme(
        panel.border = element_rect(colour = "darkgrey", fill=NA, size=0.5),
        axis.ticks = element_line(colour = "darkgrey",
                                  size = 0.3),
        legend.title = element_blank(),
        legend.position="bottom") +
  geom_density(alpha=0.6, position = 'identity', bins = 50) +
  scale_fill_uchicago(labels = c("Treated", "Controls")) +
  labs(y = "Density", x = "Propensity Score")

pdf_file = file.path(figure_path, paste0("plot_overlap_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 7,
    height = 4)

plot_overlap

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

## Ensemble bar chart
# propensity score
ps_plot = ensemble_plot(ps_ensemble)
oc_plot = ensemble_plot(oc_ensemble)

# save plots as pdf
pdf_file = file.path(figure_path, paste0("ps_plot_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 7,
    height = 2)

ps_plot

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

pdf_file = file.path(figure_path, paste0("oc_plot_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 7,
    height = 2)

oc_plot

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

# rmse of ate
rmseATE_ens = round(rmse_calc(ate$ate_t, ate$ens), 2)
rmseATE_lasso = round(rmse_calc(ate$ate_t, ate$lasso), 2)
rmseATE_xgb = round(rmse_calc(ate$ate_t, ate$xgb), 2)
rmseATE_nn = round(rmse_calc(ate$ate_t, ate$nn), 2)


rmseATE = as.data.frame(cbind(rbind(rmseATE_ens, rmseATE_lasso, rmseATE_xgb, rmseATE_nn), ml_names))
rmseATE$V1 = as.numeric(rmseATE$V1)
colnames(rmseATE) = c("RMSE", "Method")
rownames(rmseATE) = ml_names

rmseATE_plot = ggplot(rmseATE, aes(x = Method, y = RMSE, fill = Method)) +
  theme_tufte(base_family = "serif") +
  geom_col(aes(y = RMSE), stat = "fill") + 
  labs(x = "Methods") +
  coord_cartesian(ylim = c(0,0.3)) +
  geom_label(aes(x = Method, label = RMSE), position = position_stack(vjust = 1), size = 3) +
  scale_fill_uchicago(palette = "light") +
  geom_hline(yintercept = 0) +
  guides(fill = FALSE)



pdf_file = file.path(figure_path, paste0("rmse_ate_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 7,
    height = 4)

rmseATE_plot

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

# metrics for ate
MSE_ate_ens = sum((ate$ens - ate$ate_t)^2)/nrow(ate)
MSE_ate_lasso = sum((ate$lasso - ate$ate_t)^2)/nrow(ate)
MSE_ate_xgb = sum((ate$xgb - ate$ate_t)^2)/nrow(ate)
MSE_ate_nn = sum((ate$nn - ate$ate_t)^2)/nrow(ate)

BIAS_ate_ens = sum(abs(ate$ens - ate$ate_t))/nrow(ate)
BIAS_ate_lasso = sum(abs(ate$lasso - ate$ate_t))/nrow(ate)
BIAS_ate_xgb = sum(abs(ate$xgb - ate$ate_t))/nrow(ate)
BIAS_ate_nn = sum(abs(ate$nn - ate$ate_t))/nrow(ate)

SD_ate_ens = sqrt(sum((ate$ens - mean(ate$ens))^2)/(nrow(ate) - 1))
SD_ate_lasso = sqrt(sum((ate$lasso - mean(ate$lasso))^2)/(nrow(ate) - 1))
SD_ate_xgb = sqrt(sum((ate$xgb - mean(ate$xgb))^2)/(nrow(ate) - 1))
SD_ate_nn = sqrt(sum((ate$nn - mean(ate$nn))^2)/(nrow(ate) - 1))

# bias of ate

biasATE = as.data.frame(cbind(rbind(BIAS_ate_ens, BIAS_ate_lasso, BIAS_ate_xgb, BIAS_ate_nn), ml_names))
biasATE$V1 = as.numeric(biasATE$V1)
colnames(biasATE) = c("Bias", "Method")
rownames(biasATE) = ml_names

biasATE_plot = ggplot(biasATE, aes(x = Method, y = Bias, fill = Method)) +
  theme_tufte(base_family = "serif") +
  geom_col(aes(y = Bias), stat = "fill") + 
  labs(x = "Methods") +
  coord_cartesian(ylim = c(0,0.1)) +
  geom_label(aes(x = Method, label = Bias), position = position_stack(vjust = 1), size = 3) +
  scale_fill_uchicago(palette = "light") +
  geom_hline(yintercept = 0) +
  guides(fill = FALSE)

pdf_file = file.path(figure_path, paste0("bias_ate_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 7,
    height = 4)

biasATE_plot

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

# violin plot ATE
plot_data_ate_violin = melt(ate[,-1])

ate_violin = ggplot(plot_data_ate_violin, aes(x=variable, y=value, fill = variable)) + 
  theme_tufte(base_family = "serif") +
  theme(panel.grid.major.y = element_line(linetype = "dotted",
                                          colour = "darkgrey",
                                          size = 0.3),
        panel.grid.minor.y = element_line(linetype = "dotted",
                                          colour = "darkgrey",
                                          size = 0.3),
        panel.border = element_rect(colour = "darkgrey", fill=NA, size=0.5),
        axis.ticks = element_line(colour = "darkgrey",
                                  size = 0.3),
        legend.title = element_blank(),
        legend.position="bottom", 
        strip.text = element_blank()) +
  facet_wrap(.~variable, scales = "free", nrow = 1, ncol = 4) +
  geom_violin(trim=TRUE, show.legend = FALSE) +
  geom_boxplot(width=0.12, show.legend = FALSE) +
  scale_fill_uchicago(palette = "light",
                      labels = c("Ensemble", "Lasso", "XGBoost", "Neural Network")) +
  geom_hline(yintercept = 0.5) + 
  labs(x = NULL, y = "Average Treatment Effect") +
  scale_x_discrete(labels=c("ens" = "Ensemble", "lasso" = "Lasso",
                            "xgb" = "XGBoost", "nn" = "Neural Network"))


pdf_file = file.path(figure_path, paste0("violin_ate_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 7,
    height = 3)

ate_violin

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

# all error metrics combined 

all_errors_ate = as.data.frame(cbind(ml_names, 
                                 rbind(round(MSE_ate_ens, 4), 
                                       round(MSE_ate_lasso, 4),
                                       round(MSE_ate_xgb, 4),
                                       round(MSE_ate_nn, 4)), 
                                 rbind(round(BIAS_ate_ens, 2), 
                                       round(BIAS_ate_lasso, 2), 
                                       round(BIAS_ate_xgb, 2), 
                                       round(BIAS_ate_nn, 2)), 
                                 rbind(round(SD_ate_ens, 2), 
                                       round(SD_ate_lasso, 2), 
                                       round(SD_ate_xgb, 2), 
                                       round(SD_ate_nn, 2))))
colnames(all_errors_ate) = c("methods", "MSE", "Bias", "SD")

plot_data_all_erros_ate = melt(all_errors_ate, id.vars = "methods")
plot_data_all_erros_ate$methods = factor(plot_data_all_erros_ate$methods, levels = unique(plot_data_all_erros_ate$methods))
plot_data_all_erros_ate$value = as.numeric(plot_data_all_erros_ate$value)


plot_all_errors_ate = ggplot(plot_data_all_erros_ate, aes(group = variable, y = value, fill = methods)) +
  theme_tufte(base_family = "serif") +
  theme(panel.grid.major.y = element_line(linetype = "dotted",
                                          colour = "darkgrey",
                                          size = 0.3),
        panel.grid.minor.y = element_line(linetype = "dotted",
                                          colour = "darkgrey",
                                          size = 0.3),
        panel.border = element_rect(colour = "darkgrey", fill=NA, size=0.5),
        axis.ticks = element_line(colour = "darkgrey",
                                  size = 0.3),
        legend.title = element_blank(),
        legend.position="bottom") +
  geom_bar(aes(x = methods, fill = variable, y = value), position = "dodge", stat = "identity") +
  labs(x = NULL, y = NULL ) +
  geom_text(aes(x = methods, y = value, fill = variable,
            label = value,
            vjust = -1), 
            position = position_dodge(width= 0.9), 
            show.legend = FALSE, 
            size = 3.5,
            family = "serif") +
  scale_fill_uchicago(palette = "light",
                      labels = c("MSE", "Absolute Bias", "Standard Deviation")) +
  scale_y_continuous(expand = expansion(mult = c(0, .2)))

pdf_file = file.path(figure_path, paste0("plot_all_errors_ate_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 7,
    height = 3)

plot_all_errors_ate

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

# average rmse of simulations
rmse_ens = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  rmse_ens[i] = sqrt(sum((data.table(te_ens[i,]) - data.table(true_te[i,]))^2)/nrow(data.table(true_te[i,])))
}
rmse_lasso = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  rmse_lasso[i] = sqrt(sum((data.table(te_lasso[i,]) - data.table(true_te[i,]))^2)/nrow(data.table(true_te[i,])))
}
rmse_xgb = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  rmse_xgb[i] = sqrt(sum((data.table(te_xgb[i,]) - data.table(true_te[i,]))^2)/nrow(data.table(true_te[i,])))
}
rmse_nn = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  rmse_nn[i] = sqrt(sum((data.table(te_nn[i,]) - data.table(true_te[i,]))^2)/nrow(data.table(true_te[i,])))
}

avgrmse_ens = round(mean(rmse_ens), 2)
avgrmse_lasso = round(mean(rmse_lasso), 2)
avgrmse_xgb = round(mean(rmse_xgb), 2)
avgrmse_nn = round(mean(rmse_nn), 2)


avgrmse = as.data.frame(cbind(rbind(avgrmse_ens, avgrmse_lasso, avgrmse_xgb, avgrmse_nn), ml_names))
avgrmse$V1 = as.numeric(avgrmse$V1)
colnames(avgrmse) = c("RMSE", "Method")
rownames(avgrmse) = ml_names

avgrmse_plot = ggplot(avgrmse, aes(x = Method, y = RMSE, fill = Method)) +
  theme_tufte(base_family = "serif") +
  geom_col(aes(y = RMSE), stat = "fill") + 
  labs(x = "Methods") +
  coord_cartesian(ylim = c(0,10)) +
  geom_label(aes(x = Method, label = RMSE), position = position_stack(vjust = 1), size = 3) +
  scale_fill_uchicago(palette = "light") +
  geom_hline(yintercept = 0) +
  guides(fill = FALSE)

pdf_file = file.path(figure_path, paste0("rmse_ate_eachsim_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 7,
    height = 4)

avgrmse_plot

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

# average mse across all repetitions ----
# average mse per replication
MSE_true_te = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  MSE_true_te[i] = sum((data.table(true_te[i,]) - data.table(true_te[i,]))^2)/nrow(data.table(true_te[i,]))
}
MSE_ens = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  MSE_ens[i] = sum((data.table(te_ens[i,]) - data.table(true_te[i,]))^2)/nrow(data.table(true_te[i,]))
}
MSE_lasso = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  MSE_lasso[i] = sum((data.table(te_lasso[i,]) - data.table(true_te[i,]))^2)/nrow(data.table(true_te[i,]))
}
MSE_xgb = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  MSE_xgb[i] = sum((data.table(te_xgb[i,]) - data.table(true_te[i,]))^2)/nrow(data.table(true_te[i,]))
}
MSE_nn = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  MSE_nn[i] = sum((data.table(te_nn[i,]) - data.table(true_te[i,]))^2)/nrow(data.table(true_te[i,]))
}

avgMSE_ens = round(mean(MSE_ens), 2)
avgMSE_lasso = round(mean(MSE_lasso), 2)
avgMSE_xgb = round(mean(MSE_xgb), 2)
avgMSE_nn = round(mean(MSE_nn), 2)


avgMSE = as.data.frame(cbind(ml_names, rbind(avgMSE_ens, avgMSE_lasso, avgMSE_xgb, avgMSE_nn)))
colnames(avgMSE) = c("Method", "MSE")
avgMSE$MSE = as.numeric(avgMSE$MSE)
rownames(avgMSE) = ml_names

avgMSE_plot = ggplot(avgMSE, aes(x = Method, y = MSE, fill = Method)) +
  theme_tufte(base_family = "serif") +
  geom_col(aes(y = MSE), stat = "fill") + 
  labs(x = "Methods") +
  coord_cartesian(ylim = c(0,10)) +
  geom_label(aes(x = Method, label = MSE), position = position_stack(vjust = 1), size = 3) +
  scale_fill_uchicago(palette = "light") +
  geom_hline(yintercept = 0) +
  guides(fill = FALSE)

pdf_file = file.path(figure_path, paste0("MSE_ate_eachsim_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 7,
    height = 4)

avgMSE_plot

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

# average bias per replication
BIAS_ens = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  temp = data.table(te_ens[i,]) - data.table(true_te[i,])
  BIAS_ens[i] = mean(abs(temp$V1))
}
BIAS_lasso = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  temp = data.table(te_lasso[i,]) - data.table(true_te[i,])
  BIAS_lasso[i] = mean(abs(temp$V1))
}
BIAS_xgb = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  temp = data.table(te_xgb[i,]) - data.table(true_te[i,])
  BIAS_xgb[i] = mean(abs(temp$V1))
}
BIAS_nn = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  temp = data.table(te_nn[i,]) - data.table(true_te[i,])
  BIAS_nn[i] = mean(abs(temp$V1))
}

avgBIAS_ens = round(mean(BIAS_ens), 2)
avgBIAS_lasso = round(mean(BIAS_lasso), 2)
avgBIAS_xgb = round(mean(BIAS_xgb), 2)
avgBIAS_nn = round(mean(BIAS_nn), 2)

avgBIAS = rbind.data.frame(avgBIAS_ens, avgBIAS_lasso, avgBIAS_xgb, avgBIAS_nn)
colnames(avgBIAS) = c("Bias")

# average SD per replication
SD_true = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  temp = data.table(true_te[i,])
  mean_temp = mean(temp$V1)
  SD_true[i] = sqrt(sum((temp - mean_temp)^2)/(nrow(temp) - 1))
}

SD_ens = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  temp = data.table(te_ens[i,])
  mean_temp = mean(temp$V1)
  SD_ens[i] = sqrt(sum((temp - mean_temp)^2)/(nrow(temp) - 1))
}
SD_lasso = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  temp = data.table(te_lasso[i,])
  mean_temp = mean(temp$V1)
  SD_lasso[i] = sqrt(sum((temp - mean_temp)^2)/(nrow(temp) - 1))
}
SD_xgb = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  temp = data.table(te_xgb[i,])
  mean_temp = mean(temp$V1)
  SD_xgb[i] = sqrt(sum((temp - mean_temp)^2)/(nrow(temp) - 1))
}
SD_nn = matrix(NA, ncol = 1, nrow = n_simulations)
for (i in 1:n_simulations){
  temp = data.table(te_nn[i,])
  mean_temp = mean(temp$V1)
  SD_nn[i] = sqrt(sum((temp - mean_temp)^2)/(nrow(temp) - 1))
}

avgSD_ens = round(mean(SD_ens), 2)
avgSD_lasso = round(mean(SD_lasso), 2)
avgSD_xgb = round(mean(SD_xgb), 2)
avgSD_nn = round(mean(SD_nn), 2)

avgSD = rbind.data.frame(avgSD_ens, avgSD_lasso, avgSD_xgb, avgSD_nn)
colnames(avgSD) = c("SD")

JB_ens = matrix(NA, ncol = 1, nrow = n_simulations)
JB_test_ens = numeric(2000)
for (i in 1:n_simulations){
  temp = data.table(te_ens[i,])
  JB_test_ens[i] = jb.norm.test(temp$V1, nrepl = 1)$p.value
  JB_ens[i] = ifelse(JB_test_ens > 0.05, 1, 0)
}

JB_lasso = matrix(NA, ncol = 1, nrow = n_simulations)
JB_test_lasso = numeric(2000)
for (i in 1:n_simulations){
  temp = data.table(te_lasso[i,])
  JB_test_lasso[i] = jb.norm.test(temp$V1, nrepl = 1)$p.value
  JB_lasso[i] = ifelse(JB_test_lasso > 0.05, 1, 0)
}

JB_xgb = matrix(NA, ncol = 1, nrow = n_simulations)
JB_test_xgb = numeric(2000)
for (i in 1:n_simulations){
  temp = data.table(te_xgb[i,])
  JB_test_xgb[i] = jb.norm.test(temp$V1, nrepl = 1)$p.value
  JB_xgb[i] = ifelse(JB_test_xgb > 0.05, 1, 0)
}

JB_nn = matrix(NA, ncol = 1, nrow = n_simulations)
JB_test_nn = numeric(2000)
for (i in 1:n_simulations){
  temp = data.table(te_nn[i,])
  JB_test_nn[i] = jb.norm.test(temp$V1, nrepl = 1)$p.value
  JB_nn[i] = ifelse(JB_test_nn > 0.05, 1, 0)
}


# put all metrics together
avgMetrics = cbind.data.frame(avgMSE, avgBIAS, avgSD)

# development of ensemble ----
## propensity score
index = seq(1, nrow(ps_ensemble), 1)
data_ps_development = melt(ps_ensemble)
data_ps_development = cbind.data.frame(index, data_ps_development)

plot_ps_development = ggplot(data_ps_development, aes(x=index, y=value, fill=variable)) +
  theme_tufte(base_family = "serif") +
  theme(legend.position = "bottom") +
  scale_fill_uchicago(palette = "light") +
  geom_area() +
  labs(x = "Simulations", y = "Share in ensemble")

pdf_file = file.path(figure_path, paste0("ps_development_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 13,
    height = 4)

plot_ps_development

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

## conditional outcome
index = seq(1, nrow(oc_ensemble), 1)
data_oc_development = melt(oc_ensemble)
data_oc_development = cbind.data.frame(index, data_oc_development)

plot_oc_development = ggplot(data_oc_development, aes(x=index, y=value, fill=variable)) +
  theme_tufte(base_family = "serif") +
  theme(legend.position = "bottom") +
  scale_fill_uchicago(palette = "light") +
  geom_area() +
  labs(x = "Simulations", y = "Share in ensemble")

pdf_file = file.path(figure_path, paste0("oc_development_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 13,
    height = 4)

plot_oc_development

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

# show correlation between confounders and true te

data_corr_map = cbind.data.frame(data.table(true_te[1,]), x_example)
colnames(data_corr_map) = c("TE", "X1", "X2", "X3", "X4", "X5", "X6", "X7",
                            "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15")


data_corr_map_norms = data_corr_map %>% 
  dplyr::select(c(TE, X1, X3, X5, X7, X9, X11, X13, X15)) %>%
  gather(-TE, key = "var", value = "value")

data_corr_map_bins = data_corr_map %>% 
  dplyr::select(c(TE, X2, X4, X6, X8, X10, X12, X14)) %>%
  gather(-TE, key = "var", value = "value")

plot_norms = ggplot(data_corr_map_norms, aes(x = value, y = TE)) +
  facet_wrap(~ factor(var, levels = unique(var)), scales = "free") +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  stat_smooth() + 
  theme_bw(base_family = "serif") +
  theme(text = element_text(size=12)) +
  labs(x = "Confounders", y = "Treatment Effect")

plot_bins = ggplot(data_corr_map_bins, aes(x = TE, y = value)) +
  facet_wrap(~ factor(var, levels = unique(var)), scales = "free") +
  geom_point(alpha = 0.3) +
  stat_smooth() + 
  theme_bw(base_family = "serif") +
  theme(text = element_text(size=12)) +
  labs(x = "Treatment Effect", y = "Confounders")

pdf_file = file.path(figure_path, paste0("corr_plot_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 13,
    height = 15)

grid.arrange(plot_norms, plot_bins, nrow = 2)

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

r <- cor(data_corr_map, use="complete.obs")
r = round(r,2)

pdf_file = file.path(figure_path, paste0("corr_map_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 7,
    height = 6)

ggcorrplot(r, 
           hc.order = FALSE, 
           type = "lower",
           lab = TRUE)

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")


# Coverage probability, construct each interval

alpha = 0.05 # quantile necessary for construction
n = 2000 # number of observations
m = 500 # number of repetitions in MCS

## ensemble
upper_limit_ens = numeric(m)
lower_limit_ens = numeric(m)
mean_ate_ens = numeric(m)

for(i in 1:m)
{ 
  x = t(t(te_ens[i,]))
  lower_limit_ens[i] = mean(x) - qt(alpha / 2, df=n-1, lower.tail = FALSE)*sd(x)/sqrt(n)
  upper_limit_ens[i] = mean(x) + qt(alpha / 2, df=n-1, lower.tail = FALSE)*sd(x)/sqrt(n)
  mean_ate_ens[i] = mean(x)
}

inCI_ens = ifelse(lower_limit_ens < 0.5 & upper_limit_ens > 0.5, 1, 0)
plot_data_CI_ens = cbind.data.frame(c(1:500), mean_ate_ens, lower_limit_ens, upper_limit_ens, inCI_ens)
colnames(plot_data_CI_ens) = c("Index", "Mean", "Lower", "Upper", "inCI")
coverage_prob_ens = mean(lower_limit_ens < 0.5 & upper_limit_ens > 0.5)

plot_coverage_ens = ggplot(plot_data_CI_ens, aes(x = Index, y = Mean)) +
  theme_tufte(base_family = "serif") +
  geom_errorbar(aes(ymax = Upper, ymin = Lower, color = as.factor(inCI)), show.legend = FALSE) +
  scale_color_manual(values = c("blue", "grey")) +
  geom_point(size = 1, aes(color = as.factor(inCI)), show.legend = FALSE) +
  scale_fill_manual(values = c("darkblue", "darkgrey")) +
  geom_hline(yintercept = 0.5, color = "red") +
  labs(y = "Estimated treatment effect", x = "Simulation repetition", title = "Ensemble")

## lasso
upper_limit_lasso = numeric(m)
lower_limit_lasso = numeric(m)
mean_ate_lasso = numeric(m)

for(i in 1:m)
{ 
  x = t(t(te_lasso[i,]))
  lower_limit_lasso[i] = mean(x) - qt(alpha / 2, df=n-1, lower.tail = FALSE)*sd(x)/sqrt(n)
  upper_limit_lasso[i] = mean(x) + qt(alpha / 2, df=n-1, lower.tail = FALSE)*sd(x)/sqrt(n)
  mean_ate_lasso[i] = mean(x)
}

inCI_lasso = ifelse(lower_limit_lasso < 0.5 & upper_limit_lasso > 0.5, 1, 0)
plot_data_CI_lasso = cbind.data.frame(c(1:500), mean_ate_lasso, lower_limit_lasso, upper_limit_lasso, inCI_lasso)
colnames(plot_data_CI_lasso) = c("Index", "Mean", "Lower", "Upper", "inCI")
coverage_prob_lasso = mean(lower_limit_lasso < 0.5 & upper_limit_lasso > 0.5)

plot_coverage_lasso = ggplot(plot_data_CI_lasso, aes(x = Index, y = Mean)) +
  theme_tufte(base_family = "serif") +
  geom_errorbar(aes(ymax = Upper, ymin = Lower, color = as.factor(inCI)), show.legend = FALSE) +
  scale_color_manual(values = c("blue", "grey")) +
  geom_point(size = 1, aes(color = as.factor(inCI)), show.legend = FALSE) +
  scale_fill_manual(values = c("darkblue", "darkgrey")) +
  geom_hline(yintercept = 0.5, color = "red") +
  labs(y = "Estimated treatment effect", x = "Simulation repetition", title = "Lasso")

## XGBoost
upper_limit_xgb = numeric(m)
lower_limit_xgb = numeric(m)
mean_ate_xgb = numeric(m)

for(i in 1:m)
{ 
  x = t(t(te_xgb[i,]))
  lower_limit_xgb[i] = mean(x) - qt(alpha / 2, df=n-1, lower.tail = FALSE)*sd(x)/sqrt(n)
  upper_limit_xgb[i] = mean(x) + qt(alpha / 2, df=n-1, lower.tail = FALSE)*sd(x)/sqrt(n)
  mean_ate_xgb[i] = mean(x)
}

inCI_xgb = ifelse(lower_limit_xgb < 0.5 & upper_limit_xgb > 0.5, 1, 0)
plot_data_CI_xgb = cbind.data.frame(c(1:500), mean_ate_xgb, lower_limit_xgb, upper_limit_xgb, inCI_xgb)
colnames(plot_data_CI_xgb) = c("Index", "Mean", "Lower", "Upper", "inCI")
coverage_prob_xgb = mean(lower_limit_xgb < 0.5 & upper_limit_xgb > 0.5)

plot_coverage_xgb = ggplot(plot_data_CI_xgb, aes(x = Index, y = Mean)) +
  theme_tufte(base_family = "serif") +
  geom_errorbar(aes(ymax = Upper, ymin = Lower, color = as.factor(inCI)), show.legend = FALSE) +
  scale_color_manual(values = c("blue", "grey")) +
  geom_point(size = 1, aes(color = as.factor(inCI)), show.legend = FALSE) +
  scale_fill_manual(values = c("darkblue", "darkgrey")) +
  geom_hline(yintercept = 0.5, color = "red") +
  labs(y = "Estimated treatment effect", x = "Simulation repetition", title = "XGBoost")

## Neural network
upper_limit_nn = numeric(m)
lower_limit_nn = numeric(m)
mean_ate_nn = numeric(m)

for(i in 1:m)
{ 
  x = t(t(te_nn[i,]))
  lower_limit_nn[i] = mean(x) - qt(alpha / 2, df=n-1, lower.tail = FALSE)*sd(x)/sqrt(n)
  upper_limit_nn[i] = mean(x) + qt(alpha / 2, df=n-1, lower.tail = FALSE)*sd(x)/sqrt(n)
  mean_ate_nn[i] = mean(x)
}

inCI_nn = ifelse(lower_limit_nn < 0.5 & upper_limit_nn > 0.5, 1, 0)
plot_data_CI_nn = cbind.data.frame(c(1:500), mean_ate_nn, lower_limit_nn, upper_limit_nn, inCI_nn)
colnames(plot_data_CI_nn) = c("Index", "Mean", "Lower", "Upper", "inCI")
coverage_prob_nn = mean(lower_limit_nn < 0.5 & upper_limit_nn > 0.5)

plot_coverage_nn = ggplot(plot_data_CI_nn, aes(x = Index, y = Mean)) +
  theme_tufte(base_family = "serif") +
  geom_errorbar(aes(ymax = Upper, ymin = Lower, color = as.factor(inCI)), show.legend = FALSE) +
  scale_color_manual(values = c("blue", "grey")) +
  geom_point(size = 1, aes(color = as.factor(inCI)), show.legend = FALSE) +
  scale_fill_manual(values = c("darkblue", "darkgrey")) +
  geom_hline(yintercept = 0.5, color = "red") +
  labs(y = "Estimated treatment effect", x = "Simulation repetition", title = "Neural Network")


pdf_file = file.path(figure_path, paste0("coverage_probs_1_50", ".pdf"))

# create pdf file with the given input
pdf(file = pdf_file,   # The directory you want to save the file in
    family = "serif",
    fonts = "serif",
    width = 13,
    height = 16)

grid.arrange(plot_coverage_ens, plot_coverage_lasso, plot_coverage_xgb, plot_coverage_nn, ncol = 1, nrow = 4, newpage = FALSE)

dev.off()
extrafont::embed_fonts(pdf_file, outfile=pdf_file)

rm(pdf_file)
print("Plot has been saved as pdf.")

coverage_probs_all = cbind.data.frame(coverage_prob_ens, coverage_prob_lasso, coverage_prob_xgb, coverage_prob_nn)

# jarque bera test for normality
jb_ate_ens = jb.norm.test(ate$ens, nrepl = 500)
jb_ate_lasso = jb.norm.test(ate$lasso, nrepl = 500)
jb_ate_xgb = jb.norm.test(ate$xgb, nrepl = 500)
jb_ate_nn = jb.norm.test(ate$nn, nrepl = 500)

jb_ate_values_all = cbind.data.frame(ml_names, rbind.data.frame(jb_ate_ens, jb_ate_lasso, jb_ate_xgb, jb_ate_nn))

# illustration of true_te to true_p
plot_data_te_x = cbind.data.frame(data.table(true_te[1,]), x_example$V1)
colnames(plot_data_te_x) = c("TE", "X1")
plot_data_te_x = plot_data_te_x[order(plot_data_te_x$X1),]

plot_te_x = ggplot(plot_data_te_x, aes(x = X1, y = TE))+
  geom_line() +
  stat_smooth()
