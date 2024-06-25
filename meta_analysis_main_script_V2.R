# Meta survival data analysis
# Author: Canan Karakoc
# Last update: May 29 2024

# Install/load libraries

# library(styler)
library(tidyverse)
library(minpack.lm)
library(survival)
library(stats)
library(loo)
library(ggforce)
library(ggpubr)
library(nlme)
library(MASS)
# packages for models
library(lme4)
library(sjPlot)

library(dplyr) # masked from another package
set.seed(1234)
setwd("~/GitHub/metaAnalysis")

library(patchwork)
library(brms) # If using bayesian statistics 

# Loading data
main.dir <- "~/GitHub/metaAnalysis/metaData"
allFiles <- list.files(path = main.dir, recursive = TRUE, full.names = TRUE)
csv_files <- allFiles[grep("\\.csv$", allFiles)]

# Meta data for later
meta <- read.csv("~/GitHub/metaAnalysis/000_combined_library.csv",
  sep = ",", header = T, na.strings = ""
)

# Read all CSV files into a list of data frames
# dfs <- map(csv_files, ~read_csv(.x, col_names = FALSE))
# dfs <- map(dfs, ~mutate_if(.x, is.numeric, as.numeric))

dataAll <- csv_files %>%
  setNames(basename(csv_files)) %>%
  map_df(read_csv, .id = "Full_key", col_names = c("Time", "Survival")) %>%
  mutate(Time = round(Time, digits = 1)) %>%
  mutate(Survival = round(Survival, digits = 1)) %>% # I digitized too precise that digits are long
  mutate(across(Time:Survival, ~ ifelse(.x < 0, 0, .x))) %>%
  mutate(across(Time:Survival, ~ ifelse(.x > 100, 100, .x)))

dataClean <- dataAll %>%
  mutate(key = gsub("\\..*", "", Full_key)) %>%
  mutate(Key = gsub("\\_.*", "", key)) %>%
  dplyr::select(-key) %>%
  left_join(meta, by = "Key")

# Find unique datasets to reverse using Reverse_Flag column
datasets_to_reverse <- dataClean %>%
  filter(Graph == "mortality") %>%
  pull(Full_key) %>%
  unique()

# Reverse the order of the Value column for selected datasets
# Remove the data sets with more than 20% survivors
dataCleanSurvival <- dataClean %>%
  group_by(Full_key) %>%
  mutate(Survival_corrected = if_else(Full_key %in% datasets_to_reverse, rev(Survival), Survival)) %>%
  filter(last(Survival_corrected) <= 20) %>%
  ungroup()

# Find datasets with time unit in days
datasets_in_days <- dataClean %>%
  filter(Time_units == "days") %>%
  pull(Full_key) %>%
  unique()

# Convert time unit from days to hours for selected datasets
# Standardize time between 1-0
dataCleanSurvival_hours <- dataCleanSurvival %>%
  mutate(Time_hours = if_else(Full_key %in% datasets_in_days, Time * 24, Time)) %>%
  group_by(Full_key) %>%
  mutate(Time_std = Time_hours / max(Time_hours)) %>%
  ungroup()

# Standardize survival values between 0 and 1 with final survival as 0
standardized_df <- dataCleanSurvival_hours %>%
  group_by(Full_key) %>%
  mutate(
    Standardized_Survival = (Survival_corrected - min(Survival_corrected)) / (max(Survival_corrected) - min(Survival_corrected)) # Scale between 0 and 1
  ) %>%
  ungroup()

# Final data to use
standardized_final_df <- standardized_df %>%
  dplyr::select(
    Full_key, Key, Time, Survival, Time_hours, Time_std, Survival_corrected, Standardized_Survival,
    Host, Pathogen, Host_taxa, Host_other, Pathogen_taxa, Pathogen_other, Max_num_host, Data
  )

write.csv(standardized_final_df, "standardized_final_df.csv", row.names = FALSE)

######################## MODELS #############################

# Create folders to store the plots
dir.create("exponential_model_plots", showWarnings = FALSE)
dir.create("weibull_model_plots", showWarnings = FALSE)
dir.create("gompertz_hazard_model_plots", showWarnings = FALSE)
dir.create("gompertz_survival_model_plots", showWarnings = FALSE)
dir.create("loglogistic_model_plots", showWarnings = FALSE)
dir.create("generalizedgamma_model_plots", showWarnings = FALSE)

# Define model fitting functions
fit_exponential <- function(data) {
  fit <- try(nlsLM(Standardized_Survival ~ exp(-lambda * Time_std), data = data, start = list(lambda = 0.1), control = nls.lm.control(maxiter = 500)), silent = TRUE)
  if (inherits(fit, "try-error") || any(is.na(coef(fit)))) {
    return(NULL)
  } else {
    return(fit)
  }
}

fit_weibull <- function(data) {
  fit <- try(nlsLM(Standardized_Survival ~ exp(-(Time_std / lambda)^k), data = data, start = list(lambda = 1, k = 0.1), control = nls.lm.control(maxiter = 500)), silent = TRUE)
  if (inherits(fit, "try-error") || any(is.na(coef(fit)))) {
    return(NULL)
  } else {
    return(fit)
  }
}

fit_gompertz_hazard <- function(data) {
  fit <- try(nlsLM(Standardized_Survival ~ exp(-a * exp(b * Time_std)), data = data, start = list(a = 0.1, b = 0.1), control = nls.lm.control(maxiter = 500)), silent = TRUE)
  if (inherits(fit, "try-error") || any(is.na(coef(fit)))) {
    return(NULL)
  } else {
    return(fit)
  }
}

fit_gompertz_survival <- function(data) {
  fit <- try(nlsLM(Standardized_Survival ~ exp(-(a / b) * (exp(b * Time_std) - 1)), data = data, start = list(a = 0.1, b = 0.1), control = nls.lm.control(maxiter = 500)), silent = TRUE)
  if (inherits(fit, "try-error") || any(is.na(coef(fit)))) {
    return(NULL)
  } else {
    return(fit)
  }
}

fit_loglogistic <- function(data) {
  fit <- try(nlsLM(Standardized_Survival ~ 1 / (1 + (Time_std / alpha)^beta), data = data, start = list(alpha = 1, beta = 1), control = nls.lm.control(maxiter = 500)), silent = TRUE)
  if (inherits(fit, "try-error") || any(is.na(coef(fit)))) {
    return(NULL)
  } else {
    return(fit)
  }
}

fit_generalizedgamma <- function(data) {
  fit <- try(nlsLM(Standardized_Survival ~ exp(-(Time_std / beta)^gamma * (1 + (alpha - 1) * (Time_std / beta)^gamma)), data = data, start = list(beta = 1, gamma = 1, alpha = 1), control = nls.lm.control(maxiter = 500)), silent = TRUE)
  if (inherits(fit, "try-error") || any(is.na(coef(fit)))) {
    return(NULL)
  } else {
    return(fit)
  }
}

# Function to calculate AICc
calculate_aicc <- function(model, data) {
  n <- nrow(data)
  k <- length(coef(model)) + 1 # number of parameters including the variance
  if (n <= k + 1) {
    return(Inf) # Return infinity if sample size is too small
  }
  aic <- AIC(model)
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  return(aicc)
}

# Function to calculate WAIC and its components
calculate_waic <- function(model, data) {
  log_lik <- logLik(model)
  if (is.null(log_lik)) {
    return(list(mean_waic = NA))
  }
  
  log_lik_values <- as.vector(log_lik)
  n <- length(log_lik_values)
  
  # Compute the log pointwise predictive density (lppd)
  lppd <- sum(log_lik_values)
  
  # Compute the effective number of parameters (pWAIC)
  pwaic <- sum((log_lik_values - mean(log_lik_values))^2)
  
  # Compute WAIC
  waic <- -2 * (lppd - pwaic)
  
  return(list(mean_waic = waic))
}

# Function to create and save model fit plot
save_model_fit_plot <- function(data, predictions, dataset_name, model_name, folder_name) {
  prediction_df <- data.frame(Time_std = seq(0, 1, length.out = 100), Predictions = predictions)
  
  plot <- ggplot() +
    geom_point(data = data, aes(x = Time_std, y = Standardized_Survival), color = "black", size = 3, shape = 1) +
    geom_line(data = prediction_df, aes(x = Time_std, y = Predictions), color = "firebrick") +
    labs(title = paste(model_name, "Model Fit for Dataset:", dataset_name), x = "Time", y = "Standardized_Survival") +
    theme_bw()
  
  dir.create(folder_name, showWarnings = FALSE)
  plot_filename <- paste0(folder_name, "/model_fit_", dataset_name, ".png")
  ggsave(filename = plot_filename, plot = plot, width = 8, height = 6)
}

# Initialize lists to store results
all_results <- list()
all_predictions <- data.frame()

# Iterate over each dataset
unique_datasets <- unique(standardized_final_df$Full_key)

for (dataset_name in unique_datasets) {
  dataset <- standardized_final_df %>%
    filter(Full_key == dataset_name) %>%
    dplyr::select(Time_std, Standardized_Survival)
  
  time_seq <- seq(0, 1, length.out = 100)
  
  # Fit models and generate predictions
  model_fits <- list(
    exponential = fit_exponential(dataset),
    weibull = fit_weibull(dataset),
    gompertz_hazard = fit_gompertz_hazard(dataset),
    gompertz_survival = fit_gompertz_survival(dataset),
    loglogistic = fit_loglogistic(dataset),
    generalizedgamma = fit_generalizedgamma(dataset)
  )
  
  for (model_name in names(model_fits)) {
    fit <- model_fits[[model_name]]
    if (!is.null(fit)) {
      if (model_name == "gompertz_hazard") {
        # Generate predictions for Gompertz hazard model
        a <- coef(fit)[1]
        b <- coef(fit)[2]
        preds <- exp(-a * exp(b * time_seq))
      } else if (model_name == "gompertz_survival") {
        # Generate predictions for Gompertz survival model
        a <- coef(fit)[1]
        b <- coef(fit)[2]
        preds <- exp(-(a / b) * (exp(b * time_seq) - 1))
      } else {
        preds <- predict(fit, newdata = data.frame(Time_std = time_seq))
      }
      waic_values <- calculate_waic(fit, dataset)
      aicc_value <- calculate_aicc(fit, dataset)
      
      # Save results
      result <- data.frame(
        Dataset = dataset_name,
        Model = model_name,
        Mean_WAIC = waic_values$mean_waic,
        AICc = aicc_value,
        Num_Obs = nrow(dataset)
      )
      all_results <- rbind(all_results, result)
      
      # Save predictions
      predictions <- data.frame(Full_key = dataset_name, Model = model_name, Time_std = time_seq, Predictions = preds)
      all_predictions <- rbind(all_predictions, predictions)
      
      # Save plot
      save_model_fit_plot(dataset, preds, dataset_name, model_name, paste0(model_name, "_model_plots"))
    }
  }
}

# Save data frames to files
write.csv(all_results, "results_df.csv", row.names = FALSE)
write.csv(all_predictions, "all_predictions.csv", row.names = FALSE)

#####################################################################################
# READ SAVED DATA
standardized_final_df <- read.table("standardized_final_df.csv", header = T, sep = ",", dec = ".")
results_df <- read.table("results_df.csv", header = T, sep = ",", dec = ".")
all_predictions <- read.table("all_predictions.csv", header = T, sep = ",", dec = ".")

# Data
colnames(all_results)[1] <- "Full_key"

df2_first_match <- standardized_final_df %>%
  group_by(Full_key) %>%
  slice(1) %>%
  ungroup()

allData_WAIC <- all_results[-1, ] %>%
  left_join(df2_first_match, by = "Full_key") %>%
  dplyr::select(Full_key, Key, Host_taxa, Host_other, Pathogen_taxa, Pathogen_other, Num_Obs, Max_num_host, Model, Mean_WAIC, AICc, Data)

# Number of host
# Some papers (marked in the meta data) have different host numbers for different trials
# I am changing them manually here.

# ASXFW6GM_1 - ASXFW6GM_4
# c(24, 18, 8, 9)
# For (SLELFBN3_1 - SLELFBN3_30)
# c(24, 28, 23, 21, 23, 27, 25, 24, 29, 29, 27, 22, 26, 27, 25, 20, 21, 26, 23,  9,  9, 10, 25, 23, 20, 15, 14, 14, 11, 12)
# FRAQNTF6_1 - FRAQNTF6_4
# c(20, 20, 19, 19)

# NSH3NLMI_1 10
# NSH3NLMI_2 12
# 72SE5ARY_1 10
# 72SE5ARY_2 20

# Sample datasets and num_hosts vectors
datasets <- c(
  "ASXFW6GM_1", "ASXFW6GM_2", "ASXFW6GM_3", "ASXFW6GM_4", "SLELFBN3_1", "SLELFBN3_2", "SLELFBN3_3", "SLELFBN3_4", "SLELFBN3_5",
  "SLELFBN3_6", "SLELFBN3_7", "SLELFBN3_8", "SLELFBN3_9", "SLELFBN3_10", "SLELFBN3_11", "SLELFBN3_12", "SLELFBN3_13", "SLELFBN3_14", "SLELFBN3_15",
  "SLELFBN3_16", "SLELFBN3_17", "SLELFBN3_18", "SLELFBN3_19", "SLELFBN3_20", "SLELFBN3_21", "SLELFBN3_22", "SLELFBN3_23", "SLELFBN3_24", "SLELFBN3_25",
  "SLELFBN3_26", "SLELFBN3_27", "SLELFBN3_28", "SLELFBN3_29", "SLELFBN3_30", "FRAQNTF6_1", "FRAQNTF6_2", "FRAQNTF6_3", "FRAQNTF6_4",
  "NSH3NLMI_1", "NSH3NLMI_2", "72SE5ARY_1", "72SE5ARY_2"
)

num_hosts <- c(
  24, 18, 8, 9, 24, 28, 23, 21, 23, 27, 25, 24, 29, 29, 27, 22, 26, 27, 25, 20, 21, 26, 23,  9,  9, 10, 
  25, 23, 20, 15, 14, 14, 11, 12,  20, 20, 19, 19,  10, 12, 10, 20
)

# Create mergedata data frame
mergedata <- data.frame(key = datasets, Max_num_hosts = num_hosts)

# Join and fill NAs
allData_WAIC_filled <- allData_WAIC %>%
  mutate(key = gsub("\\..*", "", Full_key)) %>%
  left_join(mergedata, by = "key") %>%
  mutate(Max_num_host_filled = coalesce(Max_num_host, Max_num_hosts))

# Function to calculate standard error
standard_error <- function(x) {
  n <- sum(!is.na(x)) # Number of non-NA observations
  if (n == 0) {
    return(NA)
  } else {
    return(sd(x, na.rm = TRUE) / sqrt(n))
  }
}

###################### FIGURES ###########################
# panel.grid.major = element_blank(), panel.grid.minor = element_blank(),

mytheme <- theme_bw() +
  theme(axis.ticks.length = unit(.25, "cm")) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14, color = "black"), axis.title = element_text(size = 16)) +
  theme(panel.border = element_rect(
    fill = NA, colour = "black",
    size = 1
  )) +
  theme(strip.text.x = element_text(size = 14), strip.background = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(panel.border = element_rect(
    fill = NA, colour = "black",
    linewidth = 1
  )) +
  theme(
    axis.text.x.top = element_blank(), axis.title.x.top = element_blank(),
    axis.text.y.right = element_blank(), axis.title.y.right = element_blank()
  ) +
  theme(
    axis.title.x = element_text(margin = margin(10, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
    axis.text.x = element_text(margin = margin(10, 0, 0, 0)),
    axis.text.y = element_text(margin = margin(0, 10, 0, 0))
  )

# Color blind palette
cbpalette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#999999", "#F0E442", "#000000")


# HOST
allData_WAIC_sum <- allData_WAIC_filled %>%
  group_by(Model, Host_other) %>%
  summarize(
    meanWAIC = mean(as.numeric(Mean_WAIC), na.rm = T),
    seWAIC = standard_error(as.numeric(Mean_WAIC)),
    meanObs = round(mean(as.numeric(Num_Obs), na.rm = T), 0),
    meanHost = round(mean(as.numeric(Max_num_host_filled), na.rm = T), 0),
    numData = length(unique(Full_key))
  ) %>%
  ungroup() %>%
  filter(!Model == "gompertz_hazard") 

# Order factor levels
allData_WAIC_sum$Host_other <- factor(allData_WAIC_sum$Host_other, 
                                      levels = c("Seedlings", "Drosophila sp.", "Other insects", 
                                                 "Nematodes", "Moth larvae", "Other invertebrates", 
                                                 "Fish", "Avian", "Mice", "Other mammals"))


# Create the plot and add the annotations
host_waic <- ggplot(allData_WAIC_sum, aes(x = Model, y = meanWAIC, color = Model)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = meanWAIC - seWAIC, ymax = meanWAIC + seWAIC), width = 0.2) +
  facet_wrap(~Host_other, ncol = 2, scales = "free") +
  mytheme +
  theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = cbpalette, , labels = c("exponential", "generalized-gamma", "gompertz", "log-logistic", "weibull")) +
  ylab("WAIC") +
  geom_text(aes(label = paste0("num. datasets = ", numData)), x = Inf, y = Inf, hjust = 1.1, vjust = 1.2, size = 4, check_overlap = TRUE, show.legend = FALSE) +
  geom_text(aes(label = paste0("num. time obs. = ", meanObs)), x = Inf, y = Inf, hjust = 1.1, vjust = 2.3, size = 4, check_overlap = TRUE, show.legend = FALSE) +
  geom_text(aes(label = paste0("avg. host reps = ", meanHost)), x = Inf, y = Inf, hjust = 1.1, vjust = 3.4, size = 4, check_overlap = TRUE, show.legend = FALSE) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave("figures/host_waic.png", plot = host_waic, width = 6, height = 12, units = "in", dpi = 300)
ggsave("figures/host_waic.pdf", plot = host_waic, width = 6, height = 12, units = "in", dpi = 300)

# PATGOGEN
# Order factor levels
PallData_WAIC_sum <- allData_WAIC_filled %>%
  group_by(Model, Pathogen_other) %>%
  summarize(
    meanWAIC = mean(as.numeric(Mean_WAIC)),
    seWAIC = standard_error(as.numeric(Mean_WAIC)),
    meanObs = round(mean(as.numeric(Num_Obs)), 0),
    meanHost = round(mean(as.numeric(Max_num_host_filled), na.rm = T), 0),
    numData = length(unique(Full_key))
  ) %>%
  filter(!Model == "gompertz_hazard") 

PallData_WAIC_sum$Pathogen_other <- factor(PallData_WAIC_sum$Pathogen_other, 
                                          levels = c("Gram-positive bacteria", "Gram-negative bacteria",
                                                     "DNA virus", "RNA virus", "Fungi", "Protozoan parasite"))

# Create the plot and add the annotations
pathogen_waic <- ggplot(PallData_WAIC_sum, aes(x = Model, y = meanWAIC, color = Model)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = meanWAIC - seWAIC, ymax = meanWAIC + seWAIC), width = 0.2) +
  facet_wrap(~Pathogen_other, ncol = 2, scales = "free") +
  mytheme +
  theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = cbpalette, labels = c("exponential", "generalized-gamma", "gompertz", "log-logistic", "weibull")) +
  ylab("WAIC") +
  geom_text(aes(label = paste0("num. datasets = ", numData)), x = Inf, y = Inf, hjust = 1.1, vjust = 1.2, size = 4, check_overlap = TRUE, show.legend = FALSE) +
  geom_text(aes(label = paste0("num. time obs. = ", meanObs)), x = Inf, y = Inf, hjust = 1.1, vjust = 2.3, size = 4, check_overlap = TRUE, show.legend = FALSE) +
  geom_text(aes(label = paste0("avg. host reps = ", meanHost)), x = Inf, y = Inf, hjust = 1.1, vjust = 3.4, size = 4, check_overlap = TRUE, show.legend = FALSE) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave("figures/pathogen_waic.png", plot = pathogen_waic, width = 6, height = 8, units = "in", , dpi = 300)
ggsave("figures/pathogen_waic.pdf", plot = pathogen_waic, width = 6, height = 8, units = "in", , dpi = 300)

# Data_type
Data_WAIC <- allData_WAIC_filled %>%
  group_by(Model, Data) %>%
  summarize(meanWAIC = mean(Mean_WAIC, na.rm = TRUE), seWAIC = standard_error(Mean_WAIC)) %>%
  filter(!Model == "gompertz_hazard")

Datatype <- ggplot(Data_WAIC, aes(x = Data, y = meanWAIC, color = Model)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = meanWAIC - seWAIC, ymax = meanWAIC + seWAIC),
    width = 0.2,
    position = position_dodge(width = 0.5, preserve = "single")
  ) +
  mytheme +
  scale_color_manual(values = cbpalette) +
  theme(axis.title.x = element_blank()) +
  ylab("WAIC")

ggsave("figures/datatype.pdf", plot = Datatype, width = 6, height = 3, units = "in", dpi = 300)

Data_num <- allData_WAIC_filled %>%
  filter(Model == "gompertz_hazard")%>%
  group_by(Data) %>%
  summarize(Data = n())

# Num. of Observations
allData_WAIC_sum_obs <- allData_WAIC_filled %>%
  group_by(Full_key, Model) %>%
  summarize(
    meanObs = mean(as.numeric(Num_Obs), na.rm = T),
    meanWAIC = mean(Mean_WAIC, na.rm = T)
  ) %>%
  ungroup()

allData_WAIC_sum_num <- allData_WAIC_filled %>%
  group_by(Full_key, Model) %>%
  summarize(
    meanNum = mean(as.numeric(Max_num_host_filled), na.rm = T),
    meanWAIC = mean(Mean_WAIC, na.rm = T)
  ) %>%
  ungroup()

# Define the first plot
Obs <- ggplot(allData_WAIC_sum_obs, aes(x = meanObs, y = meanWAIC, color = Model)) +
  geom_point(size = 3, position = position_dodge(width = 0.5), alpha = 0.5) +
  mytheme +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = cbpalette) +
  ylab("WAIC") +
  xlab("Number of observations") +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left",
    label.y.npc = "bottom",
    size = 4
  ) +
  theme(strip.background = element_blank()) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# Define the second plot
Num_host_filtered <- ggplot(allData_WAIC_sum_num %>% filter(meanNum != 500), aes(x = as.numeric(meanNum), y = meanWAIC, color = Model)) +
  geom_point(size = 3, position = position_dodge(width = 0.5), alpha = 0.5) +
  mytheme +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = cbpalette) +
  ylab("WAIC") +
  xlab("Number of hosts") +
  theme(strip.background = element_blank()) +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left",
    label.y.npc = "bottom",
    size = 4
  ) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# Combine the plots using patchwork and remove the duplicate legend
combined_plot <- (Obs + Num_host_filtered) + plot_layout(guides = 'collect') &
  theme(legend.position = "bottom")

ggsave("figures/Num_host_Obs.pdf", plot = combined_plot, width = 9.5, height = 5.5, units = "in", dpi = 300)

# Combine the plots using patchwork and ensure one common legend
combined_plot <- (Obs + Num_host_filtered) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

# Supplementary figure
# Create combined_factor and plotting_data
standardized_final_df$combined_factor <- interaction(standardized_final_df$Host_other, standardized_final_df$Pathogen_taxa, sep = "\n")
plotting_data <- standardized_final_df %>%
  dplyr::select(Full_key, combined_factor, Standardized_Survival, Time_std) %>%
  mutate(key = gsub("\\..*", "", Full_key))

# Create factors
factors <- standardized_final_df %>%
  dplyr::select(Full_key, combined_factor)

# Merge and arrange prediction data
prediction_plotting_data <- all_predictions %>%
  mutate(key = gsub("\\..*", "", Full_key)) %>%
  left_join(factors, by = "Full_key") %>%
  filter(!Model == "gompertz_hazard") %>%
  arrange(combined_factor)
  
# Example figure
fig1_dat1 <- plotting_data %>%
  filter(Full_key == "7CJPGPP9_1.csv")
fig1_dat2 <- prediction_plotting_data %>%
  filter(Full_key == "7CJPGPP9_1.csv") %>%
  filter(!Model == "gompertz_hazard")

fig1 <- ggplot(fig1_dat1, aes(x = Time_std, y = Standardized_Survival)) +
  geom_point(shape = 1, size = 3, color = "grey25") + # Raw data
  geom_line(data = fig1_dat2, aes(y = Predictions, color = Model), size = 1, alpha = 0.5) + # Predictions
  labs(x = "Standardized time", y = "Standardized survival") +
  mytheme +
  scale_color_manual(values = cbpalette)

ggsave("figures/fig1.pdf", plot = fig1, width = 7, height = 4, units = "in", dpi = 300)

# Another figure 
fig1_dat1_2 <- plotting_data %>%
  filter(Full_key == "7BHF4LIC_1.csv")
fig1_dat2_2 <- prediction_plotting_data %>%
  filter(Full_key == "7BHF4LIC_1.csv") %>%
  filter(!Model == "gompertz_hazard")

aic <- allData_WAIC_filled %>% filter(Full_key == "7BHF4LIC_1.csv")

fig1_2 <- ggplot(fig1_dat1_2, aes(x = Time_std, y = Standardized_Survival)) +
  geom_point(shape = 1, size = 3, color = "grey25") + # Raw data
  geom_line(data = fig1_dat2_2, aes(y = Predictions, color = Model), size = 1, alpha = 0.5) + # Predictions
  labs(x = "Standardized time", y = "Standardized survival") +
  mytheme +
  scale_color_manual(values = cbpalette)

ggsave("figures/fig1_2.pdf", plot = fig1_2, width = 7, height = 4, units = "in", dpi = 300)
  
# Plotting
ALL <- ggplot(plotting_data, aes(x = Time_std, y = Standardized_Survival)) +
  geom_point(shape = 1, size = 3, color = "grey25") +  # Raw data
  geom_line(data = prediction_plotting_data, aes(y = Predictions, x = Time_std, color = Model), size = 0.8, alpha = 0.5) +  # Predictions
  labs(x = "Standardized time", y = "Standardized survival") +
  facet_wrap(~ combined_factor + key) +
  scale_color_manual(values = cbpalette)+
  mytheme+
  theme(legend.position = "bottom")

# Pagination parameters
ncol <- 5
nrow <- 7
n_pages <- ceiling(length(unique(standardized_final_df$Full_key)) / (ncol * nrow))

# Loop to create paginated plots
for (i in 1:n_pages) {
  paginated_plot <- ALL + facet_wrap_paginate(~ combined_factor + key, ncol = ncol, nrow = nrow, page = i)
  
  # Save each page to a file
  ggsave(paste0("combined_plot_page_", i, ".png"), paginated_plot, width = 15, height = 17)
  
  # Print each paginated plot
  print(paginated_plot)
}

######################################################################
# Evidence ratios

# Define a small constant
epsilon <- 1e-10

# Define winsorize function
winsorize <- function(x, probs = c(0.01, 0.99)) {
  quants <- quantile(x, probs = probs, na.rm = TRUE)
  x[x < quants[1]] <- quants[1]
  x[x > quants[2]] <- quants[2]
  return(x)
}

# winsorizing the AICc values before calculating the Evidence Ratios can help 
# stabilize the calculations by reducing the impact of extreme values early in the process. 
# This can prevent extremely high or low Evidence Ratios caused by outliers in the AICc values.
# We'll continue using WAIC. Winsorizing is not needed for WAIC, because it performs pretty good. 

# Winsorize AICc values
allData_WAIC_filled_win <- allData_WAIC_filled %>%
  group_by(Model) %>%
  mutate(Mean_WAIC = winsorize(Mean_WAIC)) %>%
  ungroup()

# Reshape the data to have AIC/WAIC values of both models in one row for each Full_key
allData_AIC_likelihood <- allData_WAIC_filled_win %>%
  filter(!Model %in% c("exponential", "gompertz_hazard")) %>%
  dplyr::select(Full_key, Model, Mean_WAIC) %>%
  ungroup() %>%
  pivot_wider(names_from = Model, values_from = Mean_WAIC) %>%
  rowwise() %>%
  mutate(
    AIC_min_weibull = min(gompertz_survival + epsilon, weibull + epsilon, na.rm = TRUE),
    AIC_min_loglogistic = min(gompertz_survival + epsilon, loglogistic + epsilon, na.rm = TRUE),
    AIC_min_generalizedgamma = min(gompertz_survival + epsilon, generalizedgamma + epsilon, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    Delta_gompertz_weibull = (gompertz_survival + epsilon) - AIC_min_weibull,
    Delta_weibull_weibull = (weibull + epsilon) - AIC_min_weibull,
    Delta_gompertz_loglogistic = (gompertz_survival + epsilon) - AIC_min_loglogistic,
    Delta_loglogistic_loglogistic = (loglogistic + epsilon) - AIC_min_loglogistic,
    Delta_gompertz_generalizedgamma = (gompertz_survival + epsilon) - AIC_min_generalizedgamma,
    Delta_generalizedgamma_generalizedgamma = (generalizedgamma + epsilon) - AIC_min_generalizedgamma
  ) %>%
  mutate(
    Relative_Likelihood_gompertz_weibull = exp(-0.5 * Delta_gompertz_weibull),
    Relative_Likelihood_weibull_weibull = exp(-0.5 * Delta_weibull_weibull),
    Relative_Likelihood_gompertz_loglogistic = exp(-0.5 * Delta_gompertz_loglogistic),
    Relative_Likelihood_loglogistic_loglogistic = exp(-0.5 * Delta_loglogistic_loglogistic),
    Relative_Likelihood_gompertz_generalizedgamma = exp(-0.5 * Delta_gompertz_generalizedgamma),
    Relative_Likelihood_generalizedgamma_generalizedgamma = exp(-0.5 * Delta_generalizedgamma_generalizedgamma)
  ) %>%
  mutate(
    Evidence_Ratio_gompertz_weibull = (Relative_Likelihood_gompertz_weibull + epsilon) / (Relative_Likelihood_weibull_weibull + epsilon),
    Evidence_Ratio_gompertz_loglogistic = (Relative_Likelihood_gompertz_loglogistic + epsilon) / (Relative_Likelihood_loglogistic_loglogistic + epsilon),
    Evidence_Ratio_gompertz_generalizedgamma = (Relative_Likelihood_gompertz_generalizedgamma + epsilon) / (Relative_Likelihood_generalizedgamma_generalizedgamma + epsilon)
  ) %>%
  pivot_longer(cols = starts_with("Evidence_Ratio"), names_to = "Comparison", values_to = "Evidence_Ratio")

# Inspect the resulting data for extreme values
summary(allData_AIC_likelihood)

# Plot Evidence Ratios from AICc
ER_AIC <- ggplot(allData_AIC_likelihood, aes(x = Comparison, y = log10(Evidence_Ratio+epsilon))) +
  geom_boxplot() +
  geom_jitter(shape = 21, alpha = 0.5) +
  scale_x_discrete(labels = c("Gen.-gamma", "Log-logistic", "Weibull")) +
  labs(y = expression(italic(log[10])~"(Evidence ratio"[WAIC]*")"),
       x = "Comparison with Gompertz") +
  mytheme
ER_AIC

ggsave("figures/evidence_ratio.png", plot = ER_AIC, width = 6, height = 4, units = "in", dpi = 300)
ggsave("figures/evidence_ratio.pdf", plot = ER_AIC, width = 6, height = 4, units = "in", dpi = 300)

# Test for significance 
# filtered_data <- allData_AIC_likelihood %>%
# filter_all(all_vars(!is.infinite(.))) %>%
# drop_na()

# We are reporting lm model in the manuscript
model1_evidence1 <- lm(Evidence_Ratio ~ Comparison, data = allData_AIC_likelihood)
summary(model1_evidence1)

# Fit the first Bayesian model
# Set contrasts for Comparison to contr.sum
#allData_AIC_likelihood$log_Evidence_Ratio <- log(allData_AIC_likelihood$Evidence_Ratio)
#allData_AIC_likelihood$Comparison <- as.factor(allData_AIC_likelihood$Comparison)
#contrasts(allData_AIC_likelihood$Comparison) <- contr.sum(length(levels(allData_AIC_likelihood$Comparison)))
#model <- brm(log(Evidence_Ratio) ~ Comparison, data = allData_AIC_likelihood)
#summary(model)

# Delta WAIC
allData_AICc_delta <- allData_WAIC_filled %>%
  filter(Model %in% c("gompertz_survival", "exponential")) %>%
  spread(Model, Mean_WAIC) %>%
  group_by(Full_key, Key, Host_taxa, Pathogen_taxa, Host_other, Pathogen_other) %>%
  summarise(exponential = mean(exponential, na.rm = TRUE), gompertz_survival = mean(gompertz_survival, na.rm =T), 
            Num_Obs = mean(Num_Obs), Num_Host = mean(Max_num_host_filled)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(DeltaSub = exponential - gompertz_survival)

allData_AICc_delta_sum <- allData_AICc_delta %>%
  group_by(Host_taxa) %>%
  summarize(
    delta_WAIC = mean(DeltaSub, na.rm = T),
    se = standard_error(DeltaSub)) 

allData_AICc_delta_sum2 <- allData_AICc_delta %>%
  group_by(Host_other) %>%
  summarize(
    delta_WAIC = mean(DeltaSub, na.rm = T),
    se = standard_error(DeltaSub))  

allData_AICc_delta_sum_P <- allData_AICc_delta %>%
  group_by(Pathogen_taxa) %>%
  summarize(
    delta_WAIC = mean(DeltaSub, na.rm = T),
    se = standard_error(DeltaSub))  

allData_AICc_delta_sum_P2 <- allData_AICc_delta %>%
  group_by(Pathogen_other) %>%
  summarize(
    delta_WAIC = mean(DeltaSub, na.rm = T),
    se = standard_error(DeltaSub))  

delta_host <- ggplot(allData_AICc_delta_sum, aes(x = Host_taxa, y = delta_WAIC)) +
  geom_point(size=5, shape=21)+
  geom_errorbar(aes(ymin = delta_WAIC - se, ymax = delta_WAIC + se), width = 0.2) +
  mytheme+
  labs(y = expression(Delta[WAIC]), x = NULL, title = "Host taxa") +
  theme(plot.title = element_text(hjust = 0.5))

allData_AICc_delta_sum2$Host_other <- factor(allData_AICc_delta_sum2$Host_other, 
                                             levels = rev(c("Seedlings", "Drosophila sp.", "Other insects", 
                                                        "Nematodes", "Moth larvae", "Other invertebrates", 
                                                        "Fish", "Avian", "Mice", "Other mammals")))
delta_host2 <- ggplot(allData_AICc_delta_sum2, aes(y = Host_other, x = delta_WAIC)) +
  geom_point(size=5, shape=21)+
  geom_errorbar(aes(xmin = delta_WAIC - se, xmax = delta_WAIC + se), width = 0.2) +
  mytheme+
  labs(x = expression(Delta[WAIC]), y = NULL, title = "Host taxa") +
  theme(plot.title = element_text(hjust = 0.5))

delta_path <- ggplot(allData_AICc_delta_sum_P, aes(x = Pathogen_taxa, y = delta_WAIC)) +
  geom_point(size = 5, shape = 21)+
  geom_errorbar(aes(ymin = delta_WAIC - se, ymax = delta_WAIC + se), width = 0.2) +
  mytheme+
  ylab(expression(Delta[WAIC])) +
  scale_color_manual(values = cbpalette)+
  labs(y = NULL, x = NULL, title = "Pathogen taxa") +
  theme(plot.title = element_text(hjust = 0.5))

allData_AICc_delta_sum_P2$Pathogen_other <- factor(allData_AICc_delta_sum_P2$Pathogen_other, 
                                                              levels = rev(c("Gram-positive bacteria", "Gram-negative bacteria",
                                                                         "DNA virus", "RNA virus", "Fungi", "Protozoan parasite")))

delta_path2 <- ggplot(allData_AICc_delta_sum_P2, aes(y = Pathogen_other, x = delta_WAIC)) +
  geom_point(size = 5, shape = 21)+
  geom_errorbar(aes(xmin = delta_WAIC - se, xmax = delta_WAIC + se), width = 0.2) +
  mytheme+
  xlab(expression(Delta[WAIC])) +
  scale_color_manual(values = cbpalette)+
  labs(y = NULL, title = "Pathogen taxa") +
  theme(plot.title = element_text(hjust = 0.5))

comb_delta <- ggarrange(delta_host, delta_path)
comb_delta2 <- ggarrange(delta_host2, delta_path2)

ggsave("figures/hostpathWAIC.pdf", plot = comb_delta, width = 8, height = 4, units = "in", dpi = 300)
ggsave("figures/hostpath2WAIC.pdf", plot = comb_delta2, width = 10, height = 4, units = "in", dpi = 300)
ggsave("figures/hostpath2WAIC.png", plot = comb_delta2, width = 10, height = 4, units = "in", dpi = 300)

# Observations and number of host 
delta_sum_obs <- allData_AICc_delta %>%
  group_by(Full_key) %>%
  summarize(
    meanObs = mean(as.numeric(Num_Obs), na.rm = T),
    meanEvidence = mean(DeltaSub, na.rm = T)
  ) %>%
  ungroup()

delta_sum_num <- allData_AICc_delta %>%
  group_by(Full_key) %>%
  summarize(
    meanNum = mean(as.numeric(Num_Host), na.rm = T),
    meanEvidence = mean(DeltaSub, na.rm = T)
  ) %>%
  ungroup()

Obs_delta <- ggplot(delta_sum_obs, aes(x = meanObs, y = meanEvidence)) +
  geom_point(size = 3, position = position_dodge(width = 0.5), alpha = 0.5) +
  mytheme +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = cbpalette) +
  ylab(expression(Delta[WAIC])) +
  xlab("Number of observations") +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  theme(strip.background = element_blank()) 

Num_host_filtered_delta <- ggplot(delta_sum_num %>% filter(meanNum != 500), aes(x = as.numeric(meanNum), y = meanEvidence)) +
  geom_point(size = 3, position = position_dodge(width = 0.5), alpha = 0.5) +
  mytheme +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = cbpalette) +
  ylab(expression(Delta[WAIC])) +
  xlab("Number of hosts") +
  theme(strip.background = element_blank()) +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))+
  labs(y = NULL)

comb_numb <- ggarrange(Obs_delta, Num_host_filtered_delta)
ggsave("figures/obs_numWAIC.pdf", plot = comb_numb, width = 8, height = 4, units = "in", dpi = 300)
ggsave("figures/obs_numWAIC.png", plot = comb_numb, width = 8, height = 4, units = "in", dpi = 300)


# Zoom where the most data are 
Obs_delta_zoom <- ggplot(delta_sum_obs %>% filter(meanObs < 21), aes(x = meanObs, y = meanEvidence)) +
  geom_point(size = 3, position = position_dodge(width = 0.5), alpha = 0.5) +
  mytheme +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = cbpalette) +
  ylab(expression(Delta[WAIC])) +
  xlab("Number of observations") +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  theme(strip.background = element_blank()) 

ggsave("figures/obs_WAIC_zoom.pdf", plot = Obs_delta_zoom, width = 4, height = 4, units = "in", dpi = 300)
ggsave("figures/obs_WAIC_zoom.png", plot = Obs_delta_zoom, width = 4, height = 4, units = "in", dpi = 300)


Num_host_filtered_delta_zoom <- ggplot(delta_sum_num %>% filter(meanNum < 21), aes(x = as.numeric(meanNum), y = meanEvidence)) +
  geom_point(size = 3, position = position_dodge(width = 0.5), alpha = 0.5) +
  mytheme +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = cbpalette) +
  ylab(expression(Delta[WAIC])) +
  xlab("Number of hosts") +
  theme(strip.background = element_blank()) +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))+
  labs(y = NULL)

ggsave("figures/num_WAIC_zoom.pdf", plot = Num_host_filtered_delta_zoom, width = 4, height = 4, units = "in", dpi = 300)
ggsave("figures/onum_WAIC_zoom.png", plot = Num_host_filtered_delta_zoom, width = 4, height = 4, units = "in", dpi = 300)


####################### MODELS################################
#Add observations to the model
#Phylogenetic distance

allData_AICc_delta <- allData_WAIC_filled %>%
  filter(Model %in% c("gompertz_survival", "exponential")) %>%
  spread(Model, Mean_WAIC) %>%
  group_by(Full_key, Key, Host_taxa, Pathogen_taxa, Host_other, Pathogen_other) %>%
  summarise(exponential = mean(exponential, na.rm = TRUE), gompertz_survival = mean(gompertz_survival, na.rm =T), 
            Num_Obs = mean(Num_Obs), Num_Host = mean(Max_num_host_filled)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(DeltaSub = exponential - gompertz_survival)

# Ensure that the columns are factors
allData_AICc_delta$Host_taxa <- as.factor(allData_AICc_delta$Host_taxa)
allData_AICc_delta$Pathogen_taxa <- as.factor(allData_AICc_delta$Pathogen_taxa)

allData_AICc_delta$Host_other <- as.factor(allData_AICc_delta$Host_other)
allData_AICc_delta$Pathogen_other <- as.factor(allData_AICc_delta$Pathogen_other)

# Convert Num_Obs and Num_Host to numeric if necessary
allData_AICc_delta$Num_Obs <- as.numeric(as.character(allData_AICc_delta$Num_Obs))
allData_AICc_delta$Num_Host <- as.numeric(as.character(allData_AICc_delta$Num_Host))

# Fit the mixed-effects model with centered predictors
allData_AICc_delta$Host_other <- factor(allData_AICc_delta$Host_other, 
                                      levels = c("Seedlings", "Drosophila sp.", "Other insects", 
                                                 "Nematodes", "Moth larvae", "Other invertebrates", 
                                                 "Fish", "Avian", "Mice", "Other mammals"))

allData_AICc_delta$Pathogen_other <- factor(allData_AICc_delta$Pathogen_other, 
                                           levels = c("Gram-positive bacteria", "Gram-negative bacteria",
                                                      "DNA virus", "RNA virus", "Fungi", "Protozoan parasite"))



lme_model_AIC2 <- lmer(DeltaSub ~ Host_other + Pathogen_other  + (1|Num_Obs) +  (1|Num_Host), 
                      data = allData_AICc_delta)
summary(lme_model_AIC2)

lm_model_AIC3 <- lm(DeltaSub ~ Num_Obs + Num_Host + Host_other + Pathogen_other, 
                       data = allData_AICc_delta)
summary(lm_model_AIC3)

lm_model_AIC4 <- lm(DeltaSub ~ Num_Obs + Num_Host + Host_taxa + Pathogen_taxa, 
                    data = allData_AICc_delta)
summary(lm_model_AIC4)

# Plot fixed effects
sjPlot::plot_model(lme_model_AIC, type = "est", show.values = T, show.p = TRUE)
# Plot random effects
sjPlot::plot_model(lme_model_AIC, type = "re", show.values = T, show.p = TRUE)

# Plot fixed effects
sjPlot::plot_model(lme_model_AIC2, type = "est", show.values = F, show.p = TRUE)
# Plot random effects
sjPlot::plot_model(lme_model_AIC2, type = "re", show.values = T, show.p = TRUE)

# Plot fixed effects
sjPlot::plot_model(lm_model_AIC3, type = "est", show.values = T, show.p = TRUE)
sjPlot::plot_model(lm_model_AIC4, type = "est", show.values = T, show.p = TRUE)

# Plot for the manuscript 

lmCoef <- sjPlot::plot_model(lm_model_AIC3, type = "est", show.values = T, show.p = TRUE)+
  scale_color_manual(values = cbpalette)+
  labs(y = expression(Delta[WAIC]~ estimates), title = NULL)+
  mytheme
