# Meta survival data analysis
# Author: Canan Karakoc
# Last update: Jul 03 2024

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

set.seed(1234)
setwd("~/GitHub/metaAnalysis")

library(patchwork)
library(dplyr) # masked from another package

# Manuscript figures
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

# Loading data
# Data will be uploaded to Figshare
main.dir <- "~/GitHub/metaAnalysis_largeFiles/metaData"
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
  #mutate(Time = round(Time, digits = 1)) %>%
  #mutate(Survival = round(Survival, digits = 1)) %>% # I digitized too precise that digits are long
  mutate(across(Time:Survival, ~ ifelse(.x < 0, 0, .x))) %>%
  mutate(across(Survival, ~ ifelse(.x > 100, 100, .x)))

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
    Host_taxa, Pathogen_taxa, Max_num_host, Data
  )

#check the problematic data 
check1 <- standardized_final_df %>%
  filter(Key == "59XSYFHS")

check2 <- dataCleanSurvival_hours %>%
  filter(Key == "59XSYFHS")

check3 <- dataClean %>%
  filter(Key == "59XSYFHS")

check4 <- dataAll %>%
  filter(Full_key == "59XSYFHS.csv")

write.csv(standardized_final_df, "data/standardized_final_df.csv", row.names = FALSE)

# Number of datasets after filtering
unique(standardized_final_df$Full_key)
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
write.csv(all_results, "data/results_df.csv", row.names = FALSE)
write.csv(all_predictions, "data/all_predictions.csv", row.names = FALSE)

#####################################################################################
# READ SAVED DATA
standardized_final_df <- read.table("data/standardized_final_df.csv", header = T, sep = ",", dec = ".")
all_results <- read.table("data/results_df.csv", header = T, sep = ",", dec = ".")
all_predictions <- read.table("data/all_predictions.csv", header = T, sep = ",", dec = ".")

names(all_results)[1] <- "Full_key"

freq <- all_results %>%
  group_by(Full_key) %>%
  filter(AICc == min(AICc)) %>%
  group_by(Model) %>%
  summarise(n = n()) 

# Rank all models for each group by Mean_WAIC
ranked_results <- all_results %>%
  group_by(Full_key) %>%
  arrange(AICc) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Extract where 'gompertz_hazard' is not the best or if it is the best, take the second best
adjusted_results <- ranked_results %>%
  group_by(Full_key) %>%
  summarise(
    best_model = if_else(any(Model == "gompertz_hazard" & rank == 1), 
                         Model[rank == 2], # Take second best if gompertz_hazard is the best
                         Model[rank == 1]  # Else take the best as usual
    ),
    .groups = 'drop'
  )

#find <- adjusted_results %>% filter(best_model == "gompertz_survival")

# Count how many times each model ends up as the best after adjustments
model_counts <- adjusted_results %>%
  group_by(best_model) %>%
  summarise(n = n(), .groups = 'drop')

# Density plot for Mean_WAIC by Model
all_results$Model <- factor(all_results$Model, 
                            levels = c("exponential", "gompertz_hazard", "gompertz_survival", "weibull", "loglogistic", "generalizedgamma"))

# Plot with ggplot2
density2 <- ggplot(data = all_results %>% filter(Model != "gompertz_hazard"), aes(x = AICc, color = Model)) +
  geom_density(alpha = 0.5, size = 1) +
  labs(
    x = "AICc",
    y = "Density") +
  mytheme +
  scale_color_manual(values = cbpalette, labels = c("Constant (14)", "Gompertz (77)", "Weibull (41)", "Log-logistic (67)", "Gen.-gamma (10)")) +
  scale_x_continuous(limits = c(-200, 100)) +
  theme(legend.position = c(0.1, 0.5), 
        legend.justification = c(0, 0),
        legend.key = element_blank()) +
  guides(color = guide_legend(override.aes = list(linetype = 1)))

ggsave("figures/densityAICc.png", plot = density2, width = 5.5, height = 4.5, units = "in", dpi = 300)
ggsave("figures/densityAICc.pdf", plot = density2, width = 5.5, height = 4.5, units = "in", dpi = 300)

# Data
colnames(all_results)[1] <- "Full_key"

df2_first_match <- standardized_final_df %>%
  group_by(Full_key) %>%
  slice(1) %>%
  ungroup()

allData_WAIC <- all_results[-1, ] %>%
  left_join(df2_first_match, by = "Full_key") %>%
  dplyr::select(Full_key, Key, Host_taxa, Pathogen_taxa, Num_Obs, Max_num_host, Model, Mean_WAIC, AICc, Data)

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

#species_SLELFBN3 <-  c("Drosophila ananassae", "Drosophila ananassae", " Drosophila ananassae", "Drosophila euronotus", "Drosophila flavomontana",
#"Drosophila mauritiana", "Drosophila mauritiana", "Drosophila mauritiana", "Drosophila melanogaster", "Drosophila melanogaster",
#"Drosophila melanogaster", "Drosophila prosaltans", "Drosophila prosaltans", "Drosophila santomea", "Drosophila santomea",    
#"Drosophila santomea", "Drosophila simulans", "Drosophila simulans", "Drosophila simulans", "Drosophila subobscura",  
#"Drosophila subobscura", "Drosophila subobscura", "Drosophila tropicalis", "Drosophila tropicalis", "Drosophila tropicalis",  
#"Hirtodrosophila duncani", "Hirtodrosophila duncani", "Hirtodrosophila duncani", "Zaprionus taronus", "Zaprionus taronus")

# Create mergedata data frame
mergedata <- data.frame(key = datasets, Max_num_hosts = num_hosts)

# Join and fill NAs
allData_WAIC_filled <- allData_WAIC %>%
  mutate(key = gsub("\\..*", "", Full_key)) %>%
  left_join(mergedata, by = "key") %>%
  mutate(Max_num_host_filled = coalesce(as.numeric(Max_num_host), as.numeric(Max_num_hosts)))

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

# HOST
allData_WAIC_sum <- allData_WAIC_filled %>%
  filter(
    is.finite(as.numeric(AICc)),
    is.finite(as.numeric(Num_Obs)),
    is.finite(as.numeric(Max_num_host_filled))
  ) %>%
  group_by(Model, Host_taxa) %>%
  summarize(
    meanAIC = mean(as.numeric(AICc), na.rm = TRUE),
    seAIC = standard_error(as.numeric(AICc)),
    meanObs = round(mean(as.numeric(Num_Obs), na.rm = TRUE), 0),
    meanHost = round(mean(as.numeric(Max_num_host_filled), na.rm = TRUE), 0),
    numData = length(unique(Full_key))
  ) %>%
  ungroup() %>%
  filter(!Model == "gompertz_hazard") 

# Order factor levels
allData_WAIC_sum$Host_taxa <- factor(allData_WAIC_sum$Host_taxa, 
                                      levels = c("Seedlings", "Drosophila sp.", "Other insects", 
                                                 "Nematodes", "Moth larvae", "Other invertebrates", 
                                                 "Fish", "Avian", "Mice", "Other mammals"))

allData_WAIC_sum$Model <- factor(allData_WAIC_sum$Model, 
                                     levels = c("exponential", "gompertz_survival", "weibull", 
                                                "loglogistic", "generalizedgamma"))

# Create the plot and add the annotations
host_waic <- ggplot(allData_WAIC_sum, aes(x = Model, y = meanAIC, color = Model)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = meanAIC - seAIC, ymax = meanAIC + seAIC), width = 0.2) +
  facet_wrap(~Host_taxa, ncol = 2, scales = "free") +
  mytheme +
  theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.title.x = element_blank()) +
  ylab("AICc") +
  geom_text(aes(label = paste0("num. datasets = ", numData)), x = Inf, y = Inf, hjust = 1.1, vjust = 1.2, size = 4, check_overlap = TRUE, show.legend = FALSE) +
  geom_text(aes(label = paste0("num. time obs. = ", meanObs)), x = Inf, y = Inf, hjust = 1.1, vjust = 2.3, size = 4, check_overlap = TRUE, show.legend = FALSE) +
  geom_text(aes(label = paste0("avg. host reps = ", meanHost)), x = Inf, y = Inf, hjust = 1.1, vjust = 3.4, size = 4, check_overlap = TRUE, show.legend = FALSE) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))+
  scale_color_manual(values = cbpalette, labels = c("Constant", "Gompertz", "Weibull", "Log-logistic", "Gen.-gamma")) 

ggsave("figures/host_aic.png", plot = host_waic, width = 6, height = 12, units = "in", dpi = 300)
ggsave("figures/host_aic.pdf", plot = host_waic, width = 6, height = 12, units = "in", dpi = 300)

# PATGOGEN
# Order factor levels
PallData_WAIC_sum <- allData_WAIC_filled %>%
  filter(
    is.finite(as.numeric(AICc)),
    is.finite(as.numeric(Num_Obs)),
    is.finite(as.numeric(Max_num_host_filled))
  ) %>%
  group_by(Model, Pathogen_taxa) %>%
  summarize(
    meanAIC = mean(as.numeric(AICc)),
    seAIC = standard_error(as.numeric(AICc)),
    meanObs = round(mean(as.numeric(Num_Obs)), 0),
    meanHost = round(mean(as.numeric(Max_num_host_filled), na.rm = T), 0),
    numData = length(unique(Full_key))
  ) %>%
  filter(!Model == "gompertz_hazard") 

PallData_WAIC_sum$Pathogen_taxa <- factor(PallData_WAIC_sum$Pathogen_taxa, 
                                          levels = c("Gram-positive bacteria", "Gram-negative bacteria",
                                                     "DNA virus", "RNA virus", "Fungi", "Protozoan parasite"))


PallData_WAIC_sum$Model <- factor(PallData_WAIC_sum$Model, 
                                 levels = c("exponential", "gompertz_survival", "weibull", 
                                            "loglogistic", "generalizedgamma"))


# Create the plot and add the annotations
pathogen_waic <- ggplot(PallData_WAIC_sum, aes(x = Model, y = meanAIC, color = Model)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = meanAIC - seAIC, ymax = meanAIC + seAIC), width = 0.2) +
  facet_wrap(~Pathogen_taxa, ncol = 2, scales = "free") +
  mytheme +
  theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = cbpalette, labels = c("Constant", "Gompertz", "Weibull", "Log-logistic", "Gen.-gamma")) +
  ylab("AICc") +
  geom_text(aes(label = paste0("num. datasets = ", numData)), x = Inf, y = Inf, hjust = 1.1, vjust = 1.2, size = 4, check_overlap = TRUE, show.legend = FALSE) +
  geom_text(aes(label = paste0("num. time obs. = ", meanObs)), x = Inf, y = Inf, hjust = 1.1, vjust = 2.3, size = 4, check_overlap = TRUE, show.legend = FALSE) +
  geom_text(aes(label = paste0("avg. host reps = ", meanHost)), x = Inf, y = Inf, hjust = 1.1, vjust = 3.4, size = 4, check_overlap = TRUE, show.legend = FALSE) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave("figures/pathogen_aic.png", plot = pathogen_waic, width = 6, height = 8, units = "in", , dpi = 300)
ggsave("figures/pathogen_aic.pdf", plot = pathogen_waic, width = 6, height = 8, units = "in", , dpi = 300)

# Data_type
Data_WAIC <- allData_WAIC_filled %>%
  filter(
    is.finite(as.numeric(AICc)),
  ) %>%
  group_by(Model, Data) %>%
  summarize(meanAIC = mean(AICc, na.rm = TRUE), seAIC = standard_error(AICc)) %>%
  filter(!Model == "gompertz_hazard")

Data_WAIC$Model <- factor(Data_WAIC$Model, 
                                 levels = c("exponential", "gompertz_survival", "weibull", 
                                            "loglogistic", "generalizedgamma"))

Datatype <- ggplot(Data_WAIC, aes(x = Data, y = meanAIC, color = Model)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = meanAIC - seAIC, ymax = meanAIC + seAIC),
    width = 0.2,
    position = position_dodge(width = 0.5, preserve = "single")
  ) +
  mytheme +
  theme(axis.title.x = element_blank()) +
  ylab("AICc")+
  scale_color_manual(values = cbpalette, labels = c("Constant","Gompertz", "Weibull", "Log-logistic", "Generalized-gamma"))+
theme(
  legend.position = c(0.1, 0.1), # Positioning the legend inside the plot
  legend.justification = c(0, 0) # Aligning the legend to the bottom left
)

ggsave("figures/datatypeAICc.pdf", plot = Datatype, width = 5.5, height = 4.5, units = "in", dpi = 300)
ggsave("figures/datatypeAICc.png", plot = Datatype, width = 5.5, height = 4.5, units = "in", dpi = 300)

Data_num <- allData_WAIC_filled %>%
  filter(Model == "gompertz_hazard")%>%
  group_by(Data) %>%
  summarize(Data = n())

# Num. of Observations
allData_WAIC_sum_obs <- allData_WAIC_filled %>%
  group_by(Full_key, Model) %>%
  summarize(
    meanObs = mean(as.numeric(Num_Obs), na.rm = T),
    meanAIC = mean(AICc, na.rm = T)
  ) %>%
  ungroup()

allData_WAIC_sum_num <- allData_WAIC_filled %>%
  group_by(Full_key, Model) %>%
  summarize(
    meanNum = mean(as.numeric(Max_num_host_filled), na.rm = T),
    meanAIC = mean(AICc, na.rm = T)
  ) %>%
  ungroup()

# Supplementary figure
# Create combined_factor and plotting_data
standardized_final_df$combined_factor <- interaction(standardized_final_df$Host_taxa, standardized_final_df$Pathogen_taxa, sep = "\n")
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
  filter(Full_key == "8C5XBCN2A.csv")
fig1_dat2 <- prediction_plotting_data %>%
  filter(Full_key == "8C5XBCN2A.csv") %>%
  filter(!Model == "gompertz_hazard")

fig1_dat2$Model <- factor(fig1_dat2$Model, 
                                 levels = c("exponential", "gompertz_survival", "weibull", 
                                            "loglogistic", "generalizedgamma"))


fig1 <- ggplot(fig1_dat1, aes(x = Time_std, y = Standardized_Survival)) +
  geom_point(size = 4, color = "grey25") + # Raw data
  geom_line(data = fig1_dat2, aes(y = Predictions, color = Model), size = 1) + # Predictions
  labs(x = "Standardized time", y = "Standardized survival") +
  mytheme +
  scale_color_manual(values = cbpalette, labels = c("Constant", "Gompertz", "Weibull", "Log-logistic", "Gen.-gamma"))+
  theme(
    legend.position = c(0.1, 0.1), # Positioning the legend inside the plot
    legend.justification = c(0, 0) # Aligning the legend to the bottom left
)

ggsave("figures/fig1.pdf", plot = fig1, width = 6, height = 5.5, units = "in", dpi = 300)
aic <- allData_WAIC_filled %>% filter(Full_key == "8C5XBCN2A.csv")

# Plotting

prediction_plotting_data <- prediction_plotting_data %>%
  filter(!Model == "gompertz_hazard")
prediction_plotting_data$Model <- factor(prediction_plotting_data$Model, 
                          levels = c("exponential", "gompertz_survival", "weibull", 
                                     "loglogistic", "generalizedgamma"))

ALL <- ggplot(plotting_data, aes(x = Time_std, y = Standardized_Survival)) +
  geom_point(size = 3, color = "grey25") +  # Raw data
  geom_line(data = prediction_plotting_data, aes(y = Predictions, x = Time_std, color = Model), size = 0.7) +  # Predictions
  labs(x = "Standardized time", y = "Standardized survival") +
  facet_wrap(~ combined_factor + key) +
  scale_color_manual(values = cbpalette, labels = c("Constant","Gompertz", "Weibull", "Log-logistic", "Generalized-gamma"))+
  mytheme+
  theme(legend.position = "bottom")

# Pagination parameters
ncol <- 5
nrow <- 6
n_pages <- ceiling(length(unique(standardized_final_df$Full_key)) / (ncol * nrow))

# Loop to create paginated plots
for (i in 1:n_pages) {
  paginated_plot <- ALL + facet_wrap_paginate(~ combined_factor + key, ncol = ncol, nrow = nrow, page = i)
  
  # Save each page to a file
  ggsave(paste0("combined_plot_page_", i, ".png"), paginated_plot, width = 16, height = 17)
  
  # Print each paginated plot
  print(paginated_plot)
}

######################################################################
# Delta WAIC
allData_AICc_delta <- allData_WAIC_filled %>%
  filter(Model %in% c("gompertz_survival", "exponential")) %>%
  spread(Model, AICc) %>%
  group_by(Full_key, Key, Host_taxa, Pathogen_taxa) %>%
  summarise(exponential = mean(exponential, na.rm = TRUE), gompertz_survival = mean(gompertz_survival, na.rm =T), 
            Num_Obs = mean(Num_Obs), Num_Host = mean(Max_num_host_filled)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(DeltaSub = exponential - gompertz_survival)

allData_AICc_delta_sum <- allData_AICc_delta %>%
  group_by(Host_taxa) %>%
  summarize(
    delta_AIC = mean(DeltaSub, na.rm = T),
    se = standard_error(DeltaSub)) 

allData_AICc_delta_sum_P <- allData_AICc_delta %>%
  group_by(Pathogen_taxa) %>%
  summarize(
    delta_AIC = mean(DeltaSub, na.rm = T),
    se = standard_error(DeltaSub))  

allData_AICc_delta_sum$Host_taxa <- factor(allData_AICc_delta_sum2$Host_taxa, 
                                             levels = rev(c("Seedlings", "Drosophila sp.", "Other insects", 
                                                        "Nematodes", "Moth larvae", "Other invertebrates", 
                                                        "Fish", "Avian", "Mice", "Other mammals")))

# Reorder the factor levels based on the values
allData_AICc_delta_sum$Host_taxa <- reorder(allData_AICc_delta_sum$Host_taxa , allData_AICc_delta_sum$delta_AIC)

delta_host <- ggplot(allData_AICc_delta_sum, aes(y = Host_taxa, x = delta_AIC)) +
  geom_point(size=5, shape=21)+
  geom_errorbar(aes(xmin = delta_AIC - se, xmax = delta_AIC + se), width = 0.2) +
  mytheme+
  labs(x = expression(Delta[AICc]), y = NULL) +
  theme(plot.title = element_text(hjust = 0.5))

allData_AICc_delta_sum_P$Pathogen_taxa <- reorder(allData_AICc_delta_sum_P$Pathogen_taxa , allData_AICc_delta_sum_P$delta_AIC)

delta_path <- ggplot(allData_AICc_delta_sum_P, aes(y = Pathogen_taxa, x = delta_AIC)) +
  geom_point(size = 5, shape = 21)+
  geom_errorbar(aes(xmin = delta_AIC - se, xmax = delta_AIC + se), width = 0.2) +
  mytheme+
  xlab(expression(Delta[AICc])) +
  scale_color_manual(values = cbpalette)+
  labs(y = NULL) +
  theme(plot.title = element_text(hjust = 0.5))

comb_delta <- delta_host+delta_path+
plot_annotation(tag_levels = 'A')

ggsave("figures/hostpathAIC.pdf", plot = comb_delta, width = 10, height = 4, units = "in", dpi = 300)
ggsave("figures/hostpathAIC.png", plot = comb_delta, width = 10, height = 4, units = "in", dpi = 300)

# Delta Model i - Exponential for other models 
allModels_delta <- allData_WAIC_filled %>%
  filter(Model != "gompertz_hazard") %>%
  spread(Model, AICc) %>%
  group_by(Full_key, Key, Host_taxa, Pathogen_taxa) %>%
  summarise(exponential = mean(exponential, na.rm = TRUE), gompertz_survival = mean(gompertz_survival, na.rm =T), 
            generalizedgamma = mean(generalizedgamma, na.rm = TRUE), loglogistic = mean(loglogistic, na.rm = TRUE),
            weibull = mean(weibull, na.rm = TRUE), Num_Obs = mean(Num_Obs), Num_Host = mean(Max_num_host_filled)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(DeltaGompertz = exponential - gompertz_survival, 
         DeltaWeibull  = exponential - weibull, 
         DeltaLoglogistic = exponential - loglogistic, 
         DeltaGamma = exponential - generalizedgamma) %>%
  pivot_longer(DeltaGompertz:DeltaGamma, names_to = "comparison", values_to = "DeltaAIC")%>%
filter_all(all_vars(!is.na(.))) %>%
  filter_all(all_vars(!is.nan(.))) %>%
  filter_all(all_vars(!is.infinite(.)))

allModels_delta$comparison <- factor(allModels_delta$comparison, 
        levels = c("DeltaGompertz", "DeltaWeibull", 
        "DeltaLoglogistic", "DeltaGamma"))


op <- ggplot(allModels_delta, aes(y = DeltaAIC, x = comparison)) +
  geom_hline (yintercept = 0, linetype = "dashed")+
  geom_jitter(size = 3, shape = 21, color = "grey50")+
  geom_boxplot(fill = "white", , alpha  = 0.7)+
  mytheme+
  xlab(expression(Delta[AICc])) +
  labs(y = expression(Delta[AICc]), x = NULL)+
  scale_x_discrete(labels = c("Gompertz", "Weibull", "Log-logistic", "Gen.-gamma"))+
  scale_y_continuous(limits = c(-30,30))

ggsave("figures/opAIC.pdf", plot = op, width = 6, height = 4, units = "in", dpi = 300)
ggsave("figures/opAIC.png", plot = op, width = 6, height = 4, units = "in", dpi = 300)

allData_WAIC_filled$Model <- factor(allData_WAIC_filled$Model, 
                                         levels = c("exponential", "gompertz_survival", "gompertz_hazard", "weibull", 
                                                    "loglogistic", "generalizedgamma"))

Comp_gomp <- ggplot(allData_WAIC_filled %>% filter(Model != "gompertz_hazard"), 
                    aes(y = AICc, x = Model))+
  geom_jitter(size = 3, shape = 21, color = "grey50")+
  geom_boxplot(fill = "white", , alpha  = 0.7)+
  mytheme+
  labs(y = "AICc", x = NULL)+
  scale_x_discrete(labels = c("Constant", "Gompertz", "Weibull", "Log-logistic", "Gen.-gamma"))

ggsave("figures/delta_all_aic.pdf", plot = Comp_gomp, width = 6, height = 4, units = "in", dpi = 300)
ggsave("figures/delta_all_aic.png", plot = Comp_gomp, width = 6, height = 4, units = "in", dpi = 300)

#Comparison with Gompertz 
allModels_delta_g <- allData_WAIC_filled %>%
  filter(Model != "gompertz_hazard") %>%
  spread(Model, AICc) %>%
  group_by(Full_key, Key, Host_taxa, Pathogen_taxa) %>%
  summarise(exponential = mean(exponential, na.rm =T),
    gompertz_survival = mean(gompertz_survival, na.rm =T), 
            generalizedgamma = mean(generalizedgamma, na.rm = TRUE), loglogistic = mean(loglogistic, na.rm = TRUE),
            weibull = mean(weibull, na.rm = TRUE), Num_Obs = mean(Num_Obs), Num_Host = mean(Max_num_host_filled)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(Exponential = gompertz_survival - exponential,
         Weibull  = gompertz_survival - weibull, 
         Loglogistic = gompertz_survival - loglogistic, 
         Gen.gamma = gompertz_survival -generalizedgamma) %>%
  pivot_longer(Exponential:Gen.gamma, names_to = "comparison", values_to = "DeltaAIC")%>%
  filter_all(all_vars(!is.na(.))) %>%
  filter_all(all_vars(!is.nan(.))) %>%
  filter_all(all_vars(!is.infinite(.)))

delta_all_g2 <- ggplot(allModels_delta_g, aes(y = DeltaAIC, x = comparison)) +
  geom_hline (yintercept = 0, linetype = "dashed")+
  geom_jitter(size = 3, shape = 21, color = "grey50")+
  geom_boxplot(fill = "white", , alpha  = 0.7)+
  mytheme+
  labs(y = expression(AICc[Gompertz] - AICc[Model]), x = NULL)+
  scale_x_discrete(labels = c("Gompertz", "Weibull", "Log-logistic", "Gen.-gamma"))+
  scale_y_continuous(limits = c(-30, 30))

ggsave("figures/delta_all_g2-limited.pdf", plot = delta_all_g2, width = 6, height = 4, units = "in", dpi = 300)
ggsave("figures/delta_all_g2-limited.png", plot = delta_all_g2, width = 6, height = 4, units = "in", dpi = 300)

# Figure 2 

# Combine the plots and add labels
op <- op + theme(axis.text.x = element_text(angle = 45, hjust = 1))
delta_all_g2 <- delta_all_g2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot <- density2 + delta_all_g2 + op +
  plot_annotation(tag_levels = 'A')

ggsave("figures/combined_plot_fig2.pdf", plot = combined_plot, width = 15, height = 5, units = "in", dpi = 300)
ggsave("figures/combined_plot_fig2.png", plot = combined_plot, width = 15, height = 5, units = "in", dpi = 300)

summary(lm(DeltaAIC~comparison, data = allModels_delta))
anova(lm(DeltaAIC~comparison, data = allModels_delta))

summary(lm(DeltaAIC~comparison, data = allModels_delta_g))
anova(lm(DeltaAIC~comparison, data = allModels_delta_g))

mod_all <- all_results %>% filter(Model != "gompertz_hazard")%>%
  filter_all(all_vars(!is.na(.))) %>%
  filter_all(all_vars(!is.nan(.))) %>%
  filter_all(all_vars(!is.infinite(.)))

summary(lm(AICc~Model, data = mod_all))
anova(lm(AICc~Model, data = mod_all))

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
  geom_point(size = 3, shape = 1, position = position_dodge(width = 0.5), alpha = 0.5) +
  mytheme +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = cbpalette) +
  ylab(expression(Delta[AICc])) +
  xlab("Number of observations") +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  theme(strip.background = element_blank()) +
  coord_cartesian(xlim = c(5, 50))

Num_host_filtered_delta <- ggplot(delta_sum_num, aes(x = as.numeric(meanNum), y = meanEvidence)) +
  geom_point(size = 3, shape = 1, position = position_dodge(width = 0.5), alpha = 0.5) +
  mytheme +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = cbpalette) +
  ylab(expression(Delta[AICc])) +
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
  coord_cartesian(xlim = c(5, 100))


comb_numb <- Obs_delta + Num_host_filtered_delta +
  plot_annotation(tag_levels = 'A', theme = theme(plot.tag = element_text(face = "bold")))

ggsave("figures/obs_numAIC.pdf", plot = comb_numb, width = 8, height = 3.5, units = "in", dpi = 300)
ggsave("figures/obs_numAIC.png", plot = comb_numb, width = 8, height = 3.5, units = "in", dpi = 300)


####################### MODELS################################
#Add observations to the model
#Phylogenetic distance

allData_AICc_delta <- allData_WAIC_filled %>%
  filter(Model %in% c("gompertz_survival", "exponential")) %>%
  spread(Model, AICc) %>%
  group_by(Full_key, Key, Host_taxa, Pathogen_taxa) %>%
  summarise(exponential = mean(exponential, na.rm = TRUE), gompertz_survival = mean(gompertz_survival, na.rm =T), 
            Num_Obs = mean(Num_Obs), Num_Host = mean(Max_num_host_filled)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(DeltaSub = exponential - gompertz_survival)

# Ensure that the columns are factors
allData_AICc_delta$Host_taxa <- as.factor(allData_AICc_delta$Host_taxa)
allData_AICc_delta$Pathogen_taxa <- as.factor(allData_AICc_delta$Pathogen_taxa)

# Convert Num_Obs and Num_Host to numeric if necessary
allData_AICc_delta$Num_Obs <- as.numeric(as.character(allData_AICc_delta$Num_Obs))
allData_AICc_delta$Num_Host <- as.numeric(as.character(allData_AICc_delta$Num_Host))

# Fit the mixed-effects model with centered predictors
allData_AICc_delta$Host_taxa <- factor(allData_AICc_delta$Host_taxa , 
                                      levels = c("Seedlings", "Drosophila sp.", "Other insects", 
                                                 "Nematodes", "Moth larvae", "Other invertebrates", 
                                                 "Fish", "Avian", "Mice", "Other mammals"))

allData_AICc_delta$Pathogen_taxa <- factor(allData_AICc_delta$Pathogen_taxa, 
                                           levels = c("Gram-positive bacteria", "Gram-negative bacteria",
                                                      "DNA virus", "RNA virus", "Fungi", "Protozoan parasite"))



lme_model_AIC2 <- lmer(DeltaSub ~ Host_taxa + Pathogen_taxa  + (1|Num_Obs) +  (1|Num_Host), 
                      data = allData_AICc_delta)
summary(lme_model_AIC2)

lm_model_AIC3 <- lm(DeltaSub ~ Num_Obs + Num_Host + Host_taxa + Pathogen_taxa, 
                       data = allData_AICc_delta)
summary(lm_model_AIC3)
anova(lm_model_AIC3)

library(car)
vif(lm_model_AIC3)

# Plot fixed effects
sjPlot::plot_model(lme_model_AIC2, type = "est", show.values = F, show.p = TRUE)
# Plot random effects
sjPlot::plot_model(lme_model_AIC2, type = "re", show.values = T, show.p = TRUE)
# Plot fixed effects
sjPlot::plot_model(lm_model_AIC3, type = "est", show.values = T, show.p = TRUE)

# Plot for the manuscript 
# I would not use this. Use the distance as a covariate. 
lmCoef <- sjPlot::plot_model(lm_model_AIC3, type = "est", show.values = T, show.p = TRUE)+
  scale_color_manual(values = cbpalette)+
  labs(y = expression(Delta[AICc]~ estimates), title = NULL)+
  mytheme

model_data <- lm_model_AIC3$model

# Modify the factor levels
model_data$Host_taxa <- factor(model_data$Host_taxa, 
                        levels = (c("Seedlings",  
                        "Nematodes", "Other invertebrates", "Drosophila sp.", "Moth larvae", "Other insects",
                        "Avian", "Mice", "Other mammals", "Fish")))

model_data$Pathogen_taxa <- factor(model_data$Pathogen_taxa, 
                        levels = (c("Protozoan parasite",  "Fungi", 
                                             "Gram-positive bacteria", "Gram-negative bacteria", "DNA virus", "RNA virus")))
# Re-fit the model with modified data if necessary
lm_model_AIC3 <- update(lm_model_AIC3, data = model_data)

# Plot with modified factor levels
lmCoef <- sjPlot::plot_model(lm_model_AIC3, type = "est", show.values = TRUE, show.p = TRUE) +
  scale_color_manual(values = cbpalette) +
  labs(y = expression(Delta[AICc]~ estimates), title = NULL) +
  mytheme+
  scale_x_discrete(labels = function(x) gsub("_taxa", "", sub(".*_taxa", "", x))) 

# Models with deviation contrasts 

# Set the contrast options for categorical variables
options(contrasts = c("contr.sum", "contr.poly"))

# Fit the linear model with deviation coding
model_deviation <- lm(DeltaSub ~ Num_Obs + Num_Host + Host_taxa + Pathogen_taxa, data = allData_AICc_delta)
summary(model_deviation)
anova(model_deviation)

# Plot the model estimates with actual level names
plot_model(model, type = "est", show.values = TRUE, value.offset = 0.4) +
  labs(title = "Model with Deviation Coding",
       y = "Effect Sizes (ΔAICc estimates)",
       x = "")

contrasts(allData_AICc_delta$Host_taxa)
contrasts(allData_AICc_delta$Pathogen_taxa)

library(broom)

# Fit the model
# Fit the model with deviation coding
allData_AICc_delta$Host_taxa <- factor(allData_AICc_delta$Host_taxa, 
                                       levels = c("Seedlings", "Nematodes", "Other invertebrates", "Drosophila sp.", "Moth larvae", "Other insects", "Avian", "Mice", "Other mammals", "Fish"))

allData_AICc_delta$Pathogen_taxa <- factor(allData_AICc_delta$Pathogen_taxa, 
                                           levels = c("Protozoan parasite", "Fungi", "Gram-positive bacteria", "Gram-negative bacteria", "DNA virus", "RNA virus"))

# Set contrasts to deviation coding
contrasts(allData_AICc_delta$Host_taxa) <- contr.sum(length(levels(allData_AICc_delta$Host_taxa)))
contrasts(allData_AICc_delta$Pathogen_taxa) <- contr.sum(length(levels(allData_AICc_delta$Pathogen_taxa)))

# Fit the model
model <- lm(DeltaSub ~ Num_Obs + Num_Host + Host_taxa + Pathogen_taxa, data = allData_AICc_delta)

#Extract model estimates
model_tidy <- tidy(model)

# Calculate the reference level effect size for Host_taxa and Pathogen_taxa
host_taxa_ref <- -sum(model_tidy %>% filter(str_detect(term, "Host_taxa")) %>% pull(estimate))
pathogen_taxa_ref <- -sum(model_tidy %>% filter(str_detect(term, "Pathogen_taxa")) %>% pull(estimate))

# Create data frame for reference levels
ref_levels <- tibble(
  term = c("Host_taxaRef", "Pathogen_taxaRef"),
  estimate = c(host_taxa_ref, pathogen_taxa_ref),
  std.error = c(NA, NA),  # Standard errors are not available for manually calculated reference levels
  statistic = c(NA, NA),
  p.value = c(NA, NA)
)

# Combine with original model estimates
model_tidy <- bind_rows(model_tidy, ref_levels)

# Manually replace term names using case_when
model_tidy <- model_tidy %>%
  mutate(term = case_when(
    term == "Num_Obs" ~ "Num Obs",
    term == "Num_Host" ~ "Num Host",
    term == "Host_taxa1" ~ "Nematodes",
    term == "Host_taxa2" ~ "Other invertebrates",
    term == "Host_taxa3" ~ "Drosophila sp.",
    term == "Host_taxa4" ~ "Moth larvae",
    term == "Host_taxa5" ~ "Other insects",
    term == "Host_taxa6" ~ "Avian",
    term == "Host_taxa7" ~ "Mice",
    term == "Host_taxa8" ~ "Other mammals",
    term == "Host_taxa9" ~ "Fish",
    term == "Host_taxaRef" ~ "Seedlings",
    term == "Pathogen_taxa1" ~ "Fungi",
    term == "Pathogen_taxa2" ~ "Gram-positive bacteria",
    term == "Pathogen_taxa3" ~ "Gram-negative bacteria",
    term == "Pathogen_taxa4" ~ "DNA virus",
    term == "Pathogen_taxa5" ~ "RNA virus",
    term == "Pathogen_taxaRef" ~ "Protozoan parasite",
    TRUE ~ term  # keep the original term if no match
  ))

# Reorder terms based on factor levels
term_order <- rev(c("Seedlings", "Nematodes", "Other invertebrates", "Drosophila sp.", "Moth larvae", "Other insects", "Avian", "Mice", "Other mammals", "Fish", 
                "Protozoan parasite", "Fungi", "Gram-positive bacteria", "Gram-negative bacteria", "DNA virus", "RNA virus", 
                "Num Obs", "Num Host"))

model_tidy <- model_tidy %>%
  mutate(term = factor(term, levels = term_order)) %>%
  filter(term != "(Intercept)") %>%
  mutate(sign = ifelse(estimate > 0, "Positive", "Negative"))

# Plot using ggplot2
ggplot(model_tidy, aes(x = estimate, y = term, color = sign)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
  labs(
       x = "Effect Sizes (ΔAICc estimates)",
       y = "") +
  scale_color_manual(values = c("Positive" = "#0072B2", "Negative" = "#D55E00")) +
  geom_text(aes(label = round(estimate, 2)), vjust = -0.5, size = 4) +
  mytheme+
  theme(legend.position = "none")


### PHYLO ANALYSIS ###
# Distance/tree 
library(ape)

newick_tree <- "(((Pisaster_ochraceus:634.70000000,((((Oncorhynchus_mykiss:15.06984000,Oncorhynchus_tshawytscha:15.06984000)'14':190.03834000,(Micropterus_salmoides:112.23406000,(Paralichthys_olivaceus:103.76124000,(Oreochromis_niloticus:91.74830000,Poecilia_reticulata:91.74830000)'13':12.01294000)'25':8.47282000)'37':92.87412000)'36':18.90948000,(Silurus_glanis:142.11161000,Danio_rerio:142.11161000)'35':81.90605000)'34':204.98234000,((Bos_taurus:94.00000000,((Oryctolagus_cuniculus:78.97061000,(Mus_musculus:69.91125000,Cavia_porcellus:69.91125000)'43':9.05936000)'33':8.22939000,(Macaca_mulatta:3.31741000,Macaca_fascicularis:3.31741000)'51':83.88259000)'50':6.80000000)'49':224.95000000,Gallus_gallus:318.95000000)'57':110.05000000)'60':205.70000000)'56':72.90436000,((((Melanoplus_borealis:383.15000000,((Apis_mellifera_macedonica:343.05330000,((((((Drosophila_flavomontana:29.59186000,((Drosophila_prosaltans:21.56067000,Drosophila_tropicalis:21.56067000)'48':5.99910000,((Drosophila_ananassae:25.57846000,(Drosophila_suzukii:23.00787000,(((Drosophila_simulans:4.62390000,Drosophila_mauritiana:4.62390000)'63':0.00000000,Drosophila_melanogaster:4.62390000)'47':6.94583000,Drosophila_santomea:11.56973000)'68':11.43814000)'73':2.57059000)'72':0.62224000,Drosophila_subobscura:26.20070000)'84':1.35907000)'87':2.03209000)'83':14.22660000,Drosophila_euronotus:43.81846000)'82':0.62949000,Zaprionus_taronus:44.44795000)'80':20.67730000,Penicillidia_conspicua:65.12525000)'79':175.79475000,Anopheles_arabiensis:240.92000000)'94':58.47011000,Galleria_mellonella:299.39011000)'92':44.60989000)'78':17.53224000,Meccus_pallidipennis:361.53224000)'98':21.61776000)'97':121.40000000,Daphnia_magna:504.55000000)'77':10.82000000,(Penaeus_vannamei:364.66874000,(Macrobrachium_nipponense:331.97658000,Procambarus_clarkii:331.97658000)'71':32.69216000)'67':150.70126000)'66':56.74064000,Caenorhabditis_elegans:572.11064000)'46':135.49372000)'32':890.84039000,(Arctostaphylos_glauca:125.11960000,Cucurbita_pepo:125.11960000)'106':1473.32515000);"
tree <- read.tree(text = newick_tree)
plot(tree, show.tip.label = TRUE)

dist_matrix <- cophenetic(tree)

keys <- read.csv("~/GitHub/metaAnalysis/Phylogeny_list.csv",
                 sep = ",", header = T, na.strings = ""
)

# Load the gplots package for the heatmap.2 function
if (!requireNamespace("gplots", quietly = TRUE)) {
  install.packages("gplots")
}

# Compute average pairwise distance for each species
avg_dist <- rowMeans(dist_matrix)
avg_dist_df <- data.frame(Host = names(avg_dist), AvgDist = avg_dist)

# Perform PCoA
pcoa_res <- ape::pcoa(dist_matrix)

# Extract the first few principal coordinates
pcoa_coords <- pcoa_res$vectors[, 1:2]  # Using the first two coordinates as an example
pcoa_df <- data.frame(Host = rownames(pcoa_coords), pcoa_coords)

# Data is from the main_script
modelData2 <- keys %>%
  full_join(pcoa_df, by = "Host")%>%
  full_join(avg_dist_df, by = "Host")

distanceModel <- allData_AICc_delta %>%
  full_join(modelData2, by = "Full_key")

distanceModel2 <- distanceModel %>%
  filter_all(all_vars(!is.na(.))) %>%
  filter_all(all_vars(!is.nan(.))) %>%
  filter_all(all_vars(!is.infinite(.)))

lm_model_AIC4 <- lm(DeltaSub ~ Num_Obs + Num_Host + Axis.1 + Axis.2 + Pathogen_taxa, 
                    data = distanceModel)
summary(lm_model_AIC4)
sjPlot::plot_model(lm_model_AIC4, type = "est", show.values = T, show.p = TRUE)

anova(lm_model_AIC4)
lm_model_AIC5 <- lm(DeltaSub ~ Num_Obs + Num_Host + Axis.1 + Axis.2, 
                   data = distanceModel)
summary(lm_model_AIC5)

sjPlot::plot_model(lm_model_AIC5, type = "est", show.values = T, show.p = TRUE)

lmCoef_Dist <- sjPlot::plot_model(lm_model_AIC5, type = "est", show.values = T, show.p = TRUE)+
  scale_color_manual(values = cbpalette)+
  labs(y = expression(Delta[AICc]~ estimates), title = NULL)+
  mytheme

lm_model_AIC4 <- lmer(DeltaSub ~ (1|Num_Obs) + (1|Num_Host) + Axis.1 + Axis.2 + Pathogen_taxa, 
                    data = distanceModel)
summary(lm_model_AIC4)
sjPlot::plot_model(lm_model_AIC4, type = "est", show.values = T, show.p = TRUE)
anova(lm_model_AIC4)

ggplot(distanceModel, aes( x = Axis.2, y = DeltaSub))+
  geom_point()+
  stat_smooth(method = "lm")
  
library(reshape2)

# Convert the distance matrix to long format
dist_long <- melt(as.matrix(dist_matrix))
colnames(dist_long) <- c("Species1", "Species2", "Distance")

# Create the heatmap using ggplot2
ggplot(dist_long, aes(x = Species1, y = Species2, fill = Distance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Phylogenetic Distance Heatmap", x = NULL, y = NULL, fill = "Distance")

# Mapped tree for the manuscript 

category_data <- data.frame(
  tip_label = distanceModel$Host,
  category = allData_AICc_delta$Host_taxa
)

# Ensure categories are factors for consistent coloring
category_data$category <- factor(category_data$category)

# Extract tip labels and coordinates from the tree
tip_labels <- data.frame(tip_label = tree$tip.label)
tip_coords <- as.data.frame(tree$edge[tree$edge[, 2] <= Ntip(tree), 2])
names(tip_coords) <- "node"
tip_coords$tip_label <- tree$tip.label

# Merge coordinates with category data
merged_data <- tip_coords %>%
  left_join(category_data, by = "tip_label") %>%
  filter_all(all_vars(!is.na(.)))


# Plot the tree
plot(tree, show.tip.label = FALSE)
tiplabels(tree$tip.label, frame = "none", adj = c(0.5, -0.5))

cbpalette1 <- c("#44AA99","#DDCC77","#E69F00","#CC79A7","#F0E442","#882255",
                "#6699CC", "#332288", "#000000","#888888")


category_colors <- setNames(cbpalette1[1:length(unique(category_data$category))], unique(category_data$category))

category_data$color <- category_colors[category_data$category]

plot_colored_tree <- function(tree, category_data) {
  # Extract tip labels and coordinates from the tree
  tip_labels <- data.frame(tip_label = tree$tip.label)
  tip_coords <- as.data.frame(tree$edge[tree$edge[, 2] <= Ntip(tree), 2])
  names(tip_coords) <- "node"
  tip_coords$tip_label <- tree$tip.label
  
  # Merge coordinates with category data
  merged_data <- tip_coords %>%
    left_join(category_data, by = "tip_label")
  
  # Plot the tree without tip labels
  plot(tree, show.tip.label = FALSE)
  
  # Add colored dots at the tip labels
  for (i in 1:nrow(merged_data)) {
    tip_color <- merged_data$color[i]
    tip_node <- merged_data$node[i]
    tip_coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    points(tip_coords$xx[tip_node], tip_coords$yy[tip_node], col = tip_color, pch = 16, cex = 1.5)
  }
}

plot_colored_tree(tree, category_data)

# Manually define the categories some species, because there are discrepancies between the 
# Synonyms in the phylo list and the tree 
#Poecilia_reticulata <- Fish
#Oreochromis_niloticus (Nile_tilapia) <- Fish
#Meccus_pallidipennis <- Other_insects (G9FYXDG9)

allData_AICc_delta_sum$Host_taxa <- factor(allData_AICc_delta_sum2$Host_taxa, 
                                           levels = rev(c("Seedlings",  
                                                          "Nematodes", "Other invertebrates", "Drosophila sp.", "Moth larvae", "Other insects",
                                                          "Avian", "Mice", "Other mammals", "Fish")))


cbpalette2 <- c("#6699CC","#888888","#000000","#332288","#E69F00","#F0E442","#DDCC77","#882255", 
                "#CC79A7","#44AA99")


delta_host <- ggplot(allData_AICc_delta_sum, aes(y = Host_taxa, x = delta_AIC, fill =Host_taxa)) +
  geom_point(size=5, shape=21, color = "grey25")+
  geom_errorbar(aes(xmin = delta_AIC - se, xmax = delta_AIC + se), width = 0.2, color ="grey25") +
  mytheme+
  labs(x = expression(Delta[AICc]), y = NULL) +
  scale_fill_manual(values = cbpalette2)+
  theme(legend.position = "none")

allData_AICc_delta_bars <- allData_AICc_delta %>%
  group_by(Host_taxa) %>%
  summarise (MeanObs = mean(Num_Obs), MeanHost = mean(Num_Host))

allData_AICc_delta_bars$Host_taxa <- factor(allData_AICc_delta_bars$Host_taxa, 
                                           levels = rev(c("Seedlings",  
                                                          "Nematodes", "Other invertebrates", "Drosophila sp.", "Moth larvae", "Other insects",
                                                          "Avian", "Mice", "Other mammals", "Fish")))


delta_NUM <- ggplot(allData_AICc_delta_bars, aes(y = Host_taxa, x = MeanObs, fill =Host_taxa)) +
  geom_bar(stat = "identity", width = 0.25)+
  mytheme+
  labs(x = "Num. Observations", y = NULL) +
  scale_fill_manual(values = cbpalette2)+
  theme(legend.position = "none")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

library(cowplot)
library(ggplotify)
tree_gg <- as.ggplot(~plot_colored_tree(tree, category_data))

# Combine the plots with adjusted heights and widths
combined_plot <- plot_grid(delta_host, delta_NUM, ncol = 2, align = "h", axis = "tb", rel_heights = c(1, 1), rel_widths = c(2, 1))
print(combined_plot)

# Create the delta_HOST plot with x-axis on top
delta_HOST <- ggplot(allData_AICc_delta_bars, aes(y = Host_taxa, x = MeanHost, color = Host_taxa)) +
  geom_bar(stat = "identity", width = 0.25, fill = "white") +
  mytheme +
  labs(x = "Average Num. Host", y = NULL) +
  scale_color_manual(values = cbpalette2) +
  scale_x_continuous(position = "top") +  # Place the x-axis on top
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x.top = element_line(colour = "black"),
        axis.ticks.x.top = element_line(colour = "black"),
        axis.title.x.top = element_text(margin = margin(t = 10)),
        axis.text.x.top = element_text(margin = margin(t = 10))) +  # Ensure axis text is set
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# Combine the plots with adjusted heights and widths
combined_plot <- plot_grid(delta_host, delta_NUM, delta_HOST, ncol = 3, align = "h", axis = "tb", rel_heights = c(1, 1, 1), rel_widths = c(2, 1, 1))
print(combined_plot)

###########################################################################
# Pathogen 
newick_path_tree <- "((((Pseudomonas_aeruginosa:1481.02000000,(Klebsiella_aerogenes:38.91576000,Serratia_marcescens:38.91576000)'14':1442.10424000)'13':511.86500000,((Burkholderia_cepacia:0.11759000,Burkholderia_cenocepacia:0.11759000)'25':11.46866000,Burkholderia_thailandensis:11.58625000)'37':1981.29875000)'36':1201.28614000,((((Bacillus_anthracis:918.13597000,(Staphylococcus_coagulans:80.86734000,Staphylococcus_aureus:80.86734000)'35':837.26863000)'34':853.87025000,Streptococcus_pyogenes:1772.00622000)'43':1032.87564000,(Clostridioides_difficile:19.38809000,Paeniclostridium_sordellii:19.38809000)'33':2785.49377000)'51':132.41914000,((Corynebacterium_striatum:233.80831000,Renibacterium_salmoninarum:233.80831000)'50':19.73321000,Mycobacterium_tuberculosis:253.54152000)'49':2683.75948000)'57':256.87014000)'60':1055.82886000,((((((Aspergillus_flavus:62.42848000,Aspergillus_fumigatus:62.42848000)'56':289.35300000,((Beauveria_bassiana:154.50000000,Metarhizium_anisopliae:154.50000000)'48':194.54212000,Neofusicoccum_australe:349.04212000)'63':2.73936000)'47':171.66671000,(Pichia_kudriavzevii:224.69900000,(Candida_albicans:51.90226000,Candida_orthopsilosis:51.90226000)'68':172.79674000)'73':298.74919000)'72':151.74837000,Cunninghamella_bertholletiae:675.19656000)'84':0.00000000,Entomophaga_grylli:675.19656000)'87':923.24819000,((Saprolegnia_parasitica:1105.58726000,((Ichthyophthirius_multifiliis:226.62242000,Tetrahymena_thermophila:226.62242000)'83':813.70416000,Toxoplasma_gondii:1040.32658000)'82':65.26068000)'80':409.47270000,(Naegleria:1466.45600000,(Leishmania_donovani:265.93100000,Trypanosoma_cruzi:265.93100000)'79':1200.52500000)'94':48.60396000)'92':83.38479000)'78':2651.55525000);"
treeP <- read.tree(text = newick_path_tree)
plot(treeP, show.tip.label = TRUE)

hostCat <- allData_AICc_delta %>%
  left_join(keys, by = "Full_key")

category_data2 <- data.frame(
  tip_label = hostCat$Pathogen,
  category = hostCat$Pathogen_taxa
)

# Ensure categories are factors for consistent coloring
category_data2$category <- factor(category_data2$category)

# Extract tip labels and coordinates from the tree
tip_labels2 <- data.frame(tip_label = treeP$tip.label)
tip_coords2 <- as.data.frame(treeP$edge[treeP$edge[, 2] <= Ntip(treeP), 2])
names(tip_coords2) <- "node"
tip_coords2$tip_label <- treeP$tip.label

# Merge coordinates with category data
merged_data_pat <- tip_coords2 %>%
  left_join(category_data2, by = "tip_label") %>%
  filter_all(all_vars(!is.na(.)))

# Plot the tree
plot(treeP, show.tip.label = FALSE)
tiplabels(treeP$tip.label, frame = "none", adj = c(0.5, -0.5))

cbpalette3 <- rev(c("#DDCC77","#E69F00","#000000", "#888888","#6699CC","#332288"))

category_colors2 <- setNames(cbpalette3[1:length(unique(category_data2$category))], unique(category_data2$category))
category_data2$color <- category_colors2[category_data2$category]
plot_colored_tree(treeP, category_data2)

allData_AICc_delta_sum_P$Pathogen_taxa <- factor(allData_AICc_delta_sum_P$Pathogen_taxa, 
                                           levels = rev(c("Protozoan parasite",  "Fungi", 
                                                           "Gram-positive bacteria", "Gram-negative bacteria", "DNA virus", "RNA virus")))

cbpalette4 <- c("#000000","#888888","#6699CC","#332288","#E69F00","#DDCC77")

delta_pathogen <- ggplot(allData_AICc_delta_sum_P, aes(y = Pathogen_taxa, x = delta_AIC, fill = Pathogen_taxa)) +
  geom_point(size=5, shape=21, color = "grey25")+
  geom_errorbar(aes(xmin = delta_AIC - se, xmax = delta_AIC + se), width = 0.2, color ="grey25") +
  mytheme+
  labs(x = expression(Delta[AICc]), y = NULL) +
  scale_fill_manual(values = cbpalette4)+
  theme(legend.position = "none")

allData_AICc_delta_bars_P <- allData_AICc_delta %>%
  group_by(Pathogen_taxa) %>%
  summarise(MeanObs = mean(Num_Obs), MeanHost = mean(Num_Host))

allData_AICc_delta_bars_P$Pathogen_taxa <- factor(allData_AICc_delta_bars_P$Pathogen_taxa, 
                                                 levels = rev(c("Protozoan parasite",  "Fungi", 
                                                                "Gram-positive bacteria", "Gram-negative bacteria", "DNA virus", "RNA virus")))


delta_NUM_P <- ggplot(allData_AICc_delta_bars_P, aes(y = Pathogen_taxa, x = MeanObs, fill = Pathogen_taxa)) +
  geom_bar(stat = "identity", width = 0.25)+
  mytheme+
  labs(x = "Num. Observations", y = NULL) +
  scale_fill_manual(values = cbpalette4)+
  theme(legend.position = "none")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

tree_gg_P <- as.ggplot(~plot_colored_tree(treeP, category_data2))


# Create the delta_HOST plot with x-axis on top
delta_HOST_P <- ggplot(allData_AICc_delta_bars_P, aes(y = Pathogen_taxa, x = MeanHost, color = Pathogen_taxa)) +
  geom_bar(stat = "identity", width = 0.25, fill = "white") +
  mytheme +
  labs(x = "Average Num. Host", y = NULL) +
  scale_color_manual(values = cbpalette4) +
  scale_x_continuous(position = "top") +  # Place the x-axis on top
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x.top = element_line(colour = "black"),
        axis.ticks.x.top = element_line(colour = "black"),
        axis.title.x.top = element_text(margin = margin(t = 10)),
        axis.text.x.top = element_text(margin = margin(t = 10))) +  # Ensure axis text is set
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# Combine the plots with adjusted heights and widths
combined_plotP <- plot_grid(delta_pathogen, delta_NUM_P, delta_HOST_P, ncol = 3, align = "h", axis = "tb", rel_heights = c(1, 1, 1), rel_widths = c(2, 1, 1))

#Regression with host and pathogen distance  -give an arbitrary number to the virusses 

dist_matrixH <- cophenetic(tree)
dist_matrixP <- cophenetic(treeP)

# Perform PCoA
pcoa_resP <- ape::pcoa(dist_matrixP)

# Extract the first few principal coordinates
pcoa_coords <- pcoa_resP$vectors[, 1:2]  # Using the first two coordinates as an example
pcoa_dfP <- data.frame(Pathogen = rownames(pcoa_coords), pcoa_coords)

names(pcoa_df) <- c("Host", "PCoA_Axis1_Host", "PCoA_Axis2_Host")
names(pcoa_dfP) <- c("Pathogen", "PCoA_Axis1_Pathogen", "PCoA_Axis2_Pathogen")

# Data is from the main_script
modelData3 <- keys %>%
  full_join(pcoa_df, by = "Host") %>%
  full_join(pcoa_dfP, by = "Pathogen")

# Replace NA with a constant value, for example, a large value
constant_value1 <- max(pcoa_dfP$PCoA_Axis1_Pathogen) * 2
constant_value2 <- max(pcoa_dfP$PCoA_Axis2_Pathogen) * 2

distanceModelP <- allData_AICc_delta %>%
  full_join(modelData3, by = "Full_key") %>%
  mutate(is_virus = ifelse(Pathogen_taxa %in% c("DNA virus", "RNA virus"), 1, 0))

# Replace NA values in PCoA axes with constant values only if they correspond to viruses
constant_value1 <- max(distanceModelP$PCoA_Axis1_Pathogen, na.rm = TRUE) * 2
constant_value2 <- max(distanceModelP$PCoA_Axis2_Pathogen, na.rm = TRUE) * 2

distanceModelP$PCoA_Axis1_Pathogen[is.na(distanceModelP$PCoA_Axis1_Pathogen) & distanceModelP$is_virus == 1] <- constant_value1
distanceModelP$PCoA_Axis2_Pathogen[is.na(distanceModelP$PCoA_Axis2_Pathogen) & distanceModelP$is_virus == 1] <- constant_value2

distanceModel_final <- distanceModelP %>%
  filter_all(all_vars(!is.na(.))) %>%
  filter_all(all_vars(!is.nan(.))) %>%
  filter_all(all_vars(!is.infinite(.)))

distance_model <- lm(DeltaSub ~ scale(Num_Host) + scale(Num_Obs) + scale(PCoA_Axis1_Host) + scale(PCoA_Axis2_Host) + scale(PCoA_Axis1_Pathogen) + scale(PCoA_Axis2_Pathogen), data = distanceModel_final)
summary(distance_model)

lmCoef_distance <- sjPlot::plot_model(distance_model, type = "est", show.values = T, show.p = TRUE)+
  scale_color_manual(values = cbpalette)+
  labs(y = expression(Delta[AICc]~ estimates), title = NULL)+
  mytheme

anova(distance_model)

# I want to test for other models 

allData_AICc_delta2 <- allData_WAIC_filled %>%
  filter(Model %in% c("weibull", "exponential")) %>%
  spread(Model, AICc) %>%
  group_by(Full_key, Key, Host_taxa, Pathogen_taxa) %>%
  summarise(exponential = mean(exponential, na.rm = TRUE), weibull = mean(weibull, na.rm =T), 
            Num_Obs = mean(Num_Obs), Num_Host = mean(Max_num_host_filled)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(DeltaSub = exponential - weibull)

allData_AICc_delta3 <- allData_WAIC_filled %>%
  filter(Model %in% c("loglogistic", "exponential")) %>%
  spread(Model, AICc) %>%
  group_by(Full_key, Key, Host_taxa, Pathogen_taxa) %>%
  summarise(exponential = mean(exponential, na.rm = TRUE), loglogistic = mean(loglogistic, na.rm =T), 
            Num_Obs = mean(Num_Obs), Num_Host = mean(Max_num_host_filled)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(DeltaSub = exponential - loglogistic)

allData_AICc_delta4 <- allData_WAIC_filled %>%
  filter(Model %in% c("generalizedgamma", "exponential")) %>%
  spread(Model, AICc) %>%
  group_by(Full_key, Key, Host_taxa, Pathogen_taxa) %>%
  summarise(exponential = mean(exponential, na.rm = TRUE), generalizedgamma = mean(generalizedgamma, na.rm =T), 
            Num_Obs = mean(Num_Obs), Num_Host = mean(Max_num_host_filled)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(DeltaSub = exponential - generalizedgamma)

#Weibull
distanceModelP2 <- allData_AICc_delta2 %>%
  full_join(modelData3, by = "Full_key") %>%
  mutate(is_virus = ifelse(Pathogen_taxa %in% c("DNA virus", "RNA virus"), 1, 0))
distanceModelP2$PCoA_Axis1_Pathogen[is.na(distanceModelP2$PCoA_Axis1_Pathogen) & distanceModelP2$is_virus == 1] <- constant_value1
distanceModelP2$PCoA_Axis2_Pathogen[is.na(distanceModelP2$PCoA_Axis2_Pathogen) & distanceModelP2$is_virus == 1] <- constant_value2
distanceModel_final2 <- distanceModelP2 %>%
  filter_all(all_vars(!is.na(.))) %>%
  filter_all(all_vars(!is.nan(.))) %>%
  filter_all(all_vars(!is.infinite(.)))

distance_model2 <- lm(DeltaSub ~ scale(Num_Host) + scale(Num_Obs) + scale(PCoA_Axis1_Host) + scale(PCoA_Axis2_Host) + scale(PCoA_Axis1_Pathogen) + scale(PCoA_Axis2_Pathogen), data = distanceModel_final2)
summary(distance_model2)
lmCoef_distance2 <- sjPlot::plot_model(distance_model2, type = "est", show.values = T, show.p = TRUE)+
  scale_color_manual(values = cbpalette)+
  labs(y = expression(Delta[AICc]~ estimates), title = NULL)+
  mytheme

#Loglogistic
distanceModelP3 <- allData_AICc_delta3 %>%
  full_join(modelData3, by = "Full_key") %>%
  mutate(is_virus = ifelse(Pathogen_taxa %in% c("DNA virus", "RNA virus"), 1, 0))
distanceModelP3$PCoA_Axis1_Pathogen[is.na(distanceModelP3$PCoA_Axis1_Pathogen) & distanceModelP3$is_virus == 1] <- constant_value1
distanceModelP3$PCoA_Axis2_Pathogen[is.na(distanceModelP3$PCoA_Axis2_Pathogen) & distanceModelP3$is_virus == 1] <- constant_value2
distanceModel_final3 <- distanceModelP3 %>%
  filter_all(all_vars(!is.na(.))) %>%
  filter_all(all_vars(!is.nan(.))) %>%
  filter_all(all_vars(!is.infinite(.)))

distance_model3 <- lm(DeltaSub ~ scale(Num_Host) + scale(Num_Obs) + scale(PCoA_Axis1_Host) + scale(PCoA_Axis2_Host) + scale(PCoA_Axis1_Pathogen) + scale(PCoA_Axis2_Pathogen), data = distanceModel_final3)
summary(distance_model3)
lmCoef_distance3 <- sjPlot::plot_model(distance_model3, type = "est", show.values = T, show.p = TRUE)+
  scale_color_manual(values = cbpalette)+
  labs(y = expression(Delta[AICc]~ estimates), title = NULL)+
  mytheme

#Generalized gamma 
distanceModelP4 <- allData_AICc_delta4 %>%
  full_join(modelData3, by = "Full_key") %>%
  mutate(is_virus = ifelse(Pathogen_taxa %in% c("DNA virus", "RNA virus"), 1, 0))
distanceModelP4$PCoA_Axis1_Pathogen[is.na(distanceModelP4$PCoA_Axis1_Pathogen) & distanceModelP4$is_virus == 1] <- constant_value1
distanceModelP4$PCoA_Axis2_Pathogen[is.na(distanceModelP4$PCoA_Axis2_Pathogen) & distanceModelP4$is_virus == 1] <- constant_value2
distanceModel_final4 <- distanceModelP4 %>%
  filter_all(all_vars(!is.na(.))) %>%
  filter_all(all_vars(!is.nan(.))) %>%
  filter_all(all_vars(!is.infinite(.)))

distance_model4 <- lm(DeltaSub ~ scale(Num_Host) + scale(Num_Obs) + scale(PCoA_Axis1_Host) + scale(PCoA_Axis2_Host) + scale(PCoA_Axis1_Pathogen) + scale(PCoA_Axis2_Pathogen), data = distanceModel_final4)
summary(distance_model4)
lmCoef_distance4 <- sjPlot::plot_model(distance_model4, type = "est", show.values = T, show.p = TRUE)+
  scale_color_manual(values = cbpalette)+
  labs(y = expression(Delta[AICc]~ estimates), title = NULL)+
  mytheme

# Deviation coding for other models 
#allData_AICc_delta2
#allData_AICc_delta3
#allData_AICc_delta4

# Set the contrast options for categorical variables
options(contrasts = c("contr.sum", "contr.poly"))

# Fit the linear model with deviation coding
model_deviation2 <- lm(DeltaSub ~ Num_Obs + Num_Host + Host_taxa + Pathogen_taxa, data = allData_AICc_delta2)
summary(model_deviation2)
anova(model_deviation2)

# Plot the model estimates with actual level names
plot_model(model_deviation2, type = "est", show.values = TRUE, value.offset = 0.4, axis.labels = c(rev(levels(factor(allData_AICc_delta2$Pathogen_taxa))), rev(levels(factor(allData_AICc_delta2$Host_taxa))), "Num_Obs", "Num_Host")) +
  labs(title = "Model with Deviation Coding - for Weibull",
       y = "Effect Sizes (ΔAICc estimates)",
       x = "")

# Fit the linear model with deviation coding
model_deviation3 <- lm(DeltaSub ~ Num_Obs + Num_Host + Host_taxa + Pathogen_taxa, data = allData_AICc_delta3)
summary(model_deviation3)
anova(model_deviation3)

# Plot the model estimates with actual level names
plot_model(model_deviation3, type = "est", show.values = TRUE, value.offset = 0.4, axis.labels = c(rev(levels(allData_AICc_delta3$Pathogen_taxa)), rev(levels(allData_AICc_delta3$Host_taxa)), "Num_Obs", "Num_Host")) +
  labs(title = "Model with Deviation Coding - for Log-logistic",
       y = "Effect Sizes (ΔAICc estimates)",
       x = "")

# Plot the model estimates with actual level names
plot_model(model_deviation2, type = "est", show.values = TRUE, value.offset = 0.4, axis.labels = c(rev(levels(allData_AICc_delta3$Pathogen_taxa)), rev(levels(allData_AICc_delta3$Host_taxa)), "Num_Obs", "Num_Host")) +
  labs(title = "Model with Deviation Coding - for Weibull",
       y = "Effect Sizes (ΔAICc estimates)",
       x = "")

# Fit the linear model with deviation coding
model_deviation4 <- lm(DeltaSub ~ Num_Obs + Num_Host + Host_taxa + Pathogen_taxa, data = allData_AICc_delta4)
summary(model_deviation4)
anova(model_deviation4)

# Plot the model estimates with actual level names
plot_model(model_deviation4, type = "est", show.values = TRUE, value.offset = 0.4, axis.labels = c(rev(levels(allData_AICc_delta4$Pathogen_taxa)), rev(levels(allData_AICc_delta4$Host_taxa)), "Num_Obs", "Num_Host")) +
  labs(title = "Model with Deviation Coding - for Log-logistic",
       y = "Effect Sizes (ΔAICc estimates)",
       x = "")
