# Meta survival data analysis
# Author: Canan Karakoc
# Last update: August 7, 2024

# Install/load libraries
library(tidyverse)
library(minpack.lm)
library(survival)
library(stats)
library(loo)
library(ggforce)
library(ggpubr)
library(MASS)
library(sjPlot)
library(patchwork)
library(cowplot)
library(ggplotify)
library(Matrix)
library(MCMCglmm)
library(phytools)
library(ape)
library(car)
library(caper)
library(broom)
library(dplyr) # might be loaded last 

set.seed(1234)
setwd("~/GitHub/metaAnalysis")

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

### FIGURE 2A ###
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

#################################################

###FIGURE S2 & S3 ###
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

#########################################################################
### FIGURE S4 ###
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

####################################################################

### FIGURE 1 & FIGURE S1 ###
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
### FIGURE 3 A&B ###

# Delta AIC
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

###################################################################
### FIGURE 2B & C ###

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
  labs(y = expression(AICc[C] - AICc[j]), x = NULL)+
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
  mutate(Constant = gompertz_survival - exponential,
         Weibull  = gompertz_survival - weibull, 
         Loglogistic = gompertz_survival - loglogistic, 
         Gen.gamma = gompertz_survival -generalizedgamma) %>%
  pivot_longer(Constant:Gen.gamma, names_to = "comparison", values_to = "DeltaAIC")%>%
  filter_all(all_vars(!is.na(.))) %>%
  filter_all(all_vars(!is.nan(.))) %>%
  filter_all(all_vars(!is.infinite(.)))

delta_all_g2 <- ggplot(allModels_delta_g, aes(y = DeltaAIC, x = comparison)) +
  geom_hline (yintercept = 0, linetype = "dashed")+
  geom_jitter(size = 3, shape = 21, color = "grey50")+
  geom_boxplot(fill = "white", , alpha  = 0.7)+
  mytheme+
  labs(y = expression(AICc[G] - AICc[j]), x = NULL)+
  scale_x_discrete(labels = c("Constant", "Weibull", "Log-logistic", "Gen.-gamma"))+
  scale_y_continuous(limits = c(-30, 30))

ggsave("figures/delta_all_g2-limited.pdf", plot = delta_all_g2, width = 6, height = 4, units = "in", dpi = 300)
ggsave("figures/delta_all_g2-limited.png", plot = delta_all_g2, width = 6, height = 4, units = "in", dpi = 300)

# Combine the plots and add labels
op <- op + theme(axis.text.x = element_text(angle = 45, hjust = 1))
delta_all_g2 <- delta_all_g2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot <- density2 + op + delta_all_g2 +
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

# Add observations to the model
# Phylogenetic distance

allData_AICc_delta <- allData_WAIC_filled %>%
  filter(Model %in% c("gompertz_survival", "exponential")) %>%
  spread(Model, AICc) %>%
  group_by(Full_key, Key, Host_taxa, Pathogen_taxa, Data) %>%
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
complete_datalm <- allData_AICc_delta %>%
  filter(!is.na(DeltaSub) & !is.nan(DeltaSub) & is.finite(DeltaSub))


lm_model_AIC3 <- lm(DeltaSub ~ Data + Num_Obs + Num_Host + Host_taxa + Pathogen_taxa, 
                       data = complete_datalm)
summary(lm_model_AIC3)
anova(lm_model_AIC3)

vif(lm_model_AIC3)

# Plot fixed effects
sjPlot::plot_model(lm_model_AIC3, type = "est", show.values = T, show.p = TRUE)

# Plot for the manuscript 

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
model_deviation <- lm(DeltaSub ~ Data + Num_Obs + Num_Host + Host_taxa + Pathogen_taxa, data = complete_datalm)
summary(model_deviation)
anova(model_deviation)

# Plot the model estimates with actual level names
plot_model(model, type = "est", show.values = TRUE, value.offset = 0.4) +
  labs(title = "Model with Deviation Coding",
       y = "Effect Sizes (ΔAICc estimates)",
       x = "")

contrasts(complete_datalm$Host_taxa)
contrasts(complete_datalm$Pathogen_taxa)

complete_datalm$Host_taxa <- factor(complete_datalm$Host_taxa, 
                                       levels = c("Seedlings", "Nematodes", "Other invertebrates", "Drosophila sp.", "Moth larvae", "Other insects", "Avian", "Mice", "Other mammals", "Fish"))

complete_datalm$Pathogen_taxa <- factor(complete_datalm$Pathogen_taxa, 
                                           levels = c("Protozoan parasite", "Fungi", "Gram-positive bacteria", "Gram-negative bacteria", "DNA virus", "RNA virus"))

# Set contrasts to deviation coding
contrasts(complete_datalm$Host_taxa) <- contr.sum(length(levels(complete_datalm$Host_taxa)))
contrasts(complete_datalm$Pathogen_taxa) <- contr.sum(length(levels(complete_datalm$Pathogen_taxa)))

# Fit the model
model <- lm(DeltaSub ~ Data + Num_Obs + Num_Host + Host_taxa + Pathogen_taxa, data = complete_datalm)
summary(model)
anova(model)
#Extract model estimates
model_tidy <- tidy(model)

# Calculate the reference level effect size for Host_taxa and Pathogen_taxa
host_taxa_ref <- -sum(model_tidy %>% filter(str_detect(term, "Host_taxa")) %>% pull(estimate))
pathogen_taxa_ref <- -sum(model_tidy %>% filter(str_detect(term, "Pathogen_taxa")) %>% pull(estimate))
data_ref <- -sum(model_tidy %>% filter(str_detect(term, "Data")) %>% pull(estimate))

# Create data frame for reference levels
ref_levels <- tibble(
  term = c("Host_taxaRef", "Pathogen_taxaRef", "Dataref"),
  estimate = c(host_taxa_ref, pathogen_taxa_ref, data_ref),
  std.error = c(NA, NA, NA),  # Standard errors are not available for manually calculated reference levels
  statistic = c(NA, NA, NA),
  p.value = c(NA, NA, NA)
)

# Combine with original model estimates
model_tidy <- bind_rows(model_tidy, ref_levels)

# Manually replace term names using case_when
model_tidy <- model_tidy %>%
  mutate(term = case_when(
    term == "Data1" ~ "Raw data",
    term == "Dataref" ~ "Probability data",
    term == "Num_Obs" ~ "Num. obs.",
    term == "Num_Host" ~ "Num. host",
    term == "Host_taxa1" ~ "Nematodes",
    term == "Host_taxa2" ~ "Other inverteb.",
    term == "Host_taxa3" ~ "Drosophila sp.",
    term == "Host_taxa4" ~ "Moth larvae",
    term == "Host_taxa5" ~ "Other insects",
    term == "Host_taxa6" ~ "Avian",
    term == "Host_taxa7" ~ "Mice",
    term == "Host_taxa8" ~ "Other mammals",
    term == "Host_taxa9" ~ "Fish",
    term == "Host_taxaRef" ~ "Seedlings",
    term == "Pathogen_taxa1" ~ "Fungi",
    term == "Pathogen_taxa2" ~ "Gram-(+) bacteria",
    term == "Pathogen_taxa3" ~ "Gram-(-) bacteria",
    term == "Pathogen_taxa4" ~ "DNA virus",
    term == "Pathogen_taxa5" ~ "RNA virus",
    term == "Pathogen_taxaRef" ~ "Protozoan parasite",
    TRUE ~ term  # keep the original term if no match
  ))

# Reorder terms based on factor levels
term_order <- rev(c("Seedlings", "Nematodes", "Other inverteb.", "Drosophila sp.", "Moth larvae", "Other insects", "Avian", "Mice", "Other mammals", "Fish", 
                "Protozoan parasite", "Fungi", "Gram-(+) bacteria", "Gram-(-) bacteria", "DNA virus", "RNA virus", 
                "Num. host", "Num. obs.", "Probability data", "Raw data"))

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
       x = expression(Delta[AICc]~ estimates),
       y = "") +
  scale_color_manual(values = c("Positive" = "#0072B2", "Negative" = "#D55E00")) +
  geom_text(aes(label = round(estimate, 2)), vjust = -0.3, size = 4) +
  mytheme+
  theme(legend.position = "none")

############################################################################

### FIGURE 3 & FIGURE S5 ###
### PHYLO ANALYSIS ###
# Distance/tree 

newick_tree <- "(((Pisaster_ochraceus:634.70000000,((((Oncorhynchus_mykiss:15.06984000,Oncorhynchus_tshawytscha:15.06984000)'14':190.03834000,(Micropterus_salmoides:112.23406000,(Paralichthys_olivaceus:103.76124000,(Oreochromis_niloticus:91.74830000,Poecilia_reticulata:91.74830000)'13':12.01294000)'25':8.47282000)'37':92.87412000)'36':18.90948000,(Silurus_glanis:142.11161000,Danio_rerio:142.11161000)'35':81.90605000)'34':204.98234000,((Bos_taurus:94.00000000,((Oryctolagus_cuniculus:78.97061000,(Mus_musculus:69.91125000,Cavia_porcellus:69.91125000)'43':9.05936000)'33':8.22939000,(Macaca_mulatta:3.31741000,Macaca_fascicularis:3.31741000)'51':83.88259000)'50':6.80000000)'49':224.95000000,Gallus_gallus:318.95000000)'57':110.05000000)'60':205.70000000)'56':72.90436000,((((Melanoplus_borealis:383.15000000,((Apis_mellifera_macedonica:343.05330000,((((((Drosophila_flavomontana:29.59186000,((Drosophila_prosaltans:21.56067000,Drosophila_tropicalis:21.56067000)'48':5.99910000,((Drosophila_ananassae:25.57846000,(Drosophila_suzukii:23.00787000,(((Drosophila_simulans:4.62390000,Drosophila_mauritiana:4.62390000)'63':0.00000000,Drosophila_melanogaster:4.62390000)'47':6.94583000,Drosophila_santomea:11.56973000)'68':11.43814000)'73':2.57059000)'72':0.62224000,Drosophila_subobscura:26.20070000)'84':1.35907000)'87':2.03209000)'83':14.22660000,Drosophila_euronotus:43.81846000)'82':0.62949000,Zaprionus_taronus:44.44795000)'80':20.67730000,Penicillidia_conspicua:65.12525000)'79':175.79475000,Anopheles_arabiensis:240.92000000)'94':58.47011000,Galleria_mellonella:299.39011000)'92':44.60989000)'78':17.53224000,Meccus_pallidipennis:361.53224000)'98':21.61776000)'97':121.40000000,Daphnia_magna:504.55000000)'77':10.82000000,(Penaeus_vannamei:364.66874000,(Macrobrachium_nipponense:331.97658000,Procambarus_clarkii:331.97658000)'71':32.69216000)'67':150.70126000)'66':56.74064000,Caenorhabditis_elegans:572.11064000)'46':135.49372000)'32':890.84039000,(Arctostaphylos_glauca:125.11960000,Cucurbita_pepo:125.11960000)'106':1473.32515000);"
tree <- read.tree(text = newick_tree)
plot(tree, show.tip.label = TRUE)

dist_matrix <- cophenetic(tree)

keys <- read.csv("~/GitHub/metaAnalysis/Phylogeny_list.csv",
                 sep = ",", header = T, na.strings = ""
)

category_data <- data.frame(
  tip_label = keys$Host,
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

cbpalette2 <- c("#6699CC","#888888","#000000","#332288","#E69F00","#F0E442","#DDCC77","#882255", 
                "#CC79A7","#44AA99")

allData_AICc_delta_sum$Host_taxa <- factor(allData_AICc_delta_sum$Host_taxa, 
                                            levels = rev(c("Seedlings",  
                                                           "Nematodes", "Other invertebrates", "Drosophila sp.", "Moth larvae", "Other insects",
                                                           "Avian", "Mice", "Other mammals", "Fish")))

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

tree_gg <- as.ggplot(~plot_colored_tree(tree, category_data))

# Combine the plots with adjusted heights and widths
combined_plot <- plot_grid(delta_host, delta_NUM, ncol = 2, align = "h", axis = "tb", rel_heights = c(1, 1), rel_widths = c(2, 1))

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

#############################################################################

### FIGURE 4 ###
# Distance matrix as a random effect 

ultrametric_tree <- chronos(tree)

# Get the species names in your data and tree
data_species <- unique(keys$Host)
tree_species <- ultrametric_tree$tip.label

# Find species present in data but not in tree
missing_species <- setdiff(data_species, tree_species)
if(length(missing_species) > 0) {
  warning("The following species are in the data but not in the tree: ", paste(missing_species, collapse = ", "))
}

phylo_cov <- inverseA(ultrametric_tree, nodes = "TIPS")$Ainv
phylo_cov_matrix <- as.matrix(phylo_cov)

CovData <- allData_AICc_delta %>%
  full_join(keys, by = "Full_key") %>%
  mutate(species = Host)

phylo_cov_df <- as.data.frame(phylo_cov_matrix)
phylo_cov_df$species <- rownames(phylo_cov_df)

merged_data <- merge(phylo_cov_df, CovData, by = "species", all.x = TRUE)
merged_data2 <- cbind.data.frame(category = merged_data$Host_taxa, 
                                merged_data[,2:42])
  
# Calculate average values for each category
fishd <- merged_data2 %>% filter(category == "Fish")
fishdMean <- colMeans(fishd[,-1], na.rm = TRUE)

drod <- merged_data2 %>% filter(category == "Drosophila sp.")
drodMean <- colMeans(drod[,-1], na.rm = TRUE)

insd <- merged_data2 %>% filter(category == "Other insects")
insdMean <- colMeans(insd[,-1], na.rm = TRUE)

# Create rows for missing species with average values
missing_species <- data.frame(
  rbind(fishdMean, drodMean, insdMean),
  species = c("Acanthopagrus_schlegeli", "Hirtodrosophila_duncani", "Camnula_pellucida")
  
)
rownames(missing_species) <- missing_species$species
colnames(missing_species) <- colnames(phylo_cov_df)

# Bind the new species rows
phylo_cov_df <- rbind(phylo_cov_df, missing_species)

# Add columns for missing species with average values
new_cols <- t(missing_species[,-42])
colnames(new_cols) <- missing_species$species
new_cols <- as.data.frame(new_cols)

# Initialize columns for the new species with NA values
for (species_name in missing_species$species) {
  phylo_cov_df[species_name] <- NA
}

# Fill in the new columns with the mean values
for (i in 1:nrow(missing_species)) {
  species_name <- missing_species$species[i]
  mean_values <- unlist(missing_species[i, -ncol(missing_species)])
  
  # Add the mean values to the new columns
  phylo_cov_df[, species_name] <- c(mean_values, rep(NA, nrow(phylo_cov_df) - length(mean_values)))
}

# Ensure the species column is at the end for consistent binding
phylo_cov_df$species <- rownames(phylo_cov_df)

# Now we need to fill in the lower triangular part to maintain symmetry
for (i in 1:nrow(missing_species)) {
  species_name <- missing_species$species[i]
  phylo_cov_df[species_name, species_name] <- 0
  for (j in 1:ncol(new_cols)) {
    phylo_cov_df[species_name, colnames(new_cols)[j]] <- new_cols[j, i]
    phylo_cov_df[colnames(new_cols)[j], species_name] <- new_cols[j, i]
  }
}

# Ensure the species column is at the end for consistent binding
phylo_cov_df$species <- rownames(phylo_cov_df)
phylo_cov_df <- phylo_cov_df[-42]

phylo_cov_matrix <- as.matrix(phylo_cov_df)
rownames(phylo_cov_matrix) <- rownames(phylo_cov_df)
colnames(phylo_cov_matrix) <- rownames(phylo_cov_df)

# Extract species levels from data
data_species <- unique(CovData$species)

# Extract species levels from phylo_cov_matrix
tree_species <- rownames(phylo_cov_matrix)

# Find species in data not present in phylo_cov_matrix
missing_in_tree <- setdiff(data_species, tree_species)
if(length(missing_in_tree) > 0) {
  warning("The following species are in the data but not in the tree matrix: ", paste(missing_in_tree, collapse = ", "))
}

# Find species in phylo_cov_matrix not present in data
missing_in_data <- setdiff(tree_species, data_species)
if(length(missing_in_data) > 0) {
  warning("The following species are in the tree matrix but not in the data: ", paste(missing_in_data, collapse = ", "))
}

data <- CovData[CovData$species %in% tree_species, ]

# Check for missing values in the data for Mus_musculus
mus_data1 <- data[data$Full_key == "2DRM7FWM.csv", ] 
mus_data2 <- data[data$Full_key == "5NVSCRNB_1.csv", ]
mus_data3 <- data[data$Full_key == "5NVSCRNB_2.csv", ]
missing_cols1 <- colnames(mus_data1)[colSums(is.na(mus_data1)) > 0]
missing_cols2 <- colnames(mus_data2)[colSums(is.na(mus_data2)) > 0]
missing_cols3 <- colnames(mus_data3)[colSums(is.na(mus_data3)) > 0]

# Impute missing values for Mus_musculus 
for(col in missing_cols1) {
  mus_data1[, col] <- ifelse(is.na(mus_data1[, col]), mean(data[[col]], na.rm = TRUE), mus_data1[, col])
}

for(col in missing_cols2) {
  mus_data2[, col] <- ifelse(is.na(mus_data2[, col]), mean(data[[col]], na.rm = TRUE), mus_data2[, col])
}

for(col in missing_cols3) {
  mus_data3[, col] <- ifelse(is.na(mus_data3[, col]), mean(data[[col]], na.rm = TRUE), mus_data3[, col])
}

# Replace the original row for Mus_musculus with the imputed data
data[data$Full_key == "2DRM7FWM.csv", ] <- mus_data1
data[data$Full_key == "5NVSCRNB_1.csv", ] <- mus_data2
data[data$Full_key == "5NVSCRNB_2.csv", ] <- mus_data3

complete_data <- data
complete_data$species <- complete_data$Host

# Ensure the species order in the covariance matrix matches the data
phylo_cov_matrix <- phylo_cov_matrix[data_species, data_species]

# Verify alignment after reordering
aligned <- all(rownames(phylo_cov_matrix) == data_species)

# Recheck species alignment in data and matrix
data_species <- unique(data$species)
matrix_species <- rownames(phylo_cov_matrix)

# Reorder matrix to match the order in data
phylo_cov_matrix <- phylo_cov_matrix[data_species, data_species]

# Remove the groups attribute from complete_data
attr(complete_data, "groups") <- NULL

# Verify that the attribute has been removed
str(complete_data)

# Ensure all species we in complete_data are present in phylo_cov_matrix
species_in_data <- unique(complete_data$species)
species_in_matrix <- rownames(phylo_cov_matrix)

missing_in_matrix <- setdiff(species_in_data, species_in_matrix)
missing_in_data <- setdiff(species_in_matrix, species_in_data)

# Verify the dimensions of the matrix
print(dim(phylo_cov_matrix))

# Ensure species alignment in data and matrix
species_in_matrix <- rownames(phylo_cov_matrix)
species_in_data <- unique(complete_data$species)

# Ensure species levels in filtered_data match the filtered phylogenetic matrix
complete_data <- complete_data[complete_data$species %in% rownames(phylo_cov_matrix), ]

# Ensure species factor levels match
complete_data$species <- factor(complete_data$species, levels = rownames(complete_data))

# Check for any NA values in the matrix
if (any(is.na(phylo_cov_matrix))) {
  stop("There are NA values in the phylogenetic covariance matrix.")
}

complete_data$species <- complete_data$Host

# Check if the matrix is symmetric
is_symmetric <- all(phylo_cov_matrix == t(phylo_cov_matrix))

# Force the matrix to be symmetric
phylo_cov_matrix <- (phylo_cov_matrix + t(phylo_cov_matrix)) / 2

# Verify the matrix is now symmetric
is_symmetric <- all(phylo_cov_matrix == t(phylo_cov_matrix))

# Function to make matrix positive semi-definite
make_positive_semi_definite <- function(mat) {
  eig <- eigen(mat)
  eig$values[eig$values < 0] <- 0
  mat <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  return(mat)
}

phylo_cov_matrix <- make_positive_semi_definite(phylo_cov_matrix)
rownames(phylo_cov_matrix) <- rownames(phylo_cov_df)
colnames(phylo_cov_matrix) <- rownames(phylo_cov_df)

# Ensure there are no NAs or mismatched entries
phylo_cov_matrix[is.na(phylo_cov_matrix)] <- 0

data <- as.data.frame(data)
phylo_cov_matrix <- as(phylo_cov_matrix, "dgCMatrix")

phylo_cov_matrix <- phylo_cov_matrix / max(abs(phylo_cov_matrix))
diag(phylo_cov_matrix) <- diag(phylo_cov_matrix) + 1e-5

priors <- list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1)))

complete_data <- as.data.frame(complete_data)
model <- MCMCglmm(DeltaSub ~ Num_Obs,
                  random = ~ species,
                  ginverse = list(species = phylo_cov_matrix),
                  data = complete_data,
                  family = "gaussian",
                  prior = priors,
                  nitt = 13000,
                  burnin = 3000,
                  thin = 10)

summary(model)

standard_model <- lm(DeltaSub ~ Num_Obs, data = data)

# R-squared for the standard linear regression model
r_squared_lm <- summary(standard_model)$r.squared

# Extract fixed effects variance
fixed_effects_matrix <- model$Sol %*% t(model$X)
fixed_var <- apply(fixed_effects_matrix, 1, var)

# Extract random effects and residual variances
random_var <- apply(model$VCV, 2, mean)
residual_var <- mean(model$VCV[,"units"])

# Sum of variances of the random effects
total_random_var <- sum(random_var)

# Calculate Marginal and Conditional R-squared
marginal_r2_phylo <- mean(fixed_var) / (mean(fixed_var) + total_random_var + residual_var)
conditional_r2_phylo <- (mean(fixed_var) + total_random_var) / (mean(fixed_var) + total_random_var + residual_var)

# Create a data frame for plotting
r_squared_df <- data.frame(
  Model = c("Standard LM", "Phylogenetic model (Marginal)", "Phylogenetic model (Conditional)"),
  R_Squared = c(r_squared_lm, marginal_r2_phylo, conditional_r2_phylo)
)

# Extract the slopes 
phylo_slope_Num_Obs <- summary(model)$solutions["Num_Obs", "post.mean"]
standard_slope_Num_Obs <- coef(standard_model)["Num_Obs"]

#Create a data frame for plotting the regression lines
plot_data <- data.frame(
  Num_Obs = seq(min(data$Num_Obs), max(data$Num_Obs), length.out = 100)
)

plot_data$phylo_fit_Num_Obs <- phylo_slope_Num_Obs * plot_data$Num_Obs
plot_data$standard_fit_Num_Obs <- standard_slope_Num_Obs * plot_data$Num_Obs

# Plot the data and regression lines for Num_Obs
plot1 <- ggplot(data, aes(x = Num_Obs, y = DeltaSub)) +
  geom_point(shape = 21, size = 3) +
  geom_line(data = plot_data, aes(y = phylo_fit_Num_Obs, color = "Phylogenetic model"), linetype = "dashed", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data, aes(y = standard_fit_Num_Obs, color = "Standard model"), linetype = "solid", size = 1, show.legend = TRUE) +
  labs(
    x = "Number of observations",
    y = expression(Delta[AICc]),
    color = "Model"
  ) +
  scale_color_manual(values = c("Phylogenetic model" = "#0072B2", "Standard model" = "#D55E00")) +
  mytheme +
  scale_x_continuous(limits = c(5, 50)) +
  scale_y_continuous(limits = c(0, 300))+
  theme(legend.position = c(0.4, 0.9))

# Filter Bee data with 500 hosts 
filtered_data <- as.data.frame(complete_data[!(complete_data$Num_Host == 500), ])

# Remove the species "Apis_mellifera_macedonica" from the phylogenetic matrix
filtered_phylo_cov_matrix <- phylo_cov_matrix[rownames(phylo_cov_matrix) != "Apis_mellifera_macedonica", colnames(phylo_cov_matrix) != "Apis_mellifera_macedonica"]

# Fit the phylogenetic model with the filtered data
model2 <- MCMCglmm(DeltaSub ~ Num_Host,
                   random = ~ species,
                   ginverse = list(species = filtered_phylo_cov_matrix),
                   data = filtered_data,
                   family = "gaussian",
                   prior = priors,
                   nitt = 13000,
                   burnin = 3000,
                   thin = 10)

summary(model2)

standard_model2 <- lm(DeltaSub ~ Num_Host, data = filtered_data)

# R-squared for the standard linear regression model
r_squared_lm2 <- summary(standard_model2)$r.squared

# Extract fixed effects variance
fixed_effects_matrix2 <- model2$Sol %*% t(model2$X)
fixed_var2 <- apply(fixed_effects_matrix2, 1, var)

# Extract random effects and residual variances
random_var2 <- apply(model2$VCV, 2, mean)
residual_var2 <- mean(model2$VCV[,"units"])

# Sum of variances of the random effects
total_random_var2 <- sum(random_var2)

# Calculate Marginal and Conditional R-squared
marginal_r2_phylo2 <- mean(fixed_var2) / (mean(fixed_var2) + total_random_var2 + residual_var2)
conditional_r2_phylo2 <- (mean(fixed_var2) + total_random_var2) / (mean(fixed_var2) + total_random_var2 + residual_var2)

# Create a data frame for plotting
r_squared_df2 <- data.frame(
  Model = c("Standard LM", "Phylogenetic model (Marginal)", "Phylogenetic model (Conditional)"),
  R_Squared = c(r_squared_lm2, marginal_r2_phylo2, conditional_r2_phylo2)
)

# Extract the slopes 
phylo_slope_Num_Host2 <- summary(model2)$solutions["Num_Host", "post.mean"]
standard_slope_Num_Host2 <- coef(standard_model2)["Num_Host"]

#Create a data frame for plotting the regression lines
plot_data2 <- data.frame(
  Num_Host = seq(min(filtered_data$Num_Host), max(filtered_data$Num_Host), length.out = 100)
)

plot_data2$phylo_fit_Num_Host <- phylo_slope_Num_Host2 * plot_data2$Num_Host
plot_data2$standard_fit_Num_Host <- standard_slope_Num_Host2 * plot_data2$Num_Host

# Plot the data and regression lines for Num_Obs
plot2 <- ggplot(filtered_data, aes(x = Num_Host, y = DeltaSub)) +
  geom_point(shape = 21, size = 3) +
  geom_line(data = plot_data2, aes(y = phylo_fit_Num_Host), color = "#0072B2", linetype = "dashed", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data2, aes(y = standard_fit_Num_Host), color = "#D55E00", linetype = "solid", size = 1, show.legend = TRUE) +
  labs(
    x = "Number of host",
    y = expression(Delta[AICc]),
    color = "Model"
  ) +
  scale_color_manual(values = c("Phylogenetic model" = "#0072B2", "Standard model" = "#D55E00")) +
  mytheme +
  scale_x_continuous(limits = c(0, 100))+
  scale_y_continuous(limits = c(0, 250))

combined_plot_Mod <- plot1 + plot2 +
  plot_annotation(tag_levels = 'A')

#####################################################
# PATHOGEN
# Create a distance matrix from the tree

# Check for any missing nodes
if (any(is.na(treeP$node.label))) {
  # Assign unique labels to any missing nodes
  treeP$node.label[is.na(treeP$node.label)] <- paste0("Node_", seq(sum(is.na(treeP$node.label))))
}

# Ensure the tree is ultrametric
if (!is.ultrametric(treeP)) {
  ultrametric_treeP <- chronos(treeP)
} else {
  ultrametric_treeP <- treeP
}

# Ensure zero lengths are replaced in the pathogen distance matrix
pathogen_dist_matrix <- cophenetic(ultrametric_treeP)
if (any(pathogen_dist_matrix == 0, na.rm = TRUE)) {
  pathogen_dist_matrix[pathogen_dist_matrix == 0] <- 1e-10
}

# Check for any NA values and handle them
if (any(is.na(pathogen_dist_matrix))) {
  pathogen_dist_matrix[is.na(pathogen_dist_matrix)] <- max(pathogen_dist_matrix, na.rm = TRUE) + 1
}

# Convert the distance matrix to a covariance matrix
pathogen_cov_matrix <- max(pathogen_dist_matrix) - pathogen_dist_matrix

# Make the covariance matrix positive definite
pathogen_cov_matrix <- as.matrix(nearPD(pathogen_cov_matrix)$mat)

# Ensure species order in covariance matrix matches the data
species_in_matrix <- rownames(pathogen_cov_matrix)
filtered_dataP <- complete_data %>%
  filter(Pathogen %in% species_in_matrix)

# Ensure the pathogen order in the covariance matrix matches the filtered data
pathogen_order <- unique(filtered_dataP$Pathogen)
pathogen_cov_matrix <- pathogen_cov_matrix[pathogen_order, pathogen_order]

# Verify the order
stopifnot(all(rownames(pathogen_cov_matrix) == pathogen_order))

# Identify species in the distance matrix and in the data
species_in_matrix <- rownames(pathogen_cov_matrix)
species_in_data <- unique(complete_data$Pathogen)

# Identify missing species
species_not_in_matrix <- setdiff(species_in_data, species_in_matrix)

# Calculate average distances for each taxonomic group from existing species
group_meansP <- complete_data %>%
  filter(Pathogen %in% species_in_matrix) %>%
  group_by(Pathogen_taxa) %>%
  summarise(across(all_of(species_in_matrix), mean, na.rm = TRUE)) %>%
  ungroup()

# Assign average distances to missing pathogens based on taxonomic group
assigned_distancesP <- missing_species_dataP %>%
  left_join(group_meansP, by = "Pathogen_taxa") %>%
  rowwise() %>%
  mutate(assigned_distance = list(
    case_when(
      Pathogen_taxa == "Gram-positive bacteria" ~ group_meansP %>% filter(Pathogen_taxa == "Gram-positive bacteria") %>% dplyr::select(-Pathogen_taxa) %>% as.numeric(),
      Pathogen_taxa == "Gram-negative bacteria" ~ group_meansP %>% filter(Pathogen_taxa == "Gram-negative bacteria") %>% dplyr::select(-Pathogen_taxa) %>% as.numeric(),
      Pathogen_taxa == "Protozoan parasite" ~ group_meansP %>% filter(Pathogen_taxa == "Protozoan parasite") %>% dplyr::select(-Pathogen_taxa) %>% as.numeric(),
      Pathogen_taxa == "Fungi" ~ group_meansP %>% filter(Pathogen_taxa == "Fungi") %>% dplyr::select(-Pathogen_taxa) %>% as.numeric(),
      Pathogen == "DNA_virus" ~ rep(20000, length(species_in_matrix)),
      Pathogen == "RNA_virus" ~ rep(30000, length(species_in_matrix)),
      TRUE ~ NA_real_
    )
  ))

# Create a new distance matrix with the missing species
new_distance_matrixP <- matrix(NA, nrow = length(species_in_matrix) + length(species_not_in_matrix),
                               ncol = length(species_in_matrix) + length(species_not_in_matrix),
                               dimnames = list(c(species_in_matrix, species_not_in_matrix),
                                               c(species_in_matrix, species_not_in_matrix)))

# Fill in the new distance matrix with existing distances
new_distance_matrixP[species_in_matrix, species_in_matrix] <- pathogen_dist_matrix

# Add the assigned distances for the missing species
for (i in seq_along(species_not_in_matrix)) {
  species <- species_not_in_matrix[i]
  distances <- assigned_distancesP$assigned_distance[[i]]
  new_distance_matrixP[species, species_in_matrix] <- distances
  new_distance_matrixP[species_in_matrix, species] <- distances
  new_distance_matrixP[species, species] <- 0
}

# Ensure the matrix is symmetric
new_distance_matrixP[upper.tri(new_distance_matrixP)] <- t(new_distance_matrixP)[upper.tri(new_distance_matrixP)]

# Replace any remaining NAs with a high value
new_distance_matrixP[is.na(new_distance_matrixP)] <- max(new_distance_matrixP, na.rm = TRUE) + 1

# Convert the distance matrix to a covariance matrix
max_dist <- max(new_distance_matrixP, na.rm = TRUE)
pathogen_cov_matrix_complete <- max_dist - new_distance_matrixP

# Make the covariance matrix positive definite
pathogen_cov_matrix_complete <- as.matrix(nearPD(pathogen_cov_matrix_complete)$mat)

# Ensure species order in covariance matrix matches the data
filtered_dataP <- complete_data %>%
  filter(Pathogen %in% rownames(pathogen_cov_matrix_complete))

# Ensure the pathogen order in the covariance matrix matches the filtered data
pathogen_order <- unique(filtered_dataP$Pathogen)
pathogen_cov_matrix_complete <- pathogen_cov_matrix_complete[pathogen_order, pathogen_order]

# Verify the order
stopifnot(all(rownames(pathogen_cov_matrix_complete) == pathogen_order))

# Ensure Pathogen column in filtered_dataP is a factor
filtered_dataP$Pathogen <- factor(filtered_dataP$Pathogen, levels = rownames(pathogen_cov_matrix_complete))

# Check if all levels of Pathogen in filtered_dataP are in the rownames of the covariance matrix
pathogens_in_data <- levels(filtered_dataP$Pathogen)
pathogens_in_matrix <- rownames(pathogen_cov_matrix_complete)

# Verify the order of the species in the data and the covariance matrix
all(pathogens_in_data == pathogens_in_matrix)

# Ensure Pathogen column in filtered_dataP is a factor with correct levels
filtered_dataP$Pathogen <- factor(filtered_dataP$Pathogen, levels = pathogens_in_matrix)

# Extract the row names from the covariance matrix
matrix_row_names <- rownames(pathogen_cov_matrix_complete)

# Extract the levels from the Pathogen factor
pathogen_levels <- levels(filtered_dataP$Pathogen)

# Ensure Pathogen column in filtered_dataP is a factor with correct levels
filtered_dataP$Pathogen <- factor(filtered_dataP$Pathogen, levels = matrix_row_names)

# Verify the dimensions of the covariance matrix
print(dim(pathogen_cov_matrix_complete))

# Ensure Pathogen is a factor with correct levels
filtered_dataP$Pathogen <- factor(filtered_dataP$Pathogen, levels = pathogens_in_matrix)

filtered_dataP <- as.data.frame(filtered_dataP)

# Define the ginverse list
ginverse_list <- list(Pathogen = as(pathogen_cov_matrix_complete, "dgCMatrix"))

# Define the priors for a single random effect (Pathogen)
priors2 <- list(R = list(V = 1, nu = 0.002),
               G = list(
                 G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1e-4)
               ))

# Run the MCMCglmm model
model_combinedP <- MCMCglmm(DeltaSub ~ Num_Obs,
                            random = ~ Pathogen,
                            ginverse = ginverse_list,
                            data = filtered_dataP,
                            family = "gaussian",
                            prior = priors,
                            nitt = 13000,
                            burnin = 3000,
                            thin = 10)

summary(model_combinedP)

# Extract fixed effects variance
fixed_effects_matrix2P <- model_combinedP$Sol %*% t(model_combinedP$X)
fixed_var2P <- apply(fixed_effects_matrix2P, 1, var)

# Extract random effects and residual variances
random_var2P <- apply(model_combinedP$VCV, 2, mean)
residual_var2P <- mean(model_combinedP$VCV[,"units"])

# Sum of variances of the random effects
total_random_var2P <- sum(random_var2P)

# Calculate Marginal and Conditional R-squared
marginal_r2_phylo2P <- mean(fixed_var2P) / (mean(fixed_var2P) + total_random_var2P + residual_var2P)
conditional_r2_phylo2P <- (mean(fixed_var2P) + total_random_var2P) / (mean(fixed_var2P) + total_random_var2P + residual_var2P)

# Create a data frame for plotting
r_squared_df2P <- data.frame(
  Model = c("Standard LM", "Phylogenetic model (Marginal)", "Phylogenetic model (Conditional)"),
  R_Squared = c(r_squared_lm, marginal_r2_phylo2P, conditional_r2_phylo2P)
)

# Extract the slopes 
phylo_slope_Num_Obs2P <- summary(model_combinedP)$solutions["Num_Obs", "post.mean"]

#Create a data frame for plotting the regression lines
plot_data2P <- data.frame(
  Num_Obs = seq(min(filtered_dataP$Num_Obs), max(filtered_dataP$Num_Obs), length.out = 100)
)

plot_data2P$phylo_fit_Num_Obs2P  <- phylo_slope_Num_Obs2P  * plot_data2$Num_Obs 

# Plot the data and regression lines for Num_Obs
plot1_P <- ggplot(filtered_dataP, aes(x = Num_Obs, y = DeltaSub)) +
  geom_point(shape = 21, size = 3) +
  geom_line(data = plot_data2, aes(y = phylo_fit_Num_Obs, color = "Phylogenetic model"), linetype = "dashed", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data2, aes(y = standard_fit_Num_Obs, color = "Standard model"), linetype = "solid", size = 1, show.legend = TRUE) +
  labs(
    x = "Number of observations",
    y = expression(Delta[AICc]),
    color = "Model"
  ) +
  scale_color_manual(values = c("Phylogenetic model" = "#0072B2", "Standard model" = "#D55E00")) +
  mytheme +
  scale_x_continuous(limits = c(5, 50)) +
  scale_y_continuous(limits = c(0, 300))+
  theme(legend.position = c(0.4, 0.9))



# Filter Bee data with 500 hosts 
filtered_dataP_bee <- as.data.frame(filtered_dataP[!(filtered_dataP$Num_Host == 500), ])

# Fit the phylogenetic model with the filtered data
# Run the MCMCglmm model
model_combinedP2 <- MCMCglmm(DeltaSub ~ Num_Host,
                            random = ~ Pathogen,
                            ginverse = ginverse_list,
                            data = filtered_dataP_bee,
                            family = "gaussian",
                            prior = priors,
                            nitt = 13000,
                            burnin = 3000,
                            thin = 10)
summary(model_combinedP2)

# Extract fixed effects variance
fixed_effects_matrix2NH <- model_combinedP2$Sol %*% t(model_combinedP2$X)
fixed_var2NH  <- apply(fixed_effects_matrix2NH, 1, var)

# Extract random effects and residual variances
random_var2NH <- apply(model_combinedP2$VCV, 2, mean)
residual_var2NH <- mean(model_combinedP2$VCV[,"units"])

# Sum of variances of the random effects
total_random_var2NH <- sum(random_var2NH)

# Calculate Marginal and Conditional R-squared
marginal_r2_phylo2NH <- mean(fixed_var2NH) / (mean(fixed_var2NH) + total_random_var2NH + residual_var2NH)
conditional_r2_phylo2NH <- (mean(fixed_var2NH) + total_random_var2NH) / (mean(fixed_var2NH) + total_random_var2NH + residual_var2NH)

# Create a data frame for plotting
r_squared_df2NH <- data.frame(
  Model = c("Phylogenetic model (Marginal)", "Phylogenetic model (Conditional)"),
  R_Squared = c(marginal_r2_phylo2NH, conditional_r2_phylo2NH)
)

# Extract the slopes 
phylo_slope_Num_Host2NH <- summary(model_combinedP2)$solutions["Num_Host", "post.mean"]

#Create a data frame for plotting the regression lines
plot_data2NH <- data.frame(
  Num_Host = seq(min(data$Num_Host), max(data$Num_Host), length.out = 100)
)

plot_data2NH$phylo_fit_Num_Host <- phylo_slope_Num_Host2NH * plot_data2NH$Num_Host

# Plot the data and regression lines for Num_Obs
plot2NH <- ggplot(filtered_dataP_bee, aes(x = Num_Host, y = DeltaSub)) +
  geom_point(shape = 21, size = 3) +
  geom_line(data = plot_data2NH, aes(y = phylo_fit_Num_Host), color = "#0072B2", linetype = "dashed", size = 1, show.legend = TRUE) +
  labs(
    x = "Number of host",
    y = expression(Delta[AICc]),
    color = "Model"
  ) +
  scale_color_manual(values = c("Phylogenetic model" = "#0072B2", "Standard model" = "#D55E00")) +
  mytheme +
  scale_x_continuous(limits = c(0, 100))+
  scale_y_continuous(limits = c(0, 250))

############################################

# Combined plot 
# Create a data frame for plotting
r_squared_df_ALL <- data.frame(
  Model = c("Standard LM", "Phylogenetic model (host)", "Phylogenetic model (pathogen)"),
  R_Squared = c(r_squared_lm, conditional_r2_phylo, conditional_r2_phylo2NH )
)

#Create a data frame for plotting the regression lines
plot_data_ALL <- data.frame(
  Num_Obs = seq(min(data$Num_Obs), max(data$Num_Obs), length.out = 100)
)

plot_data_ALL$phylo_fit_Num_Obs <- phylo_slope_Num_Obs * plot_data_ALL$Num_Obs
plot_data_ALL$standard_fit_Num_Obs <- standard_slope_Num_Obs * plot_data$Num_Obs
plot_data_ALL$phylo_fit_Num_Obs2P  <- phylo_slope_Num_Obs2P  * plot_data$Num_Obs 

# Plot the data and regression lines for Num_Obs
plot_ALL_NO <- ggplot(complete_data, aes(x = Num_Obs, y = DeltaSub)) +
  geom_point(shape = 21, size = 3) +
  geom_line(data = plot_data_ALL, aes(y = standard_fit_Num_Obs, color = "Standard model"), linetype = "solid", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data_ALL, aes(y = phylo_fit_Num_Obs, color = "Phylogenetic model (host)"), linetype = "dashed", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data_ALL, aes(y = phylo_fit_Num_Obs2P, color = "Phylogenetic model (pathogen)"), linetype = "dashed", size = 1, show.legend = TRUE)+
  labs(
    x = "Number of observations",
    y = expression(Delta[AICc]),
    color = "Model"
  ) +
  scale_color_manual(values = c("Phylogenetic model (host)" = "#0072B2", "Phylogenetic model (pathogen)" = "#D55E00", "Standard model" = "#CC79A7")) +
  mytheme +
  scale_x_continuous(limits = c(5, 50)) +
  scale_y_continuous(limits = c(0, 300))+
  theme(legend.position = c(0.4, 0.9))



# Create a data frame for plotting
r_squared_df_ALL_NH <- data.frame(
  Model = c("Standard LM", "Phylogenetic model (host)", "Phylogenetic model (pathogen)"),
  R_Squared = c(r_squared_lm2, conditional_r2_phylo2, conditional_r2_phylo2NH )
)

#Create a data frame for plotting the regression lines
plot_data2_ALL <- data.frame(
  Num_Host = seq(min(filtered_data$Num_Host), max(filtered_data$Num_Host), length.out = 100)
)

plot_data2_ALL$phylo_fit_Num_Host_H  <- phylo_slope_Num_Host2 * plot_data2$Num_Host
plot_data2_ALL$standard_fit_Num_Host <- standard_slope_Num_Host2 * plot_data2$Num_Host
plot_data2_ALL$phylo_fit_Num_Host_P  <- phylo_slope_Num_Host2NH * plot_data2$Num_Host


plot2NH_ALL <- plot_ALL_NO <- ggplot(complete_data, aes(x = Num_Host, y = DeltaSub)) +
  geom_point(shape = 21, size = 3) +
  geom_line(data = plot_data2_ALL, aes(y = standard_fit_Num_Host, color = "Standard model"), linetype = "solid", size = 1, show.legend = F) +
  geom_line(data = plot_data2_ALL, aes(y = phylo_fit_Num_Host_P, color = "Phylogenetic model (pathogen)"), linetype = "dashed", size = 1, show.legend = F)+
  geom_line(data = plot_data2_ALL, aes(y = phylo_fit_Num_Host_H, color = "Phylogenetic model (host)"), linetype = "dashed", size = 1, show.legend = F)+
  labs(
    x = "Number of observations",
    y = expression(Delta[AICc]),
    color = "Model"
  ) +
  scale_color_manual(values = c("Phylogenetic model (host)" = "#0072B2", "Phylogenetic model (pathogen)" = "#D55E00", "Standard model" = "#CC79A7")) +
  mytheme +
  #theme(legend.position = c(0.4, 0.9))+
  scale_x_continuous(limits = c(0, 100))+
  scale_y_continuous(limits = c(0, 250))

combined_plot_Mod_ALL <- plot_ALL_NO + plot2NH_ALL +
  plot_annotation(tag_levels = 'A')

#################################################################################
### ADDING DATA TYPE TO THE MODELS ###

# Ensure 'Category' is a factor
data$Category <- as.factor(data$Category)
filtered_data$Category <- as.factor(filtered_data$Category)
filtered_dataP$Category <- as.factor(filtered_dataP)
filtered_dataP_bee$Category <- as.factor(filtered_dataP_bee)


# Update the models to include the interaction term with Category
# Ensure 'Category' is a factor
data$Data <- as.factor(data$Data)
filtered_data$Data <- as.factor(filtered_data$Data)
filtered_dataP$Data <- as.factor(filtered_dataP)
filtered_dataP_bee$Data <- as.factor(filtered_dataP_bee)

# Subset data by Data
data_cat1 <- subset(data, Data == "raw")
data_cat2 <- subset(data, Data == "probability")

filtered_data_cat1 <- subset(filtered_data, Data == "raw")
filtered_data_cat2 <- subset(filtered_data, Data == "probability")

filtered_dataP_cat1 <- subset(filtered_dataP, Data == "raw")
filtered_dataP_cat2 <- subset(filtered_dataP, Data == "probability")

filtered_dataP_bee_cat1 <- subset(filtered_dataP_bee, Data == "raw")
filtered_dataP_bee_cat2 <- subset(filtered_dataP_bee, Data == "probability")

# Fit models for each category
model_cat1 <- MCMCglmm(DeltaSub ~ Num_Obs,
                       random = ~ species,
                       ginverse = list(species = phylo_cov_matrix),
                       data = data_cat1,
                       family = "gaussian",
                       prior = priors,
                       nitt = 13000,
                       burnin = 3000,
                       thin = 10)

model_cat2 <- MCMCglmm(DeltaSub ~ Num_Obs,
                       random = ~ species,
                       ginverse = list(species = phylo_cov_matrix),
                       data = data_cat2,
                       family = "gaussian",
                       prior = priors,
                       nitt = 13000,
                       burnin = 3000,
                       thin = 10)

model2_cat1 <- MCMCglmm(DeltaSub ~ Num_Host,
                        random = ~ species,
                        ginverse = list(species = filtered_phylo_cov_matrix),
                        data = filtered_data_cat1,
                        family = "gaussian",
                        prior = priors,
                        nitt = 13000,
                        burnin = 3000,
                        thin = 10)

model2_cat2 <- MCMCglmm(DeltaSub ~ Num_Host,
                        random = ~ species,
                        ginverse = list(species = filtered_phylo_cov_matrix),
                        data = filtered_data_cat2,
                        family = "gaussian",
                        prior = priors,
                        nitt = 13000,
                        burnin = 3000,
                        thin = 10)

model_combinedP_cat1 <- MCMCglmm(DeltaSub ~ Num_Obs,
                                 random = ~ Pathogen,
                                 ginverse = ginverse_list,
                                 data = filtered_dataP_cat1,
                                 family = "gaussian",
                                 prior = priors,
                                 nitt = 13000,
                                 burnin = 3000,
                                 thin = 10)

model_combinedP_cat2 <- MCMCglmm(DeltaSub ~ Num_Obs,
                                 random = ~ Pathogen,
                                 ginverse = ginverse_list,
                                 data = filtered_dataP_cat2,
                                 family = "gaussian",
                                 prior = priors,
                                 nitt = 13000,
                                 burnin = 3000,
                                 thin = 10)

model_combinedP2_cat1 <- MCMCglmm(DeltaSub ~ Num_Host,
                                  random = ~ Pathogen,
                                  ginverse = ginverse_list,
                                  data = filtered_dataP_bee_cat1,
                                  family = "gaussian",
                                  prior = priors,
                                  nitt = 13000,
                                  burnin = 3000,
                                  thin = 10)

model_combinedP2_cat2 <- MCMCglmm(DeltaSub ~ Num_Host,
                                  random = ~ Pathogen,
                                  ginverse = ginverse_list,
                                  data = filtered_dataP_bee_cat2,
                                  family = "gaussian",
                                  prior = priors,
                                  nitt = 13000,
                                  burnin = 3000,
                                  thin = 10)

# Summarize the models
summary(model_cat1)
summary(model_cat2)
summary(model2_cat1)
summary(model2_cat2)
summary(model_combinedP_cat1)
summary(model_combinedP_cat2)
summary(model_combinedP2_cat1)
summary(model_combinedP2_cat2)

# Fit standard linear models for each category
standard_model_cat1 <- lm(DeltaSub ~ Num_Obs, data = data_cat1)
standard_model_cat2 <- lm(DeltaSub ~ Num_Obs, data = data_cat2)

standard_model2_cat1 <- lm(DeltaSub ~ Num_Host, data = filtered_data_cat1)
standard_model2_cat2 <- lm(DeltaSub ~ Num_Host, data = filtered_data_cat2)

# Extract the slopes for plotting
standard_slope_Num_Obs_cat1 <- coef(standard_model_cat1)["Num_Obs"]
standard_slope_Num_Obs_cat2 <- coef(standard_model_cat2)["Num_Obs"]

standard_slope_Num_Host_cat1 <- coef(standard_model2_cat1)["Num_Host"]
standard_slope_Num_Host_cat2 <- coef(standard_model2_cat2)["Num_Host"]

# Create data frames for plotting
plot_data_cat1 <- data.frame(
  Num_Obs = seq(min(data_cat1$Num_Obs), max(data_cat1$Num_Obs), length.out = 100)
)

plot_data_cat2 <- data.frame(
  Num_Obs = seq(min(data_cat2$Num_Obs), max(data_cat2$Num_Obs), length.out = 100)
)

plot_data2_cat1 <- data.frame(
  Num_Host = seq(min(filtered_data_cat1$Num_Host), max(filtered_data_cat1$Num_Host), length.out = 100)
)

plot_data2_cat2 <- data.frame(
  Num_Host = seq(min(filtered_data_cat2$Num_Host), max(filtered_data_cat2$Num_Host), length.out = 100)
)

plot_data2P_cat1 <- data.frame(
  Num_Obs = seq(min(filtered_dataP_cat1$Num_Obs), max(filtered_dataP_cat1$Num_Obs), length.out = 100)
)

plot_data2P_cat2 <- data.frame(
  Num_Obs = seq(min(filtered_dataP_cat2$Num_Obs), max(filtered_dataP_cat2$Num_Obs), length.out = 100)
)

plot_data2NH_cat1 <- data.frame(
  Num_Host = seq(min(filtered_dataP_bee_cat1$Num_Host), max(filtered_dataP_bee_cat1$Num_Host), length.out = 100)
)

plot_data2NH_cat2 <- data.frame(
  Num_Host = seq(min(filtered_dataP_bee_cat2$Num_Host), max(filtered_dataP_bee_cat2$Num_Host), length.out = 100)
)

# Add the fitted values to the plot data
plot_data_cat1$phylo_fit_Num_Obs <- phylo_slope_Num_Obs_cat1 * plot_data_cat1$Num_Obs
plot_data_cat2$phylo_fit_Num_Obs <- phylo_slope_Num_Obs_cat2 * plot_data_cat2$Num_Obs
plot_data_cat1$standard_fit_Num_Obs <- standard_slope_Num_Obs_cat1 * plot_data_cat1$Num_Obs
plot_data_cat2$standard_fit_Num_Obs <- standard_slope_Num_Obs_cat2 * plot_data_cat2$Num_Obs

plot_data2_cat1$phylo_fit_Num_Host <- phylo_slope_Num_Host_cat1 * plot_data2_cat1$Num_Host
plot_data2_cat2$phylo_fit_Num_Host <- phylo_slope_Num_Host_cat2 * plot_data2_cat2$Num_Host
plot_data2_cat1$standard_fit_Num_Host <- standard_slope_Num_Host_cat1 * plot_data2_cat1$Num_Host
plot_data2_cat2$standard_fit_Num_Host <- standard_slope_Num_Host_cat2 * plot_data2_cat2$Num_Host

plot_data2P_cat1$phylo_fit_Num_Obs2P <- phylo_slope_Num_Obs2P_cat1 * plot_data2P_cat1$Num_Obs
plot_data2P_cat2$phylo_fit_Num_Obs2P <- phylo_slope_Num_Obs2P_cat2 * plot_data2P_cat2$Num_Obs

plot_data2NH_cat1$phylo_fit_Num_Host2NH <- phylo_slope_Num_Host2NH_cat1 * plot_data2NH_cat1$Num_Host
plot_data2NH_cat2$phylo_fit_Num_Host2NH <- phylo_slope_Num_Host2NH_cat2 * plot_data2NH_cat2$Num_Host

# Plot the data and regression lines for each category


plot1.1 <- ggplot(data_cat1, aes(x = Num_Obs, y = DeltaSub)) +
  geom_point(shape = 21, size = 3) +
  geom_line(data = plot_data_cat1, aes(y = standard_fit_Num_Obs, color = "Standard model"), linetype = "solid", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data_cat1, aes(y = phylo_fit_Num_Obs, color = "Phylo. model (H)"), linetype = "dashed", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data2P_cat1, aes(y = phylo_fit_Num_Obs2P, color = "Phylo. model (P)"), linetype = "dashed", size = 1, show.legend = TRUE) +
  labs(
    x = "Number of observations",
    y = expression(Delta[AICc]),
    color = "Model"
  ) +
  scale_color_manual(values = c("Standard model" = "#CC79A7", "Phylo. model (H)" = "#0072B2",
                                "Phylo. model (P)" = "#E69F00")) +
  mytheme +
  scale_x_continuous(limits = c(5, 50)) +
  scale_y_continuous(limits = c(0, 300)) +
  theme(legend.position = c(0.35, 0.8)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())


plot1.2 <- ggplot(data_cat2, aes(x = Num_Obs, y = DeltaSub)) +
  geom_point(shape = 21, size = 3) +
  geom_line(data = plot_data_cat2, aes(y = standard_fit_Num_Obs, color = "Standard model -probability"), linetype = "solid", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data_cat2, aes(y = phylo_fit_Num_Obs, color = "Phylo. model -probability (H)"), linetype = "dashed", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data2P_cat2, aes(y = phylo_fit_Num_Obs2P, color = "Phylo. model -probability (P)"), linetype = "dashed", size = 1, show.legend = TRUE) +
  
  labs(
    x = "Number of observations",
    y = expression(Delta[AICc]),
    color = "Model"
  ) +
  scale_color_manual(values = c("Standard model -probability" = "#CC79A7", "Phylo. model -probability (H)" = "#0072B2",
                                "Phylo. model -probability (P)" = "#E69F00")) +
  mytheme +
  scale_x_continuous(limits = c(5, 50)) +
  scale_y_continuous(limits = c(0, 300)) +
  #theme(legend.position = c(0.25, 0.85))
  theme(legend.position = "none")

plot2.1 <- ggplot(data_cat1, aes(x = Num_Host, y = DeltaSub)) +
  geom_point(shape = 21, size = 3) +
  geom_line(data = plot_data2_cat1, aes(y = standard_fit_Num_Host, color = "Standard model -raw"), linetype = "solid", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data2_cat1, aes(y = phylo_fit_Num_Host, color = "Phylo. model -raw (H)"), linetype = "dashed", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data2NH_cat1, aes(y = phylo_fit_Num_Host2NH, color = "Phylo. model -raw (P)"), linetype = "dashed", size = 1, show.legend = TRUE) +
  labs(
    x = "Number of hosts",
    y = expression(Delta[AICc]),
    color = "Model"
  ) +
  scale_color_manual(values = c("Standard model -raw" = "#CC79A7", "Phylo. model -raw (H)" = "#0072B2",
                                "Phylo. model -raw (P)" = "#E69F00")) +
  mytheme +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 250)) +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank())


plot2.2 <- ggplot(data_cat2, aes(x = Num_Host, y = DeltaSub)) +
  geom_point(shape = 21, size = 3) +
  geom_line(data = plot_data2_cat2, aes(y = standard_fit_Num_Host, color = "Standard model -probability"), linetype = "solid", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data2_cat2, aes(y = standard_fit_Num_Host, color = "Phylo. model -probability (H)"), linetype = "solid", size = 1, show.legend = TRUE) +
  geom_line(data = plot_data2NH_cat2, aes(y = phylo_fit_Num_Host2NH, color = "Phylo. model -probability (P)"), linetype = "dashed", size = 1, show.legend = TRUE) +
  labs(
    x = "Number of hosts",
    y = expression(Delta[AICc]),
    color = "Model"
  ) +
  scale_color_manual(values =c("Standard model -probability" = "#CC79A7", "Phylo. model -probability (H)" = "#0072B2",
    "Phylo. model -probability (P)" = "#E69F00"))+
  mytheme +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 250)) +
  theme(legend.position = "none")+
  theme(axis.title.y = element_blank())

combined_plot_Mod <- plot1.1 +  plot2.1 + plot1.2 + plot2.2 + plot_annotation(tag_levels = 'A')
combined_plot_Mod


# Calculate R-squared for standard linear models
r_squared_lm_cat1 <- summary(standard_model_cat1)$r.squared
r_squared_lm_cat2 <- summary(standard_model_cat2)$r.squared

r_squared_lm2_cat1 <- summary(standard_model2_cat1)$r.squared
r_squared_lm2_cat2 <- summary(standard_model2_cat2)$r.squared

# Calculate R-squared for phylogenetic models
# Marginal and Conditional R-squared for model_cat1
fixed_effects_matrix_cat1 <- model_cat1$Sol %*% t(model_cat1$X)
fixed_var_cat1 <- apply(fixed_effects_matrix_cat1, 1, var)
random_var_cat1 <- apply(model_cat1$VCV, 2, mean)
residual_var_cat1 <- mean(model_cat1$VCV[,"units"])
total_random_var_cat1 <- sum(random_var_cat1)
marginal_r2_phylo_cat1 <- mean(fixed_var_cat1) / (mean(fixed_var_cat1) + total_random_var_cat1 + residual_var_cat1)
conditional_r2_phylo_cat1 <- (mean(fixed_var_cat1) + total_random_var_cat1) / (mean(fixed_var_cat1) + total_random_var_cat1 + residual_var_cat1)

# Marginal and Conditional R-squared for model_cat2
fixed_effects_matrix_cat2 <- model_cat2$Sol %*% t(model_cat2$X)
fixed_var_cat2 <- apply(fixed_effects_matrix_cat2, 1, var)
random_var_cat2 <- apply(model_cat2$VCV, 2, mean)
residual_var_cat2 <- mean(model_cat2$VCV[,"units"])
total_random_var_cat2 <- sum(random_var_cat2)
marginal_r2_phylo_cat2 <- mean(fixed_var_cat2) / (mean(fixed_var_cat2) + total_random_var_cat2 + residual_var_cat2)
conditional_r2_phylo_cat2 <- (mean(fixed_var_cat2) + total_random_var_cat2) / (mean(fixed_var_cat2) + total_random_var_cat2 + residual_var_cat2)

# Marginal and Conditional R-squared for model2_cat1
fixed_effects_matrix2_cat1 <- model2_cat1$Sol %*% t(model2_cat1$X)
fixed_var2_cat1 <- apply(fixed_effects_matrix2_cat1, 1, var)
random_var2_cat1 <- apply(model2_cat1$VCV, 2, mean)
residual_var2_cat1 <- mean(model2_cat1$VCV[,"units"])
total_random_var2_cat1 <- sum(random_var2_cat1)
marginal_r2_phylo2_cat1 <- mean(fixed_var2_cat1) / (mean(fixed_var2_cat1) + total_random_var2_cat1 + residual_var2_cat1)
conditional_r2_phylo2_cat1 <- (mean(fixed_var2_cat1) + total_random_var2_cat1) / (mean(fixed_var2_cat1) + total_random_var2_cat1 + residual_var2_cat1)

# Marginal and Conditional R-squared for model2_cat2
fixed_effects_matrix2_cat2 <- model2_cat2$Sol %*% t(model2_cat2$X)
fixed_var2_cat2 <- apply(fixed_effects_matrix2_cat2, 1, var)
random_var2_cat2 <- apply(model2_cat2$VCV, 2, mean)
residual_var2_cat2 <- mean(model2_cat2$VCV[,"units"])
total_random_var2_cat2 <- sum(random_var2_cat2)
marginal_r2_phylo2_cat2 <- mean(fixed_var2_cat2) / (mean(fixed_var2_cat2) + total_random_var2_cat2 + residual_var2_cat2)
conditional_r2_phylo2_cat2 <- (mean(fixed_var2_cat2) + total_random_var2_cat2) / (mean(fixed_var2_cat2) + total_random_var2_cat2 + residual_var2_cat2)

# Marginal and Conditional R-squared for model_combinedP_cat1
fixed_effects_matrix2P_cat1 <- model_combinedP_cat1$Sol %*% t(model_combinedP_cat1$X)
fixed_var2P_cat1 <- apply(fixed_effects_matrix2P_cat1, 1, var)
random_var2P_cat1 <- apply(model_combinedP_cat1$VCV, 2, mean)
residual_var2P_cat1 <- mean(model_combinedP_cat1$VCV[,"units"])
total_random_var2P_cat1 <- sum(random_var2P_cat1)
marginal_r2_phylo2P_cat1 <- mean(fixed_var2P_cat1) / (mean(fixed_var2P_cat1) + total_random_var2P_cat1 + residual_var2P_cat1)
conditional_r2_phylo2P_cat1 <- (mean(fixed_var2P_cat1) + total_random_var2P_cat1) / (mean(fixed_var2P_cat1) + total_random_var2P_cat1 + residual_var2P_cat1)

# Marginal and Conditional R-squared for model_combinedP_cat2
fixed_effects_matrix2P_cat2 <- model_combinedP_cat2$Sol %*% t(model_combinedP_cat2$X)
fixed_var2P_cat2 <- apply(fixed_effects_matrix2P_cat2, 1, var)
random_var2P_cat2 <- apply(model_combinedP_cat2$VCV, 2, mean)
residual_var2P_cat2 <- mean(model_combinedP_cat2$VCV[,"units"])
total_random_var2P_cat2 <- sum(random_var2P_cat2)
marginal_r2_phylo2P_cat2 <- mean(fixed_var2P_cat2) / (mean(fixed_var2P_cat2) + total_random_var2P_cat2 + residual_var2P_cat2)
conditional_r2_phylo2P_cat2 <- (mean(fixed_var2P_cat2) + total_random_var2P_cat2) / (mean(fixed_var2P_cat2) + total_random_var2P_cat2 + residual_var2P_cat2)

# Marginal and Conditional R-squared for model_combinedP2_cat1
fixed_effects_matrix2NH_cat1 <- model_combinedP2_cat1$Sol %*% t(model_combinedP2_cat1$X)
fixed_var2NH_cat1 <- apply(fixed_effects_matrix2NH_cat1, 1, var)
random_var2NH_cat1 <- apply(model_combinedP2_cat1$VCV, 2, mean)
residual_var2NH_cat1 <- mean(model_combinedP2_cat1$VCV[,"units"])
total_random_var2NH_cat1 <- sum(random_var2NH_cat1)
marginal_r2_phylo2NH_cat1 <- mean(fixed_var2NH_cat1) / (mean(fixed_var2NH_cat1) + total_random_var2NH_cat1 + residual_var2NH_cat1)
conditional_r2_phylo2NH_cat1 <- (mean(fixed_var2NH_cat1) + total_random_var2NH_cat1) / (mean(fixed_var2NH_cat1) + total_random_var2NH_cat1 + residual_var2NH_cat1)

# Marginal and Conditional R-squared for model_combinedP2_cat2
fixed_effects_matrix2NH_cat2 <- model_combinedP2_cat2$Sol %*% t(model_combinedP2_cat2$X)
fixed_var2NH_cat2 <- apply(fixed_effects_matrix2NH_cat2, 1, var)
random_var2NH_cat2 <- apply(model_combinedP2_cat2$VCV, 2, mean)
residual_var2NH_cat2 <- mean(model_combinedP2_cat2$VCV[,"units"])
total_random_var2NH_cat2 <- sum(random_var2NH_cat2)
marginal_r2_phylo2NH_cat2 <- mean(fixed_var2NH_cat2) / (mean(fixed_var2NH_cat2) + total_random_var2NH_cat2 + residual_var2NH_cat2)
conditional_r2_phylo2NH_cat2 <- (mean(fixed_var2NH_cat2) + total_random_var2NH_cat2) / (mean(fixed_var2NH_cat2) + total_random_var2NH_cat2 + residual_var2NH_cat2)

# Combine all R-squared values into a data frame for plotting
r_squared_df <- data.frame(
  Model = c(
    "Standard LM/Num. obs.-raw", "Standard LM/Num. obs.-probability", 
    "Phylo. model-Host/Num. obs.-raw (Marginal)", "Phylo. model-Host/Num. obs.-raw (Conditional)", 
    "Phylo. model-Host/Num. obs.-probability (Marginal)", "Phylo. model-Host/Num. obs.-probability (Conditional)", 
    "Standard LM/Num. host-raw", "Standard LM/Num. host-probability", 
    "Phylo. model-Pathogen/Num. host.-raw (Marginal)", "Phylo. model-Pathogen/Num. host.-raw (Conditional)", 
    "Phylo. model-Pathogen/Num. host.-probability (Marginal)", "Phylo. model-Pathogen/Num. host.-probability (Conditional)",
    "Phylo. model-Pathoghen/Num. obs.-raw (Marginal)", "Phylo. model-Pathoghen/Num. obs.-raw (Conditional)",
    "Phylo. model-Pathoghen/Num. obs.-probability (Marginal)", "Phylo. model-Pathoghen/Num. obs.-probability (Conditional)",
    "Phylo. model-Pathoghen/Num. host-raw (Marginal)", "Phylo. model-Pathoghen/Num. host-raw (Conditional)",
    "Phylo. model-Pathoghen/Num. host-probability (Marginal)", "Phylo. model-Pathoghen/Num. host-probability  (Conditional)"
  ),
  R_Squared = c(
    r_squared_lm_cat1, r_squared_lm_cat2,
    marginal_r2_phylo_cat1, conditional_r2_phylo_cat1,
    marginal_r2_phylo_cat2, conditional_r2_phylo_cat2,
    r_squared_lm2_cat1, r_squared_lm2_cat2,
    marginal_r2_phylo2_cat1, conditional_r2_phylo2_cat1,
    marginal_r2_phylo2_cat2, conditional_r2_phylo2_cat2,
    marginal_r2_phylo2P_cat1, conditional_r2_phylo2P_cat1,
    marginal_r2_phylo2P_cat2, conditional_r2_phylo2P_cat2,
    marginal_r2_phylo2NH_cat1, conditional_r2_phylo2NH_cat1,
    marginal_r2_phylo2NH_cat2, conditional_r2_phylo2NH_cat2
  )
)

write.table(r_squared_df, file = "r_squared_values.csv", sep = ",", row.names = FALSE, col.names = TRUE)

# Boxplot for figure 4
ggplot(data, aes(x = Data, y = DeltaSub)) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_jitter(shape = 21, aes(size = Num_Host), color = "grey25") +
  geom_boxplot(alpha = 0.5)+
    labs(x = "Data type",
    y = expression(Delta[AICc])
  ) +
  #scale_color_manual(values =c("Standard model -probability" = "#CC79A7", "Phylo. model -probability (H)" = "#0072B2",
                               #"Phylo. model -probability (P)" = "#E69F00")+
  mytheme +
  scale_y_continuous(limits = c(0, 275)) +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())

