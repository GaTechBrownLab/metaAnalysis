# Meta survival data analysis
# Author: Canan Karakoc
# Last update: October 30, 2024

################################################################################################
# SETUP #
################################################################################################

# Install/load libraries
# Some packages are loded in-code 

library(tidyverse)
library(minpack.lm) # for regressions
library(ggforce) # for the supplementary plot 
library(patchwork) # to combine plots

# Function to calculate standard error
standard_error <- function(x) {
  n <- sum(!is.na(x)) # Number of non-NA observations
  if (n == 0) {
    return(NA)
  } else {
    return(sd(x, na.rm = TRUE) / sqrt(n))
  }
}

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
# Data will be uploaded to Figshare - Too large for GitHub
main.dir <- "~/GitHub/metaAnalysis_largeFiles/metaData"
allFiles <- list.files(path = main.dir, recursive = TRUE, full.names = TRUE)
csv_files <- allFiles[grep("\\.csv$", allFiles)]

# Meta data for later
meta <- read.csv("~/GitHub/metaAnalysis/data/000_combined_library.csv",
  sep = ",", header = T, na.strings = ""
)

###############################################################################################
# DATA CLEANING # 
###############################################################################################

# Read all CSV files into a list of data frames
dataAll <- csv_files %>%
  setNames(basename(csv_files)) %>%
  map_df(read_csv, .id = "Full_key", col_names = c("Time", "Survival")) 

# Clean the raw data before standardization
dataClean <- dataAll %>%
  mutate(key = gsub("\\..*", "", Full_key)) %>%
  mutate(Key = gsub("\\_.*", "", key)) %>%
  dplyr::select(-key) %>%
  left_join(meta, by = "Key") %>%
  mutate(across(Time:Survival, ~ ifelse(.x < 0, 0, .x))) %>% # Replace negative values with 0
  mutate(across(Survival, ~ ifelse(.x > 100, 100, .x)))%>%   # Cap survival at 100%
  mutate(Survival = round(Survival, 0)) # Counts can't be floating number 

# Function to remove repeating zeros after the first occurrence of zero
remove_repeating_zeros <- function(data) {
  # Check if there are any zeros in the Survival_corrected column
  if (any(data$Survival_corrected == 0)) {
    # Identify the first occurrence of zero in the Survival_corrected column
    first_zero_index <- min(which(data$Survival_corrected == 0))
    
    # Keep rows only up to the first occurrence of zero
    return(data[1:first_zero_index, ])
  } else {
    # If there are no zeros, return the entire dataset
    return(data)
  }
}

# Function to ensure survival is always decreasing or the same over time
ensure_decreasing_survival <- function(Survival_corrected) {
  # Loop through the survival values and ensure non-increasing values
  for (i in 2:length(Survival_corrected)) {
    if (Survival_corrected[i] > Survival_corrected[i - 1]) {
      Survival_corrected[i] <- Survival_corrected[i - 1]  # Set current value to the previous if it's higher
    }
  }
  return(Survival_corrected)
}

# Reverse datasets and ensure survival decreases over time
# Filter out datasets where the last survival value is greater than 20
# Remove trailing zeros from Survival_corrected

# Find unique datasets to reverse using Reverse_Flag column
datasets_to_reverse <- dataClean %>%
  filter(Graph == "mortality") %>%
  pull(Full_key) %>%
  unique()

dataCleanSurvival <- dataClean %>%
  group_by(Full_key) %>%
  group_modify(~ {
    .x <- .x %>%
      arrange(Time)  # Always ensure the data is sorted by Time
    
    # Reverse survival if the dataset is marked for reversal
    if (.y$Full_key %in% datasets_to_reverse) {
      .x <- .x %>%
        mutate(Survival_corrected = 100 - Survival)  # Reverse Survival for datasets needing reversal
    } else {
      .x <- .x %>%
        mutate(Survival_corrected = Survival)  # Otherwise, keep Survival unchanged
    }
    
    # Ensure the Survival values are non-increasing
    .x <- .x %>%
      mutate(Survival_corrected = ensure_decreasing_survival(Survival_corrected))  # Apply the non-increasing rule
    
    return(.x)
  }) %>%
  ungroup() %>%  # Ungroup after the group_modify operation
  group_by(Full_key) %>%
  filter(last(Survival_corrected) <= 20) %>%
  group_modify(~ remove_repeating_zeros(.x)) %>%  # Pass the whole dataframe to the function
  ungroup()

# Find datasets with time unit in days
datasets_in_days <- dataCleanSurvival %>%
  filter(Time_units == "days") %>%
  pull(Full_key) %>%
  unique()

datasets_in_weeks <- dataCleanSurvival %>%
  filter(Time_units == "weeks") %>%
  pull(Full_key) %>%
  unique()

# Continue with the rest of your pipeline (time conversion, standardization, etc.)
dataCleanSurvival_hours <- dataCleanSurvival %>%
  mutate(Time_hours = case_when(
    Full_key %in% datasets_in_days ~ Time * 24,    # Convert days to hours
    Full_key %in% datasets_in_weeks ~ Time * 7 * 24, # Convert weeks to hours
    TRUE ~ Time  # If already in hours, keep as is
  )) 

# Standardize both Survival and Time values between 0 and 1
standardized_df <- dataCleanSurvival_hours %>%
  group_by(Full_key) %>%
  
  # Check if all survival values are constant (i.e., max == min)
  mutate(is_constant = max(Survival_corrected) == min(Survival_corrected)) %>%
  
  # Apply standardization for non-constant datasets
  mutate(
    Standardized_Survival = ifelse(
      !is_constant,  # Only apply standardization when survival varies
      (Survival_corrected - min(Survival_corrected)) / (max(Survival_corrected) - min(Survival_corrected)),
      1  # For constant survival datasets, assign 1
    ),
    
    # Standardize Time_hours between 0 and 1
    Time_std = (Time_hours - min(Time_hours)) / (max(Time_hours) - min(Time_hours))
  ) %>%
  
  # Ensure the last survival value is 0
  mutate(Standardized_Survival = ifelse(row_number() == n(), 0, Standardized_Survival)) %>%
  
  # Replace very small values with 0.001 to avoid zero issues (after setting last value to 0)
  mutate(Standardized_Survival = ifelse(Standardized_Survival < 0.001 & Standardized_Survival != 0, 0.001, Standardized_Survival)) %>%
  
  ungroup() %>%
  dplyr::select(-is_constant)  # Remove the helper column

# Final data to use
standardized_final_df <- standardized_df %>%
  dplyr::select(
    Full_key, Key, Time, Survival, Time_hours, Time_std, Survival_corrected, Standardized_Survival,
    Host_taxa, Pathogen_taxa, Max_num_host, Data
  )

# Number of datasets after filtering
unique(standardized_final_df$Full_key)

ggplot(dataCleanSurvival, aes(x = Time, y = Survival_corrected, group = Full_key)) +
  geom_line() + geom_point() +
  labs(title = "Survival Corrected Before Standardization")

ggplot(standardized_final_df, aes(x = Time, y = Standardized_Survival, group = Full_key)) +
  geom_line() + geom_point() +
  labs(title = "Standardized Survival After Cleaning")

write.csv(standardized_final_df, "data/standardized_final_df.csv", row.names = FALSE)

################################################################################################
# MODELS # 
################################################################################################

# Create folders to store the plots
dir.create("exponential_model_plots", showWarnings = FALSE)
dir.create("weibull_model_plots", showWarnings = FALSE)
dir.create("gompertz_survival_model_plots", showWarnings = FALSE)
dir.create("loglogistic_model_plots", showWarnings = FALSE)
dir.create("generalizedgamma_model_plots", showWarnings = FALSE)

# Model fitting functions with explicit parameter extraction for recording

# Exponential model fitting
fit_exponential <- function(data, start_lambda = 0.1) {
  message("Fitting Exponential model...")
  fit <- try(nlsLM(Standardized_Survival ~ exp(-lambda * Time_std), data = data, 
                   start = list(lambda = start_lambda), control = nls.lm.control(maxiter = 1000)), silent = FALSE)
  
  if (inherits(fit, "try-error") || any(is.na(coef(fit)))) {
    message("Exponential model fitting failed.")
    return(NULL)
  } else {
    lambda <- coef(fit)["lambda"]
    message("Exponential model fit successful.")
    return(list(fit = fit, params = list(lambda = as.numeric(lambda)), rss = sum(residuals(fit)^2)))
  }
}

# Weibull model fitting
fit_weibull <- function(data, start_lambda = 1, start_k = 0.1) {
  message("Fitting Weibull model...")
  fit <- try(nlsLM(Standardized_Survival ~ exp(-(Time_std / lambda)^k), data = data, 
                   start = list(lambda = start_lambda, k = start_k), control = nls.lm.control(maxiter = 1000)), silent = FALSE)
  
  if (inherits(fit, "try-error") || any(is.na(coef(fit)))) {
    message("Weibull model fitting failed.")
    return(NULL)
  } else {
    lambda <- coef(fit)["lambda"]
    k <- coef(fit)["k"]
    message("Weibull model fit successful.")
    return(list(fit = fit, params = list(lambda = as.numeric(lambda), k = as.numeric(k)), rss = sum(residuals(fit)^2)))
  }
}

# Gompertz survival model with Standardized Time
fit_gompertz_survival <- function(data, start_a = 1, start_b = 1) {
  message("Fitting Gompertz Survival model (Standardized Time)...")
  fit <- try(nlsLM(Standardized_Survival ~ exp(-(a / b) * (exp(b * Time_std) - 1)), data = data, 
                   start = list(a = start_a, b = start_b), 
                   #lower = c(a = 0, b = 0), upper = c(a = Inf, b = Inf), 
                   control = nls.lm.control(maxiter = 1000)), silent = FALSE)
  
  if (inherits(fit, "try-error") || any(is.na(coef(fit)))) {
    message("Gompertz model fitting failed.")
    return(NULL)
  } else {
    a <- coef(fit)["a"]
    b <- coef(fit)["b"]
    lt50 <- (1 / b) * log((log(2) * b / a) + 1)  # LT50 calculation
    message("Gompertz model fit successful.")
    return(list(fit = fit, params = list(a = as.numeric(a), b = as.numeric(b)), LT50 = lt50, rss = sum(residuals(fit)^2)))
  }
}

# Gompertz survival model with Raw Time (Time_hours)
fit_gompertz_survival_rawtime <- function(data, start_a = 0.005, start_b = 0.0005) {
  message("Fitting Gompertz Survival model (Raw Time)...")
  fit <- try(nlsLM(Standardized_Survival ~ exp(-(a / b) * (exp(b * Time_hours) - 1)), data = data, 
                   start = list(a = start_a, b = start_b), 
                   #lower = c(a = 0, b = 0), upper = c(a = 1, b = 1), 
                   control = nls.lm.control(maxiter = 1000)), silent = FALSE)
  
  if (inherits(fit, "try-error") || any(is.na(coef(fit)))) {
    message("Gompertz Raw Time model fitting failed.")
    return(NULL)
  } else {
    a <- coef(fit)["a"]
    b <- coef(fit)["b"]
    lt50 <- (1 / b) * log((log(2) * b / a) + 1)  # LT50 calculation
    message("Gompertz Raw Time model fit successful.")
    return(list(fit = fit, params = list(a = as.numeric(a), b = as.numeric(b)), LT50 = lt50, rss = sum(residuals(fit)^2)))
  }
}

# Log-Logistic model fitting
fit_loglogistic <- function(data, start_alpha = 1, start_beta = 1) {
  message("Fitting Log-Logistic model...")
  fit <- try(nlsLM(Standardized_Survival ~ 1 / (1 + (Time_std / alpha)^beta), data = data, 
                   start = list(alpha = start_alpha, beta = start_beta), 
                   #lower = c(alpha = 0, beta = 0), upper = c(alpha = Inf, beta = Inf), 
                   control = nls.lm.control(maxiter = 1000)), silent = FALSE)
  
  if (inherits(fit, "try-error") || any(is.na(coef(fit)))) {
    message("Log-Logistic model fitting failed.")
    return(NULL)
  } else {
    alpha <- coef(fit)["alpha"]
    beta <- coef(fit)["beta"]
    message("Log-Logistic model fit successful.")
    return(list(fit = fit, params = list(alpha = as.numeric(alpha), beta = as.numeric(beta)), rss = sum(residuals(fit)^2)))
  }
}

# Generalized Gamma model fitting
fit_generalizedgamma <- function(data, start_beta = 1, start_gamma = 1, start_alpha = 1) {
  message("Fitting Generalized Gamma model...")
  fit <- try(nlsLM(Standardized_Survival ~ exp(-(Time_std / beta)^gamma * (1 + (alpha - 1) * (Time_std / beta)^gamma)), 
                   data = data, start = list(beta = start_beta, gamma = start_gamma, alpha = start_alpha), 
                   #lower = c(beta = 0, gamma = 0, alpha = 0), upper = c(beta = Inf, gamma = Inf, alpha = Inf), 
                   control = nls.lm.control(maxiter = 1000)), silent = FALSE)
  
  if (inherits(fit, "try-error") || any(is.na(coef(fit)))) { 
    message("Generalized Gamma model fitting failed.") 
    return(NULL) 
    } else { 
      beta <- coef(fit)["beta"] 
      gamma <- coef(fit)["gamma"] 
      alpha <- coef(fit)["alpha"] 
    message("Generalized Gamma model fit successful.") 
    return(list(fit = fit, params = list(beta = as.numeric(beta), gamma = as.numeric(gamma), alpha = as.numeric(alpha)), rss = sum(residuals(fit)^2))) } }
      

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

# Function to create and save model fit plot
save_model_fit_plot <- function(data, predictions, dataset_name, model_name, folder_name) {
  # Set the x variable based on the model type
  if (model_name == "gompertz_survival_rawtime") {
    x_var <- data$Time_hours
    x_label <- "Time (Hours)"
  } else {
    x_var <- data$Time_std
    x_label <- "Standardized Time"
  }
  
  # Create a prediction data frame using the correct x-axis variable
  prediction_df <- data.frame(x_var = seq(min(x_var), max(x_var), length.out = 100), Predictions = predictions)
  
  # Generate the plot
  plot <- ggplot() +
    geom_point(data = data, aes(x = x_var, y = Standardized_Survival), color = "black", size = 3, shape = 1) +
    geom_line(data = prediction_df, aes(x = x_var, y = Predictions), color = "firebrick") +
    labs(title = paste(model_name, "Model Fit for Dataset:", dataset_name), x = x_label, y = "Survival") +
    theme_bw()
  
  # Create the folder and save the plot
  dir.create(folder_name, showWarnings = FALSE)
  plot_filename <- paste0(folder_name, "/model_fit_", dataset_name, ".png")
  ggsave(filename = plot_filename, plot = plot, width = 8, height = 6)
}

# Initialize storage for results and predictions
all_results <- data.frame(Dataset = character(), Model = character(), AICc = numeric(),
                          lambda = numeric(), k = numeric(), a = numeric(), b = numeric(),
                          alpha = numeric(), beta = numeric(), gamma = numeric(), LT50 = numeric(),
                          Log_MT_Linearity = numeric(), Num_Obs = integer(), stringsAsFactors = FALSE)
all_predictions <- data.frame(Full_key = character(), Model = character(), Time = numeric(), Predictions = numeric())

unique_datasets <- standardized_final_df$Full_key

# Main loop to process each dataset with updated checks
for (dataset_name in unique_datasets) {
  message(paste("Processing dataset:", dataset_name))
  dataset <- standardized_final_df %>%
    filter(Full_key == dataset_name) %>%
    dplyr::select(Time_std, Standardized_Survival, Time_hours)
  
  # Fit models and store in a list
  model_fits <- list(
    exponential = fit_exponential(dataset),
    weibull = fit_weibull(dataset),
    gompertz_survival = fit_gompertz_survival(dataset),
    gompertz_survival_rawtime = fit_gompertz_survival_rawtime(dataset),
    loglogistic = fit_loglogistic(dataset),
    generalizedgamma = fit_generalizedgamma(dataset)
  )
  
  # Process each model
  for (model_name in names(model_fits)) {
    fit <- model_fits[[model_name]]
    
    # Skip if model fit is NULL
    if (is.null(fit) || is.null(fit$params)) {
      message(paste("Skipping model", model_name, "for dataset", dataset_name, "due to missing parameters."))
      next
    }
    
    # Extract parameters for each model to save for post-hoc analysis
    params <- fit$params
    lambda <- if ("lambda" %in% names(params)) as.numeric(params["lambda"]) else NA
    k <- if ("k" %in% names(params)) as.numeric(params["k"]) else NA
    a <- if ("a" %in% names(params)) as.numeric(params["a"]) else NA
    b <- if ("b" %in% names(params)) as.numeric(params["b"]) else NA
    alpha <- if ("alpha" %in% names(params)) as.numeric(params["alpha"]) else NA
    beta <- if ("beta" %in% names(params)) as.numeric(params["beta"]) else NA
    gamma <- if ("gamma" %in% names(params)) as.numeric(params["gamma"]) else NA
    lt50 <- if (!is.null(fit$LT50)) fit$LT50 else NA
    
    # Calculate AICc
    aicc_value <- if (!is.null(fit$fit)) calculate_aicc(fit$fit, dataset) else NA
    
    # Generate predictions for each model based on the fitted parameters
    time_points <- if (model_name == "gompertz_survival_rawtime") dataset$Time_hours else dataset$Time_std
    prediction_times <- seq(min(time_points), max(time_points), length.out = 100)
    
    # Generate predictions based on the model type
    if (model_name %in% c("gompertz_survival", "gompertz_survival_rawtime")) {
      # Check if a and b are numeric and non-NA before calculating predictions
      if (is.numeric(a) && is.numeric(b) && !is.na(a) && !is.na(b)) {
        preds <- exp(-(a / b) * (exp(b * prediction_times) - 1))
      } else {
        message(paste("Skipping predictions for model", model_name, "in dataset", dataset_name, "due to non-numeric a or b."))
        next
      }
    } else if (!is.null(fit$fit)) {
      preds <- try(predict(fit$fit, newdata = data.frame(Time_std = prediction_times)), silent = TRUE)
      if (inherits(preds, "try-error")) {
        message(paste("Prediction failed for model", model_name, "in dataset", dataset_name))
        next
      }
    } else {
      next
    }
    
    # Save predictions to all_predictions dataframe
    predictions <- data.frame(Full_key = dataset_name, Model = model_name, Time = prediction_times, Predictions = preds)
    all_predictions <- rbind(all_predictions, predictions)
    
    # Save plot for each model
    save_model_fit_plot(dataset, preds, dataset_name, model_name, paste0(model_name, "_model_plots"))
    
    # Save model parameters to all_results dataframe
    result <- data.frame(
      Dataset = dataset_name,
      Model = model_name,
      AICc = aicc_value,
      lambda = lambda,
      k = k,
      a = a,
      b = b,
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      LT50 = lt50,
      Log_MT_Linearity = NA,  # Placeholder for linearity check
      Num_Obs = nrow(dataset)
    )
    
    all_results <- rbind(all_results, result)
  }
}

# Save all results and predictions to files
# These large files won't be pushed to GitHub
write.csv(all_results, "data/results_df.csv", row.names = FALSE)
write.csv(all_predictions, "data/all_predictions.csv", row.names = FALSE)

##############################################################################################
# READ SAVED DATA #
standardized_final_df <- read.table("data/standardized_final_df.csv", header = T, sep = ",", dec = ".")
all_results <- read.table("data/results_df.csv", header = T, sep = ",", dec = ".")
all_predictions <- read.table("data/all_predictions.csv", header = T, sep = ",", dec = ".")

#all_results <- read.table("~/GitHub/metaAnalysis_largeFiles/results_predictions/results_df.csv", header = T, sep = ",", dec = ".")
#all_predictions <- read.table("~/GitHub/metaAnalysis_largeFiles/results_predictions/all_predictions.csv", header = T, sep = ",", dec = ".")

###############################################################################################
# FIGURE 1 & FIGURE S1 #
###############################################################################################
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
  filter(!Model == "gompertz_survival_rawtime") %>%
  arrange(combined_factor)

fig1_dat1 <- plotting_data %>%
  filter(Full_key == "8C5XBCN2A.csv")
fig1_dat2 <- prediction_plotting_data %>%
  filter(Full_key == "8C5XBCN2A.csv") %>%
  filter(!Model == "gompertz_survival_rawtime")

fig1_dat2$Model <- factor(fig1_dat2$Model, 
                          levels = c("exponential", "gompertz_survival", "weibull", 
                                     "loglogistic", "generalizedgamma"))

aic <- all_results %>% filter(Dataset == "8C5XBCN2A.csv")

fig1 <- ggplot(fig1_dat1, aes(x = Time_std, y = Standardized_Survival)) +
  geom_point(size = 4, color = "grey25") + # Raw data
  geom_line(data = fig1_dat2, aes(x = Time, y = Predictions, color = Model), size = 1) + # Predictions
  labs(x = "Standardized time", y = "Standardized survival") +
  mytheme +
  theme(legend.title = element_text(size = 16))+
  scale_color_manual(name = "AICc", values = cbpalette, labels = c("Constant: -1.03", "Gompertz: -27.31", "Weibull: -21.61", "Log-logistic: -12.41", "Gen.-gamma: -16.09"))+
  theme(
    legend.position = c(0.045, 0.045), # Positioning the legend inside the plot
    legend.justification = c(0, 0) # Aligning the legend to the bottom left
  )+coord_fixed(ratio = 1)
  

ggsave("figures/fig1.pdf", plot = fig1, width = 6.5, height = 5.5, units = "in", dpi = 300)

# Plotting - Save to model_plots folder
prediction_plotting_data <- prediction_plotting_data %>%
  filter(!Model == "gompertz_survival_rawtime")
prediction_plotting_data$Model <- factor(prediction_plotting_data$Model, 
                                         levels = c("exponential", "gompertz_survival", "weibull", 
                                                    "loglogistic", "generalizedgamma"))

ALL <- ggplot(plotting_data, aes(x = Time_std, y = Standardized_Survival)) +
  geom_point(size = 3, color = "grey25") +  # Raw data
  geom_line(data = prediction_plotting_data, aes(y = Predictions, x = Time, color = Model), size = 0.7) +  # Predictions
  labs(x = "Standardized time", y = "Standardized survival") +
  facet_wrap(~ combined_factor + key) +
  scale_color_manual(values = cbpalette, labels = c("Constant","Gompertz", "Weibull", "Log-logistic", "Generalized-gamma"))+
  mytheme+
  theme(legend.position = "bottom")

# Pagination parameters
ncol <- 5
nrow <- 5
n_pages <- ceiling(length(unique(standardized_final_df$Full_key)) / (ncol * nrow))

# Loop to create paginated plots
for (i in 1:n_pages) {
  paginated_plot <- ALL + facet_wrap_paginate(~ combined_factor + key, ncol = ncol, nrow = nrow, page = i)
  
  # Save each page to a file
  ggsave(paste0("combined_plot_page_", i, ".png"), paginated_plot, width = 16, height = 17)
  
  # Print each paginated plot
  print(paginated_plot)
}

# FIGURE 2A #
freq <- all_results %>%
  distinct(Dataset, Model, .keep_all = T) %>%
  group_by(Dataset) %>%
  filter(AICc == min(AICc)) %>%
  group_by(Model) %>%
  summarise(n = n()) 

# Rank all models for each group by Mean_WAIC
ranked_results <- all_results %>%
 filter(!Model == "<NA>") %>%
  distinct(Dataset, Model, .keep_all = T) %>%
  group_by(Dataset) %>%
  arrange(AICc) %>%
  mutate(rank = row_number()) %>%
  mutate(best_model = Model[rank == 1])
  
# Count how many times each model ends up as the best after adjustments
model_counts <- ranked_results %>%
  distinct(Dataset, .keep_all = T) %>%
  group_by(best_model) %>%
  summarise(n = n(), .groups = 'drop')

# Density plot for Mean_AIC by Model
# Plot with ggplot2
density2 <- ggplot(data = all_results %>% filter(Model != "gompertz_survival_rawtime") %>% 
                     distinct(Dataset, Model, .keep_all = T), aes(x = AICc, color = Model)) +
  geom_density(alpha = 0.5, size = 1) +
  labs(
    x = "AICc",
    y = "Density") +
  mytheme +
  scale_color_manual(values = cbpalette, labels = c("Constant (25)", "Gompertz (70)", "Weibull (45)", "Log-logistic (61)", "Gen.-gamma (8)")) +
  scale_x_continuous(limits = c(-200, 100)) +
  theme(legend.position = c(0.1, 0.5), 
        legend.justification = c(0, 0),
        legend.key = element_blank()) +
  guides(color = guide_legend(override.aes = list(linetype = 1)))

# Data
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
  "ASXFW6GM_1.csv", "ASXFW6GM_2.csv", "ASXFW6GM_3.csv", "ASXFW6GM_4.csv", "SLELFBN3_1.csv", "SLELFBN3_2.csv", "SLELFBN3_3.csv", "SLELFBN3_4.csv", "SLELFBN3_5.csv",
  "SLELFBN3_6.csv", "SLELFBN3_7.csv", "SLELFBN3_8.csv", "SLELFBN3_9.csv", "SLELFBN3_10.csv", "SLELFBN3_11.csv", "SLELFBN3_12.csv", "SLELFBN3_13.csv", "SLELFBN3_14.csv", "SLELFBN3_15.csv",
  "SLELFBN3_16.csv", "SLELFBN3_17.csv", "SLELFBN3_18.csv", "SLELFBN3_19.csv", "SLELFBN3_20.csv", "SLELFBN3_21.csv", "SLELFBN3_22.csv", "SLELFBN3_23.csv", "SLELFBN3_24.csv", "SLELFBN3_25.csv",
  "SLELFBN3_26.csv", "SLELFBN3_27.csv", "SLELFBN3_28.csv", "SLELFBN3_29.csv", "SLELFBN3_30.csv", "FRAQNTF6_1.csv", "FRAQNTF6_2.csv", "FRAQNTF6_3.csv", "FRAQNTF6_4.csv",
  "NSH3NLMI_1.csv", "NSH3NLMI_2.csv", "72SE5ARY_1.csv", "72SE5ARY_2.csv"
)

num_hosts <- c(
  24, 18, 8, 9, 24, 28, 23, 21, 23, 27, 25, 24, 29, 29, 27, 22, 26, 27, 25, 20, 21, 26, 23,  9,  9, 10, 
  25, 23, 20, 15, 14, 14, 11, 12,  20, 20, 19, 19,  10, 12, 10, 20
)

# Create mergedata data frame
mergedata <- data.frame(Dataset = datasets, Max_num_hosts = num_hosts)
 
allData_AIC_filled <- all_results %>%
  left_join(standardized_final_df, by = c("Dataset" = "Full_key")) %>%
  left_join(mergedata, by = "Dataset") %>%  
  mutate(Max_num_host_filled = coalesce(as.numeric(Max_num_host), as.numeric(Max_num_hosts))) %>%
  dplyr::select(Dataset, Model, AICc, Num_Obs, Host_taxa, Pathogen_taxa, Data, Max_num_host_filled)

rows_with_na <- apply(allData_AIC_filled, 1, function(x) any(is.na(x)))
allData_AIC_filled[rows_with_na, ]
  
###############################################################################################
# FIGURE S2 & S3 #
###############################################################################################

# HOST
allData_AIC_sum <- allData_AIC_filled %>%
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
    numData = length(unique(Dataset))
  ) %>%
  ungroup() %>%
  filter(!Model == "gompertz_survival_rawtime") 

# Order factor levels
allData_AIC_sum$Host_taxa <- factor(allData_AIC_sum$Host_taxa, 
                                     levels = c("Seedlings", "Drosophila sp.", "Other insects", 
                                                "Nematodes", "Moth larvae", "Other invertebrates", 
                                                "Fish", "Avian", "Mice", "Other mammals"))

allData_AIC_sum$Model <- factor(allData_AIC_sum$Model, 
                                 levels = c("exponential", "gompertz_survival", "weibull", 
                                            "loglogistic", "generalizedgamma"))

# Create the plot and add the annotations
host_aic <- ggplot(allData_AIC_sum, aes(x = Model, y = meanAIC, color = Model)) +
  geom_point(size = 5, shape = 21) +
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

ggsave("figures/host_aic.pdf", plot = host_aic, width = 6, height = 12, units = "in", dpi = 300)

# PATGOGEN
# Order factor levels
PallData_AIC_sum <- allData_AIC_filled %>%
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
    numData = length(unique(Dataset))
  ) %>%
  filter(!Model == "gompertz_survival_rawtime") 

PallData_AIC_sum$Pathogen_taxa <- factor(PallData_AIC_sum$Pathogen_taxa, 
                                          levels = c("Gram-positive bacteria", "Gram-negative bacteria",
                                                     "DNA virus", "RNA virus", "Fungi", "Protozoan parasite"))


PallData_AIC_sum$Model <- factor(PallData_AIC_sum$Model, 
                                  levels = c("exponential", "gompertz_survival", "weibull", 
                                             "loglogistic", "generalizedgamma"))



# Create the plot and add the annotations
pathogen_aic <- ggplot(PallData_AIC_sum, aes(x = Model, y = meanAIC, color = Model)) +
  geom_point(size = 5, shape = 21) +
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

ggsave("figures/pathogen_aic.pdf", plot = pathogen_aic, width = 6, height = 8, units = "in", , dpi = 300)

#################################################################################################
# DATA TYPE #
#################################################################################################

Data_AIC <- allData_AIC_filled %>%
  filter(
    is.finite(as.numeric(AICc)),
  ) %>%
  group_by(Model, Data) %>%
  summarize(meanAIC = mean(AICc, na.rm = TRUE), seAIC = standard_error(AICc)) %>%
  filter(!Model == "gompertz_survival_rawtime")

Data_AIC$Model <- factor(Data_AIC$Model, 
                          levels = c("exponential", "gompertz_survival", "weibull", 
                                     "loglogistic", "generalizedgamma"))

Datatype <- ggplot(Data_AIC, aes(x = Data, y = meanAIC, color = Model)) +
  geom_point(size = 5, shape = 21, position = position_dodge(width = 0.5)) +
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

# Num. of Observations
allData_AIC_sum_obs <- allData_AIC_filled %>%
  group_by(Dataset, Model) %>%
  summarize(
    meanObs = mean(as.numeric(Num_Obs), na.rm = T),
    meanAIC = mean(AICc, na.rm = T)
  ) %>%
  ungroup()

allData_AIC_sum_num <- allData_AIC_filled %>%
  group_by(Dataset, Model) %>%
  summarize(
    meanNum = mean(as.numeric(Max_num_host_filled), na.rm = T),
    meanAIC = mean(AICc, na.rm = T)
  ) %>%
  ungroup()

################################################################################################
# FIGURE 2B #
################################################################################################

# Delta Model i - Exponential for other models 
allModels_delta <- all_results %>%
  filter(Model != "gompertz_survival_rawtime") %>%
  distinct(Dataset, Model, .keep_all = T) %>%
  select(Dataset, Model, AICc) %>%
  spread(Model, AICc) %>%
  group_by(Dataset) %>%
  summarise(exponential = mean(exponential, na.rm = TRUE), gompertz_survival = mean(gompertz_survival, na.rm =T), 
            generalizedgamma = mean(generalizedgamma, na.rm = TRUE), loglogistic = mean(loglogistic, na.rm = TRUE),
            weibull = mean(weibull, na.rm = TRUE)) %>%
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

allData_AIC_filled$Model <- factor(allData_AIC_filled$Model, 
                                    levels = c("exponential", "gompertz_survival", "gompertz_hazard", "weibull", 
                                               "loglogistic", "generalizedgamma"))

Comp_gomp <- ggplot(all_results %>% filter(Model != "<NA>"), 
                    aes(y = AICc, x = Model))+
  geom_jitter(size = 3, shape = 21, color = "grey50")+
  geom_boxplot(fill = "white", , alpha  = 0.7)+
  mytheme+
  labs(y = "AICc", x = NULL)+
  scale_x_discrete(labels = c("Constant", "Gompertz", "Weibull", "Log-logistic", "Gen.-gamma"))

# Combine the plots and add labels
op <- op + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Manuscript Figure 2
combined_plot_1 <- density2 + op +
  plot_annotation(tag_levels = 'A')

ggsave("figures/combined_plot_fig2.pdf", plot = combined_plot_1, width = 10, height = 5, units = "in", dpi = 300)

summary(lm(DeltaAIC~comparison, data = allModels_delta)) 
anova(lm(DeltaAIC~comparison, data = allModels_delta)) #F(3/736)=1.1332 p=0.3348

#summary(lm(DeltaAIC~comparison, data = allModels_delta_g)) #F(3/736)=94.578 < 2.2e-16 ***
#anova(lm(DeltaAIC~comparison, data = allModels_delta_g))

mod_all <- all_results %>%
  filter_all(all_vars(!is.infinite(.)))%>%
  distinct(Dataset, Model, .keep_all = T)

summary(lm(AICc~Model, data = mod_all)) 
anova(lm(AICc~Model, data = mod_all)) #F(4/1018) = 22.55 p = < 2.2e-16 ***

###############################################################################################

###############################################
# Percentage of accelerating mortality # 
# Won't be included in the manuscript figure #
###############################################

results_with_classification <- all_results %>%
  filter(!Model == "gompertz_survival_rawtime") %>%
  group_by(Dataset) %>%
  filter(AICc == min(AICc)) %>%
  mutate(
    Mortality_Type = case_when(
      Model == "exponential" ~ "Constant",
      Model == "gompertz_survival" & b > 0 ~ "Accelerating",
      Model == "weibull" & k > 1 ~ "Accelerating",
      Model == "weibull" & k < 1 ~ "Decelerating",
      Model == "weibull" & k == 1 ~ "Constant",
      Model == "loglogistic" & beta > 1 ~ "Accelerating",
      Model == "loglogistic" & beta < 1 ~ "Decelerating",
      Model == "loglogistic" & beta == 1 ~ "Constant",
      Model == "generalizedgamma" & gamma > 1 ~ "Accelerating",
      Model == "generalizedgamma" & gamma < 1 ~ "Decelerating",
      Model == "generalizedgamma" & gamma == 1 ~ "Constant",
      TRUE ~ "Undetermined"
    )
  ) %>%
  ungroup() %>% 
  distinct(Dataset, .keep_all = T)

results_with_classification %>%
  group_by(Mortality_Type) %>%
  summarise(count = n()) 
#86% accelerating %12 constant %2 decelerating

altC <- results_with_classification %>%
  group_by(Mortality_Type, Model) %>%
  summarise(count = n()) %>%
  filter(Model == "gompertz_survival" | Model == "weibull")
  
altCplot <- ggplot(altC, aes(y = Mortality_Type, x = count, fill = Model))+
  geom_bar(stat = "identity", color = "grey25")+
  mytheme+
  xlab("Number of datasets")+
  ylab("m(t)")+
  scale_fill_manual(values = c("#0072B2", "#D55E00"),labels = c("Gompertz", "Weibull"))+
  theme(legend.position = c(0.8, 0.8))

results_with_classification_all <- all_results %>%
  distinct(Dataset, Model, .keep_all = T) %>%
  filter(!Model == "gompertz_survival_rawtime") %>%
  group_by(Dataset) %>%
  filter(AICc == min(AICc)) %>%
  ungroup() %>%
  mutate(
    Mortality_Type = case_when(
      Model == "exponential" ~ "Constant",
      Model == "gompertz_survival" ~ "Other",
      Model == "weibull" ~ "Other",
      Model == "loglogistic" ~ "Other",
      Model == "generalizedgamma" ~ "Other",
    )
  ) %>%
  select(Dataset, Mortality_Type) 

# Inspection 

results_with_classification_gom <- all_results %>%
  filter(!Model == "gompertz_survival_rawtime") %>%
  group_by(Dataset) %>%
  mutate(
    Mortality_Type = case_when(
      Model == "exponential" ~ "Constant",
      Model == "gompertz_survival" & b > 0 ~ "Accelerating",
      Model == "gompertz_survival" & b < 0 ~ "Decelerating",
      Model == "gompertz_survival" & b == 0 ~ "Constant",
      Model == "weibull" & k > 1 ~ "Accelerating",
      Model == "weibull" & k < 1 ~ "Decelerating",
      Model == "weibull" & k == 1 ~ "Constant",
      Model == "loglogistic" & beta > 1 ~ "Accelerating",
      Model == "loglogistic" & beta < 1 ~ "Decelerating",
      Model == "loglogistic" & beta == 1 ~ "Constant",
      Model == "generalizedgamma" & gamma > 1 ~ "Accelerating",
      Model == "generalizedgamma" & gamma < 1 ~ "Decelerating",
      Model == "generalizedgamma" & gamma == 1 ~ "Constant",
      TRUE ~ "Undetermined"
    )
  ) %>%
  ungroup() %>%
  distinct(Dataset, Model, .keep_all = T) %>%
  filter(Model == "gompertz_survival" | Model == "weibull") %>%
  filter(Mortality_Type == "Decelerating") %>%
  select(Dataset, Model, Mortality_Type, b, k) 


results_with_classification_constant <- all_results %>%
  filter(!Model == "gompertz_survival_rawtime") %>%
  group_by(Dataset) %>%
  filter(AICc == min(AICc)) %>%
  mutate(
    Mortality_Type = case_when(
      Model == "exponential" ~ "Constant",
      Model == "gompertz_survival" & b > 0 ~ "Accelerating",
      Model == "weibull" & k > 1 ~ "Accelerating",
      Model == "weibull" & k < 1 ~ "Decelerating",
      Model == "weibull" & k == 1 ~ "Constant",
      Model == "loglogistic" & beta > 1 ~ "Accelerating",
      Model == "loglogistic" & beta < 1 ~ "Decelerating",
      Model == "loglogistic" & beta == 1 ~ "Constant",
      Model == "generalizedgamma" & gamma > 1 ~ "Accelerating",
      Model == "generalizedgamma" & gamma < 1 ~ "Decelerating",
      Model == "generalizedgamma" & gamma == 1 ~ "Constant",
      TRUE ~ "Undetermined"
    )
  ) %>%
  ungroup() %>% 
  distinct(Dataset, .keep_all = T) %>% 
  filter(Mortality_Type == "Constant")

################################################################################################
# Gompertz Model Accuracy across taxa groups #
# FIGURE 3c#
################################################################################################

palette_10 <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", 
                "#ffd92f", "#e5c494", "#b3b3b3", "#a6761d", "#1f78b4")

palette_6 <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")

allModels_delta_gomp <- allModels_delta  %>%
  filter(comparison == "DeltaGompertz") %>%
  select(Dataset, DeltaAIC)

allModels_delta_gomp_for_fig <- allData_AIC_filled %>%
  filter(Model == "gompertz_survival") %>%
  distinct(Dataset, .keep_all = T) %>%
  left_join(allModels_delta_gomp, by = "Dataset") %>%
  filter_if(is.numeric, all_vars(is.finite(.)))

sum_host_AIC <- allData_AIC_filled %>%
  filter(Model == "gompertz_survival") %>%
  distinct(Dataset, .keep_all = T) %>%
  left_join(allModels_delta_gomp, by = "Dataset") %>%
  group_by(Host_taxa) %>%
  summarize(Mean_AIC = mean(AICc, na.rm = T), 
            se_AIC = standard_error(AICc), 
            Mean_delta_AIC = mean(DeltaAIC, na.rm = T),
            se_delta_AIC = standard_error(DeltaAIC),
            Mean_Max_num_host_filled  = mean(Max_num_host_filled, na.rm = T),
            se_Max_num_host_filled = standard_error(Max_num_host_filled),
            Mean_Num_Obs = mean(Num_Obs, na.rm = T),
            se_Mean_Num_Obs = standard_error(Num_Obs),
            Count_Dataset = n(), 
            Percent_Raw = sum(Data == "raw") / n() * 100,  # Percentage of 'raw' data
            Percent_Probability = sum(Data == "probability") / n() * 100  # Percentage of 'probability' data
  ) 

taxa_type <- allData_AIC_filled %>% 
  left_join(results_with_classification[,c(1,14)], by = "Dataset") %>%
  distinct(Dataset, .keep_all = T) %>%
  group_by(Host_taxa, Mortality_Type) %>%
  summarise(count = n()) %>%
mutate(total_count = sum(count)) %>%  # Calculate total count for each Host_taxa
  mutate(percentage = (count / total_count) * 100) %>%  # Calculate percentage
  select(-total_count)  # Optionally remove the total_count column

taxa_type$Host_taxa <- factor(taxa_type$Host_taxa, 
                                 levels = rev(c("Seedlings",  
                                                "Nematodes", "Other invertebrates", "Drosophila sp.", "Moth larvae", "Other insects",
                                                "Avian", "Mice", "Other mammals", "Fish")))

# Now, create the stacked bar plot
ttype <- ggplot(taxa_type, aes(y = Host_taxa, x = percentage, fill = Mortality_Type)) +
  geom_bar(stat = "identity", color = "grey25") +
  labs(x = NULL, y = NULL, fill = "Mortality Type") +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = cbpalette)+
  theme(legend.position = "bottom")


sum_host_AIC$Host_taxa <- factor(sum_host_AIC$Host_taxa, 
                                           levels = rev(c("Seedlings",  
                                                          "Nematodes", "Other invertebrates", "Drosophila sp.", "Moth larvae", "Other insects",
                                                          "Avian", "Mice", "Other mammals", "Fish")))


aicplot <- ggplot(sum_host_AIC, aes(x = Mean_delta_AIC, y = Host_taxa))+
  geom_errorbar(aes(xmin = Mean_delta_AIC - se_delta_AIC, xmax = Mean_delta_AIC + se_delta_AIC), width = 0.2)+
  geom_point(shape = 21, size = 5) +
  mytheme+
  theme(legend.position = "none")+
  labs(y = NULL, x = expression(Delta[CG]))+
  scale_x_continuous(limits = c(0,75))

standardized_final_df %>%
  distinct(Full_key, .keep_all = T) %>%
  group_by(Host_taxa)%>%
  summarise(count = n())


allModels_delta_gomp$Host_taxa <- factor(allModels_delta_gomp$Host_taxa, 
                                 levels = rev(c("Seedlings",  
                                                "Nematodes", "Other invertebrates", "Drosophila sp.", "Moth larvae", "Other insects",
                                                "Avian", "Mice", "Other mammals", "Fish")))
  
otherplots1 <- ggplot(sum_host_AIC, aes(x = Count_Dataset, y = Host_taxa, fill = Host_taxa)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  mytheme +
  theme(legend.position = "none") +
  labs(y = NULL, x = "Num.Data") +
  scale_fill_manual(values = palette_10)+
  theme(legend.position = "none")

otherplots2 <- ggplot(sum_host_AIC, aes(x = Mean_Num_Obs, y = Host_taxa, fill = Host_taxa)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_errorbar(aes(xmin = Mean_Num_Obs - se_Mean_Num_Obs, xmax = Mean_Num_Obs + se_Mean_Num_Obs), width = 0.2)+
  mytheme +
  theme(legend.position = "none") +
  labs(y = NULL, x = "Num.Obs.") +
  scale_fill_manual(values = palette_10)+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())

otherplots3 <- ggplot(sum_host_AIC, aes(x = Mean_Max_num_host_filled, y = Host_taxa, fill = Host_taxa)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_errorbar(aes(xmin = Mean_Max_num_host_filled - se_Max_num_host_filled, xmax = Mean_Max_num_host_filled + se_Max_num_host_filled), width = 0.2)+
  mytheme +
  theme(legend.position = "none") +
  labs(y = NULL, x = "Num.Host") +
  scale_fill_manual(values = palette_10)+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())

# Reshape the data to long format for Percent_Raw and Percent_Probability
sum_host_AIC_long_percent <- sum_host_AIC %>%
  select(Host_taxa, Percent_Raw, Percent_Probability) %>%
  pivot_longer(cols = c(Percent_Raw, Percent_Probability), names_to = "Data_Type", values_to = "Percentage")

# Plot Mean_AIC and the stacked bar plot for percentages
datatype <- ggplot(sum_host_AIC, aes(x = Mean_AIC, y = Host_taxa, fill = Host_taxa)) +
  # Stacked bar plot for Percent_Raw and Percent_Probability
  geom_bar(data = sum_host_AIC_long_percent, aes(x = Percentage, y = Host_taxa, fill = Data_Type), 
           stat = "identity", position = "stack", width = 0.5, inherit.aes = FALSE, color = "black") +
  
  scale_fill_manual(values = c("Percent_Raw" = "grey20", "Percent_Probability" = "grey90", "Host_taxa" = palette_10)) +
  mytheme +
  theme(legend.position = "right") +
  labs(x = "Data Type", y = NULL)+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())

# Combine the plots, aligning them with a shared y-axis on the first plot
final_plot <-  otherplots1 + otherplots2+ otherplots3 + datatype + plot_layout(ncol = 4, widths = c(1, 1, 1, 1))

# All for pathogen #

sum_host_AIC_p <- allData_AIC_filled %>%
  filter(Model == "gompertz_survival") %>%
  left_join(allModels_delta_gomp, by = "Dataset") %>%
  distinct(Dataset, .keep_all = T) %>%
  group_by(Pathogen_taxa) %>%
  summarize(Mean_AIC = mean(AICc, na.rm = T), 
            se_AIC = standard_error(AICc), 
            Mean_delta_AIC = mean(DeltaAIC, na.rm = T),
            se_delta_AIC = standard_error(DeltaAIC),
            Mean_Max_num_host_filled  = mean(Max_num_host_filled, na.rm = T),
            se_Max_num_host_filled = standard_error(Max_num_host_filled),
            Mean_Num_Obs = mean(Num_Obs, na.rm = T),
            se_Mean_Num_Obs = standard_error(Num_Obs),
            Count_Dataset = n(), 
            Percent_Raw = sum(Data == "raw") / n() * 100,  # Percentage of 'raw' data
            Percent_Probability = sum(Data == "probability") / n() * 100  # Percentage of 'probability' data
  )

taxa_typep <- allData_AIC_filled %>% 
  left_join(results_with_classification[,c(1,14)], by = "Dataset") %>%
  distinct(Dataset, .keep_all = T) %>%
  group_by(Pathogen_taxa, Mortality_Type) %>%
  summarise(count = n()) %>%
  mutate(total_count = sum(count)) %>%  # Calculate total count for each Host_taxa
  mutate(percentage = (count / total_count) * 100) %>%  # Calculate percentage
  select(-total_count)  # Optionally remove the total_count column


taxa_typep$Pathogen_taxa <- factor(taxa_typep$Pathogen_taxa, 
                                       levels = rev(c("Protozoan parasite",  
                                                      "Fungi", "Gram-negative bacteria", "Gram-positive bacteria", 
                                                      "DNA virus", "RNA virus")))
# Now, create the stacked bar plot
ttype_p <- ggplot(taxa_typep, aes(y = Pathogen_taxa, x = percentage, fill = Mortality_Type)) +
  geom_bar(stat = "identity", color = "grey25") +
  labs(y = NULL, x = NULL, fill = "m(t)") +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = cbpalette)+
  theme(legend.position = "bottom")
  

sum_host_AIC_p$Pathogen_taxa <- factor(sum_host_AIC_p$Pathogen_taxa, 
                                       levels = rev(c("Protozoan parasite",  
                                                      "Fungi", "Gram-negative bacteria", "Gram-positive bacteria", 
                                                      "DNA virus", "RNA virus")))

aicplot_p <- ggplot(sum_host_AIC_p, aes(x = Mean_delta_AIC, y = Pathogen_taxa))+
  geom_errorbar(aes(xmin = Mean_delta_AIC - se_delta_AIC, xmax = Mean_delta_AIC + se_delta_AIC), width = 0.2)+
  geom_point(shape = 21, size = 5) +
  mytheme+
  theme(legend.position = "none")+
  labs(y = NULL, x = expression(Delta[CG]))+
  scale_x_continuous(limits = c(0,75))

standardized_final_df %>%
  distinct(Full_key, .keep_all = T) %>%
  group_by(Pathogen_taxa)%>%
  summarise(count = n())

otherplots1_p <- ggplot(sum_host_AIC_p, aes(x = Count_Dataset, y = Pathogen_taxa, fill = Pathogen_taxa)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  mytheme +
  theme(legend.position = "none") +
  labs(y = NULL, x = "Num.Data") +
  scale_fill_manual(values = palette_6)+
  theme(legend.position = "none")

otherplots2_p <- ggplot(sum_host_AIC_p, aes(x = Mean_Num_Obs, y = Pathogen_taxa, fill = Pathogen_taxa)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_errorbar(aes(xmin = Mean_Num_Obs - se_Mean_Num_Obs, xmax = Mean_Num_Obs + se_Mean_Num_Obs), width = 0.2)+
  mytheme +
  theme(legend.position = "none") +
  labs(y = NULL, x = "Num.Obs.") +
  scale_fill_manual(values = palette_6)+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())

otherplots3_p <- ggplot(sum_host_AIC_p, aes(x = Mean_Max_num_host_filled, y = Pathogen_taxa, fill = Pathogen_taxa)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_errorbar(aes(xmin = Mean_Max_num_host_filled - se_Max_num_host_filled, xmax = Mean_Max_num_host_filled + se_Max_num_host_filled), width = 0.2)+
  mytheme +
  theme(legend.position = "none") +
  labs(y = NULL, x = "Num.Host") +
  scale_fill_manual(values = palette_6)+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())

# Reshape the data to long format for Percent_Raw and Percent_Probability
sum_pat_AIC_long_percent <- sum_host_AIC_p %>%
  select(Pathogen_taxa, Percent_Raw, Percent_Probability) %>%
  pivot_longer(cols = c(Percent_Raw, Percent_Probability), names_to = "Data_Type", values_to = "Percentage")

# Plot Mean_AIC and the stacked bar plot for percentages
datatype_p <- ggplot(sum_pat_AIC_long_percent, aes(x = Mean_AIC, y = Pathogen_taxa, fill = Pathogen_taxa)) +
  # Stacked bar plot for Percent_Raw and Percent_Probability
  geom_bar(data = sum_pat_AIC_long_percent, aes(x = Percentage, y = Pathogen_taxa, fill = Data_Type), 
           stat = "identity", position = "stack", width = 0.5, inherit.aes = FALSE, color = "black") +
  
  scale_fill_manual(values = c("Percent_Raw" = "grey20", "Percent_Probability" = "grey90", "Host_taxa" = palette_10)) +
  mytheme +
  theme(legend.position = "right") +
  labs(x = "Data Type", y = NULL)+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())

# Combine the plots, aligning them with a shared y-axis on the first plot

# Figure S4
final_plot_p <- otherplots1 + otherplots2 + otherplots3 + datatype + otherplots1_p + otherplots2_p + otherplots3_p + datatype_p + plot_layout(nrow = 2, ncol = 4, widths = c(2, 2, 2, 2))

# Figure 3
alt_final_plot1 <- aicplot + aicplot_p + plot_layout(ncol = 2, widths = c(1, 1))+ plot_annotation(tag_levels = 'A')

#Stats for Figure 3
dataForStats <- allData_AIC_filled %>%
  filter(Model == "gompertz_survival") %>%
  left_join(allModels_delta_gomp, by = "Dataset") %>%
  filter_if(is.numeric, all_vars(is.finite(.))) %>%
  distinct(Dataset, .keep_all = T) %>%
  select(Host_taxa, Pathogen_taxa, DeltaAIC)


results <- dataForStats %>%
  group_by(Host_taxa) %>%
  summarize(
    t_test_result = list(t.test(DeltaAIC, mu = 0)),
    .groups = 'drop'
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test_result, ~ .x$statistic),
    p_value = map_dbl(t_test_result, ~ .x$p.value),
    mean_deltaAIC = map_dbl(t_test_result, ~ .x$estimate),
    significance = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      p_value <= 0.1   ~ ".",
      TRUE             ~ ""
    )
  ) %>%
  select(-t_test_result)

results_pat <- dataForStats %>%
  group_by(Pathogen_taxa) %>%
  summarize(
    t_test_result = list(t.test(DeltaAIC, mu = 0)),
    .groups = 'drop'
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test_result, ~ .x$statistic),
    p_value = map_dbl(t_test_result, ~ .x$p.value),
    mean_deltaAIC = map_dbl(t_test_result, ~ .x$estimate),
    significance = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      p_value <= 0.1   ~ ".",
      TRUE             ~ ""
    )
  ) %>%
  select(-t_test_result)

###########################################################
# Constant mortality exploration #
# FIGURE 4#
###########################################################

sum_AIC_constant <- allData_AIC_filled %>%
  left_join(allModels_delta_gomp, by = "Dataset") %>%
  distinct(Dataset, .keep_all = T) %>%
  left_join(results_with_classification_all, by = "Dataset")  %>%
  group_by(Mortality_Type) %>%
  summarize(Mean_delta_AIC = mean(DeltaAIC, na.rm = T),
            se_delta_AIC = standard_error(DeltaAIC),
            Mean_Max_num_host_filled  = mean(Max_num_host_filled, na.rm = T),
            Mean_Num_Obs = mean(Num_Obs, na.rm = T),
            Count_Dataset = n(), 
            Percent_Raw = sum(Data == "raw") / n() * 100,  # Percentage of 'raw' data
            Percent_Probability = sum(Data == "probability") / n() * 100,  # Percentage of 'probability' data
            se_Mean_Num_Obs = standard_error(Num_Obs),
            se_Max_num_host_filled = standard_error(Max_num_host_filled)) 

sum_AIC_constant_numhost <- allData_AIC_filled %>%
left_join(allModels_delta_gomp, by = "Dataset") %>%
  distinct(Dataset, .keep_all = T) %>%
  left_join(results_with_classification_all, by = "Dataset")  %>%
  filter(!Max_num_host_filled == 500) %>%
  group_by(Mortality_Type) %>%
  summarize(Mean_Max_num_host_filled  = mean(Max_num_host_filled, na.rm = T),
            se_Max_num_host_filled = standard_error(Max_num_host_filled)) 

otherplots1_c <- ggplot(sum_AIC_constant, aes(y = Mean_delta_AIC, x = Mortality_Type)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_errorbar(aes(ymin = Mean_delta_AIC - se_delta_AIC, ymax = Mean_delta_AIC + se_delta_AIC), width = 0.2)+
  mytheme +
  theme(legend.position = "none") +
  labs(x = NULL, y = expression(Delta[AICc]))

otherplots2_c <- ggplot(sum_AIC_constant, aes(y = Mean_Num_Obs, x = Mortality_Type)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_errorbar(aes(ymin = Mean_Num_Obs - se_Mean_Num_Obs, ymax = Mean_Num_Obs + se_Mean_Num_Obs), width = 0.2)+
  mytheme +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Num. Observations")

otherplots3_c <- ggplot(sum_AIC_constant, aes(y = Mean_Max_num_host_filled, x = Mortality_Type)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_errorbar(aes(ymin = Mean_Max_num_host_filled - se_Max_num_host_filled, ymax = Mean_Max_num_host_filled + se_Max_num_host_filled), width = 0.2)+
  mytheme +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Num. Host") 


otherplots3_c_alt <- ggplot(sum_AIC_constant_numhost, aes(y = Mean_Max_num_host_filled, x = Mortality_Type)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_errorbar(aes(ymin = Mean_Max_num_host_filled - se_Max_num_host_filled, ymax = Mean_Max_num_host_filled + se_Max_num_host_filled), width = 0.2)+
  mytheme +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Num. Host") 


# Reshape the data to long format for Percent_Raw and Percent_Probability
sum_pat_AIC_long_percent_c <- sum_AIC_constant %>%
  select(Mortality_Type, Percent_Raw, Percent_Probability) %>%
  pivot_longer(cols = c(Percent_Raw, Percent_Probability), names_to = "Data_Type", values_to = "Percentage")

# Plot Mean_AIC and the stacked bar plot for percentages
datatype_p_c <- ggplot(sum_pat_AIC_long_percent_c, aes(y = Mean_AIC, x = Mortality_Type)) +
  # Stacked bar plot for Percent_Raw and Percent_Probability
  geom_bar(data = sum_pat_AIC_long_percent_c, aes(y = Percentage, x = Mortality_Type, fill = Data_Type), 
           stat = "identity", position = "stack", width = 0.5, inherit.aes = FALSE, color = "black") +
  
  scale_fill_manual(values = c("Percent_Raw" = "grey20", "Percent_Probability" = "grey90", "Host_taxa" = palette_10)) +
  mytheme +
  theme(legend.position = "none") +
  labs(y = "Data Type", x = NULL)

# Combine the plots, aligning them with a shared y-axis on the first plot

# Manuscript Figure 4
constant_plot     <- otherplots2_c + otherplots3_c + datatype_p_c + plot_annotation(tag_levels = 'A')
constant_plot_alt <- otherplots2_c + otherplots3_c_alt + datatype_p_c + plot_annotation(tag_levels = 'A')

# Stats 
sum_AIC_constant_stats <- allData_AIC_filled %>%
  left_join(allModels_delta_gomp, by = "Dataset") %>%
  distinct(Dataset, .keep_all = T) %>%
  left_join(results_with_classification_all, by = "Dataset") %>%
  select(Mortality_Type, Max_num_host_filled, Num_Obs)

sum_AIC_constant_stats_alt <- allData_AIC_filled %>%
  left_join(allModels_delta_gomp, by = "Dataset") %>%
  filter(!Max_num_host_filled == 500) %>%
  distinct(Dataset, .keep_all = T) %>%
  left_join(results_with_classification_all, by = "Dataset") %>%
  select(Mortality_Type, Max_num_host_filled, Num_Obs)

sum_AIC_constant_stats$Mortality_Type <- as.factor(sum_AIC_constant_stats$Mortality_Type)

t.test(Num_Obs ~ Mortality_Type, var.equal = FALSE, data=sum_AIC_constant_stats)
t.test(Max_num_host_filled ~ Mortality_Type, var.equal = FALSE, data=sum_AIC_constant_stats)
t.test(Max_num_host_filled ~ Mortality_Type, var.equal = FALSE, data=sum_AIC_constant_stats_alt)

# Stat for percentages
sum_AIC_constant_per <- allData_AIC_filled %>%
  left_join(allModels_delta_gomp, by = "Dataset") %>%
  distinct(Dataset, .keep_all = T) %>%
  left_join(results_with_classification_all, by = "Dataset") %>%
  select(Mortality_Type, Data)

# First, create a contingency table of Mortality_Type and Data
contingency_table <- table(sum_AIC_constant_per$Mortality_Type, sum_AIC_constant_per$Data)

# Perform the chi-square test of independence
chi_square_test <- chisq.test(contingency_table)