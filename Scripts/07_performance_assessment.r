###########################################################################
# Performance Assessment
#
# This script evaluates the spatio-temporal transferability of a cane toad
# classifier using logistic regression. Validated predictions (correct vs.
# incorrect) and BirdNET confidence scores are used along with site, month and
# season as covariates.
#
#
# References:
# - Wood et al., 2023b; Wood and Kahl, 2024.
#
# Dependencies: dplyr, ggplot2, patchwork, emmeans
#
# Usage: Set your working directory and update file paths as needed.
###########################################################################

## ----- Part 1: Data Combination ----- ##
setwd("path/to/your/working/directory")  

library(dplyr)

# Combine all CSV files (one per site) from a designated folder into one dataset
file_directory <- "path/to/your/data/directory"  # Folder containing CSV files
files <- list.files(path = file_directory, pattern = "\\.csv$", full.names = TRUE)
data_list <- lapply(files, read.csv)
combined_data <- bind_rows(data_list)
write.csv(combined_data, file = "Combined_Site_Data.csv", row.names = FALSE)

## ----- Part 2: Spatial Differences ----- ##
library(ggplot2)

# Read the spatial dataset (update file name as needed)
data <- read.csv("path/to/your/combined_data_file.csv", stringsAsFactors = FALSE)

# Fit a logistic regression model using Confidence and SITE as covariates
site.conf.model <- glm(MANUAL.ID ~ Confidence + SITE, data = data, family = "binomial")
summary(site.conf.model)

# Generate predictions from Confidence models for each SITE
confidence_model_list <- lapply(unique(data$SITE), function(site) {
  subset_data <- subset(data, SITE == site)
  model <- glm(MANUAL.ID ~ Confidence, data = subset_data, family = binomial)
  list(site = site, model = model)
})
confidence_predictions_list <- lapply(confidence_model_list, function(entry) {
  site <- entry$site
  model <- entry$model
  conf_vals <- seq(min(data$Confidence), max(data$Confidence), length.out = 100)
  data.frame(Confidence = conf_vals,
             Probability = predict(model, newdata = data.frame(Confidence = conf_vals), type = "response"),
             Site = site)
})
confidence_predictions_df <- do.call(rbind, confidence_predictions_list)

# Calculate overall mean predictions for Confidence across sites
overall_confidence_mean <- aggregate(Probability ~ Confidence, confidence_predictions_df, mean)

# Define custom line styles used in the Confidence plot
line_styles <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")

# Create spatial plot for Confidence predictions using custom line styles
base_size <- 12
confidence_sites_plot <- ggplot() +
  geom_line(data = confidence_predictions_df, aes(x = Confidence, y = Probability, 
                                                  color = Site, linetype = Site), size = 0.7) +
  geom_line(data = overall_confidence_mean, aes(x = Confidence, y = Probability), 
            color = "black", size = 1.0) +
  labs(x = "Confidence Score", y = "Accuracy") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_linetype_manual(values = rep(line_styles, length.out = length(unique(confidence_predictions_df$Site)))) +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 12)),
        axis.title.y = element_text(size = 16, margin = margin(r = 12)))

## ----- Part 3: Seasonal Differences- ##
data_season <- read.csv("path/to/your/seasonal_data_file.csv", stringsAsFactors = FALSE)

season.conf.model <- glm(MANUAL.ID ~ Confidence + SEASON, data = data_season, family = "binomial")
summary(season.conf.model)

# Seasonal predictions using Confidence Score
confidence_model_list_season <- lapply(unique(data_season$SEASON), function(s) {
  subset_data <- subset(data_season, SEASON == s)
  model <- glm(MANUAL.ID ~ Confidence, data = subset_data, family = binomial)
  list(season = s, model = model)
})
confidence_predictions_list_season <- lapply(confidence_model_list_season, function(entry) {
  s <- entry$season
  model <- entry$model
  conf_vals <- seq(min(data_season$Confidence), max(data_season$Confidence), length.out = 100)
  data.frame(Confidence = conf_vals,
             Probability = predict(model, newdata = data.frame(Confidence = conf_vals), type = "response"),
             Season = s)
})
confidence_predictions_df_season <- do.call(rbind, confidence_predictions_list_season)
overall_confidence_mean_season <- aggregate(Probability ~ Confidence, confidence_predictions_df_season, mean)

custom_palette_season <- c("#0072B2", "#D55E00")
confidence_season_plot <- ggplot(confidence_predictions_df_season, aes(x = Confidence, y = Probability, color = SEASON)) +
  geom_line(size = 0.65) +
  labs(x = "Confidence Score", y = "Accuracy") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = custom_palette_season) +
  theme_minimal(base_size = base_size)

## ----- Part 4: Pairwise Comparisons for Interaction Effects ----- ##
# Fit a model with interaction between SEASON and SITE with Confidence Score
season.conf.model.interaction <- glm(MANUAL.ID ~ Confidence + SEASON * SITE, data = data_season, family = "binomial")
summary(season.conf.model.interaction)

# Extract coefficients for the interaction terms (if present)
interaction_effects <- coef(season.conf.model.interaction)[grep("SEASONWet Season:SITE", rownames(coef(season.conf.model.interaction))), ]
print(interaction_effects)

library(emmeans)
emm <- emmeans(season.conf.model.interaction, ~ SEASON | SITE)
print(summary(emm))
pairwise_comparisons <- contrast(emm, method = "pairwise")
print(summary(pairwise_comparisons))

## ----- Save Plots ----- ##
library(patchwork)

# Save the spatial plot
ggsave("Accuracy_By_Sites.jpg", confidence_sites_plot, dpi = 800)
# Save the seasonal plot
ggsave("Accuracy_By_Seasons.jpg", confidence_season_plot, dpi = 800)