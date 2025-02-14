###########################################################################
# Classifier Preliminary Performance Assessment
#
# This script demonstrates a preliminary performance assessment of the 
# custom-trained cane toad classifier. For each paired gold standard and 
# BirdNET recogniser results file, the script:
# 1. Loads manually labeled (gold standard) and recogniser data.
# 2. Aggregates data into 3-second groupings over a range of temporal resolutions.
# 3. Evaluates classifier performance (using eventEval from monitoR) at different 
#    confidence levels.
# 4. Computes precision and recall metrics.
# 5. Visualizes performance using Recall and Precision plots across temporal resolutions.
#
# References:
# - Hafner and Katz (2018) for the eventEval function in monitoR.
#
# Dependencies: dplyr, ggplot2, gridExtra, monitoR, cowplot, ggpubr
#
# Note: Customize the file paths and names below with your actual data.
###########################################################################

library(dplyr)
library(ggplot2)
library(gridExtra)
library(monitoR)
library(cowplot)
library(ggpubr)

#############################
# 1. Set Working Directory & Define File Paths
#############################

# Set your working directory (customize this path as needed)
setwd("your/working/directory")

# Define file paths for gold standard (manually labeled) data.
# Update these file paths with your actual gold standard CSV files.
labelledDataFiles <- c(
  "your/working/directory/gold_standard_file1.csv",
  "your/working/directory/gold_standard_file2.csv",
  "your/working/directory/gold_standard_file3.csv",
  "your/working/directory/gold_standard_file4.csv",
  "your/working/directory/gold_standard_file5.csv",
  "your/working/directory/gold_standard_file6.csv",
  "your/working/directory/gold_standard_file7.csv",
  "your/working/directory/gold_standard_file8.csv",
  "your/working/directory/gold_standard_file9.csv"
)

# Define file paths for BirdNET recogniser results.
# Update these file paths with your actual BirdNET results CSV files.
recogniserDataFiles <- c(
  "your/working/directory/birdnet_results_file1.csv",
  "your/working/directory/birdnet_results_file2.csv",
  "your/working/directory/birdnet_results_file3.csv",
  "your/working/directory/birdnet_results_file4.csv",
  "your/working/directory/birdnet_results_file5.csv",
  "your/working/directory/birdnet_results_file6.csv",
  "your/working/directory/birdnet_results_file7.csv",
  "your/working/directory/birdnet_results_file8.csv",
  "your/working/directory/birdnet_results_file9.csv"
)

#############################
# 2. Define Utility Functions
#############################

# Function to determine all factors of an integer (used to define temporal groupings)
factors <- function(x) {
  x <- as.integer(x)
  div <- seq_len(abs(x))
  factors <- div[x %% div == 0L]
  return(factors)
}

#############################
# 3. Evaluate Classifier Performance Across Temporal Resolutions
#############################

# Create an empty data frame to store performance results
results <- data.frame()

# Loop through each paired gold standard and recogniser file
for (i in seq_along(labelledDataFiles)) {
  
  # Load gold standard data and assign a binary label (1 = Cane_Toad, 0 = Not Cane Toad)
  labelledData <- read.csv(labelledDataFiles[i]) %>% 
    mutate(Label = ifelse(name == 'Cane_Toad', 1, 0))
  
  # Load recogniser (BirdNET) results
  recogniserData <- read.csv(recogniserDataFiles[i])
  
  # Iterate over possible temporal grouping factors (derived from the number of rows)
  for (grouping in factors(nrow(labelledData))) {
    
    # Aggregate gold standard data into groups based on 3-second intervals
    labelledData_aggregated <- labelledData %>% 
      mutate(Group = ceiling(end.time / (grouping * 3))) %>% 
      group_by(Group) %>% 
      summarise(Label = max(Label)) %>% 
      ungroup() %>% 
      mutate(name = ifelse(Label == 1, 'Cane_toad', 'N'),
             start.time = (Group - 1) * grouping * 3,
             end.time = Group * grouping * 3,
             min.frq = 0,
             max.frq = 11025)
    
    # Aggregate recogniser data into corresponding groups
    recogniserData_aggregated <- recogniserData %>% 
      mutate(Group = ceiling(time / (grouping * 3))) %>% 
      group_by(Group) %>% 
      summarise(score = max(score, na.rm = TRUE)) %>% 
      ungroup() %>% 
      mutate(time = (Group * grouping * 3) - grouping * 3 / 2,
             template = "Cane_toad")
    
    # Remove groups with invalid scores (-Inf)
    recogniserData_aggregated <- recogniserData_aggregated %>% 
      filter(score != -Inf)
    
    # Evaluate classifier performance over a range of confidence thresholds
    for (confidence in seq(0.1, 0.9, 0.1)) {
      # Use eventEval from monitoR to compare detections with gold standard labels
      eval_result <- eventEval(
        detections = as.data.frame(recogniserData_aggregated),
        standard = as.data.frame(labelledData_aggregated %>% filter(name == 'Cane_toad')),
        score.cutoff = confidence,
        tol = grouping * 3 / 2
      )
      
      # If evaluation results exist, store them along with grouping and confidence info
      if (nrow(eval_result) > 0) {
        results <- bind_rows(
          results,
          data.frame(eval_result,
                     TemporalGrouping = grouping,
                     ModelConfidence = confidence)
        )
      }
    }
  }
}

#############################
# 4. Calculate Performance Metrics (Recall & Precision)
#############################

PerformanceResults <- results %>% 
  group_by(TemporalGrouping, ModelConfidence) %>% 
  summarise(Recall = sum(outcome == "TRUE +") / (sum(outcome == "TRUE +") + sum(outcome == "FALSE -")),
            Precision = sum(outcome == "TRUE +") / (sum(outcome == "TRUE +") + sum(outcome == "FALSE +")))

#############################
# 5. Visualize Performance Results
#############################

# Define a colorblind-friendly palette (Okabe-Ito palette)
cb_palette <- c("0.1" = "#CC79A7",  # Purple
                "0.5" = "#0072B2",  # Blue
                "0.9" = "#009E73")  # Green

# Plot Recall vs. Temporal Resolution
Plot_Recall <- PerformanceResults %>%
  filter(ModelConfidence %in% c(0.1, 0.5, 0.9)) %>%
  mutate(TemporalResolution = TemporalGrouping * 3) %>%
  ggplot(aes(x = TemporalResolution, y = Recall, 
             colour = as.factor(ModelConfidence), 
             group = as.factor(ModelConfidence))) +
  geom_path(size = 0.8) +
  scale_color_manual(values = cb_palette) +
  scale_x_log10(breaks = c(3, 10, 30, 100, 300, 600, 900, 1800, 3600)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(colour = "Model Confidence", x = "Temporal Resolution (s)", y = "Recall") +
  theme_bw() +
  theme(legend.position = c(0.835, 0.2),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        axis.title.y = element_text(size = 14, margin = margin(r = 8)),
        axis.title.x = element_text(size = 14, margin = margin(t = 8)),
        axis.text = element_text(size = 11))

# Plot Precision vs. Temporal Resolution
Plot_Precision <- PerformanceResults %>%
  filter(ModelConfidence %in% c(0.1, 0.5, 0.9)) %>%
  mutate(TemporalResolution = TemporalGrouping * 3) %>%
  ggplot(aes(x = TemporalResolution, y = Precision, 
             colour = as.factor(ModelConfidence), 
             group = as.factor(ModelConfidence))) +
  geom_path(size = 0.8) +
  scale_color_manual(values = cb_palette) +
  scale_x_log10(breaks = c(3, 10, 30, 100, 300, 600, 900, 1800, 3600)) +
  scale_y_continuous(limits = c(0.90, 1), breaks = seq(0.90, 1, 0.05)) +
  labs(colour = "Model Confidence", x = "Temporal Resolution (s)", y = "Precision") +
  theme_bw() +
  theme(legend.position = c(0.835, 0.2),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        axis.title.y = element_text(size = 14, margin = margin(r = 8)),
        axis.title.x = element_text(size = 14, margin = margin(t = 8)),
        axis.text = element_text(size = 11))

# Arrange the two plots side by side with labels
Plot_Precision_Recall <- ggarrange(Plot_Precision, Plot_Recall,
                                   labels = c("A", "B"),
                                   ncol = 2)

# Display the plots
print(Plot_Precision)
print(Plot_Recall)
print(Plot_Precision_Recall)

# Save the combined plot as a high-resolution JPEG file
ggsave("Recall_and_Precision.jpg", Plot_Precision_Recall, dpi = 1000)
