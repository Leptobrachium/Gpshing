###########################################################################
# Cane Toad Detection Threshold Optimization and Performance Evaluation
#
# This script implements a detection optimization approach adapted from
# Singer et al. (2024; https://doi.org/10.1002/rse2.385).
# To enhance accuracy, time-series features are aggregated over multiple time
# intervals and modeled via conditional inference trees (CIT) to identify optimal
# threshold rules. CIT models are ranked based on a weighted sum of precision and
# recall. Additionally, a set-cover optimization is applied and bootstrapping is
# used to rank the predictors. Finally, optimized thresholds are compared against
# universal thresholds (UNI10, UNI50, UNI90).
#
# Note: This script has been modified for a single species (Rhinella marina) and
# file paths and file names that should be updated accordingly.
#
# References:
# - Singer, D., Hagge, J., Kamp, J., Hondong, H. & Schuldt, A. (2024) Aggregated
# time-series features boost species-specific differentiation of true and false
# positives in passive acoustic monitoring of bird assemblages, Remote Sensing
# in Ecology and Conservation, https://doi.org/10.1002/rse2.385
#
# Dependencies: data.table, dplyr, stringr, lubridate, tidyr, foreach, doParallel,
#               partykit, lpSolve, ggplot2, ggsignif, ggpubr, patchwork, emmeans
###########################################################################

### ----------------------- Section 1: Calculate Aggregated Time‐Series Features (ATF) ----------------------- ###
library(data.table)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)

# function.r can be obtained from https://github.com/d-singer/BirdNET-thresholds
source("path/to/your/functions.R")

# Import human validated detections (update file path)
validated_detections <- fread("path/to/your/validated_detections.csv", encoding = "UTF-8") %>%
  mutate(timestamp = str_replace(timestamp, "T", " ") %>% str_remove("Z"))

# Import all BirdNET detections (update file path; e.g., an RDS extract)
birdnet_all <- read.csv("path/to/your/birdnet_all_detections.csv") %>%
  mutate(timestamp = str_replace(timestamp, "T", " ") %>% str_remove("Z")) %>%
  mutate(timestamp = ymd_hms(timestamp))
validated_detections <- validated_detections %>%
  mutate(timestamp = str_replace(timestamp, "T", " ") %>% str_remove("Z")) %>%
  mutate(timestamp = ymd_hms(timestamp))

# Define species name (single species)
species_name <- "Rhinella marina"

# Calculate ATF using a custom function (ensure calculate_atf() is defined in functions.R)
atf <- calculate_atf(birdnet = birdnet_all, species = species_name, validets = validated_detections)

# Save ATF data as CSV and RDS
write.csv(atf, "path/to/your/output/combined_canetoad_aggregated_timeseries_features.csv", row.names = FALSE)
saveRDS(atf, "path/to/your/output/combined_canetoad_aggregated_timeseries_features.rds")

### ----------------------- Section 2: Threshold Modeling via Conditional Inference Trees (CIT) ----------------------- ###
library(foreach)
library(doParallel)
library(partykit)

setwd("path/to/your/analysis/directory")

# Read ATF features (update file path)
atf_data <- read.csv("path/to/your/output/combined_canetoad_aggregated_timeseries_features.csv", encoding = "UTF-8")

# Use single species
species <- species_name

# Extract predictor variable names (assumes predictors start at column 6)
preds <- c("conf", colnames(atf_data)[6:(ncol(atf_data)-1)])

# --- Type 1 Models (tree depth = 1) ---
ctrl1 <- ctree_control(teststat = "max", testtype = "Bonferroni", maxdepth = 1, alpha = 0.05)
num_cores <- parallel::detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

ctrees_type1 <- foreach(m = preds, .combine = rbind, .verbose = TRUE) %dopar% {
  mydat <- atf_data[atf_data$german == species, ]
  mytree <- ctree(data = mydat, formula = as.formula(paste0("validation ~ ", m)), control = ctrl1)
  
  # Get rules
  rules <- as.data.frame(.list.rules.party(mytree))
  rules <- cbind(node = as.numeric(rownames(rules)), rules)
  colnames(rules) <- c("node", "rules")
  rules$german <- species
  rules$model <- m
  
  n <- c(); n_tp <- c(); prec <- c(); node <- c()
  for (x in 1:length(mytree[[1]])) {
    y <- mytree[[1]][[x]]$data
    prec <- c(prec, sum(y$validation == 1) / nrow(y))
    n <- c(n, nrow(y))
    n_tp <- c(n_tp, sum(y$validation == 1))
    node <- c(node, x)
  }
  ctrees <- data.frame(german = species, prec = prec, n = n, n_tp = n_tp, node = node, model = m, stringsAsFactors = FALSE)
  ctrees <- left_join(ctrees, rules, by = c("node", "german", "model"))
  return(ctrees)
}

# --- Type 2 Models (tree depth = 2) ---
ctrl2 <- ctree_control(teststat = "max", testtype = "Bonferroni", maxdepth = 2, alpha = 0.05)
models_type2 <- c()
for (i in 1:2) {
  allcomb <- combn(preds, i)
  for (j in 1:ncol(allcomb)) {
    models_type2 <- c(models_type2, paste(allcomb[, j], collapse = "+"))
  }
}
models_type2 <- sort(unique(models_type2))

ctrees_type2 <- foreach(m = models_type2, .combine = rbind, .verbose = TRUE) %dopar% {
  mydat <- atf_data[atf_data$german == species, ]
  mytree <- ctree(data = mydat, formula = as.formula(paste0("validation ~ ", m)), control = ctrl2)
  
  rules <- as.data.frame(.list.rules.party(mytree))
  rules <- cbind(node = as.numeric(rownames(rules)), rules)
  colnames(rules) <- c("node", "rules")
  rules$german <- species
  rules$model <- m
  
  n <- c(); n_tp <- c(); prec <- c(); node <- c()
  for (x in 1:length(mytree[[1]])) {
    y <- mytree[[1]][[x]]$data
    prec <- c(prec, sum(y$validation == 1) / nrow(y))
    n <- c(n, nrow(y))
    n_tp <- c(n_tp, sum(y$validation == 1))
    node <- c(node, x)
  }
  ctrees <- data.frame(german = species, prec = prec, n = n, n_tp = n_tp, node = node, model = m, stringsAsFactors = FALSE)
  ctrees <- left_join(ctrees, rules, by = c("node", "german", "model"))
  return(ctrees)
}
stopCluster(cl)

# Save CIT models (both RDS and CSV)
saveRDS(ctrees_type1, "path/to/your/output/canetoad_ctrees_type1.rds")
write.csv(ctrees_type1, "path/to/your/output/canetoad_ctrees_type1.csv", row.names = FALSE)
saveRDS(ctrees_type2, "path/to/your/output/canetoad_ctrees_type2.rds")
write.csv(ctrees_type2, "path/to/your/output/canetoad_ctrees_type2.csv", row.names = FALSE)

### ----------------------- Section 3: Model Selection ----------------------- ###
library(data.table)
library(tidyr)
library(ggplot2)

setwd("path/to/your/analysis/directory")

# Load manually verified sample data (update file path)
detsample <- fread("path/to/your/output/combined_canetoad_aggregated_timeseries_features.csv", encoding = "UTF-8") %>%
  select(german, conf, validation) %>%
  group_by(german) %>%
  mutate(tp_total = sum(validation)) %>%
  filter(german == species)

# Load CIT models and combine
models1 <- readRDS("path/to/your/output/canetoad_ctrees_type1.rds")
models1$ctree_ctrl <- "ctrl1"
models2 <- readRDS("path/to/your/output/canetoad_ctrees_type2.rds")
models2$ctree_ctrl <- "ctrl2"
models_all <- rbind(models1, models2)

# Filter for terminal nodes with maximum precision per model
models_all <- models_all %>% group_by(model, german, ctree_ctrl) %>%
  mutate(n_nodes = max(node),
         prec_max = max(prec),
         is_terminal = ifelse(n_nodes == 1 & is.na(rules), 1,
                              ifelse(n_nodes > 1 & is.na(rules), 0, 1)),
         is_max_prec = ifelse(prec == prec_max, 1, 0)) %>%
  filter(is_terminal == 1 & is_max_prec == 1)

# Calculate recall for each model and overall model performance (w = 0.75)
models_all <- models_all %>% left_join(unique(detsample %>% select(german, tp_total)), by = "german") %>%
  mutate(recall = n_tp / tp_total)
w <- 0.75
models_all <- models_all %>% group_by(german) %>%
  mutate(model_performance = prec * w + recall * (1 - w),
         rank_performance = dense_rank(-model_performance))

# Save candidate and basic models (RDS and CSV)
cand_mods <- models_all %>% filter(rank_performance == 1)
saveRDS(cand_mods, "path/to/your/output/canetoad_candidate_models.rds")
write.csv(cand_mods, "path/to/your/output/canetoad_candidate_models.csv", row.names = FALSE)
base_mods <- models_all %>% filter(model == "conf" & ctree_ctrl == "ctrl1")
saveRDS(base_mods, "path/to/your/output/canetoad_basic_models.rds")
write.csv(base_mods, "path/to/your/output/canetoad_basic_models.csv", row.names = FALSE)

### ----------------------- Section 4: Set Cover Optimization ----------------------- ###
library(lpSolve)

# Load candidate models (update file path)
cand_mods <- fread("path/to/your/output/canetoad_candidate_models.csv", encoding = "UTF-8")
cand_mods <- cand_mods %>% mutate(model_id = paste0(ifelse(ctree_ctrl == "ctrl1", "type1_", "type2_"), model))
species_vec <- unique(cand_mods$german)
matrix_lp <- cand_mods %>% select(german, model_id, rank_performance) %>%
  pivot_wider(names_from = model_id, values_from = rank_performance, values_fill = list(rank_performance = 0)) %>%
  select(-german) %>% as.matrix()
model_names <- colnames(matrix_lp)
lp_model <- lp(direction = "min",
               objective.in = rep(1, ncol(matrix_lp)),
               const.mat = matrix_lp,
               const.dir = rep(">=", length(species_vec)),
               const.rhs = rep(1, length(species_vec)))
opti_mods <- sort(model_names[which(lp_model$solution >= 1)])
opti_mods <- gsub(x = opti_mods, "type2_", "")
# Save optimized model set (RDS and CSV)
saveRDS(opti_mods, "path/to/your/output/canetoad_opti_mods_setcover.rds")
write.csv(as.data.frame(opti_mods), "path/to/your/output/canetoad_opti_mods_setcover.csv", row.names = FALSE)

### ----------------------- Section 5: Bootstrapping Order of ATF Variables ----------------------- ###
set.seed(3)
models_all <- readRDS("path/to/your/output/canetoad_all_models.rds")  # Combined models from selection
atf_vars <- c("ndets10", "ndets20", "ndets30", "ndets40", "ndets50",
              "ndets60", "ndets70", "ndets80", "ndets90", "ndets99",
              "avgconf", "medconf", "maxconf", "minconf")
opti_spec <- data.frame(prec = numeric(), recall = numeric(), model_performance = numeric(),
                        remove = character(), run = integer(), rank = integer(), stringsAsFactors = FALSE)
for (i in 1:999) {
  myatf <- c("none", sample(atf_vars))
  mymodels <- models_all
  for (v in myatf) {
    mymodels$remove <- grepl(pattern = v, mymodels$model)
    mymodels <- mymodels[mymodels$remove == FALSE, ]
    opti_mods <- mymodels %>% mutate(perf_rank = dense_rank(-model_performance)) %>%
      filter(perf_rank == 1) %>%
      summarise(prec = unique(prec),
                recall = unique(recall),
                model_performance = unique(model_performance),
                remove = v) %>%
      mutate(run = i, rank = which(myatf == v))
    opti_spec <- rbind(opti_spec, opti_mods)
  }
  if (i %% 100 == 0) print(i)
}
# Save bootstrapping results
saveRDS(opti_spec, "path/to/your/output/canetoad_bootstrapping_atf.rds")
write.csv(opti_spec, "path/to/your/output/canetoad_bootstrapping_atf.csv", row.names = FALSE)

### ----------------------- Section 6: Reduction of Predictors ----------------------- ###
library(ggplot2)
models_all <- readRDS("path/to/your/output/canetoad_all_models.rds")
timevars <- c("none", "int12", "int11", "int10", "int09", "int08", "int07",
              "int06", "int05", "int04", "int03", "int02", "int01")
opti_spec_redux <- data.frame()
for (v in timevars) {
  models_all$remove <- grepl(pattern = v, models_all$model)
  reduced_models <- models_all[models_all$remove == FALSE, ]
  opti_mods <- reduced_models %>% group_by(german) %>%
    mutate(perf_rank = dense_rank(-model_performance)) %>%
    filter(perf_rank == 1) %>%
    summarise(prec = unique(prec),
              recall = unique(recall),
              model_performance = unique(model_performance),
              remove = v)
  opti_spec_redux <- rbind(opti_spec_redux, opti_mods)
}
# Save reduction results (RDS and CSV)
saveRDS(opti_spec_redux, "path/to/your/output/canetoad_reduction_of_predictors.rds")
write.csv(opti_spec_redux, "path/to/your/output/canetoad_reduction_of_predictors.csv", row.names = FALSE)

# Summarize performance metrics across time window removals
perf_all <- opti_spec_redux %>% group_by(remove) %>%
  summarize(mean_perf = mean(model_performance),
            mean_prec = mean(prec),
            mean_recall = mean(recall)) %>%
  mutate(spec = species)
perf <- perf_all %>% pivot_longer(cols = c("mean_perf", "mean_prec", "mean_recall"))
# Create descriptive labels for time windows
perf$label <- perf$remove
perf$label <- gsub(perf$label, pattern="int01", replacement="conf. score")
perf$label <- gsub(perf$label, pattern="int02", replacement="±3 s")
perf$label <- gsub(perf$label, pattern="int03", replacement="±6 s")
perf$label <- gsub(perf$label, pattern="int04", replacement="±9 s")
perf$label <- gsub(perf$label, pattern="int05", replacement="±12 s")
perf$label <- gsub(perf$label, pattern="int06", replacement="±10 min")
perf$label <- gsub(perf$label, pattern="int07", replacement="±20 min")
perf$label <- gsub(perf$label, pattern="int08", replacement="±30 min")
perf$label <- gsub(perf$label, pattern="int09", replacement="±40 min")
perf$label <- gsub(perf$label, pattern="int10", replacement="±12 h")
perf$label <- gsub(perf$label, pattern="int11", replacement="±24 h")
perf$label <- gsub(perf$label, pattern="int12", replacement="±36 h")
perf$label <- gsub(perf$label, pattern="none", replacement="±48 h")
perf$variable <- gsub(perf$name, pattern="mean_perf", replacement="Model Performance")
perf$variable <- gsub(perf$variable, pattern="mean_prec", replacement="Precision")
perf$variable <- gsub(perf$variable, pattern="mean_recall", replacement="Recall")
perf$variable <- factor(perf$variable, levels = c("Model Performance", "Precision", "Recall"))
perf$label <- factor(perf$label, levels = rev(c("conf. score", "±3 s", "±6 s", "±9 s", "±12 s",
                                                "±10 min", "±20 min", "±30 min", "±40 min",
                                                "±12 h", "±24 h", "±36 h", "±48 h")))
perf <- perf %>% arrange(label)

# Create separate plots for each performance metric
plot_model_perf <- ggplot(perf %>% filter(variable == "Model Performance"), aes(x = label, y = value)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_line(group = 1, color = "darkgreen") + 
  labs(x = "Time Window Removal", y = "Mean Model Performance", title = "Reduction of Predictors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_precision <- ggplot(perf %>% filter(variable == "Precision"), aes(x = label, y = value)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_line(group = 1, color = "darkgreen") +
  labs(x = "Time Window Removal", y = "Mean Precision", title = "Reduction of Predictors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_recall <- ggplot(perf %>% filter(variable == "Recall"), aes(x = label, y = value)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_line(group = 1, color = "darkgreen") +
  labs(x = "Time Window Removal", y = "Mean Recall", title = "Reduction of Predictors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine the three plots vertically
library(patchwork)
reduction_plots <- plot_model_perf / plot_precision / plot_recall
ggsave("path/to/your/output/Reduction_of_Predictors.jpg", reduction_plots, dpi = 1000)
# Also save as RDS and CSV if needed (here we save the data summary)
saveRDS(perf, "path/to/your/output/reduction_of_predictors_summary.rds")
write.csv(perf, "path/to/your/output/reduction_of_predictors_summary.csv", row.names = FALSE)

### ----------------------- Section 7: Performance Evaluation ----------------------- ###
library(ggpubr)
validated_detections <- fread("path/to/your/combined_canetoad_aggregated_timeseries_features.csv", encoding = "UTF-8")
validated_detections <- validated_detections %>%
  filter(german == species) %>%
  select(german, conf, validation) %>%
  mutate(tp_total = sum(validation))
cand_mods <- readRDS("path/to/your/output/canetoad_candidate_models.rds")

# Evaluate universal thresholds UNI10, UNI50, UNI90 using raw Confidence
uni10 <- validated_detections %>%
  filter(conf >= 0.1) %>%
  summarise(german = species,
            prec = mean(validation),
            threshold = "uni10",
            recall = sum(validation) / tp_total) %>%
  mutate(model_performance = prec * 0.75 + recall * 0.25)
uni50 <- validated_detections %>%
  filter(conf >= 0.5) %>%
  summarise(german = species,
            prec = mean(validation),
            recall = sum(validation) / tp_total,
            threshold = "uni50") %>%
  mutate(model_performance = prec * 0.75 + recall * 0.25)
uni90 <- validated_detections %>%
  filter(conf >= 0.9) %>%
  summarise(german = species,
            prec = mean(validation),
            recall = sum(validation) / tp_total,
            threshold = "uni90") %>%
  mutate(model_performance = prec * 0.75 + recall * 0.25)

# Use candidate models as the optimized threshold ("opti")
opti_para <- cand_mods %>% summarize(german = species,
                                     prec = unique(prec),
                                     recall = unique(recall),
                                     threshold = "opti",
                                     model_performance = unique(prec) * 0.75 + unique(recall) * 0.25)

# Combine candidate and universal thresholds
mod_comp <- bind_rows(opti_para, uni10, uni50, uni90)
mod_comp <- mod_comp %>% mutate(threshold = factor(threshold, levels = c("opti", "uni10", "uni50", "uni90")))
mod_comp$german <- factor(mod_comp$german, levels = unique(mod_comp$german))
mod_comp$threshold <- factor(mod_comp$threshold, levels = c("opti", "uni10", "uni50", "uni90"))

# Plot example: Model Performance vs. Threshold
perf_plot <- ggplot(mod_comp, aes(x = model_performance, y = german, fill = threshold)) + 
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = 0.9, lty = 2) + 
  xlab("Model Performance") +
  theme_bw() +
  theme(legend.position = 'bottom')
print(perf_plot)
ggsave("path/to/your/output/Model_Comparison.jpg", perf_plot, dpi = 1000)
fwrite(mod_comp, "path/to/your/output/canetoad_model_comparison.csv", row.names = FALSE)

###########################################################################
# End of Script
###########################################################################