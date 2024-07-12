library(data.table)
library(tidyverse)
library(lubridate)
library(runner)
library(caret)

setwd("/home/weiwen/code")

# Function to read and clean results
read_and_clean_results <- function(filepath) {
  results <- read.table(filepath, fill = TRUE, sep = "\t", header = TRUE)
  results <- as.data.frame(sapply(results, function(x) str_replace_all(x, "[',()\\]\\[]", ""))) # Remove unwanted characters
  results <- as.data.frame(sapply(results, function(x) trimws(gsub("\\s+", " ", x)))) # Remove double spaces
  rownames(results) <- results$X
  return(results)
}

# Function to parse sublineages data
parse_sublineages <- function(results, colname) {
  lineages.temp <- as.data.frame(t(setDT(tstrsplit(as.character(results["lineages", colname]), " ", fixed = TRUE))[]))
  abundances.temp <- as.data.frame(t(setDT(tstrsplit(as.character(results["abundances", colname]), " ", fixed = TRUE))[]))
  sample.temp <- rep(sub(".variants.tsv$", "", colname), nrow(lineages.temp))
  sublineages.final <- cbind(sample.temp, lineages.temp, abundances.temp)
  names(sublineages.final) <- c("Sample", "Lineage", "res_abundance")
  sublineages.final <- arrange(sublineages.final, Lineage)
  return(sublineages.final)
}

# Function to read reference data
read_reference_data <- function(sample) {
  projectdir <- getwd()
  ref_filepath <- paste0(projectdir, "/results/sample_sets/", sample, ".known_lineages.tsv")
  reference_df <- read_tsv(ref_filepath)
  reference_df <- arrange(reference_df, Lineage)
  return(reference_df)
}

# Function to calculate MAE
calculate_mae <- function(true, predicted) {
  mean(abs(true - predicted))
}

# Function to calculate F1 score and MAE
calculate_metrics <- function(reference_df, sublineages.final, mae_threshold = 0.01) {
  reference_df$Lineage <- as.factor(reference_df$Lineage)
  sublineages.final$Lineage <- as.factor(sublineages.final$Lineage)
  
  if (setequal(levels(reference_df$Lineage), levels(sublineages.final$Lineage))) {
    comparison_df <- merge(reference_df, sublineages.final, by = "Lineage", suffixes = c("_true", "_pred"))
    comparison_df$Frequency <- as.numeric(comparison_df$Frequency)
    comparison_df$res_abundance <- as.numeric(comparison_df$res_abundance)
    
    comparison_df <- comparison_df %>%
      mutate(MAE = calculate_mae(Frequency, res_abundance),
             Prediction = ifelse(MAE < mae_threshold, 1, 0),
             Truth = 1)
    overall_mae <- round(mean(comparison_df$MAE))
    
    if (length(unique(comparison_df$Prediction)) == 1) {
      f1_score <- 1
    } else {
      confusion_mat <- confusionMatrix(as.factor(comparison_df$Prediction), as.factor(comparison_df$Truth))
      f1_score <- confusion_mat$byClass["F1"]
    }
    
  } else {
    all_lineages <- union(levels(reference_df$Lineage), levels(sublineages.final$Lineage))
    extended_reference_df <- merge(data.frame(Lineage = all_lineages), reference_df, by = "Lineage", all.x = TRUE)
    extended_sublineages_final <- merge(data.frame(Lineage = all_lineages), sublineages.final, by = "Lineage", all.x = TRUE)
    
    extended_reference_df$Frequency[is.na(extended_reference_df$Frequency)] <- 0
    extended_sublineages_final$res_abundance[is.na(extended_sublineages_final$res_abundance)] <- 0
    
    comparison_df <- merge(extended_reference_df, extended_sublineages_final, by = "Lineage", suffixes = c("_true", "_pred"))
    comparison_df$Frequency <- as.numeric(comparison_df$Frequency)
    comparison_df$res_abundance <- as.numeric(comparison_df$res_abundance)
    
    comparison_df <- comparison_df %>%
      mutate(MAE = calculate_mae(Frequency, res_abundance),
             Prediction = ifelse(MAE < mae_threshold, 1, 0),
             Truth = 1)
    overall_mae <- round(mean(comparison_df$MAE))
    
    if (length(unique(comparison_df$Prediction)) == 1) {
      f1_score <- 1
    } else {
      confusion_mat <- confusionMatrix(as.factor(comparison_df$Prediction), as.factor(comparison_df$Truth))
      f1_score <- confusion_mat$byClass["F1"]
    }
  }
  
  return(c(comparison_df$Sample[1], overall_mae, f1_score))
}

# Main function to run the analysis
run_analysis <- function(results_filepath, allstat_df) {
  results <- read_and_clean_results(results_filepath)
  sublineages.final <- parse_sublineages(results, colnames(results)[2])
  reference_df <- read_reference_data(sublineages.final$Sample[1])
  
  stat_results <- calculate_metrics(reference_df, sublineages.final)
  
  # Append results to the summary dataframe
  allstat_df <- rbind(allstat_df, data.frame(sample = stat_results[1],
                                             mae_mean = as.numeric(stat_results[2]),
                                             F1_score = as.numeric(stat_results[3]),
                                             stringsAsFactors = FALSE))
  rm(list = setdiff(ls(), "allstat_df"))
  return(allstat_df)
}

# Initialize an empty dataframe for storing results
allstat_df <- data.frame(sample = character(),
                         mae_mean = numeric(),
                         F1_score = numeric(),
                         stringsAsFactors = FALSE)

# Run the analysis
file_list <- list.files("results/demix_results/", pattern = "\\.tsv$")
for (file in file_list) {
  # Check if the file name contains "empty"
  if (grepl("empty", file)) {
    next  # Skip to the next iteration of the loop
  }
  
  # Run the analysis and update the allstat_df
  allstat_df <- run_analysis(paste0("results/demix_results/",file), allstat_df)
}
write.csv(allstat_df, file = "results/allstat_df.csv", row.names = FALSE)
