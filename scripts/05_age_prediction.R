# Project: DNA degradation study
# Script: 05_age_prediction.R
# Description: This script uses the preprocessed methylation data to predict
# epigenetic age for each sample. It then calculates the Mean Absolute Error (MAE)
# against the donor's known chronological age to assess prediction accuracy
# under different DNA degradation conditions.

#### 1. Configuration & Setup ------------------------------------------------

# Use the 'here' package for portable file paths
library(here)

# --- Paths ---
processed_data_dir <- here("data", "processed")
results_age_dir <- here("results", "5_Age_prediction")
dir.create(results_age_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Libraries ---
library(methylclock)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(writexl)
library(tibble)


#### 2. Load and Prepare Data ------------------------------------------------

# Load preprocessed beta values and metadata
betas_full <- readRDS(file.path(processed_data_dir, "betas_preprocessed.rds"))
metadata <- readRDS(file.path(processed_data_dir, "metadata_preprocessed.rds"))

# --- Filter Samples that Passed QC ---
samples_to_keep_mask <- !grepl("95bp", colnames(betas_full)) &
  !colnames(betas_full) %in% c("165bp_10ng_A", "165bp_10ng_B")

betas_subset <- betas_full[, samples_to_keep_mask]
metadata_subset <- metadata[samples_to_keep_mask, ]

# --- Prepare Beta Matrix for Age Prediction ---
# Collapse CpG probes. The DNAmAge function handles NAs.
betas_cg_coll <- betasCollapseToPfx(betas_subset[grepl("^cg", rownames(betas_subset)), ])

# Format data for methylclock: needs a data.frame with probe IDs in the first column
dat_for_age <- as.data.frame(betas_cg_coll) %>%
  rownames_to_column("probeID")


#### 3. Perform Epigenetic Age Prediction ------------------------------------

message("Performing epigenetic age prediction for all samples...")

# Define the known chronological age of the sample donor
# The original script used 40.30795306
chronological_age <- 40.3 

# Predict DNAm age using various clocks
predicted_ages <- DNAmAge(dat_for_age)

# Calculate Mean Absolute Error (MAE) for each clock and sample
mae_df <- as.data.frame(predicted_ages) %>%
  # Select clocks of interest
  select(Horvath, skinHorvath, BLUP, EN) %>%
  # Calculate absolute difference from chronological age
  mutate(across(everything(), ~ abs(. - chronological_age))) %>%
  # Add back sample identifiers
  tibble::add_column(Sample_code = rownames(predicted_ages), .before = 1) %>%
  # Join with metadata to get sample conditions
  left_join(metadata_subset, by = "Sample_code")

#### 4. Summarize and Visualize MAE Results ----------------------------------

# Average MAE for technical replicates
mae_avg_df <- mae_df %>%
  mutate(Sample_combination = gsub("_A|_B", "", Sample_code)) %>%
  group_by(Sample_combination, DNA_size, DNA_input) %>%
  summarise(across(c(Horvath, skinHorvath, BLUP, EN), mean, na.rm = TRUE), .groups = "drop") %>%
  # Exclude the control sample for degradation-specific plots
  filter(!grepl("no_degraded", Sample_combination)) %>%
  # Factor variables for correct plot ordering
  mutate(
    DNA_size = factor(DNA_size, levels = c("350", "230", "165")),
    DNA_input = factor(DNA_input, levels = c("100", "50", "20", "10"))
  )

# --- Save Summary Table ---
write_xlsx(mae_avg_df, file.path(results_age_dir, "mae_summary_by_condition.xlsx"))

# --- Visualize MAE as Heatmap Tables ---
# Note: This requires a reusable plotting function, e.g., from a 'functions.R' file
# For simplicity here, we create one plot. The pattern can be repeated for other clocks.

plot_mae_table <- function(df, clock_col, title) {
  ggplot(df, aes(x = DNA_size, y = DNA_input, fill = {{clock_col}})) +
    geom_tile(color = "black") +
    geom_text(aes(label = round({{clock_col}}, 1)), color = "black", size = 4.5) +
    scale_fill_gradientn(
      colours = c("white", "#b991d5", "#592181"),
      name = "MAE (Years)",
      limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)
    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev) +
    labs(x = "DNA Fragment Size (bp)", y = "DNA Input (ng)", title = title) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      panel.grid = element_blank()
    ) +
    coord_fixed()
}

# Create plots for the main clocks
age_pred_blup_plot <- plot_mae_table(mae_avg_df, BLUP, "BLUP Clock MAE")
age_pred_en_plot <- plot_mae_table(mae_avg_df, EN, "EN Clock MAE")
age_pred_horvath_plot <- plot_mae_table(mae_avg_df, Horvath, "Horvath Clock MAE")

# Combine plots into a single figure
combined_age_plot <- ggarrange(
  age_pred_blup_plot, age_pred_en_plot, age_pred_horvath_plot,
  ncol = 3, nrow = 1, common.legend = TRUE, legend = "right"
)

ggsave(
  file.path(results_age_dir, "mae_summary_plots.png"),
  combined_age_plot, width = 14, height = 5, dpi = 300, bg = "white"
)



#### 5. Assess Missing CpGs per Clock  --------------------------------------
# This new section checks how many CpGs required by each clock are missing
# from the detected probes in each sample.

message("Checking for missing CpGs for each clock...")

# Create a list to store results for each sample
missing_cpgs_list <- list()

# Loop over each sample column in the data frame prepared for age prediction
for (sample_name in colnames(dat_for_age)[-1]) { # Exclude the 'probeID' column
  
  # Create a temporary data frame for the single sample, removing NAs
  sample_data <- dat_for_age[, c("probeID", sample_name)]
  sample_data <- sample_data[!is.na(sample_data[, 2]), ]
  
  # Check for missing clock CpGs
  missing_info <- checkClocks(sample_data)
  
  # Extract the number of missing CpGs for each clock of interest
  # The indices correspond to specific clocks in the methylclock package
  missing_counts <- c(
    Horvath = length(missing_info$Horvath),
    skinHorvath = length(missing_info$skinHorvath),
    BLUP = length(missing_info$BLUP),
    EN = length(missing_info$EN)
  )
  
  missing_cpgs_list[[sample_name]] <- missing_counts
}

# Combine the list of results into a single data frame
missing_df <- do.call(rbind, missing_cpgs_list)
missing_df <- as.data.frame(missing_df) %>%
  tibble::rownames_to_column("Sample_code")

# Join with metadata to add sample condition information
missing_summary_df <- missing_df %>%
  left_join(metadata_subset, by = "Sample_code") %>%
  select(Sample_code, DNA_size, DNA_input, Duplicate, Horvath, skinHorvath, BLUP, EN)

# Save the detailed summary to an Excel file
write_xlsx(
  missing_summary_df,
  file.path(results_age_dir, "missing_cpgs_per_clock.xlsx")
)

message("Age prediction analysis complete. All results saved to ", results_age_dir)