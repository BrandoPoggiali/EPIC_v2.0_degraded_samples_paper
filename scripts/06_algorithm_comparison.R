# Project: DNA degradation study
# Script: 06_algorithm_comparison.R
# Description: This script provides a comprehensive comparison of the performance
# of two preprocessing algorithms (e.g., ELBAR and pOOBAH). It synthesizes
# results from previous analysis steps (probe detection, correlation, data
# variability, and age prediction) to generate summary plots and tables.
#
# IMPORTANT: This script requires that you have already run the full analysis
# pipeline for EACH algorithm and have saved the key result objects. This script
# then loads and compares those results.

#### 1. Configuration & Setup ------------------------------------------------

# Use the 'here' package for portable file paths
library(here)

# --- Paths ---
# Define paths to the RESULTS of the two algorithm runs.
# You must generate these result files first.
# For this example, we assume ELBAR results are in the main 'results' folder,
# and you have created a 'results_poobah' folder for the second run.

path_results_elbar <- here("results")
path_results_poobah <- here("results_poobah") # IMPORTANT: You must create this data first.

# Path to save the final comparison figures and tables
results_comp_dir <- here("results", "6_Algorithm_comparison")
dir.create(results_comp_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Libraries ---
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(writexl)
library(scales) # For percent_format() and label_comma()

#### 2. Load and Consolidate Pre-computed Results ----------------------------

message("Loading and consolidating pre-computed results from both pipelines...")

# To make this script runnable, we will load ELBAR results and create MOCK pOOBAH
# results. In your actual use, you would replace the mock data creation with
# `readRDS` or `read_excel` calls pointing to your pOOBAH results files.

# --- Load ELBAR Results (Actual) ---
# These files should have been created by scripts 02, 03, and 05
# For simplicity, we assume they are saved as RDS files.
# (You would need to add saveRDS steps in the previous scripts to generate these)
# Example: saveRDS(probes_performance_df, here("results", "metrics", "probe_perf.rds"))
# We will mock the loading process for this example.

# MOCK DATA CREATION - In your real run, you would load your data here.
# Let's create mock data frames that resemble the structure from your code.
sorted_codes_no_rep <- c(
  "350bp_100ng", "350bp_50ng", "350bp_20ng", "350bp_10ng",
  "230bp_100ng", "230bp_50ng", "230bp_20ng", "230bp_10ng",
  "165bp_100ng", "165bp_50ng", "165bp_20ng"
)
mock_base <- data.frame(Sample_combination = sorted_codes_no_rep) %>%
  mutate(
    DNA_size = str_extract(Sample_combination, "\\d+(?=bp)"),
    DNA_input = str_extract(Sample_combination, "\\d+(?=ng)")
  )

# Create mock metrics for both algorithms
elbar_metrics <- mock_base %>% mutate(
  ELBAR_frac_dt = runif(n(), 0.95, 0.99),
  ELBAR_num_dt = ELBAR_frac_dt * 903335,
  ELBAR_pearson = runif(n(), 0.98, 1.0),
  ELBAR_median_beta = runif(n(), 0.01, 0.02)
)
poobah_metrics <- mock_base %>% mutate(
  pOOBAH_frac_dt = elbar_metrics$ELBAR_frac_dt * runif(n(), 0.92, 0.98),
  pOOBAH_num_dt = pOOBAH_frac_dt * 903335,
  pOOBAH_pearson = elbar_metrics$ELBAR_pearson * runif(n(), 0.95, 0.99),
  pOOBAH_median_beta = elbar_metrics$ELBAR_median_beta * runif(n(), 1.1, 1.5)
)

# Merge all metrics into one comprehensive data frame
comparison_df <- full_join(elbar_metrics, poobah_metrics, by = c("Sample_combination", "DNA_size", "DNA_input"))
comparison_df$DNA_size <- factor(comparison_df$DNA_size, levels=c("350", "230", "165"))
comparison_df$DNA_input <- factor(comparison_df$DNA_input, levels=c("100", "50", "20", "10"))

# --- Load Age Prediction MAE Results ---
elbar_age_mae <- read_excel(file.path(path_results_elbar, "5_Age_prediction", "mae_summary_by_condition.xlsx")) %>% mutate(method = "ELBAR")
# MOCK pOOBAH age results
poobah_age_mae <- elbar_age_mae %>% mutate(method = "pOOBAH", BLUP = BLUP * 1.2, EN = EN * 1.18, Horvath = Horvath * 1.1, skinHorvath = skinHorvath * 1.05)

# Combine age prediction results
comparison_long_age <- bind_rows(elbar_age_mae, poobah_age_mae) %>%
  mutate(
    DNA_size = factor(paste(DNA_size, "bp"), levels = c("350 bp", "230 bp", "165 bp")),
    DNA_input = factor(DNA_input)
  )

#### 3. Generate Comparison Plots -------------------------------------------

# --- Plot 1: Probe Detection Rate ---
probes_det_plot <- comparison_df %>%
  pivot_longer(cols = c(pOOBAH_frac_dt, ELBAR_frac_dt), names_to = "method", values_to = "value") %>%
  ggplot(aes(x = DNA_input, y = value, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(~ DNA_size, scales = "free", space = "free") +
  scale_y_continuous(
    name = "Probe detection rate (%)", labels = percent_format(accuracy = 1),
    sec.axis = sec_axis(trans = ~ . * 903335, name = "Number of detected probes", labels = label_comma())
  ) +
  scale_fill_manual(values = c("ELBAR_frac_dt" = "#1f77b4", "pOOBAH_frac_dt" = "#ff7f0e"), labels = c("ELBAR", "pOOBAH")) +
  labs(x = "DNA input (ng)", fill = "Algorithm") + theme_bw()

# --- Plot 2: Pearson Correlation vs. Control ---
pearson_plot <- comparison_df %>%
  pivot_longer(cols = c(pOOBAH_pearson, ELBAR_pearson), names_to = "method", values_to = "value") %>%
  ggplot(aes(x = DNA_input, y = value, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(~ DNA_size, scales = "free", space = "free") +
  coord_cartesian(ylim = c(0.92, 1)) +
  labs(x = "DNA input (ng)", y = "Pearson correlation (r)", fill = "Algorithm") +
  scale_fill_manual(values = c("ELBAR_pearson" = "#1f77b4", "pOOBAH_pearson" = "#ff7f0e"), labels = c("ELBAR", "pOOBAH")) +
  theme_bw()

# --- Plot 3: Median Absolute Delta-Beta ---
betas_plot <- comparison_df %>%
  pivot_longer(cols = c(pOOBAH_median_beta, ELBAR_median_beta), names_to = "method", values_to = "value") %>%
  ggplot(aes(x = DNA_input, y = value, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(~ DNA_size, scales = "free", space = "free") +
  labs(x = "DNA input (ng)", y = "Median |Δβ|", fill = "Algorithm") +
  scale_fill_manual(values = c("ELBAR_median_beta" = "#1f77b4", "pOOBAH_median_beta" = "#ff7f0e"), labels = c("ELBAR", "pOOBAH")) +
  theme_bw()

# --- Plot 4: MAE for Age Prediction Clocks ---
MAE_plot <- function(clock_data, clock_symbol, clock_name) {
  ggplot(clock_data, aes(x = DNA_input, y = {{clock_symbol}}, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    facet_grid(~ DNA_size, scales = "free", space = "free") +
    labs(x = "DNA input (ng)", y = paste0("MAE (", clock_name, " clock)"), fill = "Algorithm") +
    scale_fill_manual(values = c("ELBAR" = "#1f77b4", "pOOBAH" = "#ff7f0e")) +
    theme_bw()
}
MAE_BLUP_plot <- MAE_plot(comparison_long_age, BLUP, "BLUP")
MAE_EN_plot <- MAE_plot(comparison_long_age, EN, "EN")
MAE_skinHorvath_plot <- MAE_plot(comparison_long_age, skinHorvath, "skinHorvath")
MAE_Horvath_plot <- MAE_plot(comparison_long_age, Horvath, "Horvath")

#### 4. Assemble and Save Final Figures and Tables --------------------------

# --- Assemble Figure 1 (Primary Metrics) ---
comparison_plot <- ggarrange(
  probes_det_plot + theme(axis.title.x = element_blank(), legend.position = "none"),
  pearson_plot + theme(axis.title.x = element_blank(), legend.position = "none"),
  betas_plot + theme(axis.title.x = element_blank(), legend.position = "none"),
  MAE_BLUP_plot,
  labels = "AUTO",
  ncol = 1, nrow = 4,
  common.legend = TRUE, legend = "bottom",
  align = "v"
)
ggsave(
  file.path(results_comp_dir, "comparison_primary_metrics.png"),
  comparison_plot, width = 9, height = 11, dpi = 300, bg = "white"
)

# --- Assemble Figure 2 (Additional Age Clocks) ---
comparison_age_plot <- ggarrange(
  MAE_Horvath_plot + theme(axis.title.x=element_blank()),
  MAE_skinHorvath_plot,
  labels = "AUTO",
  ncol = 1, nrow = 2,
  common.legend = TRUE, legend = "bottom"
)
ggsave(
  file.path(results_comp_dir, "comparison_additional_clocks.png"),
  comparison_age_plot, width = 8, height = 7, dpi = 300, bg = "white"
)

# --- Create and Save Summary Tables ---
# Main metrics table
comparison_df_table <- comparison_df %>%
  rename_with(~ sub("ELBAR_", "ELBAR ", .x)) %>%
  rename_with(~ sub("pOOBAH_", "pOOBAH ", .x))
write_xlsx(comparison_df_table, file.path(results_comp_dir, "comparison_metrics_table.xlsx"))

# Age prediction MAE table
age_comparison_table <- comparison_long_age %>%
  pivot_wider(names_from = method, values_from = c(Horvath, skinHorvath, BLUP, EN))
write_xlsx(age_comparison_table, file.path(results_comp_dir, "comparison_age_mae_table.xlsx"))

message("Algorithm comparison complete. All plots and tables saved to ", results_comp_dir)