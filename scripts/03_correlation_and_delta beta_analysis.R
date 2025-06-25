# Project: DNA degradation study
# Script: 03_correlation_and_variability_analysis.R
# Description: This script performs correlation and delta-beta variability
# analysis to assess the precision and accuracy of methylation measurements
# under different degradation conditions compared to a control sample.

#### 1. Configuration & Setup ------------------------------------------------

# Use the 'here' package for portable file paths
library(here)

# --- Paths ---
processed_data_dir <- here("data", "processed")
results_corr_dir <- here("results", "2_Correlation")
results_delta_beta_dir <- here("results", "3_Delta_beta_value")
dir.create(results_corr_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_delta_beta_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Libraries ---
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2) # For melt()
library(ggpubr)
library(ggridges)
library(corrplot)
library(writexl)


#### 2. Load and Prepare Data ------------------------------------------------

# Load preprocessed beta values and metadata
betas_full <- readRDS(file.path(processed_data_dir, "betas_preprocessed.rds"))
metadata <- readRDS(file.path(processed_data_dir, "metadata_preprocessed.rds"))

# --- Filter Samples that Failed QC ---
# This step removes samples that were deemed low quality in the QC analysis.
samples_to_keep_mask <- !grepl("95bp", colnames(betas_full)) &
  !colnames(betas_full) %in% c("165bp_10ng_A", "165bp_10ng_B")

betas_subset <- betas_full[, samples_to_keep_mask]
metadata_subset <- metadata[samples_to_keep_mask, ]

# Order samples for plotting
sorted_codes <- metadata_subset$Sample_code

# --- Prepare Beta Matrix for Analysis ---
# 1. Keep only CpG probes ("cg")
# 2. Remove any probes that have NA in ANY of the remaining samples to ensure a common set of probes
# 3. Collapse probes with the same ID (e.g. from different bead types)
betas_cg <- betas_subset[grepl("^cg", rownames(betas_subset)), ]
betas_cg_no_na <- na.omit(betas_cg) # Creates a matrix of probes common to all samples
betas_cg_no_na_coll <- betasCollapseToPfx(as.matrix(betas_cg_no_na))

# Sort the columns to match the desired order
betas_sorted <- betas_cg_no_na_coll[, sorted_codes]

message(
  "Data prepared for analysis. Using ",
  format(nrow(betas_sorted), big.mark = ","),
  " common CpG probes across ",
  ncol(betas_sorted), " samples."
)

#### 3. Pearson Correlation (vs. Control) -----------------------------------

# Average the technical replicates to get one value per condition
sample_conditions <- gsub("_A|_B", "", colnames(betas_sorted))
unique_conditions <- unique(sample_conditions)

# Create an empty data frame to store averaged betas
average_data <- data.frame(matrix(
  ncol = length(unique_conditions),
  nrow = nrow(betas_sorted),
  dimnames = list(rownames(betas_sorted), unique_conditions)
))

# Loop through unique conditions and average the duplicates
for (cond in unique_conditions) {
  replicate_cols <- which(sample_conditions == cond)
  if (length(replicate_cols) > 1) {
    average_data[, cond] <- rowMeans(betas_sorted[, replicate_cols], na.rm = TRUE)
  } else {
    average_data[, cond] <- betas_sorted[, replicate_cols]
  }
}

# --- Create Correlation Matrix and Heatmap ---
cor_matrix <- cor(average_data, method = "pearson")
melted_cor_matrix <- melt(cor_matrix)

# Rename the control for clarity in the plot
colnames(cor_matrix)[1] <- "Control"
rownames(cor_matrix)[1] <- "Control"

corr_heatmap <- ggplot(data = melt(cor_matrix), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(colour = "#2b2b2b") +
  geom_text(aes(label = round(value, 3)), color = "black", size = 4.5) +
  scale_fill_gradient2(
    low = "white", high = "red", mid = "pink",
    midpoint = 0.95, limit = c(0.90, 1.00),
    name = bquote('Pearson r')
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.title = element_text(size = 12.5)
  ) +
  coord_fixed()

ggsave(
  file.path(results_corr_dir, "pearson_correlation_heatmap.png"),
  corr_heatmap, width = 9, height = 8, dpi = 300, bg = "white"
)


##pearson correlation between replicates
sorted_codes <- gsub("_A|_B", "", sorted_codes)
sorted_codes <- unique(sorted_codes[-1])

betas_cg <- betas[substr(rownames(betas),1,2) %in% c("cg"),]

correlation_values <- vector()
p_values <- vector()
n_probes <- vector()

#code="350bp_100ng"
for (code in sorted_codes){
  print(code)
  df_tmp <- betas_cg[, grepl(code, colnames(betas))]
  df_tmp <- na.omit(df_tmp)
  
  result <- cor.test(df_tmp[,1], df_tmp[,2], method="pearson")
  
  correlation_values <- c(correlation_values, result$estimate)
  p_values <- c(p_values,result$p.value)
  
  current_n_probes <- nrow(df_tmp)
  n_probes <- c(n_probes, current_n_probes)
  
}

names(correlation_values) <- sorted_codes
names(p_values) <- sorted_codes
names(n_probes) <- sorted_codes

correlation_values <- correlation_values^2
square_pearson_corr_replicates <- rbind(correlation_values, p_values, n_probes)
square_pearson_corr_replicates <- as.data.frame(square_pearson_corr_replicates)
square_pearson_corr_replicates$Rownames <- rownames(square_pearson_corr_replicates)

write_xlsx(square_pearson_corr_replicates, paste0(results_corr_dir, "/pearson_corr_replicates.xlsx"))

#### 4. Delta Beta (Δβ) Analysis (vs. Control) ------------------------------

# Transpose the averaged data for easier calculations
average_data_t <- as.data.frame(t(average_data))

# The first row is our control sample
control_betas <- average_data_t[1, ]

# Calculate delta beta: (beta_sample - beta_control) for all probes
delta_beta_t <- sweep(average_data_t, 2, STATS = unlist(control_betas), FUN = "-")
delta_beta <- as.data.frame(t(delta_beta_t))
delta_beta <- delta_beta[-1, ] # Remove the control row (which is all zeros)

# --- Plot Distribution of Absolute Delta Beta Values ---
delta_beta_long <- delta_beta %>%
  tibble::rownames_to_column("Probe_ID") %>%
  pivot_longer(
    cols = -Probe_ID,
    names_to = "Sample_combination",
    values_to = "Delta_beta"
  ) %>%
  mutate(
    DNA_size = factor(str_extract(Sample_combination, "\\d+(?=bp)")),
    Sample_combination = factor(
      gsub("_", " ", Sample_combination),
      levels = rev(gsub("_", " ", unique(sample_conditions)[-1]))
    )
  )


abs_beta_dist_plot <- ggplot(delta_beta_long, aes(x = abs(Delta_beta), y = Sample_combination, fill = DNA_size)) +
  geom_density_ridges(
    alpha = 0.6, scale = 2.5,
    quantile_lines = TRUE, quantiles = 2, rel_min_height = 0.01
  ) +
  theme_bw() +
  scale_x_continuous(
    name = expression(paste(italic("|Δβ|"))),
    expand = c(0, 0), limits = c(0, 0.35)
  ) +
  scale_y_discrete(
    name = "Sample Condition",
    expand = expansion(mult = c(0.01, 0.1))
  ) +
  scale_fill_manual(
    values = c("350" = "black", "230" = "#f2d93a", "165" = "#39a2a4"),
    name = "DNA Fragment\nSize (bp)"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(
  file.path(results_delta_beta_dir, "abs_delta_beta_distribution.png"),
  abs_beta_dist_plot, width = 10, height = 8, dpi = 300, bg = "white"
)

# --- Calculate and Save Delta Beta Statistics ---
summary_stats <- delta_beta_long %>%
  group_by(Sample_combination, DNA_size) %>%
  summarize(
    Median_abs_delta_beta = median(abs(Delta_beta), na.rm = TRUE),
    IQR_abs_delta_beta = IQR(abs(Delta_beta), na.rm = TRUE),
    Percent_probes_over_0.05 = (sum(abs(Delta_beta) >= 0.05) / n()) * 100,
    Percent_probes_over_0.10 = (sum(abs(Delta_beta) >= 0.10) / n()) * 100,
    .groups = "drop"
  )

write_xlsx(summary_stats, file.path(results_delta_beta_dir, "delta_beta_summary_statistics.xlsx"))


#### 5. Delta Beta by β-value Interval -------------------------------------

# This analysis checks if variability depends on the original methylation level
control_beta_values <- average_data[, 1, drop = FALSE] # Get control betas

delta_beta_intervals <- delta_beta %>%
  tibble::rownames_to_column("Probe_ID") %>%
  pivot_longer(
    cols = c(-Probe_ID),
    names_to = "Sample_combination",
    values_to = "Delta_beta"
  ) %>%
  left_join(
    control_beta_values %>% tibble::rownames_to_column("Probe_ID"),
    by = "Probe_ID"
  ) %>%
  # Create bins for the control beta values
  mutate(
    beta_interval = cut(
      !!as.name(unique_conditions[1]), # Use the name of the control column
      breaks = seq(0, 1, by = 0.1),
      labels = paste(seq(0, 0.9, 0.1), seq(0.1, 1, 0.1), sep = "-"),
      include.lowest = TRUE
    ),
    DNA_size = factor(str_extract(Sample_combination, "\\d+(?=bp)")),
    DNA_input = factor(str_extract(Sample_combination, "\\d+(?=ng)"), levels=c("100", "50", "20", "10"))
  )

# Plot boxplots of absolute delta-beta for each interval, faceted by DNA input
interval_boxplot <- ggplot(
  delta_beta_intervals,
  aes(x = beta_interval, y = abs(Delta_beta), fill = DNA_size)
) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ DNA_input, ncol = 1, labeller = labeller(DNA_input = function(x) paste(x, "ng"))) +
  coord_cartesian(ylim = c(0, 0.4)) +
  labs(
    y = expression(paste(italic("|Δβ|"))),
    x = expression(paste("Control Sample β-value Intervals")),
    fill = "DNA Fragment\nSize (bp)"
  ) +
  scale_fill_manual(values = c("350" = "black", "230" = "#f2d93a", "165" = "#39a2a4")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing = unit(0.9, "lines")
  )

ggsave(
  file.path(results_delta_beta_dir, "delta_beta_by_methylation_interval.png"),
  interval_boxplot, width = 12, height = 10, dpi = 300, bg = "white"
)

message("Correlation and Delta Beta analysis complete. Results saved.")



## Delta beta values between replicates
sorted_codes <- qc_df_sorted$Sample_code
sorted_codes <- gsub("_A|_B", "", sorted_codes)
sorted_codes <- unique(sorted_codes[-1])

betas_cg <- betas[substr(rownames(betas),1,2) %in% c("cg"),]

delta_beta_values <- data.frame(matrix(nrow=5, ncol=length(sorted_codes)))
rownames(delta_beta_values) <- c("n_probes","Median", "IQR", "beta_05", "beta_10")

code="350bp_100ng"
n <- 1
for (code in sorted_codes){
  print(code)
  df_tmp <- betas_cg[, grepl(code, colnames(betas))]
  df_tmp <- na.omit(df_tmp)
  
  result <- abs(df_tmp[,1] - df_tmp[,2])
  delta_beta_values["n_probes", n] <- length(result)
  delta_beta_values["Median", n] <- median(result)
  delta_beta_values["IQR", n] <- IQR(result)
  delta_beta_values["beta_05", n] <- (sum(result >= 0.05))/(length(result))
  delta_beta_values["beta_10", n] <- (sum(result >= 0.1))/(length(result))
  colnames(delta_beta_values)[n] <- code
  n <- n + 1  
}

delta_beta_values <- cbind(RowName = rownames(delta_beta_values), delta_beta_values)
write_xlsx(results_delta_beta_dir, paste0(results_path, "Delta_beta_statistics_replicates.xlsx"))


