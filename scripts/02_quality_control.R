# Project: DNA degradation study
# Script: 02_quality_control.R
# Description: Performs quality control analysis on the preprocessed methylation data.
# Generates QC plots and assesses sample quality based on various metrics.

#### 1. Configuration & Setup ------------------------------------------------

library(here)
source(here("scripts", "functions.R")) # Load custom functions

# --- Paths ---
idat_dir <- here("data", "idat")
processed_data_dir <- here("data", "processed")
results_qc_dir <- here("results", "1_Quality_control")
dir.create(results_qc_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Libraries ---
library(sesame)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(scales)
library(tidyr)
library(tibble)

#### 2. Load Preprocessed Data -------------------------------------------------

betas <- readRDS(file.path(processed_data_dir, "betas_preprocessed.rds"))
metadata <- readRDS(file.path(processed_data_dir, "metadata_preprocessed.rds"))

#### 3. Calculate QC Statistics with sesame -----------------------------------

# Calculate QC stats using the QC Sesame pipeline
qcs <- openSesame(idat_dir, prep = "", func = sesameQC_calcStats) #You can perfom QC adter preprocessing adding prep = "QCDPB"
qc_df <- do.call(rbind, lapply(qcs, as.data.frame))
qc_df$Sample_code <- metadata$Sample_code

# Sort the QC data frame to match the desired plot order (by fragment size, then input)
extract_and_sort <- function(codes) {
  bp_values <- as.numeric(sub("bp.*", "", codes))
  ng_values <- as.numeric(sub("^\\d+bp_", "", sub("ng.*", "", codes)))
  sorted_indices <- order(bp_values, ng_values, decreasing = TRUE)
  # Ensure the control sample is first
  return(codes[c(which(grepl("no_degraded", codes)), sorted_indices[-which(sorted_indices == which(grepl("no_degraded", codes)))])])
}
sorted_codes <- extract_and_sort(qc_df$Sample_code)
qc_df_sorted <- qc_df[match(sorted_codes, qc_df$Sample_code), ] %>%
  left_join(metadata, by = "Sample_code")

# Factor variables for consistent plotting
qc_df_sorted$DNA_size[is.na(qc_df_sorted$DNA_size)] <- "Control"
qc_df_sorted$DNA_size <- factor(qc_df_sorted$DNA_size, levels = c("Control", "350", "230", "165", "95"))
qc_df_sorted$DNA_input <- factor(qc_df_sorted$DNA_input, levels = c("250", "100", "50", "20", "10"))

#### 4. Generate QC Plots -----------------------------------------------------

# Plot 1: Background Signal Intensity
bg_intensities_plot <- ggplot(qc_df_sorted, aes(x = mean_oob_grn, y = mean_oob_red)) +
  geom_point(aes(color = DNA_size, shape = DNA_input), size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted') +
  labs(
    x = 'Mean green background signal intensity',
    y = 'Mean red background signal intensity',
    title = 'Background Signal Intensity',
    color = "DNA fragment size (bp)",
    shape = "DNA input (ng)"
  ) +
  scale_color_manual(values = c("Control" = "darkgrey", "350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1500), labels = label_comma()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1500), labels = label_comma()) +
  theme_bw()

ggsave(file.path(results_qc_dir, "bg_intensities_plot.png"), bg_intensities_plot, width = 6, height = 4.5, dpi = 600, bg = "white")


# Plot 2 & 3: Mean Signal Intensity and Red/Green Ratio (using the custom function)
# Re-format DNA_size for the bar_plotter function
qc_df_sorted_bp <- qc_df_sorted
qc_df_sorted_bp$DNA_size <- ifelse(qc_df_sorted_bp$DNA_size == "Control", "Control", paste(qc_df_sorted_bp$DNA_size, "bp"))
qc_df_sorted_bp$DNA_size <- factor(qc_df_sorted_bp$DNA_size, levels=c("Control", "350 bp", "230 bp", "165 bp", "95 bp"))

fig2b <- qc_df_sorted_bp %>% bar_plotter(y = mean_intensity, y_axis_title = "Mean signal intensity") +
  geom_hline(yintercept = 2000, colour = "red", linetype = "dashed", linewidth = 0.6)

fig2c <- qc_df_sorted_bp %>% bar_plotter(y = RGratio, y_axis_title = "Ratio red/green median signal intensities") +
  geom_hline(yintercept = 0.5, colour = "red", linetype = "dashed", size = 0.6) +
  geom_hline(yintercept = 2, colour = "red", linetype = "dashed", size = 0.6)

combined_qc_plot <- ggarrange(fig2b, fig2c, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave(file.path(results_qc_dir, "combined_intensity_ratio_plot.png"), combined_qc_plot, width = 12, height = 8, dpi = 300)


#### 5. Genotype (SNP) Probe Analysis ----------------------------------------
# This is a form of sample identity/quality check

rs_probes <- substr(rownames(betas), 1, 2) == "rs"
genotype_rs <- as.data.frame(betasCollapseToPfx(betas[rs_probes, ]))

genotype_rs_long <- genotype_rs %>%
  rownames_to_column(var = "Probes") %>%
  pivot_longer(cols = -Probes, names_to = "Sample", values_to = "Intensity") %>%
  mutate(Sample = factor(Sample, levels = sorted_codes))

genotype_plot <- ggplot(genotype_rs_long, aes(x = Probes, y = Sample, fill = Intensity)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(
    colours = c("#c02011", "#f2d93a", "#3fb11b"),
    name = "Beta Value"
  ) +
  theme_bw() +
  labs(x = "SNP Probes", y = "Samples", title = "Genotype Probe Heatmap") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(file.path(results_qc_dir, "genotype_heatmap.png"), genotype_plot, width = 14, height = 8.5)


#### 6. Remove specific samples from the analysis (Optional)----------------------------------------------

#Remove samples from betas
betas_subset <- betas[,!grepl("95bp", colnames(betas))]
betas_subset <- betas_subset[, !colnames(betas_subset) %in% c("165bp_10ng_A","165bp_10ng_B")]

#Remove NAs
betas_noNAs_subset <- na.omit(betas_subset)
betas_noNAs_cg_subset <- betas_noNAs_subset[grepl("cg", rownames(betas_noNAs_subset)),]
betas_noNAs_cg_subset_coll <- betasCollapseToPfx(as.matrix(betas_noNAs_cg_subset))

#With NAs
# betas_cg_subset <- betas_subset[grepl("cg", rownames(betas_subset)),]
# betas_cg_subset_coll <- betasCollapseToPfx(as.matrix(betas_cg_subset))

#Remove samples from metadata
metadata <- metadata[!grepl("95bp", metadata$Sample_code),]
metadata <- metadata[! metadata$Sample_code %in% c("165bp_10ng_A","165bp_10ng_B"),]
dim(metadata)
dim(betas_noNAs_cg_subset)


#Remove samples from vector with sorted names
sorted_codes <- sorted_codes[!grepl("95bp", sorted_codes)]
sorted_codes <- sorted_codes[! sorted_codes %in% c("165bp_10ng_A","165bp_10ng_B")]

betas_noNAs_cg_coll_sorted <- betas_noNAs_cg_subset_coll[, sorted_codes]

saveRDS(betas_noNAs_cg_subset, file.path(processed_data_dir, "betas_preprocessed_QC_samples.rds"))
saveRDS(metadata, file.path(processed_data_dir, "metadata_preprocessed_QC_samples.rds"))

#### 7. PCA on beta values ----------------------------------------

# Removal of probes with NAs (If needed all samples)
betas_noNAs <- na.omit(betas)
dim(betas_noNAs) # 6200 32
betas_noNAs_cg <- betas_noNAs[grepl("cg", rownames(betas_noNAs)),]


#Principal Component analysis, only on not missing probes which are present in all the samples
betas_noNAs_cg_coll <- betasCollapseToPfx(as.matrix(betas_noNAs_cg)) # Perform PCA with all samples or only QC passing samples
betas_noNAs_cg_subset_coll <- betasCollapseToPfx(as.matrix(betas_noNAs_cg_subset))
dim(betas_noNAs_cg_coll)

PCA_result <- prcomp(t(betas_noNAs_cg_coll), center = TRUE, scale. = TRUE) #creates the principal components
PCA_result_scores <- as.data.frame(PCA_result$x)
percentage <- round(PCA_result$sdev^2 / sum(PCA_result$sdev^2) * 100, 2)
percentage <- paste(colnames(PCA_result_scores), "(", paste( as.character(percentage), "%", ")", sep="") )

PCA_result_scores$Sample_code <- as.factor(metadata$Sample_code)
DNA_sizes <- metadata$DNA_size
DNA_sizes[1] <- "Not fragmented"
PCA_result_scores$DNA_size <- factor(DNA_sizes, levels=c("Not fragmented", "350", "230", "165")) # if you do with all samples: levels=c("Not fragmented", "350", "230", "165", "95")
PCA_result_scores$DNA_size <- factor(DNA_sizes, levels=c("Not fragmented", "350", "230", "165", "95")) # if you do with all samples: levels=c("Not fragmented", "350", "230", "165", "95")
PCA_result_scores$DNA_input <- factor(metadata$DNA_input, levels=c("250", "100", "50", "20", "10"))

PCA <- ggplot(PCA_result_scores,aes(x = PC1, y = PC2, color = DNA_size, shape = DNA_input)) +
  geom_point(size=3) + 
  scale_color_manual(values = c("Not fragmented" = "darkgrey", "350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000")) +
  theme_bw() +
  xlab(percentage[1]) + ylab(percentage[2]) +
  labs(color = "DNA size", shape = "DNA input", title= "6,160 CpGs") +
  theme(plot.title = element_text(hjust = 0.5)) #+

ggsave(paste0(results_path,"/1_Quality_control/PCA_on_subset_samples_200K_cg_24-03-2024.png"), 
       PCA, width = 7, height = 5.1)

PCA_all <- PCA
PCA_subset <- PCA

spacer <- ggplot() + theme_void()
combined_PCA_plot <- ggarrange(PCA_all, spacer, PCA_subset,
                               ncol = 3, nrow = 1,
                               legend ="right",
                               labels = c("A", "", "B"),
                               common.legend = TRUE,
                               widths = c(1, 0.2, 1))

ggsave(file.path(results_qc_dir, "PCA_all_samples_and_QC_passing_samples.png"), combined_PCA_plot, width = 6, height = 4.5, dpi = 600, bg = "white")



#### 8. Missing probes investigation ----------------------------------------

#Calculation of all missing probes number and recurrent for DNA size and DNA input condition
#get masked probes
library(VennDiagram)
masked_probes <- unique(unlist(KYCG_getDBs(sprintf("%s.Mask", "EPICv2"), 
                                           recommendedMaskNames()[["EPICv2"]], silent = TRUE, 
                                           ignore.case = TRUE)))

betas_non_masked <- betas_subset[!rownames(betas_subset) %in% masked_probes,]
nas_all_samples <- betas_non_masked[rowSums(is.na(betas_non_masked)) == ncol(betas_non_masked),]
dim(nas_all_samples)[1]
probes_name_all_samples <- rownames(nas_all_samples)

betas_non_masked_no_control <- betas_non_masked[,-1]
nas_all_samples_no_control <- betas_non_masked_no_control[rowSums(is.na(betas_non_masked_no_control)) == ncol(betas_non_masked_no_control),]
dim(nas_all_samples_no_control)[1]
probes_name_all_samples_no_control <- rownames(nas_all_samples_no_control)

betas_non_masked_350bp <- betas_non_masked_no_control[,grepl("350", colnames(betas_non_masked_no_control))]
nas_all_samples_350bp <- betas_non_masked_350bp[rowSums(is.na(betas_non_masked_350bp)) == ncol(betas_non_masked_350bp),]
dim(nas_all_samples_350bp)[1]
probes_name_all_samples_350 <- rownames(nas_all_samples_350bp)

betas_non_masked_230bp <- betas_non_masked_no_control[,grepl("230", colnames(betas_non_masked_no_control))]
nas_all_samples_230bp <- betas_non_masked_230bp[rowSums(is.na(betas_non_masked_230bp)) == ncol(betas_non_masked_230bp),]
dim(nas_all_samples_230bp)[1]
probes_name_all_samples_230 <- rownames(nas_all_samples_230bp)

betas_non_masked_165bp <- betas_non_masked_no_control[,grepl("165", colnames(betas_non_masked_no_control))]
nas_all_samples_165bp <- betas_non_masked_165bp[rowSums(is.na(betas_non_masked_165bp)) == ncol(betas_non_masked_165bp),]
dim(nas_all_samples_165bp)[1]
probes_name_all_samples_165 <- rownames(nas_all_samples_165bp)


venn.plot_bp <- venn.diagram(
  x = list("350 bp" = probes_name_all_samples_350, "230 bp" = probes_name_all_samples_230,
           "165 bp" = probes_name_all_samples_165),
  filename = NULL,  # Plot in RStudio viewer
  fill = c("#93C03F", "#3F93C0", "#C03F93"),
  alpha = 0.5,
  cat.col = c("black", "black", "black"),
  cat.cex = 1
)

idx <- sapply(venn.plot_bp, function(i) grepl("text", i$name))

for(i in 1:7){
  venn.plot_bp[idx][[i]]$label <- 
    format(as.numeric(venn.plot_bp[idx][[i]]$label), big.mark=",", scientific=FALSE)
}
grid.newpage()
grid.draw(venn.plot_bp) 


#Per concentrations
betas_non_masked_100ng <- betas_non_masked_no_control[,grepl("100ng", colnames(betas_non_masked_no_control))]
nas_all_samples_100ng <- betas_non_masked_100ng[rowSums(is.na(betas_non_masked_100ng)) == ncol(betas_non_masked_100ng),]
dim(nas_all_samples_100ng)[1]
probes_name_all_samples_100 <- rownames(nas_all_samples_100ng)

betas_non_masked_50ng <- betas_non_masked_no_control[,grepl("50ng", colnames(betas_non_masked_no_control))]
nas_all_samples_50ng <- betas_non_masked_50ng[rowSums(is.na(betas_non_masked_50ng)) == ncol(betas_non_masked_50ng),]
dim(nas_all_samples_50ng)[1]
probes_name_all_samples_50 <- rownames(nas_all_samples_50ng)

betas_non_masked_20ng <- betas_non_masked_no_control[,grepl("20ng", colnames(betas_non_masked_no_control))]
nas_all_samples_20ng <- betas_non_masked_20ng[rowSums(is.na(betas_non_masked_20ng)) == ncol(betas_non_masked_20ng),]
dim(nas_all_samples_20ng)[1]
probes_name_all_samples_20 <- rownames(nas_all_samples_20ng)

betas_non_masked_10ng <- betas_non_masked_no_control[,grepl("10ng", colnames(betas_non_masked_no_control))]
nas_all_samples_10ng <- betas_non_masked_10ng[rowSums(is.na(betas_non_masked_10ng)) == ncol(betas_non_masked_10ng),]
dim(nas_all_samples_10ng)[1]
probes_name_all_samples_10 <- rownames(nas_all_samples_10ng)

venn.plot_ng <- venn.diagram(
  x = list("100 ng" = probes_name_all_samples_100, "50 ng" = probes_name_all_samples_50,
           "20 ng" = probes_name_all_samples_20, "10 ng" = probes_name_all_samples_10),
  filename = NULL,  # Plot in RStudio viewer
  fill = c("#7B53AC", "#AC5357", "#84AC53", "#53ACA8"),
  alpha = 0.5,
  cat.col = c("black", "black", "black", "black"),
  cat.cex = 0.5)

idx <- sapply(venn.plot_ng, function(i) grepl("text", i$name))

for(i in 1:15){
  venn.plot_ng[idx][[i]]$label <- 
    format(as.numeric(venn.plot_ng[idx][[i]]$label), big.mark=",", scientific=FALSE)
}
grid.newpage()
grid.draw(venn.plot_ng) 

#All plots toghether
spacer <- ggplot() + theme_void()
combined_plot_nas <- ggarrange(spacer, venn.plot_bp , spacer, venn.plot_ng, spacer,
                               labels = c("A","","B", "", ""),
                               ncol = 5, nrow = 1,
                               widths = c(0.1, 1, 0.1, 1, 0.1))
combined_plot_nas

ggsave(file.path(results_qc_dir, "Missing_data_per_bp_and_ng.png"), combined_plot_nas, width = 6, height = 4.5, dpi = 600, bg = "white")

#Create excel file with missing probes
# Combine vectors into a data frame with different lengths
max_length <- max(length(probes_name_all_samples_no_control), length(probes_name_all_samples),
                  length(probes_name_all_samples_350), length(probes_name_all_samples_230), 
                  length(probes_name_all_samples_165),length(probes_name_all_samples_100), 
                  length(probes_name_all_samples_50), length(probes_name_all_samples_20),
                  length(probes_name_all_samples_10))



# Fill shorter vectors with NA
data <- data.frame(
  "Control and degraded samples" = c(probes_name_all_samples, rep(NA, max_length - length(probes_name_all_samples))),
  "Degraded samples" = c(probes_name_all_samples_no_control, rep(NA, max_length - length(probes_name_all_samples_no_control))),
  "350_bp" = c(probes_name_all_samples_350, rep(NA, max_length - length(probes_name_all_samples_350))),
  "230_bp" = c(probes_name_all_samples_230, rep(NA, max_length - length(probes_name_all_samples_230))),
  "165_bp" = c(probes_name_all_samples_165, rep(NA, max_length - length(probes_name_all_samples_165))),
  "100_ng" = c(probes_name_all_samples_100, rep(NA, max_length - length(probes_name_all_samples_100))),
  "50_ng" = c(probes_name_all_samples_50, rep(NA, max_length - length(probes_name_all_samples_50))),
  "20_ng" = c(probes_name_all_samples_20, rep(NA, max_length - length(probes_name_all_samples_20))),
  "10_ng" = c(probes_name_all_samples_10, rep(NA, max_length - length(probes_name_all_samples_10)))
)

# Export to Excel

write.xlsx(data, file = paste0(results_qc_dir,"/Missing_probes.xlsx"),
           sheetName = "Missing data", rowNames = FALSE)


message("QC analysis complete. Plots saved to ", results_qc_dir)


#### 9. Check Beta value distribution (This step require high amount of memory). -----------------------------
beta_vals_long <- betas %>% as.tibble() %>% #If you want to test type I or II probes just change the initial dataset
  add_column(CpG = rownames(betas)) %>%
  gather(Sample_code, beta_value, 1:33) 

colnames(beta_vals_long)[1] <- "Probe_ID"

beta_vals_long$Sample_code[beta_vals_long$Sample_code == "no_degraded_250ng"] <- "Control"

#Crete column for plotting
beta_vals_long <- beta_vals_long %>%
  mutate(Sample_group = ifelse(substr(beta_vals_long$Sample_code, nchar(beta_vals_long$Sample_code) - 1, nchar(beta_vals_long$Sample_code)) %in% c("_A", "_B"),
                               substr(Sample_code, 1, nchar(Sample_code) - 2),
                               Sample_code))

sorted_codes[1] <- "Control"

beta_vals_long$Sample_code <- factor(beta_vals_long$Sample_code, level=sorted_codes)
beta_vals_long$Sample_group <- factor(beta_vals_long$Sample_group, level=c(sorted_codes[1],
                                                                           unique(substring(sorted_codes[-1], 1, nchar(sorted_codes[-1])-2))))

# Plot DNA sizes of 350 bp
unique_samples <- sort(unique(beta_vals_long$Sample_group))
unique_samples <- as.character(unique_samples)

samples_distribution_350bp <- ggplot(data = beta_vals_long[substr(beta_vals_long$Sample_code,1,3) %in% c("350","Con"),]) +
  geom_density(aes(x = beta_value, color = Sample_code), size = 1) +
  labs(x = "Beta value", y = "Density", title="350 bp") +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0)) + theme_bw() +
  scale_color_manual(values = setNames(c("darkgrey", "black","black", "#f2d93a","#f2d93a", 
                                         "#39a2a4","#39a2a4", "#f20000", "#f20000"),
                                       c("Control", sorted_codes[grepl("350", sorted_codes)])),
                     labels = c("Control", substring(sorted_codes[grepl("350", sorted_codes)],7,nchar(sorted_codes[grepl("350", sorted_codes)])-2))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color="Sample type") 


# Plot DNA sizes of 230 bp
samples_distribution_230bp <- ggplot(data = beta_vals_long[substr(beta_vals_long$Sample_code,1,3) %in% c("230","Con"),]) +
  geom_density(aes(x = beta_value, color = Sample_code), size = 1) +
  labs(x = "Beta value", y = "Density", title="230 bp") +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0)) + theme_bw() +
  scale_color_manual(values = setNames(c( "darkgrey", "black","black", "#f2d93a","#f2d93a", 
                                          "#39a2a4","#39a2a4", "#f20000", "#f20000"),
                                       c("Control", sorted_codes[grepl("230", sorted_codes)])),
                     labels = c("Control", substring(sorted_codes[grepl("230", sorted_codes)],7,nchar(sorted_codes[grepl("230", sorted_codes)])-2))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color="Sample type") 

# Plot DNA sizes of 165 bp
samples_distribution_165bp <- ggplot(data = beta_vals_long[substr(beta_vals_long$Sample_code,1,3) %in% c("165","Con"),]) +
  geom_density(aes(x = beta_value, color = Sample_code), size = 1) +
  labs(x = "Beta value", y = "Density", title="165 bp") +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0)) + theme_bw() +
  scale_color_manual(values = setNames(c( "darkgrey", "black","black", "#f2d93a","#f2d93a", 
                                          "#39a2a4","#39a2a4", "#f20000", "#f20000"),
                                       c("Control", sorted_codes[grepl("165", sorted_codes)])),
                     labels = c("Control", substring(sorted_codes[grepl("165", sorted_codes)],7,nchar(sorted_codes[grepl("165", sorted_codes)])-2))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color="Sample type") 

# Plot DNA sizes of 95 bp
samples_distribution_95bp <- ggplot(data = beta_vals_long[substr(beta_vals_long$Sample_code,1,2) %in% c("95","Co"),]) +
  geom_density(aes(x = beta_value, color = Sample_code), size = 1) +
  labs(x = "Beta value", y = "Density", title="95 bp") +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0)) + theme_bw() +
  scale_color_manual(values = setNames(c( "darkgrey", "black","black", "#f2d93a","#f2d93a", 
                                          "#39a2a4","#39a2a4", "#f20000", "#f20000"),
                                       c("Control", sorted_codes[grepl("95", sorted_codes)])),
                     labels = c("Control", substring(sorted_codes[grepl("95", sorted_codes)],6,nchar(sorted_codes[grepl("95", sorted_codes)])-2))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color="Sample type") 

spacer <- ggplot() + theme_void()
combined_plot <- ggarrange(samples_distribution_350bp, spacer, samples_distribution_230bp, spacer,
                           samples_distribution_165bp, spacer, samples_distribution_95bp, spacer,
                           labels = c("A", "", "B","", "C", "", "D", ""),
                           ncol = 4, nrow = 2,
                           common.legend = TRUE,
                           legend = "bottom",
                           widths = c(1, 0.01, 1, 0.01))

ggsave(file.path(results_qc_dir, "Beta_value_distributions.png"), combined_plot_nas, width = 6, height = 4.5, dpi = 600, bg = "white")



#### 10. Probe detection rate ----------------------------------------
unique(substring(row.names(betas_subset),1 ,2))
dim(betas_subset)
betas_subset <- betas_subset[substring(row.names(betas_subset),1 ,2) %in% c("cg", "ch", "rs"),]
dim(betas_subset)

frac_df <- colSums(!is.na(betas_subset))/(nrow(betas_subset)-32896 ) #Total number of probes 903335
frac_df <- frac_df[match(sorted_codes[1:23], names(frac_df))]
num_df <- colSums(!is.na(betas_subset))
num_df <- num_df[match(sorted_codes[1:23], names(num_df))]

probes_performance_df <- data.frame(
  Sample_code = names(frac_df),
  frac_dt = as.numeric(frac_df),
  num_dt = as.numeric(num_df)
)

probes_performance_df <- probes_performance_df[-1,] # remove degraded samples, %fraction: 0.9998627, and 903211 probes 
probes_performance_df$Sample_combination <- substr(probes_performance_df$Sample_code, 1, nchar(probes_performance_df$Sample_code) - 2)

probes_performance_df <- probes_performance_df %>%
  group_by(Sample_combination) %>%
  summarize(across(.cols = everything(), .fns = mean, na.rm = TRUE), .groups = "drop")

probes_performance_df$DNA_size <- str_extract(probes_performance_df$Sample_combination, "\\d+(?=bp)")
probes_performance_df$DNA_input <- str_extract(probes_performance_df$Sample_combination, "\\d+(?=ng)")

probes_performance_df$DNA_size <- factor(probes_performance_df$DNA_size, level=c("350", "230", "165"))
probes_performance_df$DNA_input <- factor(probes_performance_df$DNA_input, level=c("100", "50", "20", "10"))

probe_success_rate_plot <- plot_metric_table(probes_performance_df, x_variable=probes_performance_df$DNA_size, 
                                      y_variable=probes_performance_df$DNA_input, frac_succ_probes=probes_performance_df$frac_dt,
                                      number_of_probes=probes_performance_df$num_dt,
                                      x_axis_title = "DNA fragment size (bp)", y_axis_title = "DNA input (ng)")


ggsave(file.path(results_qc_dir, "Plot_success_rate_and_n_probes.png"), probe_success_rate_plot, width = 6, height = 4.5, dpi = 600, bg = "white")

