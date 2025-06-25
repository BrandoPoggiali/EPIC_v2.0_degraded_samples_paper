# Project: DNA degradation study
# Script: 04_genomic_feature_analysis.R
# Description: This script analyzes the genomic distribution of CpG probes that are
# successfully detected in each sample. It assesses potential biases related to
# CpG island context, probe type (I vs. II), and gene feature location.

#### 1. Configuration & Setup ------------------------------------------------

# Use the 'here' package for portable file paths
library(here)

# --- Paths ---
processed_data_dir <- here("data", "processed")
results_dist_dir <- here("results", "4_CpG_proportion")
dir.create(results_dist_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Libraries ---
# Data manipulation and plotting
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)

# Annotation libraries
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(AnnotationHub)
library(ensembldb)


#### 2. Load and Prepare Data ------------------------------------------------

# Load preprocessed beta values and metadata
betas_full <- readRDS(file.path(processed_data_dir, "betas_preprocessed.rds"))
metadata <- readRDS(file.path(processed_data_dir, "metadata_preprocessed.rds"))

# --- Filter Samples that Passed QC ---
samples_to_keep_mask <- !grepl("95bp", colnames(betas_full)) &
  !colnames(betas_full) %in% c("165bp_10ng_A", "165bp_10ng_B")

betas_subset <- betas_full[, samples_to_keep_mask]
metadata_subset <- metadata[samples_to_keep_mask, ]
sorted_codes <- metadata_subset$Sample_code

# --- Prepare Beta Matrix for Analysis ---
# Collapse CpG probes, keeping NAs for now
betas_cg_coll <- betasCollapseToPfx(betas_subset[grepl("^cg", rownames(betas_subset)), ])
betas_cg_coll_sorted <- betas_cg_coll[, sorted_codes]

# Rename control for clarity
colnames(betas_cg_coll_sorted)[1] <- "Control"
sorted_codes[1] <- "Control"


#### 3. Load and Prepare Full Annotation Data -----------------------------

message("Loading and preparing comprehensive probe annotations...")

# --- Basic Annotations (CpG Island, Probe Type) ---
ann_epic_v2 <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
# Create a collapsed CpG name for matching (e.g., cg00000029)
ann_epic_v2$Name_cpg <- ifelse(
  grepl("cg", ann_epic_v2$Name),
  substring(ann_epic_v2$Name, 1, 10),
  NA
)

# --- Gene Feature Annotations (Promoter, Gene Body, Intergenic) ---
# Use AnnotationHub to get the latest human gene database
ah <- AnnotationHub()
# Query for the Ensembl database for Homo sapiens, release 111
ensdb_query <- query(ah, c("EnsDb", "Homo sapiens", "111"))
edb <- ensdb_query[[1]] # Load the newest database

# Get gene and promoter locations as GRanges objects
tx <- transcripts(edb)
prom <- promoters(edb, upstream = 2000, downstream = 200)

# Create a GRanges object for the EPICv2 probes
probes_gr <- GRanges(
  seqnames = ann_epic_v2$chr,
  ranges = IRanges(start = ann_epic_v2$pos, end = ann_epic_v2$pos),
  strand = ann_epic_v2$strand,
  Name = ann_epic_v2$Name_cpg
)

# Find overlaps between probes and promoters/gene bodies
overlaps_prom <- findOverlaps(probes_gr, prom)
probes_gr$is_promoter <- FALSE
probes_gr$is_promoter[queryHits(overlaps_prom)] <- TRUE

overlaps_tx <- findOverlaps(probes_gr, tx)
probes_gr$is_gene_body <- FALSE
probes_gr$is_gene_body[queryHits(overlaps_tx)] <- TRUE

# Define gene location based on overlaps
gene_location_df <- as.data.frame(mcols(probes_gr)) %>%
  mutate(
    Gene_location = case_when(
      is_promoter ~ "Promoter",
      is_gene_body ~ "Gene Body",
      TRUE ~ "Intergenic"
    )
  ) %>%
  select(Name, Gene_location) %>%
  distinct() # Keep unique mappings

# Combine all annotations into one master data frame
ann_complete <- ann_epic_v2 %>%
  left_join(gene_location_df, by = c("Name_cpg" = "Name")) %>%
  filter(!is.na(Name_cpg))

message("Annotation preparation complete.")

#### 4. Analyze Probe Distribution per Sample -------------------------------

# Prepare data frames to store results
island_results <- data.frame()
type_results <- data.frame()
gene_loc_results <- data.frame()

# Loop through each sample
for (sample_id in colnames(betas_cg_coll_sorted)) {
  # Get names of probes detected in this sample (not NA)
  detected_probes <- names(which(!is.na(betas_cg_coll_sorted[, sample_id])))
  
  # Subset the annotation to only the detected probes
  ann_subset <- ann_complete %>%
    filter(Name_cpg %in% detected_probes) %>%
    distinct(Name_cpg, .keep_all = TRUE) # Ensure one entry per CpG
  
  # --- Tally proportions ---
  # CpG Island Relation
  island_counts <- as.data.frame(table(ann_subset$Relation_to_Island))
  island_counts$Sample <- sample_id
  island_results <- rbind(island_results, island_counts)
  
  # Probe Type
  type_counts <- as.data.frame(table(ann_subset$Type))
  type_counts$Sample <- sample_id
  type_results <- rbind(type_results, type_counts)
  
  # Gene Location
  gene_loc_counts <- as.data.frame(table(ann_subset$Gene_location))
  gene_loc_counts$Sample <- sample_id
  gene_loc_results <- rbind(gene_loc_results, gene_loc_counts)
}

#### 5. Visualize Genomic Feature Distributions -----------------------------

# --- Plot 1: CpG Island Relation ---
island_plot_data <- island_results %>%
  group_by(Sample) %>%
  mutate(Proportion = Freq / sum(Freq), PercentLabel = paste0(round(Proportion * 100), "%")) %>%
  ungroup() %>%
  mutate(Sample = factor(Sample, levels = sorted_codes))

island_proportion_plot <- ggplot(island_plot_data, aes(x = Sample, y = Proportion, fill = Var1)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = PercentLabel), position = position_fill(vjust = 0.5), size = 3, color = "white") +
  labs(
    x = "Sample", y = "Relative Proportion of Detected Probes",
    title = "Distribution by CpG Island Context", fill = "Genomic Location"
  ) +
  scale_fill_manual(values = c("Island" = "#3D69C2", "N_Shore" = "#3DC253", "S_Shore" = "#64C23D", "N_Shelf" = "#C2963D", "S_Shelf" = "#C2783D", "OpenSea" = "#C23DAC")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(file.path(results_dist_dir, "cpg_island_distribution.png"), island_proportion_plot, width = 8, height = 6)


# --- Plot 2: Probe Type (I vs II) ---
type_plot_data <- type_results %>%
  group_by(Sample) %>%
  mutate(Proportion = Freq / sum(Freq), PercentLabel = paste0(round(Proportion * 100), "%")) %>%
  ungroup() %>%
  mutate(Sample = factor(Sample, levels = sorted_codes))

type_proportion_plot <- ggplot(type_plot_data, aes(x = Sample, y = Proportion, fill = Var1)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = PercentLabel), position = position_fill(vjust = 0.5), size = 3, color = "white") +
  labs(
    x = "Sample", y = "Relative Proportion of Detected Probes",
    title = "Distribution by Probe Type", fill = "Probe Type"
  ) +
  scale_fill_manual(values = c("I" = "#664DB2", "II" = "#99B24D")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(file.path(results_dist_dir, "probe_type_distribution.png"), type_proportion_plot, width = 8, height = 6)


# --- Plot 3: Gene Feature Location ---
gene_loc_plot_data <- gene_loc_results %>%
  group_by(Sample) %>%
  mutate(Proportion = Freq / sum(Freq), PercentLabel = paste0(round(Proportion * 100), "%")) %>%
  ungroup() %>%
  mutate(Sample = factor(Sample, levels = sorted_codes))

gene_location_plot <- ggplot(gene_loc_plot_data, aes(x = Sample, y = Proportion, fill = Var1)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = PercentLabel), position = position_fill(vjust = 0.5), size = 3, color = "white") +
  labs(
    x = "Sample", y = "Relative Proportion of Detected Probes",
    title = "Distribution by Gene Feature Location", fill = "Gene Location"
  ) +
  scale_fill_manual(values = c("Promoter" = "#38C781", "Gene Body" = "#8138C7", "Intergenic" = "#C78138")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(file.path(results_dist_dir, "gene_location_distribution.png"), gene_location_plot, width = 8, height = 6)


# --- Combine Plots into a Single Figure ---
combined_distribution_plot <- ggarrange(
  island_proportion_plot, type_proportion_plot, gene_location_plot,
  ncol = 1, nrow = 3,
  common.legend = FALSE, # Each plot has its own legend
  labels = "AUTO"
)

ggsave(
  file.path(results_dist_dir, "combined_genomic_distributions.png"),
  combined_distribution_plot, width = 10, height = 15, dpi = 300, bg = "white"
)

message("Genomic feature analysis complete. Plots saved to ", results_dist_dir)