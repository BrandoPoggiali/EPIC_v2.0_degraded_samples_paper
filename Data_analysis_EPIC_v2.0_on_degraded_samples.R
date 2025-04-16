# Project: DNA degradation study
# Script Title: Data analysis of EPIC data on DNA with different combination of DNA sizes and input. 
# Description: The script was done to analyzed the EPIC data obtained from the DNA degradation study.
# The study aims to investigate the impact of different combinations of DNA amounts and fragmentations 
# on DNA methylation measurements with EPIC array.
# DNA fragments of 350, 230, 165, and 95 bp
# DNA amounts: 100 ng, 50 ng, 20 ng, and 10 ng (in duplicates)
# This script contain the analysis also of the ***** REFERENCE SAMPLE (250 ng)*****
# Author: Brando Poggiali
# Date: 22-03-2024

#### 1.Setting libraries and paths. ------------------------------------------------
BiocManager::install("")
install.packages("ggpointdensity") 

#Load package
library(sesame)
library(sesameData)
library(SummarizedExperiment)
library(readxl)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(factoextra)
library(methylclock)
library(writexl)
library(ggpubr)
library(grid)
library(gridExtra) 
library(reshape2)
library(plotly)
library(stringr)
library(ggpointdensity)
library(openxlsx)
library(scales)
library(ggridges)
library(ggplot2)

sesameDataCache() #In case you install or update SeSAMe

idat_dir <- "N:/projects/degraded_samples_epic_array_7644/local_data/"
project_dir <- "N:/projects/degraded_samples_epic_array_7644/"
brando_path <- "N:/projects/degraded_samples_epic_array_7644/users/lfw156/"
results_path <- paste0(brando_path, "Results/Results_with_reference_and_swapped_samples")
annotation_files_path <- "G:/FAELLES/Dokumenter/BRP/5_Projects/EPIC_annotations_files/"


#### 2. Upload data & metadata. ----------------------------------------------------------------------------------
betas <-  openSesame(idat_dir, prep = "QCDPB", func = getBetas) # dim: 866553    144 
metadata <- read_xlsx(path = paste0(project_dir,"doc/20240322_SampleSheet_with_reference_swapped_sample.xlsx"), sheet = 2)

#Test upload data with different ELBAR algorithm, as requested by reviewer 1.
betas <-  openSesame(idat_dir, prep = "QCDIB", func = getBetas) # dim: 866553    144 


Zhou_probe_annotation_EPIC_v2 <- read_tsv(paste0(annotation_files_path, "EPICv2.hg38.manifest.tsv.gz")) #Two different annotation files, this one contains annotation of probes with probes specification.
#Zhou_probe_annotation_EPIC_v2 <- read_tsv(paste0(annotation_files_path, "EPICv2.hg38.manifest.gencode.v41.tsv.gz"))

#### 3. Tidy up dataset. ----------------------------------------------------------------------------------
metadata$Chip_pos <- paste0(metadata$Sentrix_ID, "_", metadata$Sentrix_position)
metadata <- metadata[order(metadata$Chip_pos),]
metadata$Chip_pos == colnames(betas) #check if order of name is the same

metadata$DNA_size <- str_extract(metadata$Sample_name_after_bioanalyser_measurement, "\\d+(?=bp)")
metadata$DNA_input <- str_extract(metadata$Sample_name_after_bioanalyser_measurement, "\\d+(?=ng)")
metadata$Duplicate <- str_extract(metadata$Sample_name_after_bioanalyser_measurement, "[A-Za-z]$")
metadata$Sample_code <- sub("^[^_]*_", "", metadata$Sample_name_after_bioanalyser_measurement)
metadata$DNA_input[1] <- 250
metadata$Sample_code[1] <- "no_degraded_250ng"
metadata$DNA_input <- as.numeric(metadata$DNA_input)
metadata$DNA_size <- as.numeric(metadata$DNA_size)

colnames(betas) <- metadata$Sample_code

#Save data
saveRDS(betas, paste0(brando_path, "Data/beta_values_DNA_degradation_with_ref_swapped_sample_24-03-2024.rds"))
write.table(metadata, file.path(brando_path, "Metadata/metadata_with_ref_swapped_sample_24-03-2024.tsv"), sep="\t")

#### 4. Remove specific samples from the analysis (Optional)----------------------------------------------

#Remove samples from betas
betas_subset <- betas[,!grepl("95bp", colnames(betas))]
betas_subset <- betas_subset[, !colnames(betas_subset) %in% c("165bp_10ng_A","165bp_10ng_B")]

#Remove NAs
betas_noNAs_subset <- na.omit(betas_subset)
dim(betas_noNAs_subset) # All samples: 6,200 32; samples after subset of original samples: 204,000 23
betas_noNAs_cg_subset <- betas_noNAs_subset[grepl("cg", rownames(betas_noNAs_subset)),]
betas_noNAs_cg_subset_coll <- betasCollapseToPfx(as.matrix(betas_noNAs_cg_subset))

#With NAs
betas_cg_subset <- betas_subset[grepl("cg", rownames(betas_subset)),]
betas_cg_subset_coll <- betasCollapseToPfx(as.matrix(betas_cg_subset))

#Remove samples from metadata
metadata <- metadata[!grepl("95bp", metadata$Sample_code),]
metadata <- metadata[! metadata$Sample_code %in% c("165bp_10ng_A","165bp_10ng_B"),]
dim(metadata)
dim(betas_noNAs_cg_subset)


#Remove samples from vector with sorted names
sorted_codes <- sorted_codes[!grepl("95bp", sorted_codes)]
sorted_codes <- sorted_codes[! sorted_codes %in% c("165bp_10ng_A","165bp_10ng_B")]

betas_noNAs_cg_coll_sorted <- betas_noNAs_cg_subset_coll[, sorted_codes]

saveRDS(betas_cg_subset_coll, paste0(brando_path, "Data/beta_values_DNA_degradation_with_ref_swapped_sample_only_QC_pass_12-07-2024.rds"))

#### 5. Quality control. ----------------------------------------------------------------------------------

## 5.1 Check Sesame quality control statistics
#The quality control steps are performed after applying the pre-processing steps "QCDPB" on the sesame package
qcs <- openSesame(idat_dir, prep="QCDPB", func=sesameQC_calcStats) 
qc_df <- do.call(rbind, lapply(qcs, as.data.frame))
qc_df$Sample_code <- metadata$Sample_code

#Order codes, some of the plot could be more useful if the Sample_codes are ordered
extract_and_sort <- function(codes) { # Function to extract and convert the numeric parts for sorting
  # Extract the bp values
  bp_values <- as.numeric(sub("bp.*", "", codes))
  # Extract the ng values
  ng_values <- as.numeric(sub("^\\d+bp_", "", sub("ng.*", "", codes)))
  
  # Combine the extracted values for sorting
  sorted_indices <- order(bp_values, ng_values, decreasing = TRUE)
  
  # Return the sorted codes based on the extracted values
  return(codes[c(1,sorted_indices[-length(sorted_indices)])])
}

sorted_codes <- extract_and_sort(qc_df$Sample_code)

qc_df_sorted <- qc_df[match(sorted_codes, qc_df$Sample_code),]

qc_df_sorted$DNA_size <- str_extract(qc_df_sorted$Sample_code, "\\d+(?=bp)")
qc_df_sorted$DNA_size[1] <- "Control"
qc_df_sorted$DNA_size <- factor(qc_df_sorted$DNA_size, levels=c("Control", "350", "230", "165", "95"))

qc_df_sorted$DNA_input <- str_extract(qc_df_sorted$Sample_code, "\\d+(?=ng)")
qc_df_sorted$DNA_input[1] <- 250
qc_df_sorted$DNA_input <- factor(qc_df_sorted$DNA_input, levels=c("250", "100", "50", "20", "10"))

qc_df_sorted$Duplicate <- str_extract(qc_df_sorted$Sample_code, "[A-Za-z]$")
qc_df_sorted$Duplicate[1] <- NA


saveRDS(qc_df_sorted, file.path(brando_path, "Data/qc_df_sorted_with_ref_swapped_sample_24-03-2024.rds"))

# QC Plots
qc_df_sorted$Sample_code[1] <- "Control"
qc_df_sorted$DNA_size[1] <- "Control"
background_intensities_plot <- ggplot(qc_df_sorted, aes(x = mean_oob_grn, y= mean_oob_red, label = Sample_code)) + #Background
  geom_point(aes(color=DNA_size, shape=DNA_input), size=2.75) + #geom_text(hjust = -0.3, vjust = 0.3, position=position_jitter()) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted') +
  xlab('Mean green background signal intesity') + ylab('Mean red background signal intensity') + 
  scale_color_manual(values = c("Control" = "darkgrey", "350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1500)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1500)) + theme_bw() 


ggsave(paste0(results_path,"/1_Quality_control/Green_and_red_background_intensities_05-11-2024.png"), 
       background_intensities_plot, width = 4.8, height = 3.7)

Mean_intensities_plot <- ggplot(qc_df) + #Mean intensity
  geom_bar(aes(Sample_code, mean_intensity), stat='identity') +
  xlab('Sample Name') + ylab('Mean Intensity') +
  ylim(0,18000) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df$Sample_code) 

ggsave(paste0(results_path,"/1_Quality_control/Mean_intensities_plot_23-03-2024.png"), 
       Mean_intensities_plot, width = 6, height = 4.5)

x_position <- which(qc_df_sorted$Sample_code == "165bp_20ng_A")

Mean_intensities_plot <- ggplot(qc_df_sorted) + #Mean intensity ordered by DNA sizes and input
  geom_bar(aes(Sample_code, mean_intensity, fill=DNA_size), stat='identity') +
  xlab('Sample Name') + ylab('Mean Intensity') +
  scale_fill_manual(values = c("no_degraded" = "darkgrey", "350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000")) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df_sorted$Sample_code) +
  #geom_vline(xintercept = 23.5, linetype = "dashed", color = "red", size = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,10000))

ggsave(paste0(results_path,"/1_Quality_control/Mean_intensities_plot_ordered_23-03-2024.png"), 
       Mean_intensities_plot, width = 6, height = 4.5)

ratio_mean_intensity_red_background_plot <- ggplot(qc_df_sorted) + #Mean intensity ordered by DNA sizes and input
  geom_bar(aes(Sample_code, mean_intensity / mean_oob_red) , stat='identity') +
  xlab('Sample Name') + ylab('Mean Intensity / Red background ') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df_sorted$Sample_code) 

ggsave(paste0(results_path,"/1_Quality_control/ratio_mean_intensity_red_background_plot_23-03-2024.png"), 
       ratio_mean_intensity_red_background_plot, width = 6, height = 4.5)

ratio_mean_intensity_green_background_plot <- ggplot(qc_df_sorted) + #Mean intensity ordered by DNA sizes and input
  geom_bar(aes(Sample_code, mean_intensity / mean_oob_grn) , stat='identity') +
  xlab('Sample Name') + ylab('Mean Intensity / Red background ') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df_sorted$Sample_code) 

ggsave(paste0(results_path,"/1_Quality_control/ratio_mean_intensity_green_background_plot_23-03-2024.png"), 
       ratio_mean_intensity_green_background_plot, width = 6, height = 4.5)

Ratio_red_to_green_median_int_plot <- ggplot(qc_df_sorted) + #Ratio red to green background
  geom_bar(aes(Sample_code, medR / medG), stat='identity') +
  xlab('Sample Name') + ylab('Ratio of Red to Green median Intens.') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df_sorted$Sample_code)

ggsave(paste0(results_path,"/1_Quality_control/Ratio_red_to_green_median_int_plot_02-02-2024.png"), 
       Ratio_red_to_green_median_int_plot, width = 6, height = 4.5)

Number_detected_probes_plot <- ggplot(qc_df_sorted) + #Number detected Probes
  geom_bar(aes(Sample_code, num_dt), stat='identity') +
  xlab('Sample Name') + ylab('Number of detected probes') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df_sorted$Sample_code) +  #ylim(c(0,1000000)) +
  geom_hline(yintercept = 937690, linetype = "dashed", color = "black", size = 0.4) +
  scale_y_continuous(breaks = c(0, 250000, 500000, 750000, 937690)) +
  labs(title="Number of detected probes")

ggsave(paste0(results_path,"/1_Quality_control/Number_detected_probes_plot_23-03-2024.png"), 
       Number_detected_probes_plot, width = 6, height = 4.5)

Fraction_detected_probes_plot <- ggplot(qc_df_sorted) + #Fraction of detected Probes
  geom_bar(aes(Sample_code, frac_dt), stat='identity') +
  xlab('Sample Name') + ylab('Fraction of detected probes') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df_sorted$Sample_code) +  ylim(c(0,1)) +
  labs(title="Fraction of detected probes") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0(results_path,"/1_Quality_control/Fraction_detected_probes_plot_23-03-2024.png"), 
       Fraction_detected_probes_plot, width = 6, height = 4.5)

Number_of_NAs_plot <- ggplot(qc_df_sorted) + #Number of NAs
  geom_bar(aes(Sample_code, num_na_cg), stat='identity') +
  xlab('Sample Name') + ylab('Number of NAs') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df_sorted$Sample_code) +
  labs(title="Number of NAs") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0(results_path,"/1_Quality_control/Number_of_NAs_plot_25-03-2024.png"), 
       Number_of_NAs_plot, width = 6, height = 4.5)


bisulfite_conv_efficiency_plot <- ggplot(qc_df_sorted) + #Bisulfite conversion efficiency
  geom_bar(aes(Sample_code, frac_unmeth_ch), stat='identity') +
  xlab('Sample Name') + ylab('Fraction of unmethylated ch probes') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df_sorted$Sample_code)

ggsave(paste0(results_path,"/1_Quality_control/bisulfite_conv_efficiency_plot_24-03-2024.png"), 
       bisulfite_conv_efficiency_plot, width = 6, height = 4.5)

bisulfite_conv_efficiency_2_plot <-ggplot(qc_df_sorted) + #Bisulfite conversion efficiency 2
  geom_bar(aes(Sample_code, median_beta_ch), stat='identity') +
  xlab('Sample Name') + ylab('Median beta value in ch probes') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df_sorted$Sample_code) + ylim(0,1)

ggsave(paste0(results_path,"/1_Quality_control/bisulfite_conv_efficiency_2_plot_24-03-2024.png"), 
       bisulfite_conv_efficiency_2_plot, width = 6, height = 4.5)


objects_to_remove <- grep("plot", ls(), value = TRUE) # remove plots from environment
rm(list = objects_to_remove)


# Removal of probes with NAs
betas_noNAs <- na.omit(betas)
dim(betas_noNAs) # 6200 32
betas_noNAs_cg <- betas_noNAs[grepl("cg", rownames(betas_noNAs)),]


#Principal Component analysis, only on not missing probes which are present in all the samples
betas_noNAs_cg_coll <- betasCollapseToPfx(as.matrix(betas_noNAs_cg))
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

ggsave(paste0(results_path,"/1_Quality_control/Combined_PCA_plot_17-07-2024.png"), 
       combined_PCA_plot, width = 10, height = 4, dpi = 600, bg = "white")

## 5.6 Keeping only cg or ch sites
cg_probes <- substr(rownames(betas),1,2) == "cg"
betas_cg <- betas[cg_probes,]


## 5.7 Collapse duplicates probes
betas_cg_collapsed <- betasCollapseToPfx(as.matrix(betas_cg)) 
saveRDS(betas_cg_collapsed, paste0(brando_path, "Data/betas_cg_collapsed_DNA_degradation_with_ref_swapped_sample_24-03-2024.rds"))

## 5.8 Plot genotype probes 
rs_probes <- substr(rownames(betas),1,2) == "rs"
genotype_rs <- as.data.frame(betasCollapseToPfx(betas[rs_probes,]))


genotype_rs <- genotype_rs %>% 
  rownames_to_column(var="Probes")

genotype_rs_long <- genotype_rs %>% 
  pivot_longer(cols = -Probes,
               names_to = "Sample", 
               values_to = "Intensity")


genotype_rs_long$Sample <- factor(genotype_rs_long$Sample, level=sorted_codes)


genotype_plot <- ggplot(genotype_rs_long, aes(x = Probes, y = Sample, fill = Intensity)) +
  geom_tile(color = "black") +
  #theme_void() +
  coord_fixed() +
  scale_fill_gradientn(colours=c("#c02011", "#f2d93a","#3fb11b"), 
                       breaks = c(0,0.5,1), limits = c(0,1), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                       name="Intensity") +
  theme(
    legend.title = element_text(size=14),
    legend.text  = element_text(size=11, face="bold"),
    legend.key.height = grid::unit(1,"cm"),           # height of legend key
    legend.key.width  = grid::unit(0.6,"cm"),         # width of legend key
    legend.title.align=0.20,
    axis.text.x = element_text(size=11, angle = 45, hjust = 1, vjust = 1.05),              # axis text size
    axis.text.y = element_text(vjust=0.2, size=12),            # axis text alignment
    axis.ticks = element_line(size=0.6),
    #axis.ticks.length = unit(3, "pt"),
    #axis.ticks = element_blank(),
    axis.title = element_text(size=12, face="bold"),
    legend.box.margin=margin(10,10,10,10),
    panel.background = element_rect(color = "black", size = 1)# axis title size and bold
  ) +
  scale_x_discrete(position = "bottom", expand = c(0, 0)) +  # Attach x-axis ticks to the graph
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  labs(x = "SNPs", y = "Samples")

ggsave(paste0(results_path,"/1_Quality_control/Genotypes_25-06-2024.png"), 
       genotype_plot, width = 14, height = 8.5, dpi = 600, bg = "white")

## Quality control plots for publication

background_intensities_plot <- ggplot(qc_df_sorted, aes(x = mean_oob_grn, y= mean_oob_red, label = Sample_code)) + #Background
  geom_point(aes(color=DNA_size, shape=DNA_input), size=2.75) + #geom_text(hjust = -0.3, vjust = 0.3, position=position_jitter()) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted') +
  xlab('Mean green background signal intesity') + ylab('Mean red background signal intensity') + 
  scale_color_manual(values = c("Control" = "darkgrey", "350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1500), labels = label_comma()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1500), labels = label_comma()) + theme_bw() +
  labs(color = "DNA fragment \n size (bp)", shape = "DNA input (ng)") +
  theme(axis.title = element_text(size = 13),  # Adjust x-axis title size
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),   # Adjust legend text size
        legend.title = element_text(size = 13))

ggsave(paste0(results_path,"/1_Quality_control/Plot_red_green_background_intensities_27-06-2024.png"), 
       background_intensities_plot, width = 6.1, height = 4.25, dpi = 600)

Mean_intensities_plot <- ggplot(qc_df_sorted) + #Mean intensity ordered by DNA sizes and input
  geom_bar(aes(Sample_code, mean_intensity, fill=DNA_size), stat='identity') +
  xlab('Sample name') + ylab('Mean signal intensity') + theme_bw() +
  scale_fill_manual(values = c("Control" = "darkgrey", "350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000")) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df_sorted$Sample_code) +
  #geom_vline(xintercept = 23.5, linetype = "dashed", color = "red", size = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,10000), labels = label_comma()) +
  guides(fill=guide_legend(title="DNA fragment size (bp)")) +
  theme(axis.title = element_text(size = 14.5),  # Adjust x-axis title size
        axis.text = element_text(size = 12.5))


Ratio_red_to_green_median_int_plot <- ggplot(qc_df_sorted) + #Ratio red to green background
  geom_bar(aes(Sample_code, RGratio, fill = DNA_size), stat='identity') + theme_bw() +
  scale_fill_manual(values = c("Control" = "darkgrey", "350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000")) +
  xlab('Sample name') + ylab('Ratio red/green median signal intensities') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df_sorted$Sample_code) +
  #geom_vline(xintercept = 23.5, linetype = "dashed", color = "red", size = 1) +
  guides(fill=guide_legend(title="DNA fragment size (bp)")) +
  scale_y_continuous(expand = c(0, 0.05), labels = label_comma()) +
  theme(axis.title = element_text(size = 14.5),  # Adjust x-axis title size
        axis.text = element_text(size = 12.5))


Number_of_NAs_plot <- ggplot(qc_df_sorted) + 
  geom_bar(aes(Sample_code, num_na_cg, fill = DNA_size), stat='identity') +
  xlab('Sample Name') + 
  ylab('Number of NAs') + theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_fill_manual(values = c("Control" = "darkgrey", "350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000")) +
  scale_y_continuous(breaks = c(0, 150000, 300000, 450000, 600000), expand = expansion(mult = c(0, 0.05)), labels = label_comma()) +
  scale_x_discrete(limits = qc_df_sorted$Sample_code) +
  guides(fill=guide_legend(title="DNA fragment size (bp)")) +
  geom_vline(xintercept = 23.5, linetype = "dashed", color = "red", size = 1) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(qc_df_sorted) + 
  geom_bar(aes(Sample_code, frac_unmeth, fill = DNA_size), stat='identity') +
  xlab('Sample Name') + 
  ylab('Number of NAs') + theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_fill_manual(values = c("Not fragmented" = "darkgrey", "350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000")) +
  scale_y_continuous( expand = expansion(mult = c(0, 0.05)), labels = label_comma()) +
  scale_x_discrete(limits = qc_df_sorted$Sample_code) +
  guides(fill=guide_legend(title="DNA fragment size (bp)")) +
  theme(plot.title = element_text(hjust = 0.5))


combined_plot_QC <- ggarrange(Mean_intensities_plot,
                              Ratio_red_to_green_median_int_plot, #Number_of_NAs_plot,
                              ncol = 1, nrow = 2,
                              legend ="bottom",
                              common.legend = TRUE)

ggsave(paste0(results_path,"/1_Quality_control/Quality_control_for_publication_14-04-2025.png"), 
       combined_plot_QC, width = 12.7, height = 8, dpi = 600, bg = "white")

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
ggsave(paste0(results_path,"/1_Quality_control/Missing_data_per_bp_and_ng_20-12-2024.png"), 
       combined_plot_nas, width = 8.5, height = 3.6, dpi = 600, bg = "white")

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

write.xlsx(data, file = paste0(results_path,"/1_Quality_control/Missing_probes_20-12-2024.xlsx"),
           sheetName = "Missing data", rowNames = FALSE)


#### 6. Upload of RSD data --------------------------------------------------------------------
betas <- readRDS(paste0(brando_path,"/Data/beta_values_DNA_degradation_with_ref_swapped_sample_24-03-2024.rds")) # no normalize for batch effect
metadata <- read.table(file = paste0(brando_path,'/Metadata/metadata_with_ref_swapped_sample_24-03-2024.tsv'), sep = '\t', header = TRUE)

betas_cg_collapsed <- readRDS(paste0(brando_path, "/Data/betas_cg_collapsed_DNA_degradation_with_ref_swapped_sample_24-03-2024.rds"))

qc_df_sorted <- readRDS(paste0(brando_path, "/Data/qc_df_sorted_with_ref_swapped_sample_24-03-2024.rds"))
sorted_codes <- qc_df_sorted$Sample_code
sorted_codes[1] <- "Control"
#sorted_codes <- metadata$Sample_code

sorted_codes_no_rep <- c(sorted_codes[1], gsub("_A|_B", "", sorted_codes[-1]))
sorted_codes_no_rep <- unique(sorted_codes_no_rep)
#### 7. Advance Quality control plots. ----------------------------------------------------------------------------------

## 7.1 Check Beta value distribution (This step require high amount of memory).

# betas_type_I <- betas[Zhou_probe_annotation_EPIC_v2$type != "I",]
# dim(betas_type_I)

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
# colors <- c("darkgrey", "black", "#f2d93a", "#39a2a4", "#f20000")
# colors <- c( "darkgrey", "black","black", "#f2d93a","#f2d93a", 
#              "#39a2a4","#39a2a4", "#f20000", "#f20000")
# names(colors) <- unique_samples

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


ggsave(paste0(results_path,"/1_Quality_control/beta_value_distribution/Samples_beta_distribution_350bp_03-04-2024.png"), 
       samples_distribution_350bp, width = 8, height = 6)

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

ggsave(paste0(results_path,"/1_Quality_control/beta_value_distribution/Samples_beta_distribution_230bp_03-04-2024.png"), 
       samples_distribution_230bp, width = 8, height = 6, dpi = 600, bg = "white")

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

ggsave(paste0(results_path,"/1_Quality_control/beta_value_distribution/Samples_beta_distribution_165bp_03-04-2024.png"), 
       samples_distribution_165bp, width = 8, height = 6)


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

ggsave(paste0(results_path,"/1_Quality_control/beta_value_distribution/Samples_beta_distribution_95bp_03-04-2024.png"), 
       samples_distribution_95bp, width = 8, height = 6)

spacer <- ggplot() + theme_void()
combined_plot <- ggarrange(samples_distribution_350bp, spacer, samples_distribution_230bp, spacer,
                           samples_distribution_165bp, spacer, samples_distribution_95bp, spacer,
                           labels = c("A", "", "B","", "C", "", "D", ""),
                           ncol = 4, nrow = 2,
                           common.legend = TRUE,
                           legend = "bottom",
                           widths = c(1, 0.01, 1, 0.01))

ggsave(paste0(results_path,"/1_Quality_control/beta_value_distribution/All_sample_distributions_04-11-2024.png"), 
       combined_plot, width = 10, height = 6, dpi = 600, bg = "white")



## 7.3 Probe success rate plot
qc_df_subset <- qc_df_sorted[qc_df_sorted$DNA_size != 95,]
qc_df_subset <- qc_df_subset[!c(qc_df_subset$DNA_size == 165 & qc_df_subset$DNA_input == 10),]
qc_df_subset <- qc_df_subset[-2,]

qc_df_sorted_do_degrad_samples <- qc_df_subset[-1,] #remove no degraded DNA

qc_df_sorted_do_degrad_samples$Sample_combination <- substr(qc_df_sorted_do_degrad_samples$Sample_code, 1, nchar(qc_df_sorted_do_degrad_samples$Sample_code) - 2)

qc_df_sorted_do_degrad_samples_average <- qc_df_sorted_do_degrad_samples %>%
  group_by(Sample_combination) %>%
  summarize(across(.cols = everything(), .fns = mean, na.rm = TRUE), .groups = "drop")

qc_df_sorted_do_degrad_samples_average$DNA_size <- str_extract(qc_df_sorted_do_degrad_samples_average$Sample_combination, "\\d+(?=bp)")
qc_df_sorted_do_degrad_samples_average$DNA_input <- str_extract(qc_df_sorted_do_degrad_samples_average$Sample_combination, "\\d+(?=ng)")

qc_df_sorted_do_degrad_samples_average$DNA_size <- factor(qc_df_sorted_do_degrad_samples_average$DNA_size, level=c("350", "230", "165"))
qc_df_sorted_do_degrad_samples_average$DNA_input <- factor(qc_df_sorted_do_degrad_samples_average$DNA_input, level=c("100", "50", "20", "10"))


plot_table <- function(dataframe, x_variable, y_variable, frac_succ_probes, number_of_probes,
                       x_axis_title, y_axis_title){
  graph <- ggplot(dataframe, aes(x = x_variable, y = y_variable, fill = frac_succ_probes)) +
    geom_tile(color = "black") +
    #theme_void() +
    geom_text(aes(label = paste0(round(frac_succ_probes* 100, 2), "%")), color = "black", size = 4.6, nudge_y = 0.10) +
    geom_text(aes(label = format(number_of_probes, big.mark = ",")), size = 4.6, nudge_y = -0.15) +
    coord_fixed() +
    scale_fill_gradientn(colours=c("#FFFFFF", "#d0b5e3","#792cae"), 
                         breaks = c(0,0.5,1), limits = c(0,1), 
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                         name="Probe detection rate\n   No. of probes passing filtering",
                         labels = c("0.00\n0", "50.00\n468,845", "100.00\n937,690")) +
    theme(
      legend.title = element_text(size=13, margin = margin(b = 10)),
      legend.text  = element_text(size=11, face="bold"),
      legend.key.height = grid::unit(1,"cm"),           # height of legend key
      legend.key.width  = grid::unit(0.6,"cm"),         # width of legend key
      legend.spacing = unit(10, "cm"),
      legend.title.align=0.20,
      axis.text.x = element_text(size=12),              # axis text size
      axis.text.y = element_text(vjust=0.2, size=12),            # axis text alignment
      axis.ticks = element_line(size=0.8),
      #axis.ticks.length = unit(3, "pt"),
      #axis.ticks = element_blank(),
      axis.title = element_text(size=13, face="bold"),
      legend.box.margin=margin(10,10,10,10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(color = "black", size = 0.2)# axis title size and bold
    ) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +  # Attach x-axis ticks to the graph
    scale_y_discrete(limits=rev, expand = c(0, 0)) +
    labs(x = x_axis_title, y = y_axis_title)
  
  print(graph)
}

probe_success_rate_plot <- plot_table(qc_df_sorted_do_degrad_samples_average, x_variable=qc_df_sorted_do_degrad_samples_average$DNA_size, 
           y_variable=qc_df_sorted_do_degrad_samples_average$DNA_input, frac_succ_probes=qc_df_sorted_do_degrad_samples_average$frac_dt,
           number_of_probes=qc_df_sorted_do_degrad_samples_average$num_dt,
           x_axis_title = "DNA fragment size (bp)", y_axis_title = "DNA input (ng)")


ggsave(paste0(results_path,"/1_Quality_control/Plot_success_rate_and_nprobes_14-04-2025.png"), 
       probe_success_rate_plot, width = 6.5, height = 4.3, dpi = 600, bg = "white")



#### 8. Pearson correlation. ----------------------------------------------------------------------------------
#install.packages("ggcorrplot")
library(corrplot)
sorted_codes <- qc_df_sorted$Sample_code

betas_noNAs <- na.omit(betas)
dim(betas_noNAs) # 6200 32
betas_noNAs_cg <- betas_noNAs[grepl("cg", rownames(betas_noNAs)),]

betas_noNAs_cg_sorted <- betas_noNAs_cg[, sorted_codes] # 6200 33
betas_noNAs_cg_sorted <- betas_noNAs_cg_subset[, sorted_codes] # If you want the dataset with 204K CpGs
### With duplicates

# To decide which statistics use


### With average of duplicates 
#Prepare a new dataframe to store the averages
average_data <- data.frame(matrix(ncol = length(colnames(betas_noNAs_cg_sorted))/2 +0.5, nrow = nrow(betas_noNAs_cg_sorted)))

# Setting column names for the new dataframe
names(average_data) <- gsub("_A|_B", "", colnames(betas_noNAs_cg_sorted)[seq(1, ncol(betas_noNAs_cg_sorted), 2)])

# Loop to calculate averages
i=1
average_data[,1] <- betas_noNAs_cg_sorted[,1]
for (i in seq(2, ncol(betas_noNAs_cg_sorted), 2)) {
  col_name <- gsub("_A|_B", "", colnames(betas_noNAs_cg_sorted)[i])  # General column name without _A or _B
  average_data[[col_name]] <- rowMeans(betas_noNAs_cg_sorted[, c(i, i+1)])
}


average_data_02 <- average_data[average_data[,1] < 0.2,]
average_data_02_08 <- average_data[average_data[,1] > 0.2 & average_data[,1] < 0.8,]
average_data_08 <- average_data[average_data[,1] > 0.8,]

cor_matrix <- cor(average_data_02_08) # substitute with average_data, average_data_02, average_data_02_08, or average_data_08 
melted_cor_matrix <- melt(cor_matrix)


### Calculate confidence p-value for statistic course (this analysis is not for the paper):
num_cols <- ncol(average_data_02_08)
p_value_matrix <- matrix(NA, ncol = num_cols, nrow = num_cols, dimnames = list(colnames(average_data) , colnames(average_data)))



# Loop over each pair of columns
for (i in 1:num_cols) {
  for (j in i:num_cols) {  # Start from 'i' to avoid redundant calculations and self-correlation
    test_result <- cor.test(average_data[, i], average_data[, j], method = "pearson")
    p_value_matrix[i, j] <- test_result$p.value
    p_value_matrix[j, i] <- test_result$p.value  # Mirror the matrix (since correlation matrix is symmetric)
  }
}

adjusted_p_matrix <- matrix(p.adjust(as.vector(p_value_matrix), method = "BH"),
                            nrow = nrow(p_value_matrix),
                            ncol = ncol(p_value_matrix))

colnames(adjusted_p_matrix) <- colnames(p_value_matrix)
rownames(adjusted_p_matrix) <- rownames(p_value_matrix)


# Print the p-value matrix
round(p_value_matrix,15)
p_value_matrix <- as.data.frame(p_value_matrix)
adjusted_p_matrix <- as.data.frame(adjusted_p_matrix)
write_xlsx(adjusted_p_matrix, paste0(results_path, "/2_Correlation/p_value_corrected_matrix_04_11_2024.xlsx"))



# Pearson correlation
pearson_correlation_matrix <- ggplot(data = melted_cor_matrix, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(colour = "white") + # This creates the tiles for heatmap
  geom_text(aes(label=round(value, 3)), color="black") + # This adds the text labels
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.875, limit = c(0.75,1), space = "Lab", # For the all samples use: midpoint = 0, limit = c(-1,1)
                       name="Pearson\nCorrelation") + # This defines the color scale
  theme_minimal() + # Minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank()) + # Adjust the text angle for X axis
  coord_fixed() # Ensures the tiles are square

ggsave(paste0(results_path,"/2_Correlation/Pearson_correlation_subset_samples_16-04-2024.png"), 
       pearson_correlation_matrix, width = 11, height = 9, dpi = 600, bg = "white")



# Squared Pearson correlation (it was substitute with Pearson correlation after Athina)
#squared_cor_matrix <- cor_matrix^2
squared_cor_matrix <- cor_matrix
colnames(squared_cor_matrix)[1] <- "Control"
rownames(squared_cor_matrix)[1] <- "Control"
melted_squared_cor_matrix <- melt(squared_cor_matrix)


squared_pearson_correlation_matrix <- ggplot(data = melted_squared_cor_matrix, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(colour = "#2b2b2b") + # This creates the tiles for heatmap
  geom_text(aes(label=round(value, 3)), color="black") + # This adds the text labels
  scale_fill_gradient2(low = "white", high = "red", 
                       midpoint = 0.925, limit = c(0.850,1.000), space = "Lab", # All samples midpoint=0.25, subset 0.875
                       name=bquote('Pearson correlation coefficient (r)'),
                       guide = guide_colourbar(frame.colour = "black", frame.linewidth = 0.5,
                       ticks.colour = "black", ticks.linewidth = 0.5)
  ) + # This defines the color scale
  theme_minimal() + # Minimal theme
  theme(axis.text.x = element_text(angle = 0, vjust = 1, 
                                   size = 12, hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        legend.title = element_text(size = 12.5)) + # Adjust the text angle for X axis
  coord_fixed() + 
  scale_x_discrete(labels = function(x) ifelse(nchar(x) > 7, substr(x, 7, nchar(x)), x)) +
  scale_y_discrete(labels = function(x) ifelse(nchar(x) > 7, substr(x, 7, nchar(x)), x))

squared_pearson_correlation_matrix

ggsave(paste0(results_path,"/2_Correlation/Pearson_correlation_subset_samples_15-04-2025.png"), 
       squared_pearson_correlation_matrix, width = 11, height = 9, dpi = 600, bg = "white")

ggsave(paste0(results_path,"/2_Correlation/Squared_Pearson_correlation_all_samples_16-04-2024.png"), 
       squared_pearson_correlation_matrix, width = 11, height = 9, dpi = 600, bg = "white")


#Plot pearson correlations based on DNA fragment size and input amount.
df_for_plot_correlation <- melted_squared_cor_matrix[melted_squared_cor_matrix$Var1 == "Control",]

df_for_plot_correlation$DNA_size <- substr(df_for_plot_correlation$Var2,1,3)
df_for_plot_correlation$DNA_input <- str_extract(df_for_plot_correlation$Var2, "\\d+(?=ng)")
df_for_plot_correlation <- df_for_plot_correlation[-1,]
df_for_plot_correlation$DNA_input <- factor(df_for_plot_correlation$DNA_input, levels = c(10, 20, 50, 100))
df_for_plot_correlation$DNA_size <- factor(df_for_plot_correlation$DNA_size, levels = c(350, 230, 165))


plot_correlation_02_08 <- ggplot(df_for_plot_correlation, aes(x = DNA_input, y = value, color = DNA_size, group = DNA_size)) +
  geom_line(size = 1) +  # Connect points with lines
  geom_point(size = 2) +  # Optional: Add points at each time point
  theme_bw()  +
  scale_color_manual(values = c("350" = "black", "230" = "#f2d93a", "165" = "#39a2a4")) +
  labs(color= "DNA fragment size (bp)",
       x = "DNA input (ng)",
      y = bquote('Pearson correlation coefficient (r)'),
      color = "DNA size") +
  scale_y_continuous(limits = c(0.2, 1), expand = c(0, 0))

plot_correlation_02 <- plot_correlation_02 + ggtitle("CpG Methylation 0–0.2")
plot_correlation_02_08 <- plot_correlation_02_08 + ggtitle("CpG Methylation 0.2–0.8")
plot_correlation_08 <- plot_correlation_08 + ggtitle("CpG Methylation 0.8–1")

combined_plot_pearson_corr <- ggarrange(plot_correlation_02,
                                        plot_correlation_02_08,
                                        plot_correlation_08,
                                        labels = c("A", "B", "C"),
                                        common.legend = TRUE,
                                        legend = "bottom",
                                        ncol = 3, nrow = 1)

ggsave(paste0(results_path,"/2_Correlation/Plot_correlations_15-05-2025.png"), 
       combined_plot_pearson_corr, width = 12, height = 4.8, dpi = 600, bg = "white")

#Plot of pearson correlation to visualize the values
scatter_plot_08 <- ggplot(average_data_08, aes(x = no_degraded_250ng, y = `165bp_20ng`)) +
  geom_point(color = "#3D69C2", size = 1, alpha = 0.3) +  # 👈 transparency here
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  stat_cor(method = "pearson", label.x = 0.05, label.y = 0.95, size = 5) +
  theme_minimal() +
  labs(x = "Control",
       y = "165ng 20 bp",
       title = "")




install.packages("ggpointdensity")
library(ggpointdensity)

scatter_plot_02_08 <- ggplot(average_data_02_08, aes(x = no_degraded_250ng, y = `165bp_20ng`)) +
  geom_pointdensity(adjust = 0.5, size = 1) +
  scale_color_viridis_c(option = "C") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  stat_cor(method = "pearson", label.x = 0.05, label.y = 0.95, size = 5) +
  theme_minimal() +
  labs(x = "Control",
       y = "165ng 20 bp",
       color = "Density",
       title = "")

scatter_plot_pearson_corr <- ggarrange(scatter_plot_02,
                                        scatter_plot_02_08,
                                        scatter_plot_08,
                                        labels = c("A", "B", "C"),
                                        common.legend = TRUE,
                                        legend = "bottom",
                                        ncol = 3, nrow = 1)

ggsave(paste0(results_path,"/2_Correlation/Plot_correlations_test_15-05-2025.png"), 
       scatter_plot_pearson_corr, width = 12, height = 4.8, dpi = 600, bg = "white")


#### 9. Pearson correlation between replicates --------------------------------------

sorted_codes <- qc_df_sorted$Sample_code
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
write_xlsx(square_pearson_corr_replicates, paste0(results_path, "/2_Correlation/square_pearson_corr_replicates.xlsx"))


#plot points
df_tmp <- as.data.frame(df_tmp)
ggplot(df_tmp, aes(x = df_tmp[[1]], y = df_tmp[[2]])) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_gradient(low = "white", high = "darkblue") +  # Custom colors
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position='none')

#### 10. Delta beta value. ----------------------------------------------------------------------------------
#betas_noNAs_cg_sorted <- betas_noNAs_cg[, sorted_codes]
betas_noNAs <- na.omit(betas)
dim(betas_noNAs) # 6200 32
betas_noNAs_cg <- betas_noNAs[grepl("cg", rownames(betas_noNAs)),]
betas_noNAs_cg_coll <- betasCollapseToPfx(as.matrix(betas_noNAs_cg))
betas_noNAs_cg_coll_sorted <- betas_noNAs_cg_coll[, sorted_codes] # 6200 33

#In case you use the subset dataframe, which does not contain all samples.
sorted_codes[1] <- "no_degraded_250ng"
betas_noNAs_cg_subset_coll <- betasCollapseToPfx(as.matrix(betas_noNAs_cg_subset))
betas_noNAs_cg_coll_sorted <- betas_noNAs_cg_subset_coll[, sorted_codes] # 200K 23

#betas_noNAs_cg_coll_sorted_t <- betas_noNAs_cg_coll_200K_sorted_t
betas_noNAs_cg_coll_sorted_t <- as.data.frame(t(betas_noNAs_cg_coll_sorted)) #Substitute here the 6200 or 200K datasets
betas_noNAs_cg_coll_sorted_t$Sample_combination <- substr(rownames(betas_noNAs_cg_coll_sorted_t), 1, nchar(rownames(betas_noNAs_cg_coll_sorted_t)) - 2) 

#Calculate average between replicates
betas_noNAs_cg_coll_sorted_t <- betas_noNAs_cg_coll_sorted_t %>% #Slow step for 200K
  group_by(Sample_combination) %>%
  summarize(across(.cols = everything(), .fns = mean, na.rm = TRUE), .groups = "drop")

#betas_noNAs_cg_coll_200K_sorted_t <- betas_noNAs_cg_coll_sorted_t

no_degraded_sample <- betas_noNAs_cg_coll_sorted_t[betas_noNAs_cg_coll_sorted_t$Sample_combination == "no_degraded_250",] 
delta_beta <- sweep(betas_noNAs_cg_coll_sorted_t[,-1], 2, STATS=unlist(no_degraded_sample[,-1]), FUN="-")
#delta_beta <- abs(delta_beta)
delta_beta <- cbind(betas_noNAs_cg_coll_sorted_t[,1],delta_beta)
delta_beta <- delta_beta[delta_beta$Sample_combination != "no_degraded_250",]

# Plot delta beta value distribution
delta_beta_long <- delta_beta %>% pivot_longer(!Sample_combination, 
                                               names_to = "Probes",
                                               values_to = "Delta_beta"
)
sorted_codes_no_rep <- gsub("_A|_B", "", sorted_codes[-1])
sorted_codes_no_rep <- unique(sorted_codes_no_rep)
delta_beta_long$Sample_combination <- factor(delta_beta_long$Sample_combination, level=rev(sorted_codes_no_rep))

delta_beta_long$DNA_size <- str_extract(delta_beta_long$Sample_combination, "\\d+(?=bp)")
delta_beta_long$DNA_size <- factor(delta_beta_long$DNA_size, level=c("350", "230", "165", "95"))
delta_beta_long$Sample_combination <- gsub("_", " ", delta_beta_long$Sample_combination)
delta_beta_long$Sample_combination <- factor(delta_beta_long$Sample_combination, level=rev(gsub("_", " ", sorted_codes_no_rep)))


#Plot absolute Beta value difference
abs_beta_value_distribution_plot <- ggplot(delta_beta_long, aes(x = abs(Delta_beta), y = Sample_combination, fill=DNA_size)) +
  geom_density_ridges(alpha=0.5, scale = 2.5, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 2) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.05, 0.10, 0.15, 0.2, 0.25, 0.3, 0.35),  
                     labels = c("0", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35"),limits = c(0, 0.35)) +
  theme(panel.grid.minor.x = element_blank()) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_fill_manual(values = c("350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000"), name="DNA fragment\nsize (bp)") +
  xlab(expression(paste(italic("|Δβ|")))) + ylab("Samples") +
  theme(
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 14),   # Increase x-axis text size
    axis.text.y = element_text(size = 14),
    legend.key.size = unit(2, "lines"),
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),
    axis.ticks.length.y = unit(0.21, "cm")
  ) +
  coord_cartesian(ylim = c(NA, max(as.numeric(factor(delta_beta_long$Sample_combination))) + 2.6)) +
  scale_y_discrete(labels = function(x) ifelse(nchar(x) > 7, substr(x, 7, nchar(x)), x))

abs_beta_value_distribution_plot  

ggsave(paste0(results_path,"/3_Delta_beta_value/abs_beta_value_distribution_all_samples_plot_15-04-2025.png"), 
       abs_beta_value_distribution_plot, width = 12, height = 10, dpi = 600, bg = "white")

#Plot Beta value difference
# beta_value_distribution_plot <- ggplot(delta_beta_long, aes(x = Delta_beta, y = Sample_combination, fill=DNA_size)) +
#   geom_density_ridges(alpha=0.5, scale = 3, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 2) +
#   theme_minimal() +
#   scale_x_continuous(breaks = c(-0.50, -0.25, -0.10, -0.05, 0, 0.05, 0.10, 0.25, 0.50),  
#                      labels = c("-0.50", "-0.25", "-0.10", "-0.05", "0", "0.05", "0.10", "0.25", "0.50"),limits = c(-0.5, 0.5)) +
#   scale_fill_manual(values = c("350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000"), name="DNA sizes") +
#   theme(panel.grid.minor.x = element_blank()) +
#   xlab("Δβ-value") + ylab("Samples")


ggsave(paste0(results_path,"/3_Delta_beta_value/beta_value_distribution_all_samples_plot_27-06-2024.png"), 
       beta_value_distribution_plot, width = 10, height = 7, dpi = 600, bg = "white")


#Calculate statistics for Beta value distribution and save in a file
summary_statistics <- delta_beta_long %>%
  group_by(Sample_combination) %>%
  summarize(
    Abs_Median = median(abs(Delta_beta)),
    Abs_IQR = IQR(abs(Delta_beta)),
    Median = median(Delta_beta),
    IQR = IQR(Delta_beta),
    beta_05 = (sum(abs(Delta_beta) >= 0.05))/(length(Sample_combination)),
    beta_10 = (sum(abs(Delta_beta) >= 0.1))/(length(Sample_combination)))

summary_statistics$beta_05 <- summary_statistics$beta_05 * 100
summary_statistics$beta_10 <- summary_statistics$beta_10 * 100

summary_statistics <- summary_statistics[rev(seq_len(nrow(summary_statistics))),]

colnames(summary_statistics) <- c("Sample","Median |Δβ|","IQR |Δβ|","Median Δβ","IQR Δβ","|Δβ| ≥ 0.05 (%)","|Δβ| ≥ 0.10 (%)")

write_xlsx(summary_statistics, paste0(results_path,"/3_Delta_beta_value/beta_value_distribution_statistics_all_samples_28-06-2024.xlsx"))
#### 11. Delta beta value in DNAm intervals -------------------------------------------
sorted_codes_no_rep <- gsub("_A|_B", "", sorted_codes)
sorted_codes_no_rep <- unique(sorted_codes_no_rep)
beta_value_diff_and_intervals <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(beta_value_diff_and_intervals) <- c("Delta_beta", "Range", "Conditions")

###

#Prepare a new dataframe to store the averages
average_data <- data.frame(matrix(ncol = length(colnames(betas_noNAs_cg_coll_sorted))/2 +0.5, nrow = nrow(betas_noNAs_cg_coll_sorted)))

# Setting column names for the new dataframe
names(average_data) <- gsub("_A|_B", "", colnames(betas_noNAs_cg_coll_sorted)[seq(1, ncol(betas_noNAs_cg_coll_sorted), 2)])

# Loop to calculate averages
i=1
average_data[,1] <- betas_noNAs_cg_coll_sorted[,1]
for (i in seq(2, ncol(betas_noNAs_cg_coll_sorted), 2)) {
  col_name <- gsub("_A|_B", "", colnames(betas_noNAs_cg_coll_sorted)[i])  # General column name without _A or _B
  average_data[[col_name]] <- rowMeans(betas_noNAs_cg_coll_sorted[, c(i, i+1)])
}

betas_noNAs_cg_coll_sorted_average <- average_data

#beta_values_SCD_wide <- beta_values_SCD_wide[1:80000,]

num_row <- nrow(betas_noNAs_cg_coll_sorted_average)
rownames(betas_noNAs_cg_coll_sorted_average) <- rownames(betas_noNAs_cg_coll_sorted)
sample <- sorted_codes_no_rep[2]
for (sample in sorted_codes_no_rep[-1]){ # Sometimes it could work with this sorted_codes_no_rep[-1], depends on the vector
   
  empty_df <- data.frame(matrix(nrow = num_row, ncol = 3))
  colnames(empty_df) <- c("Delta_beta", "Range", "Conditions")
  rownames(empty_df) <- rownames(betas_noNAs_cg_coll_sorted_average)
  empty_df["Delta_beta"] <- abs(betas_noNAs_cg_coll_sorted_average[,sorted_codes_no_rep[1]] - betas_noNAs_cg_coll_sorted_average[,sample])
  value <- empty_df["Delta_beta"]       
  #identify the correct beta value range      
  
  empty_df$Range <- cut(betas_noNAs_cg_coll_sorted_average[,sorted_codes_no_rep[1]], breaks =  seq(0, 1, by = 0.1), labels = FALSE)
  empty_df$Conditions <- sample
  beta_value_diff_and_intervals <- rbind(beta_value_diff_and_intervals, empty_df) 
}

#Plot aesthetics and boxplot plot of ranges

beta_value_diff_and_intervals$Range <- as.factor(beta_value_diff_and_intervals$Range)
beta_value_diff_and_intervals$Conditions <- factor(beta_value_diff_and_intervals$Conditions, levels=sorted_codes_no_rep)
beta_value_diff_and_intervals$DNA_size <- str_extract(beta_value_diff_and_intervals$Conditions, "\\d+(?=bp)")
beta_value_diff_and_intervals$DNA_size <- factor(beta_value_diff_and_intervals$DNA_size, level=c("350", "230", "165", "95"))
beta_value_diff_and_intervals$DNA_input <- str_extract(beta_value_diff_and_intervals$Conditions, "\\d+(?=ng)")
beta_value_diff_and_intervals$DNA_input <- factor(beta_value_diff_and_intervals$DNA_input, level=c("100", "50", "20", "10"))

range_labels <- c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6",
                  "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
condition_colors <- c("350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000")


boxplot_beta_value_diff <- ggplot(beta_value_diff_and_intervals,aes(x = Range, y = Delta_beta, fill = DNA_size)) +
  geom_boxplot(aes(color=DNA_size), size=0.5, outlier.shape = NA) + coord_cartesian(ylim = c(0, 0.4)) +
  facet_wrap(~ DNA_input, ncol=1, labeller=labeller(DNA_input = c("100" = "100 ng",
                                                                  "50" = "50 ng",
                                                                  "20" = "20 ng",
                                                                  "10" = "10 ng"))) +
  scale_x_discrete(labels = range_labels) +
  labs(y= expression(italic("|Δβ|")), x = expression(paste(italic("β"), "-value intervals"))) +
  scale_fill_manual(values = condition_colors) + theme_bw() + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=17), 
        legend.title=element_text(size=16), legend.text = element_text(size=15)) +
  scale_color_manual(values=c( "#808080","#808080","#808080")) +
  labs(fill="DNA fragment\nsize (bp)",color="DNA fragment\nsize (bp)") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.background = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.9, "lines"),
        text = element_text(size = 18),
        legend.key.size = unit(2, "lines"),
        legend.title = element_text(size = 16))



ggsave(paste0(results_path,"/3_Delta_beta_value/Boxplot_beta_value_diff_beta_intervals_07-11-2024.png"), 
       boxplot_beta_value_diff, width = 13, height = 11, dpi = 600, bg = "white")

#### 12. Delta beta value between replicates -----------------------------------------
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

write_xlsx(delta_beta_values, paste0(results_path, "/3_Delta_beta_value/Delta_beta_statistics_replicates.xlsx"))

#### 11. Linear models--------------------------------------


qc_df_sorted_do_degrad_samples <- qc_df_sorted[-1,] #remove no degraded DNA

qc_df_sorted_do_degrad_samples$Sample_combination <- substr(qc_df_sorted_do_degrad_samples$Sample_code, 1, nchar(qc_df_sorted_do_degrad_samples$Sample_code) - 2)

qc_df_sorted_do_degrad_samples_average <- qc_df_sorted_do_degrad_samples %>%
  group_by(Sample_combination) %>%
  summarize(across(.cols = everything(), .fns = mean, na.rm = TRUE), .groups = "drop")

qc_df_sorted_do_degrad_samples_average$DNA_size <- str_extract(qc_df_sorted_do_degrad_samples_average$Sample_combination, "\\d+(?=bp)")
qc_df_sorted_do_degrad_samples_average$DNA_input <- str_extract(qc_df_sorted_do_degrad_samples_average$Sample_combination, "\\d+(?=ng)")

#qc_df_sorted_do_degrad_samples_average$DNA_size <- factor(qc_df_sorted_do_degrad_samples_average$DNA_size, level=c("350", "230", "165", "95"))
#qc_df_sorted_do_degrad_samples_average$DNA_input <- factor(qc_df_sorted_do_degrad_samples_average$DNA_input, level=c("100", "50", "20", "10"))

qc_subset <- qc_df_sorted_do_degrad_samples_average[qc_df_sorted_do_degrad_samples_average$DNA_size != 95,]
qc_subset <- qc_subset[-2,]


# plot observations
plot(num_dtna ~ DNA_size,
     data = qc_df_sorted_do_degrad_samples_average,
     col= ifelse(irlpolwomen$country == 4,"blue","red"),
     pch=19,
     xlab = "BMI",
     ylab = expression(paste(log[10],"(vitamin D)"))
)

# We create two dataset containing the observations of covariates for which we want the
# corresponding estimated mean outcome (one for each country)
ireland <- expand.grid(bmi=irlpolwomen$bmi[which(irlpolwomen$Country == "Ireland")],
                       Country=c("Ireland"))
ireland$bmi5 <- ireland$bmi/5
poland <- expand.grid(bmi=irlpolwomen$bmi[which(irlpolwomen$Country == "Poland")],
                      Country=c("Poland"))
poland$bmi5 <- poland$bmi/5
# add the estimated means
lines(ireland$bmi, predict(lm3, ireland),col="blue",lwd=2)
lines(poland$bmi, predict(lm3, poland),col="red",lwd=2)
legend("topright",pch=19,col=c("blue","red"),legend=c("Ireland","Poland"))



ggplot(data=qc_subset, aes(x=DNA_size, y=num_dt, color=DNA_input))+
  geom_point()



#Model all samples
qc_df_sorted_do_degrad_samples_average$DNA_size <- as.numeric(qc_df_sorted_do_degrad_samples_average$DNA_size)
qc_df_sorted_do_degrad_samples_average$DNA_input <- as.numeric(qc_df_sorted_do_degrad_samples_average$DNA_input)
mult_model <- lm(num_dt ~ DNA_size*DNA_input,data=qc_df_sorted_do_degrad_samples_average)
summary(mult_model)

#Model with only good samples
qc_subset$DNA_size <- as.numeric(qc_subset$DNA_size)
qc_subset$DNA_input <- as.numeric(qc_subset$DNA_input)
mult_model_subset <- lm(num_dt ~ DNA_size+DNA_input,data=qc_subset)
summary(mult_model_subset)

str(qc_subset)

data <- data.frame(
  DNA_size =c(350) ,
  DNA_input = c(250)
)

predict(mult_model_subset, data)



#### 12. Proportion of CpGs in different genomic regions --------------------------
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
annEPIC_v2 <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
head(annEPIC_v2)
annEPIC_v2$Name_cpg <- ifelse(grepl("cg",annEPIC_v2$Name), substring(annEPIC_v2$Name,1,nchar(annEPIC_v2$Name)-5),NA)

df_island_proportion <- data.frame(matrix(nrow=ncol(betas_cg_subset_coll), ncol=4))
colnames(df_island_proportion) <- c("Island", "OpenSea", "Shelf", "Shore")
#sample <- colnames(betas_cg_collapsed)[4]

colnames(betas_cg_subset_coll)[1] <- "Control"

n <- 1
for (sample in colnames(betas_cg_subset_coll)){
  betas_cg_collapsed_subset <- betas_cg_subset_coll[,sample]
  cpg_names_subset <- names(na.omit(betas_cg_collapsed_subset))
  annEPIC_v2_subset <- annEPIC_v2[annEPIC_v2$Name_cpg %in% cpg_names_subset,]  
  annEPIC_v2_subset_no_dup <- annEPIC_v2_subset[!duplicated(annEPIC_v2_subset$Name_cpg),]
  cpg_island_location <- table(annEPIC_v2_subset_no_dup$Relation_to_Island)
  if(sum(names(cpg_island_location) != colnames(df_island_proportion)) != 0){break}
  rownames(df_island_proportion)[n] <- sample
  df_island_proportion[n,] <- cpg_island_location
  n <- n + 1
  } 


df_island_proportion <- rownames_to_column(df_island_proportion, var = "Samples")
df_island_proportion_long <- pivot_longer(df_island_proportion, 
                        cols = -Samples, 
                        names_to = "Island_location", 
                        values_to = "Number_cpg")

df_island_proportion_long$Samples <- factor(df_island_proportion_long$Samples, levels = sorted_codes) 
df_island_proportion_long$Island_location[df_island_proportion_long$Island_location == "OpenSea"] <- "Open Sea" 


#Calculate perrcentages
df_island_proportion_long <- df_island_proportion_long %>%
  group_by(Samples) %>%
  mutate(Proportion = Number_cpg / sum(Number_cpg),
         PercentLabel = paste0(round(Proportion * 100), "%"))

#Plot
# island_proportion_plot <- ggplot(df_island_proportion_long, aes(fill=Island_location, y = Number_cpg, x = Samples)) + 
#   geom_bar(position="fill", stat="identity") + 
#   theme_bw() + xlab("Samples") + ylab("Relative proportion of CpG probes") + labs(fill = "Genomic location") +  
#   scale_fill_manual(values=c("#C2963D","#C23DAC","#3D69C2", "#3DC253")) +
#   theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0)) +
#   scale_y_continuous(expand = c(0, 0)) 


island_proportion_plot <- ggplot(df_island_proportion_long, aes(fill = Island_location, y = Number_cpg, x = Samples)) + 
  geom_bar(position = "fill", stat = "identity") + 
  geom_text(aes(label = PercentLabel), 
            position = position_fill(vjust = 0.5), size = 3, color = "white") + 
  theme_bw() + 
  xlab("Samples") + 
  ylab("Relative proportion of CpG probes") + 
  labs(fill = "Target genomic location") +  
  scale_fill_manual(values = c("#C2963D","#C23DAC","#3D69C2", "#3DC253")) +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) +
  scale_y_continuous(expand = c(0, 0))


ggsave(paste0(results_path,"/4_CpG_proportion/Island_proportion_plot_07-03-2024.png"), 
       island_proportion_plot, width = 7, height = 4.5, dpi = 600)



#Plot type of probes (I and II)
type_of_probe_proportion <- data.frame(matrix(nrow=ncol(betas_cg_subset_coll), ncol=2))
colnames(type_of_probe_proportion) <- c("I", "II")
#sample <- colnames(betas_cg_collapsed)[4]
n <- 1

for (sample in colnames(betas_cg_subset_coll)){
  betas_cg_collapsed_subset <- betas_cg_subset_coll[,sample]
  cpg_names_subset <- names(na.omit(betas_cg_collapsed_subset))
  annEPIC_v2_subset <- annEPIC_v2[annEPIC_v2$Name_cpg %in% cpg_names_subset,]  
  annEPIC_v2_subset_no_dup <- annEPIC_v2_subset[!duplicated(annEPIC_v2_subset$Name_cpg),]
  type_of_probe <- table(annEPIC_v2_subset_no_dup$Type)
  if(sum(names(type_of_probe) != colnames(type_of_probe_proportion)) != 0){break}
  rownames(type_of_probe_proportion)[n] <- sample
  type_of_probe_proportion[n,] <- type_of_probe
  n <- n + 1
} 


type_of_probe_proportion <- rownames_to_column(type_of_probe_proportion, var = "Samples")
type_of_probe_proportion_long <- pivot_longer(type_of_probe_proportion, 
                                          cols = -Samples, 
                                          names_to = "Type_of_probe", 
                                          values_to = "Number_cpg")

type_of_probe_proportion_long$Samples <- factor(type_of_probe_proportion_long$Samples, levels = sorted_codes) 

#Calculate percentage
type_of_probe_proportion_long <- type_of_probe_proportion_long %>%
  group_by(Samples) %>%
  mutate(Proportion = Number_cpg / sum(Number_cpg),
         PercentLabel = paste0(round(Proportion * 100), "%"))


#Plot
type_of_probe_proportion_plot <- ggplot(type_of_probe_proportion_long, aes(fill = Type_of_probe, y = Number_cpg, x = Samples)) + 
  geom_bar(position = "fill", stat = "identity") + 
  geom_text(aes(label = PercentLabel), 
            position = position_fill(vjust = 0.5), 
            size = 3, color = "white") + 
  theme_bw() + 
  xlab("Samples") + 
  ylab("Relative proportion of CpG probes") + 
  labs(fill = "Type of probe") +  
  scale_fill_manual(values = c("#664DB2", "#99B24D")) +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) +
  scale_y_continuous(expand = c(0, 0))

# type_of_probe_proportion_plot <- ggplot(type_of_probe_proportion_long, aes(fill=Type_of_probe, y = Number_cpg, x = Samples)) + 
#   geom_bar(position="fill", stat="identity") + 
#   theme_bw() + xlab("Samples") + ylab("Percentage of CpGs") + labs(fill = "Type of probe") +  
#   scale_fill_manual(values=c("#664DB2", "#99B24D")) +
#   theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0)) +
#   scale_y_continuous(expand = c(0, 0)) 

ggsave(paste0(results_path,"/4_CpG_proportion/Type_of_probe_07-03-2024.png"), 
       type_of_probe_proportion_plot, width = 7, height = 4.5, dpi = 600)



#Plot gene location of probes (Promoter, Gene body, and Intergenic region)
library(EnsDb.Hsapiens.v111)
edb <- EnsDb.Hsapiens.v86

#For revision
library(AnnotationHub)
library(ensembldb)
ah <- AnnotationHub()

# Search for latest EnsDb database (e.g., for Homo sapiens)
query(ah, c("EnsDb", "Homo sapiens"))
ensdb <- ah[["AH116291"]]
edb <- ensdb
rm(ensdb)
#Get gene and promoters
tx <- transcripts(edb)
prom <- promoters(edb, upstream = 2000, downstream = 200)

#Create Grange object
annEPIC_v2_gr <- GRanges(
  seqnames = Rle(annEPIC_v2$chr),
  ranges = IRanges(start = annEPIC_v2$pos  , end = annEPIC_v2$pos),  # single nucleotide positions
  strand = Rle(annEPIC_v2$strand),
  metadata1 = annEPIC_v2$Name)

#Tidy up promoter and transcript Grange object
prom_no_mt <- granges(prom)
mcols(prom_no_mt) <- rep("Promoter",length(prom))
colnames(mcols(prom_no_mt)) <- "Promoter"
seqlevels(prom_no_mt) <- paste0("chr", seqlevels(prom_no_mt))

tx_no_mt <- granges(tx)
mcols(tx_no_mt) <- rep("Gene Body",length(prom))
colnames(mcols(tx_no_mt)) <- "Gene"
seqlevels(tx_no_mt) <- paste0("chr", seqlevels(tx_no_mt))

#Merge Granges object
m <- findOverlaps(prom_no_mt, annEPIC_v2_gr)
mcols(annEPIC_v2_gr)$Promoter <- NA
mcols(annEPIC_v2_gr)[subjectHits(m), "Promoter"] = mcols(prom_no_mt)[queryHits(m), "Promoter"]

m <- findOverlaps(tx_no_mt, annEPIC_v2_gr)
mcols(annEPIC_v2_gr)$Gene <- NA
mcols(annEPIC_v2_gr)[subjectHits(m), "Gene"] = mcols(tx_no_mt)[queryHits(m), "Gene"]

#Add info about Gene location in main annotation file
annEPIC_v2_gene <- as.data.frame(mcols(annEPIC_v2_gr))
annEPIC_v2_gene$Gene_location <- ifelse(is.na(annEPIC_v2_gene$Promoter), NA, annEPIC_v2_gene$Promoter) 
annEPIC_v2_gene$Gene_location <- ifelse(!is.na(annEPIC_v2_gene$Gene) & is.na(annEPIC_v2_gene$Promoter) & is.na(annEPIC_v2_gene$Gene_location), annEPIC_v2_gene$Gene, annEPIC_v2_gene$Gene_location) 
annEPIC_v2_gene$Gene_location[is.na(annEPIC_v2_gene$Gene_location)] <- "Intergenic"
annEPIC_v2_gene <- annEPIC_v2_gene[, c(1,4)]
colnames(annEPIC_v2_gene)[1] <- "Name"

annEPIC_v2_complete <- full_join(annEPIC_v2, annEPIC_v2_gene, by = "Name")


# Official plot part
Gene_location_proportion <- data.frame(matrix(nrow=ncol(betas_cg_subset_coll), ncol=3))
colnames(Gene_location_proportion) <- c("Gene Body", "Intergenic", "Promoter")
#sample <- colnames(betas_cg_subset_coll)[4]
n <- 1

for (sample in colnames(betas_cg_subset_coll)){
  betas_cg_collapsed_subset <- betas_cg_subset_coll[,sample]
  cpg_names_subset <- names(na.omit(betas_cg_collapsed_subset))
  annEPIC_v2_subset <- annEPIC_v2_complete[annEPIC_v2_complete$Name_cpg %in% cpg_names_subset,]  
  annEPIC_v2_subset_no_dup <- annEPIC_v2_subset[!duplicated(annEPIC_v2_subset$Name_cpg),]
  Gene_location <- table(annEPIC_v2_subset_no_dup$Gene_location)
  if(sum(names(Gene_location) != colnames(Gene_location_proportion)) != 0){break}
  rownames(Gene_location_proportion)[n] <- sample
  Gene_location_proportion[n,] <- Gene_location
  n <- n + 1
} 


Gene_location_proportion <- rownames_to_column(Gene_location_proportion, var = "Samples")
Gene_location_proportion_long <- pivot_longer(Gene_location_proportion, 
                                              cols = -Samples, 
                                              names_to = "Gene_location", 
                                              values_to = "Number_cpg")

Gene_location_proportion_long$Samples <- factor(Gene_location_proportion_long$Samples, levels = sorted_codes) 

#Add percentages
Gene_location_proportion_long <- Gene_location_proportion_long %>%
  group_by(Samples) %>%
  mutate(Proportion = Number_cpg / sum(Number_cpg),
         PercentLabel = paste0(round(Proportion * 100), "%"))


#Plot
Gene_location_proportion_plot <- ggplot(Gene_location_proportion_long, aes(fill = Gene_location, y = Number_cpg, x = Samples)) + 
  geom_bar(position = "fill", stat = "identity") + 
  geom_text(aes(label = PercentLabel), 
            position = position_fill(vjust = 0.5), 
            size = 3, color = "white") + 
  theme_bw() + 
  xlab("Samples") + 
  ylab("Relative proportion of CpG probes") + 
  labs(fill = "Target gene location") +  
  scale_fill_manual(values = c("#8138C7", "#C78138", "#38C781")) +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) +
  scale_y_continuous(expand = c(0, 0))

ggsave(paste0(results_path,"/4_CpG_proportion/Gene_location_07-03-2024.png"), 
       Gene_location_proportion, width = 7, height = 4.5, dpi = 600)



#All plots toghether
combined_plot_ditribtution <- ggarrange(type_of_probe_proportion_plot, island_proportion_plot, Gene_location_proportion_plot,
                                        labels = c("", "", ""),
                                        ncol = 1, nrow = 3)
ggsave(paste0(results_path,"/4_CpG_proportion/Proportion_cpg_genome_15-04-2025.png"), 
       combined_plot_ditribtution, width = 10, height = 10, dpi = 600, bg = "white")



#### 13. Genotype SNPs ---------------------------------------
library(tidyverse)
anno = read_tsv(sesameAnno_download("EPICv2.hg38.snp.tsv.gz"))
sdfs = lapply(searchIDATprefixes(idat_dir), readIDATpair)

genotypes <- data.frame(matrix(nrow=65, ncol=33))


for (n in (1:length(sdfs))){
  vcf = formatVCF(sdfs[[n]], anno)
  vcf_rs <- vcf[grepl("Probe_ID=rs",vcf$INFO),]
  
  genotypes[[n]] <- sub(".*GT=([0-9/]{3}).*", "\\1", vcf_rs$INFO)
  colnames(genotypes)[n] <- names(sdfs[n])
  
  #write.table(vcf[grepl("Probe_ID=rs",vcf$INFO),], file=paste0(results_path,"/5_SNP_probes/", names(sdfs[n]), ".vcf"), quote=FALSE, sep='\t', 
  #            row.names = FALSE)
  
}

colnames(genotypes) <- metadata$Sample_code
rownames(genotypes) <- str_extract(vcf_rs$INFO, "(?<=Probe_ID=)[^;]+")


check_genotypes <- data.frame(matrix(nrow=65, ncol=32))

for (n in 2:length(sdfs)){
  check_genotypes[,n-1] <- genotypes[,1] == genotypes[,n]
  colnames(check_genotypes)[n-1] <- colnames(genotypes)[n]
}

rownames(check_genotypes) <- rownames(genotypes)

check_genotypes <- check_genotypes[, sorted_codes[-1]]
check_genotypes <- cbind(" "=rownames(check_genotypes), check_genotypes)
write_xlsx(check_genotypes, paste0(results_path, "/5_SNP_probes/Genotypes.xlsx"))
