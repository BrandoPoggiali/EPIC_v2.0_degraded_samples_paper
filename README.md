# DNA Degradation EPIC Array Analysis

This repository contains the R scripts for the analysis of EPIC array data from a DNA degradation study. The study investigates the impact of different DNA input amounts and fragment sizes on DNA methylation measurements.

## Project Description

The analysis explores quality control metrics, inter-sample correlation, sources of data variability (delta-beta), genomic context of affected probes, and the performance of epigenetic age prediction models under various conditions of DNA degradation. It also compares the performance of the `pOOBAH` and `ELBAR` preprocessing algorithms.

**DNA Fragment Sizes**: 350, 230, 165, and 95 bp
**DNA Amounts**: 100 ng, 50 ng, 20 ng, and 10 ng (in duplicates)
**Reference Sample**: 250 ng, non-degraded

## How to Run the Analysis

### 1. Installation

You will need R and the following packages. You can install them by running the code below in your R console.

```R
# CRAN Packages
install.packages(c("dplyr", "tibble", "tidyr", "tidyverse", "pheatmap", "factoextra", "writexl", "ggpubr", "reshape2", "plotly", "stringr", "openxlsx", "scales", "ggridges", "ggplot2", "ggpointdensity", "VennDiagram", "corrplot", "here", "remotes"))

# Bioconductor Packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("sesame", "sesameData", "methylclock", "SummarizedExperiment", "ExperimentHub", "ensembldb", "AnnotationHub"))

# GitHub Packages
remotes::install_github("BrandoPoggiali/EPIC_v2.0_degraded_samples_paper")
