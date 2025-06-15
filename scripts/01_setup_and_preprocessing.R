# Project: DNA degradation study
# Script: 01_setup_and_preprocessing.R
# Description: This script loads necessary libraries, defines file paths,
# reads the raw IDAT files, preprocesses them to get beta values,
# and tidies the associated metadata.

#### 1. Configuration & Setup ------------------------------------------------

# --- EDIT THESE PATHS ---
# Use the 'here' package for portable paths relative to your project root
library(here)
idat_dir <- here("data", "idat") # Assumes .idat files are in data/idat/
metadata_path <- here("data", "SampleSheet.xlsx") # Assumes metadata sheet is in data/
results_dir <- here("results")
processed_data_dir <- here("data", "processed")

# Create directories if they don't exist
dir.create(results_dir, showWarnings = FALSE)
dir.create(processed_data_dir, showWarnings = FALSE)

# --- Load Libraries ---
library(sesame)
library(sesameData)
library(readxl)
library(dplyr)
library(stringr)
library(tibble)

sesameDataCache() # Cache sesame annotation data

#### 2. Upload and Preprocess Data ------------------------------------------

# Read IDAT files and compute beta values using the "QCDB" pipeline
# This includes quality control and dye-bias correction.
betas <- openSesame(idat_dir, prep = "QCDPB", func = getBetas)

# Load metadata from the sample sheet
metadata <- read_xlsx(path = metadata_path, sheet = 2)

#### 3. Tidy and Align Datasets ----------------------------------------------

# Create a unique chip position identifier to sort and align data
metadata$Chip_pos <- paste0(metadata$Sentrix_ID, "_", metadata$Sentrix_position)
metadata <- metadata[order(metadata$Chip_pos), ]

# Check if metadata and beta matrix columns are in the same order
if (!all(metadata$Chip_pos == colnames(betas))) {
  stop("Mismatch between metadata and beta value column names. Please check your files.")
}

# Extract DNA size, input amount, and replicate info from sample names
metadata <- metadata %>%
  mutate(
    DNA_size = as.numeric(str_extract(Sample_name_after_bioanalyser_measurement, "\\d+(?=bp)")),
    DNA_input = as.numeric(str_extract(Sample_name_after_bioanalyser_measurement, "\\d+(?=ng)")),
    Duplicate = str_extract(Sample_name_after_bioanalyser_measurement, "[A-Za-z]$"),
    Sample_code = sub("^[^_]*_", "", Sample_name_after_bioanalyser_measurement)
  )

# Manually correct the reference sample information
metadata$DNA_input[1] <- 250
metadata$DNA_size[1] <- NA # Or a descriptive string like "Control" / "Not Fragmented"
metadata$Sample_code[1] <- "no_degraded_250ng" # Or a descriptive string like "Control" / "Not Fragmented"
metadata$Duplicate[1] <- NA

# Assign the clean sample codes to the beta matrix columns
colnames(betas) <- metadata$Sample_code

#### 4. Save Processed Data -------------------------------------------------

# Save the tidy beta values and metadata for downstream analysis
saveRDS(betas, file.path(processed_data_dir, "betas_preprocessed.rds"))
saveRDS(metadata, file.path(processed_data_dir, "metadata_preprocessed.rds"))

message("Preprocessing complete. Beta values and metadata saved to ", processed_data_dir)