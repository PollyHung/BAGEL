## ---------------------------
##
## Script name: runMe.R
##
## Purpose of script: Complete BAGEL v2.5 analysis pipeline for all available
## cancer types using GISTIC2-style parameters. Performs breakpoint-based 
## arm-level copy number analysis with statistical significance testing.
##
## Author: Polly Hung
##
## Date Created: 2025-09-12
## Date Updated: 2025-09-12
##
## Copyright (c) Polly Hung, 2025
## Email: u3012149@connect.hku.hk
##
## ---------------------------
##
## Notes:
## - Uses updated GISTIC2 parameters (amp_threshold=0.25, del_threshold=-0.25)
## - Processes all cancer types with BISCUT breakpoints and BAGEL analysis
## - Skips Prologue section (BISCUT already run)
## - Generates comprehensive reports and visualizations for each cancer type
##
## ---------------------------


# Load Libraries 
library(BAGEL)
library(BISCUT)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(patchwork)


# Some Preparations ============================================================
# Set Global Control 
data_dir <- "/Users/polly_hung/Desktop/BAGEL/data/biscut+bagel"

# Discover all available cancer types (excluding .zip files)
all_cancer_dirs <- list.dirs(data_dir, recursive = FALSE)
all_cancer_types <- basename(all_cancer_dirs)
valid_cancer_types <- all_cancer_types[!grepl("\\.zip$", all_cancer_types)]

# Take Ovarian Cancer as an Example 
cancer_type <- "tcga_ovarian_serous_cystadenocarcinoma"

# Define Pathways 
output_dir <- file.path(data_dir, cancer_type, "bagel_v4_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Setting up Logging 
log_bagel(log_level = "INFO", log_file = file.path(output_dir, "analysis.log"))


# Creation of Breakpoints ======================================================
# Run BISCUT
run_biscut_pipeline(cancer_folder = cancer_type,
                    results_dir = data_dir,
                    cores = 4)


# Clean up the Breakpoints =====================================================
# Load Breakpoint Information using unified function
custom_biscut_file <- list.files(path = file.path(data_dir, cancer_type),
                                 pattern = "all_BISCUT_results.txt",
                                 recursive = TRUE,
                                 full.names = TRUE)
if (length(custom_biscut_file) == 0) { custom_biscut_file <- NULL } 
arm_definitions <- define_arm(custom_biscut_file = custom_biscut_file,
                              cancer_type = cancer_type,
                              min_statistical_support = 10,
                              min_combined_sig = 1.0) 


# Run Main BAGEL Function ======================================================
# Load Segmentation File 
seg_file <- file.path(data_dir, cancer_type, "segmentation.seg")
segments <- load_segmentation(cancer_type, data_dir)
bagel_results <- calculate_copynumber(segments = segments, 
                                      breakpoints = arm_definitions,
                                      amp_threshold = 0.25, 
                                      del_threshold = -0.25,
                                      stringent_threshold = 0.95)
create_matrix(bagel_results$arm_summaries, output_dir)


# Inspect and Plot the Result ==================================================
# Generate Chromosome Ideograms
ideogram_results <- plot_ideograms(bagel_results = bagel_results,
                                   output_dir = output_dir,
                                   save_plots = TRUE)












