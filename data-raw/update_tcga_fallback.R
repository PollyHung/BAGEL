## ---------------------------
##
## Script name: update_tcga_fallback.R
##
## Purpose of script:
##   Consolidate updated BAGEL v4 arm_definitions from 29 TCGA cancer types
##   into a single R data object for use as TCGA fallback in define_arm.R
##   Input: Individual arm_definitions.txt files from biscut+bagel analysis
##   Output: updated_tcga_breakpoints.rda for package data/
##
## Author: Polly Hung
##
## Date Created: 2025-09-30
## Date Updated: 2025-09-30
##
## Copyright (c) Polly Hung, 2025
## Email: u3012149@connect.hku.hk
##
## ---------------------------
##
## Notes:
## - Reads from /Users/polly_hung/Desktop/BAGEL/data/biscut+bagel/
## - Creates updated_tcga_breakpoints list with 29 cancer types
## - Use this to replace cancer_specific_breakpoints in define_arm.R
##
## ---------------------------

library(dplyr)
library(readr)

# Base directory
base_dir <- "/Users/polly_hung/Desktop/BAGEL/data/biscut+bagel"

# Find all TCGA arm_definitions files
bagel_v4_files <- list.files(
  path = base_dir,
  pattern = "arm_definitions\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

# Filter for TCGA only
tcga_files <- bagel_v4_files[grepl("bagel_v4_", bagel_v4_files)]
tcga_files <- tcga_files[grepl("tcga_", tcga_files)]
cat("Found", length(tcga_files), "TCGA arm definition files\n")

# Initialize list to store arm definitions by cancer type
updated_tcga_breakpoints <- list()

# Process each file
for (file_path in tcga_files) {
  # Extract cancer type from path
  # Example: "/path/tcga_breast_invasive_carcinoma/bagel_v4_analysis/arm_definitions.txt"
  cancer_type <- basename(dirname(dirname(file_path)))
  cancer_type <- gsub("^tcga_", "", cancer_type)

  cat("Processing:", cancer_type, "\n")

  # Read arm definitions
  tryCatch({
    arm_defs <- readr::read_delim(file_path, delim = "\t", show_col_types = FALSE)
    arm_defs <- as.data.frame(arm_defs)

    # Validate required columns
    required_cols <- c("chr", "chr_arm", "arm", "functional_start", "functional_end", "direction")
    missing_cols <- setdiff(required_cols, colnames(arm_defs))

    if (length(missing_cols) > 0) {
      warning("Missing columns in ", cancer_type, ": ", paste(missing_cols, collapse = ", "))
      next
    }

    # Store in list
    updated_tcga_breakpoints[[cancer_type]] <- arm_defs
    cat("  ✅ Loaded", nrow(arm_defs), "arm definitions\n")

  }, error = function(e) {
    warning("Failed to load ", cancer_type, ": ", e$message)
  })
}

usethis::use_data(updated_tcga_breakpoints)

# Summary
usethis::use_data(updated_tcga_breakpoints, )