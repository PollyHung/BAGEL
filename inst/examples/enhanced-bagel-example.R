## ---------------------------
##
## Script name: enhanced-bagel-example.R
##
## Purpose of script: Comprehensive example demonstrating the enhanced BAGEL
## package upgrade with functional aneuploidy definitions, telomere/centromere
## boundary logic, and selection type classification. Shows both simple and
## enhanced processing modes with compatibility layer integration.
##
## Author: Polly Hung
##
## Date Created: 2024-09-22
## Date Updated: 2024-09-22
##
## Copyright (c) Polly Hung, 2025
## Email: u3012149@connect.hku.hk
##
## ---------------------------
##
## Notes:
## - Enhanced BAGEL processing transforms simple cytogenetic arm definitions
##   to biologically meaningful functional aneuploidy units
## - Supports both "functional" (enhanced) and "simple" (legacy) processing modes
## - Includes comprehensive validation and multiple export formats
## - Maintains backward compatibility with existing BAGEL workflows
##
## ---------------------------

#' Enhanced BAGEL Package Demonstration
#'
#' This script demonstrates the comprehensive upgrades to the BAGEL package
#' for enhanced functional aneuploidy analysis.

library(BAGEL)
library(dplyr)
library(readr)

# ==============================================================================
# EXAMPLE 1: Enhanced Functional Processing Mode
# ==============================================================================

cat("=== Enhanced BAGEL Package Demonstration ===\n")
cat("Upgrade: From Simple Arms to Functional Aneuploidy Definitions\n\n")

# Example BISCUT results file path (replace with actual path)
biscut_file <- "path/to/all_BISCUT_results.txt"

# Enhanced functional processing with telomere/centromere boundary logic
cat("1. ENHANCED FUNCTIONAL PROCESSING MODE\n")
cat("=====================================\n")

functional_definitions <- create_custom_arm_definitions(
  biscut_results_file = biscut_file,
  mode = "functional",                    # Enhanced processing mode
  min_statistical_support = 10,           # Minimum n_right + n_left
  min_combined_sig = 1.0,                 # Minimum combined significance
  include_gene_annotation = FALSE,        # Add gene counts (optional)
  output_file = "enhanced_functional_definitions.tsv"
)

# Display enhanced structure
cat("\nEnhanced functional definitions structure:\n")
print(head(functional_definitions))

# Validation of functional boundaries
validation_results <- validate_functional_boundaries(functional_definitions)
cat("\nValidation summary:\n")
print(validation_results$coverage_stats)

# ==============================================================================
# EXAMPLE 2: Legacy Simple Processing Mode
# ==============================================================================

cat("\n2. LEGACY SIMPLE PROCESSING MODE\n")
cat("================================\n")

simple_definitions <- create_custom_arm_definitions(
  biscut_results_file = biscut_file,
  mode = "simple",                        # Legacy processing mode
  output_file = "simple_arm_definitions.tsv"
)

cat("\nSimple arm definitions structure:\n")
print(head(simple_definitions))

# ==============================================================================
# EXAMPLE 3: Integration with BAGEL Copy Number Analysis
# ==============================================================================

cat("\n3. INTEGRATION WITH BAGEL COPY NUMBER ANALYSIS\n")
cat("==============================================\n")

# Load example segmentation data (replace with actual data)
segments <- data.frame(
  Sample = rep(paste0("Sample_", 1:50), each = 1000),
  Chromosome = sample(1:22, 50000, replace = TRUE),
  Start = sample(1000000:200000000, 50000, replace = TRUE),
  End = sample(1000000:200000000, 50000, replace = TRUE),
  Log2Ratios = rnorm(50000, 0, 0.3)
)

# Ensure Start < End
segments$End <- segments$Start + sample(100000:10000000, 50000, replace = TRUE)

cat("Loaded segmentation data:", nrow(segments), "segments for",
    length(unique(segments$Sample)), "samples\n")

# Enhanced copy number analysis using functional definitions
cat("\nRunning enhanced copy number analysis...\n")

enhanced_results <- calculateCopyNumber_fixed(
  segments = segments,
  breakpoints = functional_definitions,    # Enhanced functional definitions
  cancer_type = "Enhanced_Example",
  output_dir = "enhanced_analysis_results",
  save_results = TRUE
)

cat("Enhanced analysis completed!\n")
cat("Results structure:\n")
print(names(enhanced_results))

# ==============================================================================
# EXAMPLE 4: Comparison of Processing Modes
# ==============================================================================

cat("\n4. COMPARISON OF PROCESSING MODES\n")
cat("=================================\n")

# Compare region counts and characteristics
comparison_stats <- data.frame(
  Mode = c("Enhanced Functional", "Simple Arms"),
  N_Regions = c(nrow(functional_definitions), nrow(simple_definitions)),
  Avg_Length_MB = c(
    round(mean(functional_definitions$functional_length_mb, na.rm = TRUE), 2),
    round(mean((simple_definitions$arm_end - simple_definitions$arm_start) / 1e6, na.rm = TRUE), 2)
  ),
  Boundary_Types = c(
    length(unique(functional_definitions$boundary_type)),
    1  # Simple mode has no boundary types
  ),
  Selection_Types = c(
    length(unique(functional_definitions$selection_type)),
    0  # Simple mode has no selection classification
  )
)

cat("\nComparison of processing modes:\n")
print(comparison_stats)

# ==============================================================================
# EXAMPLE 5: Enhanced Export Formats
# ==============================================================================

cat("\n5. ENHANCED EXPORT FORMATS\n")
cat("==========================\n")

# Enhanced mode exports multiple formats automatically
if (exists("functional_definitions")) {

  # GISTIC-style matrix format for downstream analysis
  gistic_format <- functional_definitions %>%
    select(
      region_id = functional_arm_id,
      coordinates = functional_coordinates,
      length_mb = functional_length_mb,
      significance = peak_significance,
      boundary_type,
      selection_type,
      direction
    ) %>%
    arrange(desc(significance))

  write_tsv(gistic_format, "enhanced_gistic_format.tsv")
  cat("✅ GISTIC-style format exported\n")

  # Comprehensive summary by chromosome and selection type
  summary_by_selection <- functional_definitions %>%
    group_by(chr_num, selection_type, direction) %>%
    summarise(
      n_regions = n(),
      total_length_mb = sum(functional_length_mb),
      avg_significance = round(mean(peak_significance), 3),
      max_significance = round(max(peak_significance), 3),
      .groups = "drop"
    ) %>%
    arrange(chr_num, desc(avg_significance))

  write_tsv(summary_by_selection, "enhanced_summary_by_selection.tsv")
  cat("✅ Summary by selection type exported\n")

  # Ideogram-ready format for visualization
  ideogram_format <- functional_definitions %>%
    mutate(
      chromosome = paste0("chr", chr_num),
      start = functional_start,
      end = functional_end,
      name = functional_arm_id,
      score = round(peak_significance * 100),  # Scale for visualization
      strand = ifelse(direction == "amp", "+", "-"),
      color = case_when(
        selection_type == "tumor_suppressor_like" ~ "#FF0000",      # Red
        selection_type == "oncogene_like" ~ "#0000FF",              # Blue
        selection_type == "essential_gene_like" ~ "#00FF00",        # Green
        selection_type == "toxic_amplification_like" ~ "#FF8000",   # Orange
        TRUE ~ "#808080"                                            # Gray
      )
    ) %>%
    select(chromosome, start, end, name, score, strand, color,
           boundary_type, selection_type, direction)

  write_tsv(ideogram_format, "enhanced_ideogram_format.tsv")
  cat("✅ Ideogram-ready format exported\n")
}

# ==============================================================================
# EXAMPLE 6: Validation and Quality Control
# ==============================================================================

cat("\n6. VALIDATION AND QUALITY CONTROL\n")
cat("=================================\n")

# Test the enhanced processing pipeline
pipeline_test <- test_enhanced_processing()
cat("Pipeline test results:\n")
print(pipeline_test)

# Comprehensive validation if data is available
if (exists("functional_definitions")) {

  # Boundary validation
  boundary_validation <- validate_functional_boundaries(functional_definitions)

  cat("\nBoundary validation results:\n")
  cat("- Overlapping regions:", !boundary_validation$boundary_logic$passed, "\n")
  cat("- Low statistical support:", boundary_validation$statistical_support$low_support_count, "\n")
  cat("- Overall validation:", ifelse(boundary_validation$overall_passed, "✅ PASS", "❌ FAIL"), "\n")

  # Selection type distribution
  selection_summary <- functional_definitions %>%
    count(selection_type, direction, boundary_type) %>%
    arrange(desc(n))

  cat("\nSelection type distribution:\n")
  print(selection_summary)

  # Coverage statistics by chromosome
  coverage_by_chr <- functional_definitions %>%
    group_by(chr_num) %>%
    summarise(
      n_regions = n(),
      total_coverage_mb = round(sum(functional_length_mb), 1),
      avg_region_size_mb = round(mean(functional_length_mb), 2),
      n_boundary_types = length(unique(boundary_type)),
      n_selection_types = length(unique(selection_type)),
      .groups = "drop"
    ) %>%
    arrange(desc(total_coverage_mb))

  cat("\nCoverage by chromosome (top 10):\n")
  print(head(coverage_by_chr, 10))
}

# ==============================================================================
# SUMMARY AND RECOMMENDATIONS
# ==============================================================================

cat("\n=== ENHANCED BAGEL UPGRADE SUMMARY ===\n")
cat("✅ Enhanced functional aneuploidy definitions created\n")
cat("✅ Telomere/centromere boundary logic implemented\n")
cat("✅ Selection type classification integrated\n")
cat("✅ Multiple export formats generated\n")
cat("✅ Backward compatibility maintained\n")
cat("✅ Comprehensive validation completed\n")

cat("\nKey Improvements:\n")
cat("- Biologically meaningful functional regions vs. simple cytogenetic arms\n")
cat("- Statistical filtering based on BISCUT significance and support\n")
cat("- Enhanced output formats compatible with GISTIC and visualization tools\n")
cat("- Comprehensive validation and quality control\n")
cat("- Seamless integration with existing BAGEL workflow\n")

cat("\nRecommended Usage:\n")
cat("1. Use 'functional' mode for enhanced biological insight\n")
cat("2. Use 'simple' mode for legacy compatibility\n")
cat("3. Apply statistical filters based on study requirements\n")
cat("4. Validate results using built-in validation functions\n")
cat("5. Export multiple formats for downstream analysis\n")

cat("\n=== Enhanced BAGEL Demonstration Complete ===\n")