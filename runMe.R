# BAGEL v2.0 Demonstration Script
# Example analysis using adrenocortical_cancer data
# This script demonstrates the complete BAGEL v2.0 workflow

library(BAGEL)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(patchwork)

# Get available cancer types
data_dir <- "/Users/polly_hung/Desktop/BAGEL/results"
cancer_types <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
cancer_types <- cancer_types[!grepl("\\.(tsv|csv)$", cancer_types)]

for(cancer in cancer_types){
  
  # Configuration
  cancer_type <- cancer
  data_dir <- "/Users/polly_hung/Desktop/BAGEL/results"
  output_dir <- file.path(data_dir, cancer_type, "bagel_v2_analysis")
  
  cat("Cancer type is ", cancer)
  
  # Setup output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Setup logging
  setup_bagel_logging(log_level = "INFO", log_file = file.path(output_dir, "analysis.log"))
  
  tryCatch({
    
    cat("=== Step 1: Environment Check ===\n")
    
    # Check if data files exist
    seg_file <- file.path(data_dir, cancer_type, "segmentation.seg")
    if (!file.exists(seg_file)) {
      stop("Segmentation file not found at: ", seg_file)
    }
    cat("✅ Segmentation file found:", seg_file, "\n")
    
    # Check file size
    file_size_mb <- round(file.size(seg_file) / 1024^2, 1)
    cat("✅ File size:", file_size_mb, "MB\n\n")
    
    cat("=== Step 2: Load and Validate Data ===\n")
    
    # Load breakpoint data
    cat("Loading breakpoint data...\n")
    breakpoint_data <- load_breakpoint_data()
    
    available_cancer_types <- get_available_cancer_types(data_dir)
    cat("Available cancer types:", length(available_cancer_types), "\n")
    cat("Target cancer type available:", cancer_type %in% available_cancer_types, "\n")
    
    # Load segmentation data
    cat("Loading segmentation data...\n")
    segments <- load_segmentation_data(cancer_type, data_dir)
    
    cat("✅ Loaded", nrow(segments), "segments from", length(unique(segments$Sample)), "samples\n")
    
    # Show sample of data
    cat("\nSegmentation data preview:\n")
    print(head(segments, 3))
    
    cat("\n=== Step 3: Get Arm Definitions ===\n")
    
    # Get arm definitions for this cancer type
    arm_definitions <- get_arm_definitions(cancer_type, breakpoint_data)
    
    cat("✅ Loaded arm definitions for", nrow(arm_definitions), "chromosome arms\n")
    cat("Arms included:", paste(head(arm_definitions$arm, 10), collapse = ", "), "...\n\n")
    
    cat("=== Step 4: Run BAGEL v2.0 Analysis ===\n")
    cat("This may take several minutes for", length(unique(segments$Sample)), "samples...\n")
    
    # Run complete BAGEL analysis with progress tracking
    start_time <- Sys.time()
    
    bagel_results <- calculateCopyNumber_fixed(
      segments = segments,
      breakpoints = arm_definitions,
      amp_threshold = log2(2.5/2),      # Amplification threshold
      del_threshold = log2(1.5/2),      # Deletion threshold
      stringent_threshold = 0.9,         # Stringent threshold
      output_dir = output_dir,
      cancer_type = cancer_type,
      use_gistic = TRUE,                 # Enable GISTIC analysis
      save_results = TRUE                # Save intermediate results
    )
    
    analysis_time <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
    cat("✅ Analysis completed in", analysis_time, "minutes\n\n")
    
    cat("=== Step 5: Create Output Matrices ===\n")
    
    # Create chromosome arm copy number matrices
    matrices <- create_arm_matrix(bagel_results$arm_summaries, output_dir)
    bagel_results$matrices <- matrices
    
    cat("✅ Created matrices:\n")
    cat("  - Copy number matrix:", nrow(matrices$cn_matrix), "arms ×", ncol(matrices$cn_matrix), "samples\n")
    cat("  - Log2 ratio matrix:", nrow(matrices$log2_matrix), "arms ×", ncol(matrices$log2_matrix), "samples\n")
    cat("  - Long format data:", nrow(matrices$long_format), "arm-sample combinations\n\n")
    
    cat("=== Step 6: Analysis Results Summary ===\n")
    
    # Show key results
    cat("Samples analyzed:", length(unique(bagel_results$segments$Sample)), "\n")
    cat("Chromosome arms:", length(unique(bagel_results$arm_summaries$Arm)), "\n")
    
    # Show significant arms if any
    if (!is.null(bagel_results$significant_arms) && nrow(bagel_results$significant_arms) > 0) {
      cat("Significant arms (q < 0.25):", nrow(bagel_results$significant_arms), "\n")
      cat("Top 5 significant arms:\n")
      top_arms <- head(bagel_results$significant_arms[order(bagel_results$significant_arms$mean_z_score, decreasing = TRUE), ], 5)
      for (i in 1:nrow(top_arms)) {
        cat(sprintf("  %s: z-score = %.2f (q = %.3f)\n", 
                    top_arms$Arm[i], top_arms$mean_z_score[i], top_arms$q_value[i]))
      }
    } else {
      cat("No significant arms found (q < 0.25)\n")
    }
    
    cat("\n=== Step 7: Copy Number Analysis ===\n")
    
    # Analyze alteration frequencies
    alteration_summary <- matrices$long_format %>%
      mutate(
        Alteration_Type = case_when(
          Copy_Number >= 2.5 ~ "Amplification",
          Copy_Number <= 1.5 ~ "Deletion",
          TRUE ~ "Normal"
        )
      ) %>%
      count(Arm, Alteration_Type) %>%
      pivot_wider(names_from = Alteration_Type, values_from = n, values_fill = 0) %>%
      mutate(
        Total_Samples = Normal + Amplification + Deletion,
        Amp_Freq = round(Amplification / Total_Samples * 100, 1),
        Del_Freq = round(Deletion / Total_Samples * 100, 1),
        Alt_Freq = Amp_Freq + Del_Freq
      ) %>%
      arrange(desc(Alt_Freq))
    
    cat("Top 10 most frequently altered arms:\n")
    print(head(alteration_summary[c("Arm", "Amp_Freq", "Del_Freq", "Alt_Freq")], 10))
    
    cat("\n=== Step 8: Generate Analysis Report ===\n")
    
    # Generate comprehensive analysis report
    generate_analysis_report(cancer_type, bagel_results, output_dir)
    
    # Save complete results
    save(bagel_results, file = file.path(output_dir, "bagel_v2_complete_results.RData"))
    
    cat("✅ Analysis report generated\n")
    cat("✅ Complete results saved\n\n")
    
    cat("=== Step 8.5: Generate Chromosome Ideograms ===\n")
    
    # Generate chromosome ideograms showing functional arm definitions
    tryCatch({
      ideogram_results <- plot_chromosome_ideograms(
        bagel_results = bagel_results,
        output_dir = output_dir,
        save_plots = TRUE
      )
      
      if (!is.null(ideogram_results)) {
        n_chromosomes <- length(ideogram_results$individual_plots)
        cat("✅ Generated ideograms for", n_chromosomes, "chromosomes\n")
        cat("✅ Ideogram plots saved to:", ideogram_results$output_directory, "\n")
      } else {
        cat("⚠️ No ideograms generated (no valid chromosome data)\n")
      }
      
    }, error = function(e) {
      cat("❌ Error generating ideograms:", e$message, "\n")
    })
    
    cat("\n")
    
    cat("=== Step 9: Output Files Summary ===\n")
    
    # List all output files
    output_files <- list.files(output_dir, recursive = TRUE, full.names = FALSE)
    cat("Generated", length(output_files), "output files:\n")
    
    # Key output files
    key_files <- c(
      "arm_copynumber_matrix.csv",
      "arm_log2ratio_matrix.csv", 
      "arm_copynumber_long.csv",
      "arm_copynumber_summary.csv",
      "BAGEL_V2_ANALYSIS_REPORT.md",
      "demo_bagel_results.RData"
    )
    
    # BAGEL results files (in BAGEL_results_* subdirectory)
    bagel_results_dir <- list.dirs(output_dir, recursive = FALSE, full.names = FALSE)
    bagel_results_subdir <- bagel_results_dir[grepl("^BAGEL_results_", bagel_results_dir)]
    if (length(bagel_results_subdir) > 0) {
      bagel_subdir_files <- list.files(file.path(output_dir, bagel_results_subdir[1]), full.names = FALSE)
      cat("\nBAGEL results files (in", bagel_results_subdir[1], "):\n")
      expected_bagel_files <- c("analysis_parameters.txt", "arm_level_summaries.txt", "arm_definitions.txt", 
                                "gistic_results.txt", "significant_arms.txt", "stringent_calls.txt")
      for (file in expected_bagel_files) {
        if (file %in% bagel_subdir_files) {
          file_path <- file.path(output_dir, bagel_results_subdir[1], file)
          file_size <- if (file.exists(file_path)) {
            paste0("(", round(file.size(file_path) / 1024, 1), " KB)")
          } else {
            "(not found)"
          }
          cat("  ✅", file, file_size, "\n")
        } else {
          cat("  ❌", file, "(missing)\n")
        }
      }
    }
    
    # Check for ideogram files
    ideogram_dir <- file.path(output_dir, "newly_defined_arms_ideogram")
    if (dir.exists(ideogram_dir)) {
      ideogram_files <- list.files(ideogram_dir, full.names = FALSE)
      cat("\nIdeogram files (in newly_defined_arms_ideogram/):\n")
      
      # Look for key ideogram files
      pdf_files <- ideogram_files[grepl("\\.pdf$", ideogram_files)]
      png_files <- ideogram_files[grepl("\\.png$", ideogram_files)]
      
      if (length(pdf_files) > 0) {
        cat("  📊 Individual chromosome PDFs:", length(pdf_files[grepl("chromosome_\\d+", pdf_files)]), "files\n")
        if (any(grepl("all_chromosomes", pdf_files))) {
          cat("  📊 Combined panel plot: all_chromosomes_ideogram_panel.pdf\n")
        }
      }
      
      if (length(png_files) > 0) {
        cat("  🖼️ PNG versions:", length(png_files), "files\n")
      }
    }
    
    cat("\nKey output files:\n")
    for (file in key_files) {
      if (file %in% output_files) {
        file_path <- file.path(output_dir, file)
        file_size <- if (file.exists(file_path)) {
          paste0("(", round(file.size(file_path) / 1024, 1), " KB)")
        } else {
          "(not found)"
        }
        cat("  ✅", file, file_size, "\n")
      } else {
        cat("  ❌", file, "(missing)\n")
      }
    }
    
    cat("\n=== Step 10: Usage Examples ===\n")
    
    cat("To use the copy number matrix in your analysis:\n\n")
    cat("# R code example:\n")
    cat("library(readr)\n")
    cat("copy_numbers <- read_csv('", file.path(output_dir, "arm_copynumber_matrix.csv"), "')\n")
    cat("log2_ratios <- read_csv('", file.path(output_dir, "arm_log2ratio_matrix.csv"), "')\n")
    cat("long_data <- read_csv('", file.path(output_dir, "arm_copynumber_long.csv"), "')\n\n")
    
    cat("Matrix dimensions:\n")
    cat("- Rows (chromosome arms):", nrow(matrices$cn_matrix), "\n")
    cat("- Columns (samples):", ncol(matrices$cn_matrix), "\n")
    cat("- Copy number range:", round(min(matrices$cn_matrix, na.rm=TRUE), 2), "to", round(max(matrices$cn_matrix, na.rm=TRUE), 2), "\n\n")
    
    cat("=== DEMONSTRATION COMPLETED SUCCESSFULLY ===\n")
    cat("Cancer type:", cancer_type, "\n")
    cat("Total runtime:", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1), "minutes\n")
    cat("Output directory:", output_dir, "\n")
    cat("Report file:", file.path(output_dir, "BAGEL_V2_ANALYSIS_REPORT.md"), "\n\n")
    
    cat("✅ BAGEL v2.0 analysis pipeline completed successfully!\n")
    cat("The chromosome arm copy number matrix is ready for downstream analysis.\n")
    
  }, error = function(e) {
    cat("❌ ERROR during demonstration:\n")
    cat("Error message:", e$message, "\n")
    cat("\nPlease check:\n")
    cat("1. Data files are present in", data_dir, "\n") 
    cat("2. BAGEL package functions are properly loaded\n")
    cat("3. Required R packages are installed\n")
    cat("4. Log file for detailed error information:", file.path(output_dir, "demo_analysis.log"), "\n")
    
    stop(e)
  })
  
}




