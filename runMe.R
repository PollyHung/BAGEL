# BAGEL v2.0 Demonstration Script
# Example analysis using ovarian_serous_cystadenocarcinoma
# This script demonstrates the complete BAGEL v2.0 workflow in single cancer 


# Load Libraries 
library(BAGEL)
library(BISCUT)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(patchwork)

# Set Global Control 
data_dir <- "/Users/polly_hung/Desktop/BAGEL/validation"
cancer_type <- "ovarian_serous_cystadenocarcinoma"
output_dir <- file.path(data_dir, cancer_type, "bagel_v2_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set Logging Information 
setup_bagel_logging(log_level = "INFO", log_file = file.path(output_dir, "analysis.log"))


# ==============================================================================
# Prologue - Creation of Breakpoints
# ==============================================================================
# Run BISCUT 
tryCatch({
  run_biscut_pipeline(cancer_folder = cancer_type, 
                      results_dir = data_dir, 
                      cores = 4)
}, error = function(e) {
  cat("❌ ERROR during BISCUT\n")
  cat("Error message:", e$message, "\n")
  stop(e)
})


# ==============================================================================
# Chapter 1 - Functionally Defined Aneuploidy 
# ==============================================================================
# Run BAGEL 
tryCatch({
  cat("=== Step 1: Load Data ===\n") # =========================================
  
  # Load Segmentation File 
  seg_file <- file.path(data_dir, cancer_type, "segmentation.seg")
  segments <- load_segmentation_data(cancer_type, data_dir)
  
  # Load Breakpoint Information 
  custom_biscut_file <- list.files(path = file.path(data_dir, cancer_type), 
                                   pattern = "all_BISCUT_results.txt", 
                                   recursive = TRUE, 
                                   full.names = TRUE)
  if (file.exists(custom_biscut_file)) {
    tryCatch({
      arm_definitions <- create_custom_arm_definitions(custom_biscut_file)
    }, error = function(e) {
      cat("⚠️ Error processing TCGA BISCUT results:", e$message, "\n")
      cat("Falling back to pre-defined TCGA consensus breakpoints\n")
      breakpoint_data <- load_breakpoint_data()
      arm_definitions <- get_arm_definitions(tcga_cancer_type, breakpoint_data)
    })
  } ## IF file does not exist, fall back to TCGA breakpoints 
  
  
  cat("=== Step 2: Run main BAGEL function ===\n") # ===========================
  bagel_results <- calculateCopyNumber_fixed(
    segments = segments,
    breakpoints = arm_definitions,
    amp_threshold = log2(2.5/2),       # Amplification threshold
    del_threshold = log2(1.5/2),       # Deletion threshold
    stringent_threshold = 0.9,         # Stringent threshold
    output_dir = output_dir,
    cancer_type = cancer_type,
    use_gistic = TRUE,                 # Enable GISTIC analysis
    save_results = TRUE                # Save intermediate results
  )
  
  
  cat("=== Step 3: Organise and Retrieve Results ===\n") # =====================
  
  # Create chromosome arm copy number matrices
  matrices <- create_arm_matrix(bagel_results$arm_summaries, output_dir)
  bagel_results$matrices <- matrices
  
  # Show significant arms if any
  if (!is.null(bagel_results$significant_arms) && nrow(bagel_results$significant_arms) > 0) {
    n = min(c(10, nrow(bagel_results$significant_arms)))
    cat("Significant arms (q < 0.25):", nrow(bagel_results$significant_arms), "\n")
    cat("Top", n, "significant arms:\n")
    top_arms <- head(bagel_results$significant_arms[order(bagel_results$significant_arms$mean_z_score, decreasing = TRUE), ], n)
  } else {
    cat("No significant arms found (q < 0.25)\n")
  }
  
  # Analyze alteration frequencies
  alteration_summary <- matrices$long_format %>%
    dplyr::mutate(Alteration_Type = case_when(
      Copy_Number >= 0.5 ~ "Amplification",
      Copy_Number <= -0.5 ~ "Deletion",
      TRUE ~ "Normal")) %>%
    count(Arm, Alteration_Type) 
  missing <- setdiff(c("Normal", "Amplification", "Deletion"), unique(alteration_summary$Alteration_Type))
  missing_df <- data.frame(Arm = unique(alteration_summary$Arm), 
                           Alteration_Type = rep(missing, length(unique(alteration_summary$Arm))), 
                           n = rep(0, length(unique(alteration_summary$Arm))))  
  alteration_summary <- rbind(alteration_summary, missing_df)
  alteration_summary <- alteration_summary %>% 
    pivot_wider(names_from = Alteration_Type, values_from = n, values_fill = 0) %>% 
    dplyr::mutate(Total_Samples = Normal + Amplification + Deletion,
                  Amp_Freq = round(Amplification / Total_Samples * 100, 1),
                  Del_Freq = round(Deletion / Total_Samples * 100, 1),
                  Alt_Freq = Amp_Freq + Del_Freq) %>%
    arrange(desc(Alt_Freq))
  
  # Generate comprehensive analysis report
  generate_analysis_report(cancer_type, bagel_results, output_dir)
  
  # Save complete results
  save(bagel_results, file = file.path(output_dir, "bagel_v2_complete_results.RData"))
  
  
  cat("=== Step 4: Generate Chromosome Ideograms ===\n") # =====================
  tryCatch({
    ideogram_results <- plot_chromosome_ideograms(
      bagel_results = bagel_results,
      output_dir = output_dir,
      save_plots = TRUE
    )
    if (!is.null(ideogram_results)) {
      n_chromosomes <- length(ideogram_results$individual_plots)
      cat("✅ Generated ideograms for", n_chromosomes, "chromosomes\n")
    } else {
      cat("⚠️ No ideograms generated (no valid chromosome data)\n")
    }
    
  }, error = function(e) {
    cat("❌ Error generating ideograms:", e$message, "\n")
  })
  
  
  cat("=== Step 5: Verify All File Exists ===\n") # ============================
  summarize_bagel_outputs(output_dir = output_dir)
  
}, error = function(e) {
  cat("❌ ERROR during BAGEL\n")
  cat("Error message:", e$message, "\n")
  stop(e)
})


# ==============================================================================
# Chapter 2 - Cross Validation Against TCGA 
# ==============================================================================
query_path <- "~/Desktop/BAGEL/validation/ovarian_serous_cystadenocarcinoma/bagel_v2_analysis/BAGEL_results_ovarian_serous_cystadenocarcinoma/arm_definitions.txt"
ref_cancer <- ""
results <- analyze_breakpoint_similarity(query_path = query_path,
                                         reference = ref_cancer,
                                         output_dir = file.path(data_dir, cancer_type, "similarity_with_"))


# ==============================================================================
# Chapter 3 - Clustering 
# ==============================================================================




















