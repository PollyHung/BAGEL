#' Calculate Copy Number Analysis with Breakpoint-based Arms (Fixed Version)
#'
#' Performs comprehensive copy number analysis using biologically meaningful chromosome 
#' arm breakpoints. Includes proper GISTIC-style statistical analysis and improved 
#' error handling.
#'
#' @param segments A data frame containing segmentation data with columns:
#'   - Sample: Sample identifier
#'   - Chromosome: Chromosome number (1-22)
#'   - Start: Segment start position
#'   - End: Segment end position
#'   - Log2Ratios: Log2 copy number ratios (or Segment_Mean)
#' @param breakpoints Breakpoint data, either:
#'   - Data frame with arm definitions (new format)
#'   - List with tel_bound element (legacy format)
#'   - Character string specifying cancer type for built-in data
#' @param amp_threshold Amplification threshold in log2 scale (default: log2(2.5/2) = 0.32)
#' @param del_threshold Deletion threshold in log2 scale (default: log2(1.5/2) = -0.415)
#' @param stringent_threshold More stringent threshold for high-confidence calls (default: 0.9)
#' @param output_dir Directory for saving results (default: current directory)
#' @param cancer_type Cancer type for analysis (used for logging and file naming)
#' @param use_gistic Logical, whether to perform GISTIC-style statistical analysis (default: TRUE)
#' @param chunk_size Number of samples to process at once for memory efficiency (default: 100)
#' @param save_results Logical, whether to save results to files (default: TRUE)
#'
#' @return A list containing:
#'   - segments: Processed segments with arm annotations
#'   - arm_summaries: Arm-level copy number summaries  
#'   - gistic_results: GISTIC statistical analysis results (if use_gistic=TRUE)
#'   - significant_arms: Arms with significant copy number alterations
#'   - parameters: Analysis parameters used
#'
#' @examples
#' \dontrun{
#' # Using built-in pan-cancer breakpoints
#' result <- calculateCopyNumber_fixed(segments, 
#'                                   breakpoints = "consensus",
#'                                   cancer_type = "BRCA")
#' 
#' # Using custom breakpoints  
#' custom_breakpoints <- get_breakpoint_data("breast_invasive_carcinoma")
#' result <- calculateCopyNumber_fixed(segments, 
#'                                   breakpoints = custom_breakpoints,
#'                                   cancer_type = "BRCA")
#' }
#'
#' @export
# =============================================================================
# GISTIC2-Style Parameter Defaults Fix
# 
# ISSUE: Original BAGEL used non-standard thresholds:
# - amp_threshold = log2(2.5/2) = 0.32
# - del_threshold = log2(1.5/2) = -0.415
#
# FIX: Use GISTIC2 parameter file values:
# - amp_threshold = 0.25 (from GISTIC2 parameter file)
# - del_threshold = -0.25 (symmetric, GISTIC2 standard)
# This matches standard GISTIC2 behavior and published analyses
# =============================================================================

calculateCopyNumber_fixed <- function(segments,
                                    breakpoints,
                                    # OLD CODE - COMMENTED OUT:
                                    # amp_threshold = log2(2.5/2),
                                    # del_threshold = log2(1.5/2),
                                    
                                    # NEW CODE - GISTIC2 Standard:
                                    amp_threshold = 0.25,    # From GISTIC2 parameter file
                                    del_threshold = -0.25,   # From GISTIC2 parameter file (symmetric)
                                    stringent_threshold = 0.9,
                                    output_dir = ".",
                                    cancer_type = "Unknown",
                                    use_gistic = TRUE,
                                    chunk_size = 100,
                                    save_results = TRUE) {
  
  bagel_log(sprintf("Starting copy number analysis for cancer type: %s", cancer_type), "INFO")
  
  # =============================================================================
  # GISTIC2 Parameter Validation - NEW
  # 
  # PURPOSE: Validate that parameters follow GISTIC2 conventions before analysis
  # =============================================================================
  
  # Validate GISTIC2 parameters
  validated_params <- validate_gistic2_parameters(amp_threshold, del_threshold, stringent_threshold)
  
  # Input validation
  segments <- validate_segments(segments)
  
  # Handle different breakpoint formats
  arm_definitions <- process_breakpoint_input(breakpoints)
  
  # Create output directory if needed
  if (save_results && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    bagel_log(sprintf("Created output directory: %s", output_dir), "INFO")
  }
  
  # Process segments using improved pipeline
  n_samples <- length(unique(segments$Sample))
  
  if (n_samples > chunk_size) {
    bagel_log(sprintf("Processing %d samples in chunks of %d", n_samples, chunk_size), "INFO")
    processed_results <- process_segments_chunked(segments, arm_definitions, chunk_size,
                                                 amp_threshold, del_threshold)
  } else {
    processed_results <- process_segments_pipeline(segments, arm_definitions, 
                                                  amp_threshold, del_threshold)
  }
  
  # Perform GISTIC-style statistical analysis if requested
  gistic_results <- NULL
  significant_arms <- NULL
  
  if (use_gistic) {
    bagel_log("Performing GISTIC-style statistical analysis", "INFO")
    gistic_results <- gistic_analysis(processed_results$segments, 
                                     processed_results$arm_definitions)
    significant_arms <- gistic_results$significant_peaks
  }
  
  # Apply stringent thresholds for high-confidence calls
  stringent_results <- apply_stringent_thresholds(processed_results$arm_summaries, 
                                                 stringent_threshold, 
                                                 amp_threshold, del_threshold)
  
  # Prepare final results
  final_results <- list(
    segments = processed_results$segments,
    arm_summaries = processed_results$arm_summaries,
    stringent_summaries = stringent_results,
    arm_definitions = processed_results$arm_definitions,
    gistic_results = gistic_results,
    significant_arms = significant_arms,
    parameters = list(
      amp_threshold = amp_threshold,
      del_threshold = del_threshold,
      stringent_threshold = stringent_threshold,
      cancer_type = cancer_type,
      use_gistic = use_gistic,
      analysis_date = Sys.Date()
    )
  )
  
  # Save results to files if requested
  if (save_results) {
    save_analysis_results(final_results, output_dir, cancer_type)
  }
  
  # Generate summary statistics
  generate_analysis_summary(final_results)
  
  bagel_log(sprintf("Copy number analysis completed for %s", cancer_type), "INFO")
  
  return(final_results)
}

#' Process Breakpoint Input
#'
#' Handle different formats of breakpoint input and return standardized format
#'
#' @param breakpoints Input breakpoints in various formats
#' @return Standardized arm definitions
#' @keywords internal
process_breakpoint_input <- function(breakpoints) {
  
  if (is.character(breakpoints)) {
    # Built-in breakpoint data
    if (breakpoints == "consensus") {
      return(get_breakpoint_data("consensus", "arm_definitions"))
    } else {
      return(get_breakpoint_data(breakpoints, "arm_definitions"))
    }
  } else if (is.data.frame(breakpoints)) {
    # Direct data frame input
    return(validate_breakpoints(breakpoints))
  } else if (is.list(breakpoints) && "tel_bound" %in% names(breakpoints)) {
    # Legacy format
    return(convert_legacy_breakpoints(breakpoints))
  } else {
    stop("Invalid breakpoints format. Must be character string, data frame, or legacy list format.")
  }
}

#' Apply Stringent Thresholds for High-Confidence Calls
#'
#' @param arm_summaries Arm-level summaries
#' @param stringent_threshold Stringent threshold value
#' @param amp_threshold Amplification threshold
#' @param del_threshold Deletion threshold
#' @return Stringent analysis results
#' @keywords internal
apply_stringent_thresholds <- function(arm_summaries, stringent_threshold, 
                                      amp_threshold, del_threshold) {
  
  stringent_amp <- amp_threshold * stringent_threshold
  stringent_del <- del_threshold * stringent_threshold
  
  arm_summaries %>%
    dplyr::mutate(
      stringent_alteration = case_when(
        mean_log2ratio >= stringent_amp ~ "High-confidence Amplification",
        mean_log2ratio <= stringent_del ~ "High-confidence Deletion",
        mean_log2ratio >= amp_threshold ~ "Moderate Amplification",
        mean_log2ratio <= del_threshold ~ "Moderate Deletion",
        TRUE ~ "Neutral"
      ),
      confidence_level = case_when(
        abs(mean_log2ratio) >= max(abs(stringent_amp), abs(stringent_del)) ~ "High",
        abs(mean_log2ratio) >= max(abs(amp_threshold), abs(del_threshold)) ~ "Moderate",
        TRUE ~ "Low"
      )
    )
}

#' Save Analysis Results to Files
#'
#' @param results Analysis results
#' @param output_dir Output directory
#' @param cancer_type Cancer type
#' @keywords internal
save_analysis_results <- function(results, output_dir, cancer_type) {
  
  # Create cancer-specific subdirectory
  cancer_dir <- file.path(output_dir, paste0("BAGEL_results_", cancer_type))
  if (!dir.exists(cancer_dir)) {
    dir.create(cancer_dir, recursive = TRUE)
  }
  
  # Save arm summaries
  arm_file <- file.path(cancer_dir, "arm_level_summaries.txt")
  readr::write_tsv(results$arm_summaries, arm_file)
  
  # Save stringent results
  stringent_file <- file.path(cancer_dir, "stringent_calls.txt")
  readr::write_tsv(results$stringent_summaries, stringent_file)
  
  # Save GISTIC results if available
  if (!is.null(results$gistic_results)) {
    gistic_file <- file.path(cancer_dir, "gistic_results.txt")
    readr::write_tsv(results$gistic_results$q_values, gistic_file)
    
    sig_file <- file.path(cancer_dir, "significant_arms.txt")
    readr::write_tsv(results$significant_arms, sig_file)
  }
  
  # Save parameters
  param_file <- file.path(cancer_dir, "analysis_parameters.txt")
  param_df <- data.frame(
    parameter = names(results$parameters),
    value = sapply(results$parameters, as.character),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(param_df, param_file)
  
  # Save arm definitions
  arm_def_file <- file.path(cancer_dir, "arm_definitions.txt")
  readr::write_tsv(results$arm_definitions, arm_def_file)
  
  bagel_log(sprintf("Results saved to: %s", cancer_dir), "INFO")
}

#' Generate Analysis Summary Statistics
#'
#' @param results Analysis results
#' @keywords internal
generate_analysis_summary <- function(results) {
  
  n_samples <- length(unique(results$segments$Sample))
  n_arms <- length(unique(results$arm_summaries$Arm))
  
  # Count alterations
  amp_count <- sum(results$arm_summaries$mean_log2ratio >= results$parameters$amp_threshold, na.rm = TRUE)
  del_count <- sum(results$arm_summaries$mean_log2ratio <= results$parameters$del_threshold, na.rm = TRUE)
  
  # Count significant arms if GISTIC was run
  sig_count <- 0
  if (!is.null(results$significant_arms)) {
    sig_count <- nrow(results$significant_arms)
  }
  
  bagel_log("=== ANALYSIS SUMMARY ===", "INFO")
  bagel_log(sprintf("Samples processed: %d", n_samples), "INFO")
  bagel_log(sprintf("Chromosome arms analyzed: %d", n_arms), "INFO")
  bagel_log(sprintf("Amplifications detected: %d", amp_count), "INFO")
  bagel_log(sprintf("Deletions detected: %d", del_count), "INFO")
  
  if (sig_count > 0) {
    bagel_log(sprintf("Statistically significant arms: %d", sig_count), "INFO")
  }
  
  bagel_log("=== END SUMMARY ===", "INFO")
}

#' Complete BAGEL Analysis Workflow
#'
#' High-level wrapper function that performs the complete BAGEL analysis workflow
#' with recommended default settings.
#'
#' @param segments Segmentation data
#' @param cancer_type Cancer type (will use appropriate built-in breakpoints)
#' @param output_dir Output directory (default: "./BAGEL_output")
#' @param ... Additional parameters passed to calculateCopyNumber_fixed
#' @return Complete analysis results
#' @export
bagel_workflow <- function(segments, cancer_type, output_dir = "./BAGEL_output", ...) {
  
  bagel_log(sprintf("Starting complete BAGEL workflow for %s", cancer_type), "INFO")
  
  # Try to use cancer-specific breakpoints, fall back to consensus
  breakpoints <- tryCatch({
    get_breakpoint_data(cancer_type, "arm_definitions")
  }, error = function(e) {
    bagel_log(sprintf("Cancer-specific breakpoints not found for %s, using consensus", cancer_type), "WARN")
    get_breakpoint_data("consensus", "arm_definitions")
  })
  
  # Run complete analysis
  results <- calculateCopyNumber_fixed(
    segments = segments,
    breakpoints = breakpoints,
    cancer_type = cancer_type,
    output_dir = output_dir,
    use_gistic = TRUE,
    save_results = TRUE,
    ...
  )
  
  bagel_log("BAGEL workflow completed successfully", "INFO")
  
  return(results)
}