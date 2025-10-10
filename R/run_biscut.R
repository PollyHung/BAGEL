#' Run BISCUT Breakpoint Detection Pipeline
#'
#' Executes the complete BISCUT analysis pipeline on cancer-type specific
#' segmentation files. This function wraps the breakpoint file creation and
#' BISCUT analysis steps into a single convenient interface.
#'
#' @param cancer_folder Character, name of the cancer folder containing segmentation data
#' @param results_dir Character, path to the directory containing cancer-specific results
#' @param cores Integer, number of CPU cores to use for parallel processing (default: 4)
#' @param skip Logical, whether to skip BISCUT processing when result file exist (default: TRUE)
#' @param cleanup Logical, whether to remove existing output directories before running (default: TRUE)
#'
#' @return List containing:
#'   \item{status}{Character indicating success, error, or no segmentation file}
#'   \item{error}{Character with error message if status is "error"}
#'   \item{result}{BISCUT analysis results if successful}
#'   \item{n_samples}{Number of unique samples in segmentation file}
#'   \item{n_segments}{Number of segments in segmentation file}
#'   \item{breakpoint_created}{Logical indicating if breakpoint files were created}
#'   \item{biscut_completed}{Logical indicating if BISCUT analysis completed}
#'
#' @details
#' This function performs the following steps:
#' 1. Checks for the existence of a segmentation.seg file in the cancer folder
#' 2. Creates breakpoint files using BISCUT::make_breakpoint_files()
#' 3. Runs BISCUT analysis using BISCUT::do_biscut()
#' 4. Returns comprehensive results and status information
#'
#' The function expects the following directory structure:
#' \code{results_dir/cancer_folder/segmentation.seg}
#'
#' @examples
#' \dontrun{
#' # Run BISCUT on a single cancer type
#' result <- run_biscut(
#'   cancer_folder = "breast_invasive_carcinoma",
#'   results_dir = "/path/to/cancer_results",
#'   cores = 8
#' )
#'
#' # Check results
#' if (result$status == "success") {
#'   print(paste("Analysis completed with", result$n_samples, "samples"))
#' } else {
#'   print(paste("Analysis failed:", result$error))
#' }
#' }
#'
#' @import dplyr
#' @importFrom readr read_tsv
#' @export
run_biscut <- function(cancer_folder, 
                       results_dir, 
                       cores = 4, 
                       skip = FALSE, 
                       cleanup = TRUE, 
                       coordinates = NULL, 
                       geneloc = NULL) {
  
  # Check if BISCUT package is available
  if (!requireNamespace("BISCUT", quietly = TRUE)) {
    stop("BISCUT package is required but not available. Please install it first.")
  }
  
  cat("Running BISCUT pipeline for:", cancer_folder, "... ")
  
  # Define paths
  cancer_dir <- file.path(results_dir, cancer_folder)
  seg_file <- file.path(cancer_dir, "segmentation.seg")
  breakpoint_dir <- file.path(cancer_dir, paste0("breakpoint_files_", cancer_folder))
  biscut_results_dir <- file.path(cancer_dir, paste0("results_", cancer_folder))
  
  # Check if target file exists 
  biscut_result <- file.exists(file.path(biscut_results_dir, "all_BISCUT_results.txt"))
  
  # Check if segmentation file exists
  if (!file.exists(seg_file)) {
    cat("SKIP - No segmentation file\n")
    return(list(
      status = "no_seg_file", 
      error = "Segmentation file not found", 
      result = NULL,
      n_samples = 0,
      n_segments = 0,
      breakpoint_created = FALSE,
      biscut_completed = FALSE
    ))
  }
  
  # Read basic info about the segmentation file
  tryCatch({
    seg_data <- readr::read_tsv(seg_file, show_col_types = FALSE, col_names = F)
    n_samples <- length(unique(seg_data$Sample))
    n_segments <- nrow(seg_data)
  }, error = function(e) {
    cat("FAILED - Cannot read segmentation file:", e$message, "\n")
    return(list(
      status = "read_error",
      error = paste("Cannot read segmentation file:", e$message),
      result = NULL,
      n_samples = 0,
      n_segments = 0,
      breakpoint_created = FALSE,
      biscut_completed = FALSE
    ))
  })
  
  result <- list(
    status = "unknown",
    error = "",
    result = NULL,
    n_samples = n_samples,
    n_segments = n_segments,
    breakpoint_created = FALSE,
    biscut_completed = FALSE
  )
  
  # Step 1: Create breakpoint files
  tryCatch({
    # Remove existing directories if cleanup requested
    if (cleanup) {
      if (dir.exists(breakpoint_dir)) {
        unlink(breakpoint_dir, recursive = TRUE)
      }
      if (dir.exists(biscut_results_dir)) {
        unlink(biscut_results_dir, recursive = TRUE)
      }
    }
    
    # Create breakpoint files
    if(is.null(coordinates)){
      cat("Running make_breakpoint_files on hg19 reference", "\n")
      BISCUT::make_breakpoint_files(
        segment_file = seg_file,
        output_dir = breakpoint_dir,
        cores = cores
      )
    } else {
      cat("Running make_breakpoint_files on hg38 reference", "\n")
      BISCUT::make_breakpoint_files(
        segment_file = seg_file,
        output_dir = breakpoint_dir,
        cores = cores, 
        chromosome_coordinates = coordinates
      )
    }
    
    
    result$breakpoint_created <- TRUE
    
    # Step 2: Run BISCUT analysis
    if(is.null(coordinates)){
      cat("Running do_biscut on hg19 reference", "\n")
      biscut_result <- BISCUT::do_biscut(
        breakpoint_file_dir = breakpoint_dir,
        results_dir = biscut_results_dir,
        cores = cores
      )
    } else {
      cat("Running do_biscut on hg38 reference", "\n")
      biscut_result <- BISCUT::do_biscut(
        breakpoint_file_dir = breakpoint_dir,
        results_dir = biscut_results_dir,
        cores = cores,
        chromosome_coordinates = coordinates, 
        gene_locations = geneloc
      )
    }
    
    temp_path <- file.path(biscut_results_dir, "all_BISCUT_results.txt")
    temp <- read.delim(temp_path, sep = "\t")
    temp$total_sample_size <- n_samples
    write.table(temp, temp_path, sep = "\t", quote = F, row.names = F, col.names = T)
    
    result$biscut_completed <- TRUE
    result$status <- "success"
    result$result <- biscut_result
    cat("SUCCESS\n")
    
  }, error = function(e) {
    result$status <<- "error"
    result$error <<- as.character(e$message)
    cat("FAILED -", e$message, "\n")
  })
  
  return(result)
  
}


