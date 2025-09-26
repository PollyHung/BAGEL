#' Create Enhanced Chromosome Arm Copy Number Matrices
#'
#' Creates comprehensive arm × sample matrices from arm-level copy number summaries,
#' following GISTIC2 conventions and generating multiple output formats for downstream
#' analysis. This function follows the BAGEL coding standard with maximum 2-level
#' function nesting for maintainability and comprehensive self-contained design.
#'
#' @section Matrix Output Formats:
#' The function generates several matrix formats to support different analysis needs:
#' \itemize{
#'   \item **Log2 Ratio Matrix**: Continuous values centered at 0 (GISTIC2 style)
#'   \item **Discrete Calls Matrix**: Integer calls (-1, 0, +1) for discrete analysis
#'   \item **Summary Statistics**: Per-arm statistical summaries with alteration frequencies
#'   \item **Long Format Data**: Tidy format for ggplot2 and statistical modeling
#' }
#'
#' @section GISTIC2 Compliance:
#' All outputs follow GISTIC2 conventions:
#' \itemize{
#'   \item Log2 ratios used directly (not absolute copy numbers)
#'   \item Standard thresholds: amp ≥ 0.25, del ≤ -0.25
#'   \item Discrete calls: -1 (deletion), 0 (neutral), +1 (amplification)
#'   \item Alteration frequencies calculated per arm across samples
#' }
#'
#' @param arm_summaries Data frame with arm-level summaries containing columns:
#'   Sample, Arm, mean_log2ratio (at minimum). Additional columns are preserved.
#' @param output_dir Character or NULL. Directory path for saving matrix files.
#'   If NULL, matrices are returned without file output (default: NULL)
#' @param amp_threshold Numeric. Amplification threshold in log2 scale (default: 0.25)
#' @param del_threshold Numeric. Deletion threshold in log2 scale (default: -0.25)
#' @param save_csv Logical. Whether to save CSV files when output_dir provided (default: TRUE)
#' @param include_metadata Logical. Whether to include analysis metadata in output (default: TRUE)
#'
#' @return List containing comprehensive matrix outputs:
#' \describe{
#'   \item{log2_matrix}{Matrix of raw log2 ratios (arms × samples)}
#'   \item{cn_matrix}{Copy number matrix (same as log2_matrix, GISTIC2 style)}
#'   \item{calls_matrix}{Discrete calls matrix (-1, 0, +1)}
#'   \item{summary_stats}{Data frame with per-arm statistics and frequencies}
#'   \item{long_format}{Long format data with calls and labels}
#'   \item{gistic_thresholds}{List with threshold parameters used}
#'   \item{metadata}{Analysis metadata (if include_metadata=TRUE)}
#' }
#'
#' @section File Outputs (when output_dir provided):
#' \describe{
#'   \item{arm_log2ratio_matrix.csv}{Raw log2 ratio values}
#'   \item{arm_copynumber_matrix.csv}{GISTIC2-style copy number matrix}
#'   \item{arm_copynumber_summary.csv}{Statistical summaries per arm}
#'   \item{arm_calls_matrix_gistic_style.csv}{Discrete calls matrix}
#'   \item{arm_copynumber_long.csv}{Long format with calls and labels}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage - return matrices only
#' matrices <- create_matrix(arm_summaries)
#'
#' # Save to files with custom thresholds
#' matrices <- create_matrix(arm_summaries,
#'                          output_dir = "results/matrices",
#'                          amp_threshold = 0.3,
#'                          del_threshold = -0.3)
#'
#' # Access different formats
#' log2_values <- matrices$log2_matrix
#' discrete_calls <- matrices$calls_matrix
#' summary_table <- matrices$summary_stats
#' }
#'
#' @seealso \code{\link{calculate_copynumber}}, \code{\link{define_arm}}
#'
#' @import dplyr
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stats median sd
#' @export
create_matrix <- function(arm_summaries,
                         output_dir = NULL,
                         amp_threshold = 0.25,      # GISTIC2 standard
                         del_threshold = -0.25,     # GISTIC2 standard
                         save_csv = TRUE,
                         include_metadata = TRUE) {

  cat("=== Creating Enhanced Arm Copy Number Matrices ===\n")
  cat("Input data:", nrow(arm_summaries), "arm-sample combinations\n")

  # ============================================================================
  # LEVEL 1 NESTED FUNCTIONS: Core matrix creation components
  # ============================================================================

  # Level 1: Validate and prepare input data
  prepare_matrix_data <- function(arm_data, params) {

    # Level 2: Validate input format and clean duplicates
    validate_and_clean_data <- function(data) {
      # Check required columns
      required_cols <- c("Sample", "Arm", "mean_log2ratio")
      missing_cols <- setdiff(required_cols, names(data))
      if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
      }

      # Handle duplicates by taking mean for each Arm-Sample combination
      clean_data <- data %>%
        dplyr::select(Sample, Arm, mean_log2ratio) %>%
        dplyr::group_by(Sample, Arm) %>%
        dplyr::summarise(mean_log2ratio = mean(mean_log2ratio, na.rm = TRUE), .groups = "drop") %>%
        dplyr::filter(!is.na(mean_log2ratio))

      # Check for data after cleaning
      if (nrow(clean_data) == 0) {
        stop("No valid data remaining after cleaning")
      }

      cat("Data validation: ", nrow(data), "->", nrow(clean_data), "valid combinations\n")
      return(clean_data)
    }

    # Level 2: Create base arm × sample matrix
    create_base_matrix <- function(clean_data) {
      # Create arm-by-sample matrix
      arm_matrix_df <- clean_data %>%
        tidyr::pivot_wider(names_from = Sample, values_from = mean_log2ratio,
                          names_sort = TRUE, values_fill = NA)

      # Convert to matrix
      arm_names <- arm_matrix_df$Arm
      arm_matrix <- as.matrix(arm_matrix_df[, -1])

      # Ensure matrix is numeric
      if (!is.numeric(arm_matrix)) {
        arm_matrix <- apply(arm_matrix, 2, as.numeric)
      }

      rownames(arm_matrix) <- arm_names

      # Validate matrix creation
      if (any(is.infinite(arm_matrix), na.rm = TRUE)) {
        warning("Infinite values detected in matrix")
      }

      cat("Matrix created:", nrow(arm_matrix), "arms ×", ncol(arm_matrix), "samples\n")
      return(arm_matrix)
    }

    # Main Level 1 logic
    validated_data <- validate_and_clean_data(arm_data)
    base_matrix <- create_base_matrix(validated_data)

    # Validate threshold parameters
    if (params$amp_threshold <= 0 || params$del_threshold >= 0) {
      warning("Non-standard thresholds: amp should be positive, del should be negative")
    }

    return(list(
      clean_data = validated_data,
      log2_matrix = base_matrix,
      n_arms = nrow(base_matrix),
      n_samples = ncol(base_matrix)
    ))
  }

  # Level 1: Generate comprehensive matrix formats
  generate_matrix_formats <- function(matrix_data, params) {

    # Level 2: Create GISTIC2-style copy number matrix
    create_gistic_matrix <- function(log2_matrix) {
      # In GISTIC2 style, copy number matrix is the same as log2 matrix
      # (centered at 0, not absolute copy numbers)
      cn_matrix <- log2_matrix
      return(cn_matrix)
    }

    # Level 2: Generate discrete calls matrix
    create_calls_matrix <- function(cn_matrix, amp_thresh, del_thresh) {
      calls_matrix <- matrix(0, nrow = nrow(cn_matrix), ncol = ncol(cn_matrix))
      rownames(calls_matrix) <- rownames(cn_matrix)
      colnames(calls_matrix) <- colnames(cn_matrix)

      # Apply thresholds
      calls_matrix[cn_matrix >= amp_thresh] <- 1   # Amplification
      calls_matrix[cn_matrix <= del_thresh] <- -1  # Deletion
      # Neutral remains 0

      return(calls_matrix)
    }

    # Main Level 1 logic
    cn_matrix <- create_gistic_matrix(matrix_data$log2_matrix)
    calls_matrix <- create_calls_matrix(cn_matrix, params$amp_threshold, params$del_threshold)

    cat("Matrix formats generated: log2, copy number, and discrete calls\n")

    return(list(
      log2_matrix = matrix_data$log2_matrix,
      cn_matrix = cn_matrix,
      calls_matrix = calls_matrix
    ))
  }

  # Level 1: Calculate comprehensive summary statistics
  calculate_matrix_statistics <- function(matrices, params) {

    # Level 2: Generate GISTIC2-style arm statistics
    compute_arm_statistics <- function(cn_matrix, calls_matrix, thresholds) {
      arm_stats <- data.frame(
        Arm = rownames(cn_matrix),
        Mean_Log2Ratio = round(rowMeans(cn_matrix, na.rm = TRUE), 4),
        Median_Log2Ratio = round(apply(cn_matrix, 1, median, na.rm = TRUE), 4),
        Min_Log2Ratio = round(apply(cn_matrix, 1, min, na.rm = TRUE), 4),
        Max_Log2Ratio = round(apply(cn_matrix, 1, max, na.rm = TRUE), 4),
        SD_Log2Ratio = round(apply(cn_matrix, 1, sd, na.rm = TRUE), 4),
        stringsAsFactors = FALSE
      )

      return(arm_stats)
    }

    # Level 2: Calculate alteration frequencies
    compute_alteration_frequencies <- function(calls_matrix, cn_matrix, thresholds) {
      n_samples <- ncol(cn_matrix)

      freq_stats <- data.frame(
        Arm = rownames(cn_matrix),
        N_Samples = n_samples,
        N_Amplified = rowSums(calls_matrix == 1, na.rm = TRUE),
        N_Deleted = rowSums(calls_matrix == -1, na.rm = TRUE),
        N_Neutral = rowSums(calls_matrix == 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::mutate(
          Pct_Amplified = round(N_Amplified / N_Samples * 100, 1),
          Pct_Deleted = round(N_Deleted / N_Samples * 100, 1),
          Pct_Altered = round((N_Amplified + N_Deleted) / N_Samples * 100, 1),
          Pct_Neutral = round(N_Neutral / N_Samples * 100, 1)
        )

      return(freq_stats)
    }

    # Main Level 1 logic
    basic_stats <- compute_arm_statistics(matrices$cn_matrix, matrices$calls_matrix, params)
    frequency_stats <- compute_alteration_frequencies(matrices$calls_matrix, matrices$cn_matrix, params)

    # Combine statistics
    combined_stats <- basic_stats %>%
      dplyr::left_join(frequency_stats, by = "Arm") %>%
      dplyr::arrange(desc(Pct_Altered), Arm)

    cat("Statistics calculated:", nrow(combined_stats), "arms with frequency data\n")

    return(combined_stats)
  }

  # Level 1: Generate long format data
  create_long_format <- function(matrices, summary_stats, params) {

    # Level 2: Create continuous data in long format
    create_continuous_long <- function(cn_matrix) {
      long_continuous <- as.data.frame(cn_matrix) %>%
        dplyr::mutate(Arm = rownames(.)) %>%
        tidyr::pivot_longer(-Arm, names_to = "Sample", values_to = "Log2Ratio") %>%
        dplyr::arrange(Arm, Sample)

      return(long_continuous)
    }

    # Level 2: Create discrete calls in long format
    create_calls_long <- function(calls_matrix) {
      calls_long <- as.data.frame(calls_matrix) %>%
        dplyr::mutate(Arm = rownames(.)) %>%
        tidyr::pivot_longer(-Arm, names_to = "Sample", values_to = "Call") %>%
        dplyr::arrange(Arm, Sample) %>%
        dplyr::mutate(
          Call_Label = case_when(
            Call == 1 ~ "Amplification",
            Call == -1 ~ "Deletion",
            Call == 0 ~ "Normal",
            TRUE ~ "Unknown"
          )
        )

      return(calls_long)
    }

    # Main Level 1 logic
    continuous_long <- create_continuous_long(matrices$cn_matrix)
    calls_long <- create_calls_long(matrices$calls_matrix)

    # Combine continuous and discrete data
    combined_long <- continuous_long %>%
      dplyr::left_join(calls_long %>% dplyr::select(Arm, Sample, Call, Call_Label),
                      by = c("Arm", "Sample")) %>%
      # Add arm-level statistics
      dplyr::left_join(summary_stats %>% dplyr::select(Arm, Pct_Altered, Pct_Amplified, Pct_Deleted),
                      by = "Arm") %>%
      dplyr::arrange(Arm, Sample)

    cat("Long format created:", nrow(combined_long), "arm-sample observations\n")

    return(combined_long)
  }

  # Level 1: Save matrices to files
  save_matrix_files <- function(all_matrices, output_path, params) {

    # Level 2: Create output directory and save CSV files
    write_matrix_files <- function(matrices, stats, long_data, dir_path) {
      # Create directory
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

      # Save all matrix formats
      write.csv(matrices$log2_matrix,
               file.path(dir_path, "arm_log2ratio_matrix.csv"))
      write.csv(matrices$cn_matrix,
               file.path(dir_path, "arm_copynumber_matrix.csv"))
      write.csv(matrices$calls_matrix,
               file.path(dir_path, "arm_calls_matrix_gistic_style.csv"))

      # Save summary and long format
      write.csv(stats,
               file.path(dir_path, "arm_copynumber_summary.csv"), row.names = FALSE)
      write.csv(long_data,
               file.path(dir_path, "arm_copynumber_long.csv"), row.names = FALSE)

      return(dir_path)
    }

    # Main Level 1 logic
    if (!is.null(output_path) && params$save_csv) {
      saved_path <- write_matrix_files(
        all_matrices$matrices,
        all_matrices$summary_stats,
        all_matrices$long_format,
        output_path
      )
      cat("Matrix files saved to:", saved_path, "\n")
      return(saved_path)
    } else {
      cat("File output skipped (output_dir=NULL or save_csv=FALSE)\n")
      return(NULL)
    }
  }

  # ============================================================================
  # MAIN MATRIX CREATION WORKFLOW
  # ============================================================================

  # Prepare analysis parameters
  analysis_params <- list(
    amp_threshold = amp_threshold,
    del_threshold = del_threshold,
    save_csv = save_csv,
    include_metadata = include_metadata,
    creation_time = Sys.time()
  )

  # Step 1: Prepare and validate matrix data
  cat("=== Step 1: Preparing Matrix Data ===\n")
  matrix_data <- prepare_matrix_data(arm_summaries, analysis_params)

  # Step 2: Generate all matrix formats
  cat("=== Step 2: Generating Matrix Formats ===\n")
  matrices <- generate_matrix_formats(matrix_data, analysis_params)

  # Step 3: Calculate comprehensive statistics
  cat("=== Step 3: Calculating Statistics ===\n")
  summary_stats <- calculate_matrix_statistics(matrices, analysis_params)

  # Step 4: Create long format data
  cat("=== Step 4: Creating Long Format Data ===\n")
  long_format <- create_long_format(matrices, summary_stats, analysis_params)

  # Step 5: Save files if requested
  cat("=== Step 5: Saving Matrix Files ===\n")
  output_path <- save_matrix_files(
    list(matrices = matrices, summary_stats = summary_stats, long_format = long_format),
    output_dir,
    analysis_params
  )

  # Step 6: Compile final results
  final_results <- list(
    log2_matrix = matrices$log2_matrix,
    cn_matrix = matrices$cn_matrix,
    calls_matrix = matrices$calls_matrix,
    summary_stats = summary_stats,
    long_format = long_format,
    gistic_thresholds = list(
      amp_threshold = amp_threshold,
      del_threshold = del_threshold
    )
  )

  # Add metadata if requested
  if (include_metadata) {
    final_results$metadata <- list(
      n_arms = matrix_data$n_arms,
      n_samples = matrix_data$n_samples,
      n_observations = nrow(long_format),
      thresholds_used = analysis_params$gistic_thresholds,
      creation_time = analysis_params$creation_time,
      output_path = output_path
    )
  }

  # Final summary
  cat("=== Matrix Creation Complete ===\n")
  cat("Matrices:", matrix_data$n_arms, "arms ×", matrix_data$n_samples, "samples\n")
  cat("Formats created: log2, copy number, discrete calls, summary, long\n")

  return(final_results)
}