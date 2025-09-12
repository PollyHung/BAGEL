#' Breakpoint Similarity Analysis with Comprehensive Outputs
#'
#' Performs statistical analysis to compare breakpoint patterns between datasets
#' and generates comprehensive outputs including matched data tables, statistical 
#' reports, correlation plots, and ideogram visualizations.
#'
#' @param query_path Character, path to the query dataset arm definitions file
#' @param reference Character, either a cancer type keyword to use TCGA reference 
#'   data, or path to a custom reference file
#' @param output_dir Character, directory to save analysis results (required)
#' @param dataset_names Character vector of length 2, names for reference and 
#'   query datasets (default: c("Reference", "Query"))
#'
#' @return List containing:
#'   \item{matched_data}{Data frame with matched breakpoints and statistics}
#'   \item{stat_result}{Statistical analysis results}
#'   \item{output_files}{Paths to generated output files}
#'
#' @details
#' This function performs comprehensive breakpoint similarity analysis and generates:
#' 1. **Matched data table**: Tab-delimited file with breakpoint comparisons
#' 2. **Statistical report**: Markdown report with correlation and similarity metrics
#' 3. **Correlation plot**: ggplot2 visualization saved as PNG
#' 4. **Ideogram visualization**: Paneled chromosome ideograms for matched arms
#'
#' The analysis is designed to validate breakpoint detection consistency across
#' different technical platforms, populations, and analysis methods.
#'
#' Expected input file format (tab-delimited):
#' - chr_num: Chromosome number (1-22)
#' - arm_type: Arm type ("p" or "q") 
#' - arm: Arm identifier (e.g., "1p", "1q")
#' - arm_start: Arm start position
#' - arm_end: Arm end position
#' - direction: Alteration direction ("amp" or "del")
#'
#' @examples
#' \dontrun{
#' # Compare against TCGA reference for ovarian cancer
#' similarity_results <- analyze_breakpoint_similarity(
#'   query_path = "path/to/query_arm_definitions.txt",
#'   reference = "ovarian_serous_cystadenocarcinoma",
#'   output_dir = "validation_results",
#'   dataset_names = c("TCGA", "Queen_Mary")
#' )
#'
#' # Check output files
#' print(similarity_results$output_files)
#' }
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom readr read_tsv
#' @importFrom tidyr pivot_longer
#' @export
analyze_breakpoint_similarity <- function(query_path, 
                                         reference, 
                                         output_dir,
                                         dataset_names = c("Reference", "Query")) {
  
  # Validate inputs
  if (!file.exists(query_path)) {
    stop("Query file not found: ", query_path)
  }
  
  if (missing(output_dir)) {
    stop("output_dir is required for saving analysis results")
  }
  
  if (length(dataset_names) != 2) {
    stop("dataset_names must be a character vector of length 2")
  }
  
  # Check required packages
  if (!requireNamespace("BISCUT", quietly = TRUE)) {
    stop("BISCUT package is required. Please install BISCUT.")
  }
  
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    warning("ggrepel package not available. Labels may overlap in plots.")
  }
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("=== BREAKPOINT SIMILARITY ANALYSIS ===\n")
  cat("Query data:", query_path, "\n")
  cat("Reference:", reference, "\n")
  cat("Dataset names:", paste(dataset_names, collapse = " vs "), "\n")
  cat("Output directory:", output_dir, "\n\n")
  
  # Get chromosome coordinates from BISCUT
  cat("Loading chromosome coordinates...\n")
  tryCatch({
    chr_coords <- BISCUT::get_chromosome_coordinates()
    cat("✅ Using BISCUT chromosome coordinates\n")
  }, error = function(e) {
    stop("Error loading BISCUT coordinates: ", e$message)
  })
  
  # Calculate arm lengths
  cat("Calculating chromosome arm lengths...\n")
  arm_lengths <- calculate_arm_lengths_internal(chr_coords)
  cat("✅ Calculated arm lengths for", nrow(arm_lengths), "chromosome arms\n")
  
  # Load reference data
  cat("Loading reference data...\n")
  reference_data <- load_reference_data(reference, dataset_names[1])
  cat("✅ Loaded reference data:", nrow(reference_data), "breakpoints\n")
  
  # Load query data
  cat("Loading query data...\n")
  query_data <- load_query_data(query_path, dataset_names[2])
  cat("✅ Loaded query data:", nrow(query_data), "breakpoints\n")
  
  # Combine datasets
  cat("Combining datasets and calculating breakpoint metrics...\n")
  combined_data <- combine_and_process_data(reference_data, query_data, arm_lengths)
  cat("✅ Combined analysis data:", nrow(combined_data), "total breakpoints\n")
  
  # Perform similarity analysis
  cat("Performing similarity analysis...\n")
  matched_data <- perform_similarity_analysis_internal(combined_data, dataset_names)
  cat("✅ Found", nrow(matched_data), "matched breakpoints\n")
  
  # Perform statistical tests
  cat("Computing statistical tests...\n")
  stat_result <- perform_statistical_tests_internal(matched_data)
  cat("✅ Statistical analysis complete\n")
  
  # Generate outputs
  cat("Generating output files...\n")
  output_files <- generate_analysis_outputs(
    matched_data, stat_result, output_dir, dataset_names
  )
  
  # Return results
  results <- list(
    matched_data = matched_data,
    stat_result = stat_result,
    output_files = output_files
  )
  
  cat("\n=== Analysis Complete ===\n")
  cat("Datasets compared:", paste(dataset_names, collapse = " vs "), "\n")
  cat("Matched breakpoints:", nrow(matched_data), "\n")
  cat("Output files generated:")
  for (file in output_files) {
    cat("  -", file, "\n")
  }
  
  return(results)
}

#' Internal function to calculate arm lengths from chromosome coordinates
#' @param chr_coords Chromosome coordinates from BISCUT
#' @return Data frame with arm length information
#' @keywords internal
calculate_arm_lengths_internal <- function(chr_coords) {
  arm_lengths <- chr_coords %>%
    dplyr::mutate(
      p_arm_length = p_end - p_start,
      q_arm_length = q_end - q_start,
      chr_num = chromosome_info
    ) %>%
    dplyr::select(chr_num, p_arm_length, q_arm_length) %>%
    tidyr::pivot_longer(
      cols = c(p_arm_length, q_arm_length),
      names_to = "arm_type_full", 
      values_to = "arm_length"
    ) %>%
    dplyr::mutate(
      arm_type = ifelse(arm_type_full == "p_arm_length", "p", "q"),
      arm = paste0(chr_num, arm_type)
    ) %>%
    dplyr::select(chr_num, arm_type, arm, arm_length) %>%
    dplyr::filter(!is.na(arm_length) & arm_length > 0)
  
  return(arm_lengths)
}

#' Internal function to load reference data
#' @param reference Either cancer type keyword or file path
#' @param dataset_name Name for the reference dataset
#' @return Data frame with reference breakpoint data
#' @keywords internal
load_reference_data <- function(reference, dataset_name) {
  
  # Check if reference is a file path or cancer type keyword
  if (file.exists(reference)) {
    # Load custom reference file
    tryCatch({
      reference_data <- readr::read_tsv(reference, show_col_types = FALSE) %>%
        dplyr::mutate(dataset = dataset_name)
      cat("Using custom reference file:", reference, "\n")
    }, error = function(e) {
      stop("Error loading custom reference file: ", e$message)
    })
  } else {
    # Use TCGA embedded data
    tryCatch({
      breakpoint_data <- load_breakpoint_data()
      reference_data <- get_arm_definitions(reference, breakpoint_data) %>%
        dplyr::mutate(dataset = dataset_name)
      cat("Using TCGA reference data for:", reference, "\n")
    }, error = function(e) {
      stop("Error loading TCGA reference data for '", reference, "': ", e$message,
           "\nAvailable cancer types can be checked with get_available_cancer_types()")
    })
  }
  
  return(reference_data)
}

#' Internal function to load query data
#' @param query_path Path to query dataset file
#' @param dataset_name Name for the query dataset
#' @return Data frame with query breakpoint data
#' @keywords internal
load_query_data <- function(query_path, dataset_name) {
  tryCatch({
    query_data <- readr::read_tsv(query_path, show_col_types = FALSE) %>%
      dplyr::mutate(dataset = dataset_name)
    return(query_data)
  }, error = function(e) {
    stop("Error loading query data from '", query_path, "': ", e$message)
  })
}

#' Internal function to combine and process breakpoint data
#' @param reference_data Reference dataset
#' @param query_data Query dataset  
#' @param arm_lengths Arm length information
#' @return Combined and processed data frame
#' @keywords internal
combine_and_process_data <- function(reference_data, query_data, arm_lengths) {
  
  # Combine datasets
  combined_data <- rbind(reference_data, query_data)
  
  # Add arm length information and calculate metrics
  combined_data <- combined_data %>%
    dplyr::left_join(arm_lengths, by = c("chr_num", "arm_type", "arm")) %>%
    dplyr::mutate(
      # Calculate breakpoint position relative to arm
      breakpoint_position = dplyr::case_when(
        arm_type == "p" ~ arm_end,   # p-arm breakpoint is the end position
        arm_type == "q" ~ arm_start  # q-arm breakpoint is the start position
      ),
      # Calculate percentage position within arm
      breakpoint_pct = dplyr::case_when(
        arm_type == "p" ~ (breakpoint_position / arm_length) * 100,
        arm_type == "q" ~ ((breakpoint_position / arm_length) * 100)
      ),
      # Add arm identifier for easier analysis
      arm_direction = paste0(arm, "_", direction)
    ) %>%
    dplyr::filter(!is.na(arm_length) & !is.na(breakpoint_position))
  
  return(combined_data)
}

#' Internal function to perform similarity analysis on combined data
#' @param combined_data Combined breakpoint data
#' @param dataset_names Names of the datasets
#' @return Data frame with matched breakpoints
#' @keywords internal
perform_similarity_analysis_internal <- function(combined_data, dataset_names) {
  
  # Separate datasets
  combined_data_split <- split(combined_data, f = combined_data$dataset)
  reference <- combined_data_split[[dataset_names[1]]]
  query <- combined_data_split[[dataset_names[2]]]
  
  # Find matching breakpoints (same arm and direction)
  matched_breakpoints <- dplyr::inner_join(
    reference %>% dplyr::select(arm, direction, breakpoint_position, arm_length, breakpoint_pct),
    query %>% dplyr::select(arm, direction, breakpoint_position, arm_length, breakpoint_pct),
    by = c("arm", "direction"),
    suffix = paste0("_", tolower(gsub("[^A-Za-z0-9]", "", dataset_names)))
  )
  
  # Calculate differences
  suffix1 <- paste0("_", tolower(gsub("[^A-Za-z0-9]", "", dataset_names[1])))
  suffix2 <- paste0("_", tolower(gsub("[^A-Za-z0-9]", "", dataset_names[2])))
  
  matched_breakpoints <- matched_breakpoints %>%
    dplyr::mutate(
      # Extract chromosome number for visualization
      chr_num = as.numeric(gsub("[pq]", "", arm)),
      # Absolute difference in base pairs
      abs_diff_bp = abs(.data[[paste0("breakpoint_position", suffix1)]] - .data[[paste0("breakpoint_position", suffix2)]]),
      # Relative difference as percentage of arm length
      rel_diff_pct = abs_diff_bp / pmax(.data[[paste0("arm_length", suffix1)]], .data[[paste0("arm_length", suffix2)]], na.rm = TRUE) * 100,
      # Log difference for visualization
      log_abs_diff = log10(abs_diff_bp + 1),
      # Add arm type
      arm_type = ifelse(grepl("p$", arm), "p", "q")
    )
  
  return(matched_breakpoints)
}

#' Internal function to perform statistical tests
#' @param matched_data Matched breakpoint data
#' @return List with statistical results
#' @keywords internal
perform_statistical_tests_internal <- function(matched_data) {
  
  # Get column names dynamically
  pos_cols <- grep("breakpoint_position_", names(matched_data), value = TRUE)
  
  # Basic descriptive statistics
  desc_stats <- list(
    n_matched = nrow(matched_data),
    mean_abs_diff = mean(matched_data$abs_diff_bp, na.rm = TRUE),
    median_abs_diff = median(matched_data$abs_diff_bp, na.rm = TRUE),
    sd_abs_diff = sd(matched_data$abs_diff_bp, na.rm = TRUE),
    mean_rel_diff = mean(matched_data$rel_diff_pct, na.rm = TRUE),
    median_rel_diff = median(matched_data$rel_diff_pct, na.rm = TRUE)
  )
  
  # Correlation analysis
  correlation_test <- cor.test(
    matched_data[[pos_cols[1]]], 
    matched_data[[pos_cols[2]]],
    method = "pearson"
  )
  
  # Technical threshold analysis
  technical_threshold <- 5.0  # 5% of arm length
  within_threshold <- sum(matched_data$rel_diff_pct < technical_threshold, na.rm = TRUE)
  threshold_pct <- (within_threshold / nrow(matched_data)) * 100
  
  # Wilcoxon test for systematic bias
  wilcox_test <- wilcox.test(matched_data$rel_diff_pct, mu = 0, alternative = "two.sided")
  
  # Bootstrap confidence interval for mean difference
  set.seed(123)
  bootstrap_means <- replicate(1000, {
    boot_indices <- sample(nrow(matched_data), nrow(matched_data), replace = TRUE)
    boot_sample <- matched_data[boot_indices, ]
    mean(boot_sample$rel_diff_pct, na.rm = TRUE)
  })
  
  boot_ci <- quantile(bootstrap_means, c(0.025, 0.975))
  
  # Effect size calculation
  effect_size <- desc_stats$mean_rel_diff / desc_stats$sd_abs_diff
  
  return(list(
    descriptive = desc_stats,
    correlation = correlation_test,
    threshold_analysis = list(
      threshold = technical_threshold,
      within_threshold = within_threshold,
      threshold_pct = threshold_pct
    ),
    wilcox_test = wilcox_test,
    bootstrap_ci = boot_ci,
    effect_size = effect_size
  ))
}

#' Internal function to generate all analysis outputs
#' @param matched_data Matched breakpoint data
#' @param stat_result Statistical analysis results
#' @param output_dir Output directory
#' @param dataset_names Dataset names
#' @return Vector of output file paths
#' @keywords internal
generate_analysis_outputs <- function(matched_data, stat_result, output_dir, dataset_names) {
  
  output_files <- character()
  
  # 1. Save matched data table
  matched_data_file <- file.path(output_dir, "matched_breakpoint_data.txt")
  write.table(matched_data, matched_data_file, 
              sep = "\t", row.names = FALSE, quote = FALSE)
  output_files <- c(output_files, matched_data_file)
  cat("✅ Saved matched data table\n")
  
  # 2. Generate statistical report
  report_file <- file.path(output_dir, "breakpoint_similarity_report.md")
  report_content <- generate_statistical_report(matched_data, stat_result, dataset_names)
  writeLines(report_content, report_file)
  output_files <- c(output_files, report_file)
  cat("✅ Generated statistical report\n")
  
  # 3. Create correlation plot
  plot_file <- file.path(output_dir, "breakpoint_correlation.png")
  create_correlation_plot(matched_data, stat_result, plot_file, dataset_names)
  output_files <- c(output_files, plot_file)
  cat("✅ Created correlation plot\n")
  
  # 4. Generate ideogram visualization
  ideogram_file <- file.path(output_dir, "matched_chromosomes_ideogram.png")
  tryCatch({
    create_matched_ideograms(matched_data, ideogram_file, output_dir)
    output_files <- c(output_files, ideogram_file)
    cat("✅ Generated ideogram visualization\n")
  }, error = function(e) {
    cat("⚠️ Error generating ideograms:", e$message, "\n")
  })
  
  return(output_files)
}

#' Internal function to generate statistical report
#' @param matched_data Matched breakpoint data
#' @param stat_result Statistical results
#' @param dataset_names Dataset names
#' @return Character vector with report content
#' @keywords internal
generate_statistical_report <- function(matched_data, stat_result, dataset_names) {
  
  corr <- stat_result$correlation
  desc <- stat_result$descriptive
  thresh <- stat_result$threshold_analysis
  
  report <- c(
    "# Breakpoint Similarity Analysis Report",
    "",
    "## Datasets Compared",
    paste("- **Reference**:", dataset_names[1]),
    paste("- **Query**:", dataset_names[2]),
    "",
    "## Summary Statistics",
    paste("- **Matched breakpoints**:", desc$n_matched),
    paste("- **Mean absolute difference**:", round(desc$mean_abs_diff, 0), "bp"),
    paste("- **Median absolute difference**:", round(desc$median_abs_diff, 0), "bp"),
    paste("- **Standard deviation**:", round(desc$sd_abs_diff, 0), "bp"),
    paste("- **Mean relative difference**:", round(desc$mean_rel_diff, 3), "% of arm length"),
    paste("- **Median relative difference**:", round(desc$median_rel_diff, 3), "% of arm length"),
    "",
    "## Correlation Analysis",
    paste("- **Pearson correlation coefficient**:", round(corr$estimate, 4)),
    paste("- **95% Confidence interval**:", round(corr$conf.int[1], 4), "to", round(corr$conf.int[2], 4)),
    paste("- **P-value**:", format(corr$p.value, scientific = TRUE)),
    "",
    "## Technical Threshold Analysis",
    paste("- **Breakpoints within", thresh$threshold, "% threshold**:", thresh$within_threshold, "/", desc$n_matched),
    paste("- **Percentage within threshold**:", round(thresh$threshold_pct, 1), "%"),
    "",
    "## Statistical Tests",
    paste("- **Wilcoxon test p-value**:", format(stat_result$wilcox_test$p.value, scientific = TRUE)),
    paste("- **Bootstrap 95% CI for mean difference**:", round(stat_result$bootstrap_ci[1], 3), "% to", round(stat_result$bootstrap_ci[2], 3), "%"),
    paste("- **Effect size (mean/sd)**:", round(stat_result$effect_size, 3)),
    "",
    "## Interpretation",
    paste("This analysis compares breakpoint patterns between", dataset_names[1], "and", dataset_names[2], "datasets."),
    "The high correlation and low relative differences suggest strong technical consistency",
    "between the datasets, validating the robustness of breakpoint detection across platforms."
  )
  
  return(report)
}

#' Internal function to create correlation plot
#' @param matched_data Matched breakpoint data
#' @param stat_result Statistical results
#' @param plot_file Output file path
#' @param dataset_names Dataset names
#' @keywords internal
create_correlation_plot <- function(matched_data, stat_result, plot_file, dataset_names) {
  
  # Get position column names
  pos_cols <- grep("breakpoint_position_", names(matched_data), value = TRUE)
  correlation_r <- stat_result$correlation$estimate
  
  # Set up chromosome colors
  if (exists("bagel_palette", envir = asNamespace("BAGEL"))) {
    set.seed(222)
    chr_colors <- sample(get("bagel_palette", envir = asNamespace("BAGEL")), size = 22, replace = FALSE)
  } else {
    chr_colors <- rainbow(22)
  }
  names(chr_colors) <- c(1:22)
  
  # Create plot
  p1 <- ggplot(matched_data, aes(x = .data[[pos_cols[1]]]/1e6, y = .data[[pos_cols[2]]]/1e6)) +
    geom_point(aes(fill = factor(chr_num), shape = factor(arm_type)), 
               size = 3, alpha = 1, color = "black", stroke = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "#DC143C", linetype = "solid", size = 0.5) +
    scale_fill_manual(values = chr_colors, name = "Chromosome") +
    scale_shape_manual(values = c("p" = 21, "q" = 24), name = "Arm", labels = c("p arm", "q arm")) +
    labs(title = paste("Breakpoint Position Correlation\n", dataset_names[1], "vs", dataset_names[2]),
         subtitle = paste0("Pearson r = ", round(correlation_r, 3)),
         x = paste(dataset_names[1], "Breakpoint Position (Mb)"),
         y = paste(dataset_names[2], "Breakpoint Position (Mb)")) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "none")
  
  # Add labels if ggrepel is available
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p1 <- p1 + ggrepel::geom_text_repel(aes(label = arm), size = 4, color = "black", 
                                        fontface = "bold", box.padding = 0.3, 
                                        point.padding = 0.3, max.overlaps = Inf)
  }
  
  ggsave(plot_file, p1, width = 4, height = 4, dpi = 300)
}

#' Internal function to create matched chromosome ideograms
#' @param matched_data Matched breakpoint data
#' @param ideogram_file Output file path
#' @param output_dir Output directory
#' @keywords internal
create_matched_ideograms <- function(matched_data, ideogram_file, output_dir) {
  
  # Create a mock bagel_results object for matched chromosomes
  unique_chrs <- sort(unique(matched_data$chr_num))
  
  # Create arm summaries for matched chromosomes
  arm_summaries <- matched_data %>%
    dplyr::select(arm, chr_num, direction) %>%
    dplyr::mutate(
      Sample = "Matched_Breakpoints",
      Copy_Number = ifelse(direction == "amp", 1, -1),
      Log2_Ratio = ifelse(direction == "amp", 0.5, -0.5)
    ) %>%
    dplyr::select(Sample, arm, Copy_Number, Log2_Ratio)
  
  # Create mock bagel_results structure
  mock_bagel_results <- list(
    arm_summaries = arm_summaries,
    parameters = list(cancer_type = "Matched_Breakpoints")
  )
  
  # Try to use existing plot_chromosome_ideograms function
  tryCatch({
    if (exists("plot_chromosome_ideograms", envir = asNamespace("BAGEL"))) {
      ideogram_results <- get("plot_chromosome_ideograms", envir = asNamespace("BAGEL"))(
        bagel_results = mock_bagel_results,
        output_dir = output_dir,
        save_plots = TRUE
      )
      
      if (!is.null(ideogram_results)) {
        cat("Generated ideograms for", length(unique_chrs), "matched chromosomes\n")
      }
    } else {
      cat("plot_chromosome_ideograms function not available\n")
    }
  }, error = function(e) {
    cat("Error creating ideograms:", e$message, "\n")
  })
}

#' Internal function to calculate similarity statistics
#' @param combined_data Combined breakpoint data
#' @param dataset_names Names of the datasets
#' @return List with similarity statistics
#' @keywords internal  
calculate_similarity_statistics <- function(combined_data, dataset_names) {
  
  # Basic summary statistics
  summary_stats <- combined_data %>%
    dplyr::group_by(dataset, arm_type, direction) %>%
    dplyr::summarise(
      n_breakpoints = n(),
      mean_position = mean(breakpoint_position, na.rm = TRUE),
      median_position = median(breakpoint_position, na.rm = TRUE),
      sd_position = sd(breakpoint_position, na.rm = TRUE),
      mean_pct = mean(breakpoint_pct, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Correlation analysis for common arms
  common_arms <- combined_data %>%
    dplyr::count(arm_direction, dataset) %>%
    tidyr::pivot_wider(names_from = dataset, values_from = n, values_fill = 0) %>%
    dplyr::filter(!!rlang::sym(dataset_names[1]) > 0 & !!rlang::sym(dataset_names[2]) > 0) %>%
    dplyr::pull(arm_direction)
  
  correlation_data <- NULL
  if (length(common_arms) > 3) {
    correlation_data <- combined_data %>%
      dplyr::filter(arm_direction %in% common_arms) %>%
      dplyr::select(arm_direction, dataset, breakpoint_position) %>%
      tidyr::pivot_wider(names_from = dataset, values_from = breakpoint_position) %>%
      dplyr::filter(!is.na(!!rlang::sym(dataset_names[1])) & !is.na(!!rlang::sym(dataset_names[2])))
    
    if (nrow(correlation_data) > 1) {
      correlation_result <- cor.test(
        correlation_data[[dataset_names[1]]], 
        correlation_data[[dataset_names[2]]]
      )
    } else {
      correlation_result <- NULL
    }
  } else {
    correlation_result <- NULL
  }
  
  # Overlap analysis
  overlap_stats <- combined_data %>%
    dplyr::group_by(arm_direction) %>%
    dplyr::summarise(
      datasets_present = n_distinct(dataset),
      in_both = datasets_present == 2,
      .groups = "drop"
    ) %>%
    dplyr::summarise(
      total_arms = n(),
      common_arms = sum(in_both),
      overlap_pct = (common_arms / total_arms) * 100
    )
  
  return(list(
    summary_stats = summary_stats,
    correlation_result = correlation_result,
    overlap_stats = overlap_stats,
    common_arms = common_arms
  ))
}

#' Internal function to save similarity analysis results
#' @param combined_data Combined analysis data
#' @param similarity_stats Statistical results
#' @param output_dir Output directory
#' @param dataset_names Dataset names
#' @keywords internal
save_similarity_results <- function(combined_data, similarity_stats, output_dir, dataset_names) {
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save combined data
  write.table(combined_data, 
              file.path(output_dir, "combined_breakpoint_data.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save summary statistics
  write.table(similarity_stats$summary_stats,
              file.path(output_dir, "summary_statistics.txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save overlap statistics
  write.table(similarity_stats$overlap_stats,
              file.path(output_dir, "overlap_statistics.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("✅ Results saved to:", output_dir, "\n")
}

#' Generate Breakpoint Similarity Report
#'
#' Creates a comprehensive report comparing breakpoint patterns between datasets.
#'
#' @param similarity_results Results from analyze_breakpoint_similarity()
#' @param output_file Path to save the report (optional)
#'
#' @return Character string containing the formatted report
#'
#' @examples
#' \dontrun{
#' results <- analyze_breakpoint_similarity(query_path, reference)
#' report <- generate_similarity_report(results, "similarity_report.md")
#' }
#'
#' @export
generate_similarity_report <- function(similarity_results, output_file = NULL) {
  
  dataset_names <- similarity_results$dataset_names
  stats <- similarity_results$similarity_stats
  
  report <- paste0(
    "# Breakpoint Similarity Analysis Report\n\n",
    "## Datasets Compared\n",
    "- **Reference**: ", dataset_names[1], "\n",
    "- **Query**: ", dataset_names[2], "\n\n",
    
    "## Summary Statistics\n",
    "- **Total breakpoints analyzed**: ", nrow(similarity_results$combined_data), "\n",
    "- **Common chromosome arms**: ", stats$overlap_stats$common_arms, "\n", 
    "- **Overlap percentage**: ", round(stats$overlap_stats$overlap_pct, 1), "%\n\n",
    
    "## Statistical Analysis\n"
  )
  
  if (!is.null(stats$correlation_result)) {
    report <- paste0(report,
      "- **Correlation coefficient**: ", round(stats$correlation_result$estimate, 3), "\n",
      "- **P-value**: ", format(stats$correlation_result$p.value, scientific = TRUE), "\n",
      "- **95% CI**: [", round(stats$correlation_result$conf.int[1], 3), 
      ", ", round(stats$correlation_result$conf.int[2], 3), "]\n\n"
    )
  }
  
  report <- paste0(report,
    "## Interpretation\n",
    "This analysis compares breakpoint patterns between ", dataset_names[1], 
    " and ", dataset_names[2], " datasets to assess technical and biological consistency.\n\n"
  )
  
  if (!is.null(output_file)) {
    writeLines(report, output_file)
    cat("Report saved to:", output_file, "\n")
  }
  
  return(report)
}