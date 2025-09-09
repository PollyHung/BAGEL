#' Process BISCUT Results to Create Custom Arm Definitions
#'
#' Processes BISCUT analysis results file to create chromosome arm definitions
#' compatible with BAGEL analysis pipeline. Uses the same logic as the original
#' TCGA pan-cancer breakpoint creation script for consistency.
#'
#' @param biscut_results_file Character, path to BISCUT all_BISCUT_results.txt file
#' @return Data frame with arm definitions compatible with BAGEL analysis
#' @import dplyr
#' @importFrom readr read_delim
#' @export
process_biscut_results <- function(biscut_results_file) {
  
  if (!file.exists(biscut_results_file)) {
    stop("BISCUT results file not found: ", biscut_results_file)
  }
  
  # Load BISCUT results using same method as original script
  tryCatch({
    biscut_data <- readr::read_delim(biscut_results_file, delim = "\t", show_col_types = FALSE)
    cat("Loaded BISCUT results:", nrow(biscut_data), "entries\n")
  }, error = function(e) {
    stop("Error reading BISCUT results file: ", e$message)
  })
  
  # Create arm-level breakpoint summary (following original script logic lines 47-55)
  arm_breakpoints <- biscut_data %>%
    group_by(arm, direction) %>%
    summarise(
      chr = unique(Chr)[1],
      n_peaks = n(),
      peak_start = min(Peak.Start),
      peak_end = max(Peak.End),
      .groups = "drop"
    )
  
  cat("Processed", nrow(arm_breakpoints), "arm-level breakpoints\n")
  
  # Get BISCUT chromosome coordinates for consistent boundaries
  tryCatch({
    biscut_coords <- BISCUT::get_chromosome_coordinates()
    cat("Using BISCUT coordinates for chromosome arm boundaries\n")
  }, error = function(e) {
    warning("Could not load BISCUT coordinates: ", e$message)
    biscut_coords <- NULL
  })
  
  # Create standardized arm definitions (following original script logic lines 74-114)
  arm_definitions <- arm_breakpoints %>%
    mutate(
      chr_num = as.numeric(chr),
      arm_type = ifelse(grepl("p$", arm), "p", "q")
    )
  
  if (!is.null(biscut_coords)) {
    # Use BISCUT coordinates for arm boundaries
    arm_definitions <- arm_definitions %>%
      left_join(
        biscut_coords %>%
          select(chromosome_info, p_start, p_end, q_start, q_end) %>%
          rename(chr_num = chromosome_info),
        by = "chr_num"
      ) %>%
      mutate(
        arm_start = case_when(
          arm_type == "p" ~ p_start,            # Use BISCUT p arm start
          arm_type == "q" ~ peak_start          # q arm starts at functional breakpoint
        ),
        arm_end = case_when(
          arm_type == "p" ~ peak_end,           # p arm ends at functional breakpoint
          arm_type == "q" ~ q_end               # Use BISCUT q arm end
        )
      ) %>%
      select(-p_start, -p_end, -q_start, -q_end)
  } else {
    # Fallback to hardcoded coordinates if BISCUT not available
    arm_definitions <- arm_definitions %>%
      mutate(
        arm_start = case_when(
          arm_type == "p" ~ 1,                  # p arm always starts at position 1
          arm_type == "q" ~ peak_start          # q arm starts at functional breakpoint
        ),
        arm_end = case_when(
          arm_type == "p" ~ peak_end,           # p arm ends at functional breakpoint
          arm_type == "q" ~ case_when(          # q arm ends at chromosome end (hardcoded)
            chr_num == 1 ~ 249250621,
            chr_num == 2 ~ 243199373,
            chr_num == 3 ~ 198022430,
            chr_num == 4 ~ 191154276,
            chr_num == 5 ~ 180915260,
            chr_num == 6 ~ 171115067,
            chr_num == 7 ~ 159138663,
            chr_num == 8 ~ 146364022,
            chr_num == 9 ~ 141213431,
            chr_num == 10 ~ 135534747,
            chr_num == 11 ~ 135006516,
            chr_num == 12 ~ 133851895,
            chr_num == 13 ~ 115169878,
            chr_num == 14 ~ 107349540,
            chr_num == 15 ~ 102531392,
            chr_num == 16 ~ 90354753,
            chr_num == 17 ~ 81195210,
            chr_num == 18 ~ 78077248,
            chr_num == 19 ~ 59128983,
            chr_num == 20 ~ 63025520,
            chr_num == 21 ~ 48129895,
            chr_num == 22 ~ 51304566,
            TRUE ~ NA_real_
          )
        )
      )
  }
  
  arm_definitions <- arm_definitions %>%
    # Filter valid arms and select final columns (following original script lines 112-113)
    filter(!is.na(arm_start) & !is.na(arm_end) & arm_start < arm_end) %>%
    select(arm, chr_num, arm_type, arm_start, arm_end, direction) %>%
    arrange(chr_num, arm_type)
  
  cat("Generated arm definitions for", length(unique(arm_definitions$chr_num)), 
      "chromosomes covering", nrow(arm_definitions), "functional regions\n")
  
  # Summary of detected alterations
  summary_stats <- arm_definitions %>%
    group_by(direction) %>%
    summarise(n_arms = n(), .groups = "drop")
  
  cat("Alteration summary:\n")
  for (i in 1:nrow(summary_stats)) {
    cat("  ", summary_stats$direction[i], ":", summary_stats$n_arms[i], "arms\n")
  }
  
  return(arm_definitions)
}

#' Process BISCUT Results with Output Validation
#'
#' Wrapper function that processes BISCUT results and validates output format
#' for compatibility with BAGEL workflow.
#'
#' @param biscut_results_file Character, path to BISCUT all_BISCUT_results.txt file
#' @param output_file Character, optional path to save arm definitions (default: NULL)
#' @return Data frame with validated arm definitions
#' @import dplyr
#' @export
create_custom_arm_definitions <- function(biscut_results_file, output_file = NULL) {
  
  cat("=== Processing Custom BISCUT Results ===\n")
  cat("Input file:", biscut_results_file, "\n")
  
  # Process BISCUT results
  arm_definitions <- process_biscut_results(biscut_results_file)
  
  # Validate format matches expected BAGEL structure
  required_cols <- c("arm", "chr_num", "arm_type", "arm_start", "arm_end", "direction")
  missing_cols <- setdiff(required_cols, colnames(arm_definitions))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns in arm definitions: ", paste(missing_cols, collapse = ", "))
  }
  
  # Validate data types and ranges
  validation_checks <- list(
    chr_num = all(arm_definitions$chr_num %in% 1:22 & !is.na(arm_definitions$chr_num)),
    arm_type = all(arm_definitions$arm_type %in% c("p", "q")),
    positions = all(arm_definitions$arm_start < arm_definitions$arm_end),
    direction = all(arm_definitions$direction %in% c("amp", "del"))
  )
  
  failed_checks <- names(validation_checks)[!unlist(validation_checks)]
  
  if (length(failed_checks) > 0) {
    stop("Validation failed for: ", paste(failed_checks, collapse = ", "))
  }
  
  cat("✅ Arm definitions validated successfully\n")
  
  # Save to file if requested
  if (!is.null(output_file)) {
    write.table(arm_definitions, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("💾 Arm definitions saved to:", output_file, "\n")
  }
  
  cat("=== Custom Arm Definitions Ready ===\n")
  
  return(arm_definitions)
}