#' Validate Segmentation Data Frame
#'
#' @param segments Data frame with segmentation data
#' @return Validated data frame or throws error
#' @keywords internal
validate_segments <- function(segments) {
  required_cols <- c("Sample", "Chromosome", "Start", "End", "Log2Ratios")
  missing_cols <- setdiff(required_cols, names(segments))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Additional validation
  if (nrow(segments) == 0) {
    stop("Segments data frame is empty")
  }
  
  # Check for numeric columns
  numeric_cols <- c("Chromosome", "Start", "End", "Log2Ratios")
  for (col in numeric_cols) {
    if (!is.numeric(segments[[col]])) {
      stop(sprintf("Column '%s' must be numeric", col))
    }
  }
  
  # Check for valid coordinate ranges
  invalid_coords <- segments$Start >= segments$End
  if (any(invalid_coords, na.rm = TRUE)) {
    stop("Invalid coordinates found: Start >= End")
  }
  
  return(segments)
}

#' Validate Breakpoint Data Structure
#'
#' @param breakpoints Breakpoint data structure
#' @return Validated breakpoints or throws error
#' @keywords internal
validate_breakpoints <- function(breakpoints) {
  if (is.null(breakpoints)) {
    stop("breakpoints cannot be NULL")
  }
  
  # Handle different breakpoint formats
  if (is.data.frame(breakpoints)) {
    required_cols <- c("arm", "chr_num", "arm_start", "arm_end")
    missing_cols <- setdiff(required_cols, names(breakpoints))
    if (length(missing_cols) > 0) {
      stop("Breakpoints data frame missing required columns: ", 
           paste(missing_cols, collapse = ", "))
    }
  } else if (is.list(breakpoints) && "tel_bound" %in% names(breakpoints)) {
    # Legacy format support
    if (!is.data.frame(breakpoints$tel_bound)) {
      stop("breakpoints$tel_bound must be a data frame")
    }
  } else {
    stop("breakpoints must be a data frame with required columns or a list containing 'tel_bound'")
  }
  
  return(breakpoints)
}

#' Clean and Standardize Segments
#'
#' @param segments Raw segmentation data
#' @param remove_outliers Logical, whether to remove extreme outliers
#' @return Cleaned and standardized segments
#' @keywords internal
validate_and_clean_segments <- function(segments, remove_outliers = TRUE) {
  # Input validation
  segments <- validate_segments(segments)
  
  # Clean and standardize
  segments_clean <- segments %>%
    dplyr::filter(
      Chromosome %in% 1:22,  # Only autosomal chromosomes
      Start < End,           # Valid segment coordinates
      !is.na(Log2Ratios)     # Valid log2 ratios
    ) %>%
    dplyr::mutate(
      Chromosome = as.numeric(Chromosome),
      Segment_Length = End - Start + 1
    )
  
  # Remove outliers (optional, configurable)
  if (remove_outliers) {
    q99 <- quantile(segments_clean$Log2Ratios, 0.99, na.rm = TRUE)
    q01 <- quantile(segments_clean$Log2Ratios, 0.01, na.rm = TRUE)
    
    segments_clean <- segments_clean %>%
      dplyr::filter(Log2Ratios >= q01 & Log2Ratios <= q99)
  }
  
  return(segments_clean)
}

#' Setup BAGEL Environment for Logging
#' 
#' @keywords internal
.bagel_env <- new.env(parent = emptyenv())

#' Setup Logging System
#'
#' @param log_level Character, logging level ("DEBUG", "INFO", "WARN", "ERROR")
#' @param log_file Character, optional file path for logging
#' @export
setup_bagel_logging <- function(log_level = "INFO", log_file = NULL) {
  .bagel_env$log_level <- log_level
  .bagel_env$log_file <- log_file
}

#' Internal Logging Function
#'
#' @param message Character, message to log
#' @param level Character, log level
#' @keywords internal
bagel_log <- function(message, level = "INFO") {
  if (!exists("log_level", envir = .bagel_env)) {
    setup_bagel_logging()
  }
  
  timestamp <- Sys.time()
  log_entry <- sprintf("[%s] %s: %s", timestamp, level, message)
  
  # Always show ERROR and WARN
  if (level %in% c("ERROR", "WARN") || 
      (.bagel_env$log_level %in% c("DEBUG", "INFO") && level == "INFO") ||
      (.bagel_env$log_level == "DEBUG" && level == "DEBUG")) {
    cat(log_entry, "\n")
  }
  
  if (!is.null(.bagel_env$log_file)) {
    cat(log_entry, "\n", file = .bagel_env$log_file, append = TRUE)
  }
}

#' Robust Function Wrapper with Error Handling
#'
#' @param func Function to execute
#' @param ... Arguments to pass to function
#' @return Function result or throws error
#' @keywords internal
with_error_handling <- function(func, ...) {
  func_name <- deparse(substitute(func))
  
  tryCatch({
    result <- func(...)
    bagel_log(sprintf("Successfully completed %s", func_name))
    return(result)
  }, error = function(e) {
    error_msg <- sprintf("Error in %s: %s", func_name, e$message)
    bagel_log(error_msg, "ERROR")
    stop(error_msg, call. = FALSE)
  }, warning = function(w) {
    warning_msg <- sprintf("Warning in %s: %s", func_name, w$message)
    bagel_log(warning_msg, "WARN")
    warning(warning_msg, call. = FALSE)
  })
}

#' Get Built-in Breakpoint Data
#'
#' @param cancer_type Character, specific cancer type or "consensus" for pan-cancer
#' @param data_type Character, type of breakpoint data to retrieve
#' @return Breakpoint data frame
#' @export
get_breakpoint_data <- function(cancer_type = "consensus", data_type = "arm_definitions") {
  if (cancer_type == "consensus") {
    if (data_type == "arm_definitions") {
      if (!exists("consensus_arm_definitions")) {
        stop("Pan-cancer consensus breakpoint data not found. Please load BAGEL package data.")
      }
      
      # Update consensus arm definitions to use BISCUT coordinates for consistency
      tryCatch({
        updated_consensus <- update_arm_definitions_with_biscut(consensus_arm_definitions)
        return(updated_consensus)
      }, error = function(e) {
        warning("Could not update with BISCUT coordinates, using original: ", e$message)
        return(consensus_arm_definitions)
      })
    } else if (data_type == "breakpoints") {
      if (!exists("consensus_breakpoints")) {
        stop("Pan-cancer consensus breakpoint data not found. Please load BAGEL package data.")
      }
      return(consensus_breakpoints)
    }
  } else {
    if (!exists("cancer_specific_breakpoints")) {
      stop("Cancer-specific breakpoint data not found. Please load BAGEL package data.")
    }
    
    if (!cancer_type %in% names(cancer_specific_breakpoints)) {
      available_types <- paste(names(cancer_specific_breakpoints), collapse = ", ")
      stop(sprintf("Cancer type '%s' not found. Available types: %s", 
                   cancer_type, available_types))
    }
    
    return(cancer_specific_breakpoints[[cancer_type]])
  }
}

#' Convert Legacy Breakpoint Format to New Format
#'
#' @param legacy_breakpoints List with tel_bound element (legacy format)
#' @return Standardized breakpoint data frame
#' @keywords internal
convert_legacy_breakpoints <- function(legacy_breakpoints) {
  if (!"tel_bound" %in% names(legacy_breakpoints)) {
    stop("Legacy breakpoints must contain 'tel_bound' element")
  }
  
  tel_bound <- legacy_breakpoints$tel_bound
  
  # Convert to new format
  standardized <- tel_bound %>%
    dplyr::mutate(
      chr_num = as.numeric(gsub("[pq]", "", arm)),
      arm_type = ifelse(grepl("p$", arm), "p", "q"),
      arm_start = case_when(
        arm_type == "p" ~ 1,
        arm_type == "q" ~ smallest_start
      ),
      arm_end = case_when(
        arm_type == "p" ~ largest_end,
        arm_type == "q" ~ case_when(
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
    ) %>%
    dplyr::select(arm, chr_num, arm_type, arm_start, arm_end, direction) %>%
    dplyr::arrange(chr_num, arm_type)
  
  return(standardized)
}