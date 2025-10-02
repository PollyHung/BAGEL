#' Load and Clean Segmentation Data for a Cancer Type
#'
#' @param cancer_type Character, cancer type name
#' @param data_dir Character, path to data directory
#' @return Cleaned segmentation data frame
#' @export
load_segments <- function(cancer_type, data_dir) {
    seg_file <- file.path(data_dir, cancer_type, "segmentation.seg")
    
    if (!file.exists(seg_file)) {
        stop(paste("Segmentation file not found for", cancer_type, "at", seg_file))
    }
    
    # Load data
    segments <- read_tsv(seg_file, show_col_types = FALSE)
    
    # Convert column name for consistency
    if ("Segment_Mean" %in% names(segments) && !"Log2Ratios" %in% names(segments)) {
        segments$Log2Ratios <- segments$Segment_Mean
    }
    
    # Clean data
    original_rows <- nrow(segments)
    segments_clean <- segments %>%
        filter(
            !is.na(Start) & !is.na(End) & !is.na(Log2Ratios),
            Start < End,
            Chromosome %in% 1:22
        ) %>%
        mutate(
            Chromosome = as.numeric(Chromosome),
            Segment_Length = End - Start + 1
        )
    
    cleaned_rows <- nrow(segments_clean)
    log_bagel(paste("Loaded", cancer_type, ":", original_rows, "->", cleaned_rows, "segments (",
        original_rows - cleaned_rows, "removed)"))
    
    return(segments_clean)
}


#' Setup Logging System
#'
#' @param log_level Character, logging level ("DEBUG", "INFO", "WARN", "ERROR")
#' @param log_file Character, optional file path for logging
#' @export
log_bagel <- function(log_level = "INFO", log_file = NULL) {
  .bagel_env <- new.env(parent = emptyenv())
  .bagel_env$log_level <- log_level
  .bagel_env$log_file <- log_file
}


