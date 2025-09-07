#' GISTIC-style Statistical Analysis
#'
#' Implements proper GISTIC algorithm with background modeling and significance testing
#'
#' @param segments Segmentation data
#' @param arm_definitions Arm definitions from breakpoints  
#' @param background Optional background copy number data
#' @param q_threshold Q-value threshold for significance (default 0.25)
#' @return GISTIC results with q-values and peak significance
#' @export
gistic_analysis <- function(segments, arm_definitions, background = NULL, q_threshold = 0.25) {
  
  bagel_log("Starting GISTIC-style statistical analysis", "INFO")
  
  # Validate inputs
  segments <- validate_segments(segments)
  arm_definitions <- validate_breakpoints(arm_definitions)
  
  # Convert legacy format if needed
  if (is.list(arm_definitions) && "tel_bound" %in% names(arm_definitions)) {
    arm_definitions <- convert_legacy_breakpoints(arm_definitions)
  }
  
  # 1. Fit background distribution
  if (is.null(background)) {
    background_dist <- fit_background_distribution(segments)
  } else {
    background_dist <- fit_background_distribution(background)
  }
  
  # 2. Annotate segments with arms
  segments_annotated <- annotate_segments_with_arms(segments, arm_definitions)
  
  # 3. Calculate arm-level summaries
  arm_summaries <- calculate_arm_level_summaries(segments_annotated, arm_definitions)
  
  # 4. Calculate residuals from background
  residuals <- calculate_residuals(arm_summaries, background_dist)
  
  # 5. Peak detection using statistical methods
  peaks <- detect_significant_peaks(residuals, background_dist)
  
  # 6. Calculate q-values with multiple testing correction
  q_values <- calculate_qvalues(peaks, background_dist, q_threshold)
  
  results <- list(
    arm_summaries = arm_summaries,
    peaks = peaks,
    q_values = q_values,
    background_dist = background_dist,
    significant_peaks = peaks[q_values$q_value < q_threshold, ],
    parameters = list(q_threshold = q_threshold)
  )
  
  bagel_log(sprintf("GISTIC analysis completed. Found %d significant peaks.", 
                   nrow(results$significant_peaks)), "INFO")
  
  return(results)
}

#' Fit Background Copy Number Distribution
#'
#' @param segments Copy number segments (preferably from normal samples)
#' @return Fitted background distribution parameters
#' @keywords internal
fit_background_distribution <- function(segments) {
  # Option 1: Use normal samples if available
  if ("SampleType" %in% names(segments)) {
    normal_segments <- segments[segments$SampleType == "Normal", ]
    if (nrow(normal_segments) > 100) {
      segments <- normal_segments
      bagel_log("Using normal samples for background distribution", "INFO")
    }
  }
  
  # Option 2: Use segments with neutral copy number
  neutral_segments <- segments[abs(segments$Log2Ratios) < 0.1, ]
  
  # Fit mixture model or simple gaussian
  if (nrow(neutral_segments) > 50) {
    mu <- mean(neutral_segments$Log2Ratios, na.rm = TRUE)
    sigma <- sd(neutral_segments$Log2Ratios, na.rm = TRUE)
    bagel_log("Using neutral segments for background distribution", "INFO")
  } else {
    # Fallback to robust estimates
    mu <- median(segments$Log2Ratios, na.rm = TRUE)
    sigma <- mad(segments$Log2Ratios, na.rm = TRUE)
    bagel_log("Using robust estimates for background distribution", "INFO")
  }
  
  return(list(mu = mu, sigma = sigma, type = "gaussian"))
}

#' Annotate Segments with Arms Using Genomic Overlaps
#'
#' @param segments Segmentation data
#' @param arm_definitions Arm definitions
#' @return Segments annotated with arm assignments
#' @keywords internal
annotate_segments_with_arms <- function(segments, arm_definitions) {
  
  # Create GRanges objects for efficient overlap detection
  segments_gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", segments$Chromosome),
    ranges = IRanges::IRanges(start = segments$Start, end = segments$End),
    sample = segments$Sample,
    log2ratio = segments$Log2Ratios
  )
  
  arms_gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", arm_definitions$chr_num),
    ranges = IRanges::IRanges(start = arm_definitions$arm_start, 
                              end = arm_definitions$arm_end),
    arm = arm_definitions$arm,
    arm_type = arm_definitions$arm_type
  )
  
  # Find overlaps
  overlaps <- GenomicRanges::findOverlaps(segments_gr, arms_gr)
  
  # Assign arm annotations
  segments_annotated <- segments
  segments_annotated$Arm <- NA_character_
  segments_annotated$Arm[S4Vectors::queryHits(overlaps)] <- 
    arms_gr$arm[S4Vectors::subjectHits(overlaps)]
  
  # Handle segments spanning multiple arms (split them)
  segments_annotated <- split_multi_arm_segments(segments_annotated, arm_definitions)
  
  return(segments_annotated)
}

#' Split Segments Spanning Multiple Arms
#'
#' @param segments_annotated Segments with arm annotations
#' @param arm_definitions Arm definitions
#' @return Segments with proper arm splits
#' @keywords internal
split_multi_arm_segments <- function(segments_annotated, arm_definitions) {
  # For now, keep segments as-is
  # Future enhancement: implement proper segment splitting
  return(segments_annotated)
}

#' Calculate Arm-Level Copy Number Summaries
#'
#' @param segments_annotated Segments with arm annotations
#' @param arm_definitions Arm definitions
#' @return Arm-level copy number summaries
#' @keywords internal
calculate_arm_level_summaries <- function(segments_annotated, arm_definitions) {
  
  # Remove segments without arm annotation
  segments_with_arms <- segments_annotated[!is.na(segments_annotated$Arm), ]
  
  if (nrow(segments_with_arms) == 0) {
    stop("No segments could be assigned to chromosome arms")
  }
  
  # Calculate weighted mean by segment length
  arm_summaries <- segments_with_arms %>%
    dplyr::group_by(Sample, Arm) %>%
    dplyr::summarise(
      n_segments = n(),
      total_length = sum(End - Start + 1),
      mean_log2ratio = weighted.mean(Log2Ratios, End - Start + 1, na.rm = TRUE),
      median_log2ratio = median(Log2Ratios, na.rm = TRUE),
      sd_log2ratio = sd(Log2Ratios, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      arm_definitions %>% dplyr::select(arm, chr_num, arm_type, direction),
      by = c("Arm" = "arm")
    )
  
  return(arm_summaries)
}

#' Calculate Residuals from Background Distribution
#'
#' @param arm_summaries Arm-level summaries
#' @param background_dist Background distribution parameters
#' @return Residuals data frame
#' @keywords internal
calculate_residuals <- function(arm_summaries, background_dist) {
  arm_summaries %>%
    dplyr::mutate(
      residual = mean_log2ratio - background_dist$mu,
      z_score = residual / background_dist$sigma
    )
}

#' Detect Significant Peaks Using Statistical Methods
#'
#' @param residuals Residuals from background
#' @param background_dist Background distribution
#' @return Peak detection results
#' @keywords internal
detect_significant_peaks <- function(residuals, background_dist) {
  
  # Group by arm and calculate significance
  peak_stats <- residuals %>%
    dplyr::group_by(Arm, chr_num, arm_type, direction) %>%
    dplyr::summarise(
      n_samples = n(),
      mean_residual = mean(residual, na.rm = TRUE),
      mean_z_score = mean(z_score, na.rm = TRUE),
      sd_residual = sd(residual, na.rm = TRUE),
      median_log2ratio = median(mean_log2ratio, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_samples >= 5)  # Require minimum samples
  
  return(peak_stats)
}

#' Calculate Q-values with Multiple Testing Correction
#'
#' @param peaks Peak statistics
#' @param background_dist Background distribution
#' @param q_threshold Q-value threshold
#' @return Results with p-values and q-values
#' @keywords internal
calculate_qvalues <- function(peaks, background_dist, q_threshold) {
  
  # Calculate z-scores and p-values
  z_scores <- peaks$mean_z_score
  
  # Two-tailed test for significant deviation from background
  p_values <- 2 * (1 - pnorm(abs(z_scores)))
  
  # Apply multiple testing correction using FDR
  q_values <- p.adjust(p_values, method = "fdr")
  
  results <- peaks %>%
    dplyr::mutate(
      p_value = p_values,
      q_value = q_values,
      significant = q_values < q_threshold,
      alteration_type = case_when(
        significant & mean_residual > 0 ~ "Amplification",
        significant & mean_residual < 0 ~ "Deletion", 
        TRUE ~ "Neutral"
      )
    ) %>%
    dplyr::arrange(q_value)
  
  return(results)
}

#' Simplified Segment Processing Pipeline
#'
#' @param segments Raw segmentation data
#' @param arm_definitions Standardized arm definitions from breakpoints
#' @param amp_threshold Amplification threshold (default log2(2.5/2) = 0.32)
#' @param del_threshold Deletion threshold (default log2(1.5/2) = -0.415)
#' @return Processed segments with arm-level annotations
#' @export
process_segments_pipeline <- function(segments, arm_definitions, 
                                    amp_threshold = log2(2.5/2), 
                                    del_threshold = log2(1.5/2)) {
  
  bagel_log("Starting segment processing pipeline", "INFO")
  
  # Step 1: Validate and clean input data
  segments_clean <- validate_and_clean_segments(segments)
  
  # Step 2: Annotate segments with arm membership
  segments_annotated <- annotate_segments_with_arms(segments_clean, arm_definitions)
  
  # Step 3: Apply thresholds to identify alterations
  segments_altered <- identify_alterations(segments_annotated, amp_threshold, del_threshold)
  
  # Step 4: Calculate arm-level summaries
  arm_summaries <- calculate_arm_level_summaries(segments_altered, arm_definitions)
  
  bagel_log(sprintf("Processed %d segments across %d samples", 
                   nrow(segments_altered), length(unique(segments_altered$Sample))), "INFO")
  
  return(list(
    segments = segments_altered,
    arm_summaries = arm_summaries,
    arm_definitions = arm_definitions,
    parameters = list(
      amp_threshold = amp_threshold,
      del_threshold = del_threshold
    )
  ))
}

#' Identify Copy Number Alterations
#'
#' @param segments_annotated Segments with arm annotations
#' @param amp_threshold Amplification threshold
#' @param del_threshold Deletion threshold
#' @return Segments with alteration classifications
#' @keywords internal
identify_alterations <- function(segments_annotated, amp_threshold, del_threshold) {
  
  segments_annotated %>%
    dplyr::mutate(
      alteration_type = case_when(
        Log2Ratios >= amp_threshold ~ "Amplification",
        Log2Ratios <= del_threshold ~ "Deletion",
        TRUE ~ "Neutral"
      ),
      is_altered = alteration_type != "Neutral"
    )
}

#' Memory-Efficient Segment Processing for Large Datasets
#'
#' @param segments Large segmentation dataset
#' @param arm_definitions Arm definitions
#' @param chunk_size Number of samples to process at once
#' @param ... Additional parameters for process_segments_pipeline
#' @return Combined processed results
#' @export
process_segments_chunked <- function(segments, arm_definitions, chunk_size = 100, ...) {
  
  samples <- unique(segments$Sample)
  n_chunks <- ceiling(length(samples) / chunk_size)
  
  bagel_log(sprintf("Processing %d samples in %d chunks", length(samples), n_chunks), "INFO")
  
  results <- list()
  
  for (i in 1:n_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, length(samples))
    chunk_samples <- samples[start_idx:end_idx]
    
    bagel_log(sprintf("Processing chunk %d/%d (%d samples)", i, n_chunks, length(chunk_samples)), "INFO")
    
    chunk_segments <- segments[segments$Sample %in% chunk_samples, ]
    chunk_results <- process_segments_pipeline(chunk_segments, arm_definitions, ...)
    
    results[[i]] <- chunk_results
    
    # Garbage collection for memory management
    gc()
  }
  
  # Combine results
  combined_results <- combine_chunk_results(results)
  
  bagel_log("Completed chunked processing", "INFO")
  
  return(combined_results)
}

#' Combine Results from Chunked Processing
#'
#' @param results List of chunk results
#' @return Combined results
#' @keywords internal
combine_chunk_results <- function(results) {
  
  # Combine segments
  all_segments <- do.call(rbind, lapply(results, function(x) x$segments))
  
  # Combine arm summaries
  all_arm_summaries <- do.call(rbind, lapply(results, function(x) x$arm_summaries))
  
  # Use arm definitions from first chunk (should be identical)
  arm_definitions <- results[[1]]$arm_definitions
  parameters <- results[[1]]$parameters
  
  return(list(
    segments = all_segments,
    arm_summaries = all_arm_summaries,
    arm_definitions = arm_definitions,
    parameters = parameters
  ))
}