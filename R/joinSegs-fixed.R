#' Join Segments Based on Alteration Type (Fixed Version)
#'
#' This function processes genomic segments based on their alteration type, joining adjacent 
#' segments if the gap between them is smaller than their combined lengths. Includes improved
#' error handling and cleaner logic.
#'
#' @param segments A dataframe containing segment information with columns for `Sample`, 
#'   `Chromosome`, `Start`, `End`, `Segment_Mean` (or `Log2Ratios`).
#' @param aneu A character string indicating the alteration type. Options are "amp", "del", 
#'   "non-del", and "non-amp".
#' @param arm_definition A data frame or list defining the arm boundaries and coordinates.
#'   If legacy format (list with coord element), will be converted.
#' @param amp_threshold Amplification threshold (default: log2(2.5/2))
#' @param del_threshold Deletion threshold (default: log2(1.5/2))  
#' @param telcent A character string indicating whether to check for 'telomere' or 'centromere' alignment.
#' @param integer A logical value indicating whether to return integer lengths or normalized percentages.
#' @param arm Character string specifying the arm (e.g., "1p", "1q")
#'
#' @return A dataframe summarizing the joined segments with columns for `Sample`, `percent`,
#'   `start`, and `end`.
#'
#' @export
joinSegs_fixed <- function(segments,
                          aneu,
                          arm_definition,
                          amp_threshold = log2(2.5/2),
                          del_threshold = log2(1.5/2),
                          telcent = "tel",
                          integer = FALSE,
                          arm) {
  
  # Input validation
  segments <- validate_segments(segments)
  
  if (missing(arm) || is.null(arm)) {
    stop("Arm parameter is required (e.g., '1p', '1q')")
  }
  
  if (!aneu %in% c("amp", "del", "non-del", "non-amp")) {
    stop("aneu must be one of: 'amp', 'del', 'non-del', 'non-amp'")
  }
  
  bagel_log(sprintf("Processing segments for arm %s, alteration type: %s", arm, aneu), "INFO")
  
  # Handle different arm_definition formats
  if (is.list(arm_definition) && !is.data.frame(arm_definition)) {
    # Legacy format - extract coordinates
    if (!"tel_bound" %in% names(arm_definition)) {
      stop("Legacy arm_definition must contain 'tel_bound' element")
    }
    
    arm_info <- arm_definition$tel_bound[arm_definition$tel_bound$arm == arm, ]
    if (nrow(arm_info) == 0) {
      stop(sprintf("Arm '%s' not found in breakpoint data", arm))
    }
    
    coord <- c(arm_info$smallest_start, arm_info$largest_end)
  } else if (is.data.frame(arm_definition)) {
    # New format
    arm_info <- arm_definition[arm_definition$arm == arm, ]
    if (nrow(arm_info) == 0) {
      stop(sprintf("Arm '%s' not found in arm_definition data", arm))
    }
    
    coord <- c(arm_info$arm_start, arm_info$arm_end)
  } else {
    stop("arm_definition must be a data frame or list with breakpoint information")
  }
  
  # Determine direction and thresholds
  direction <- switch(aneu,
                     "amp" = 1,
                     "del" = -1,
                     "non-del" = -2,
                     "non-amp" = 2,
                     0)
  
  threshold <- switch(aneu,
                     "amp" = amp_threshold,
                     "del" = del_threshold,
                     "non-del" = del_threshold,
                     "non-amp" = amp_threshold,
                     0)
  
  # Calculate arm length
  arm_length <- coord[2] - coord[1] + 1
  
  # Extract chromosome number
  chr_num <- as.numeric(gsub("[pq]$", "", arm))
  
  # Filter segments for the specific chromosome and arm region
  arm_segments <- segments %>%
    dplyr::filter(
      Chromosome == chr_num,
      End >= coord[1],
      Start <= coord[2]
    ) %>%
    dplyr::mutate(
      # Clip segments to arm boundaries
      Start = pmax(Start, coord[1]),
      End = pmin(End, coord[2]),
      # Handle both column names for segment mean
      segment_value = if ("Segment_Mean" %in% names(.)) Segment_Mean else Log2Ratios
    ) %>%
    dplyr::filter(Start < End)  # Remove invalid segments after clipping
  
  if (nrow(arm_segments) == 0) {
    bagel_log(sprintf("No segments found for arm %s", arm), "WARN")
    return(data.frame(Sample = character(0), percent = numeric(0), 
                     start = numeric(0), end = numeric(0)))
  }
  
  # Process each sample separately
  results <- list()
  samples <- unique(arm_segments$Sample)
  
  for (sample in samples) {
    sample_segments <- arm_segments %>% 
      dplyr::filter(Sample == sample) %>%
      dplyr::arrange(Start)
    
    # Classify segments based on alteration type
    classified_segments <- sample_segments %>%
      dplyr::mutate(
        is_altered = case_when(
          aneu == "del" ~ segment_value <= threshold,
          aneu == "amp" ~ segment_value >= threshold,
          aneu == "non-del" ~ segment_value > threshold,
          aneu == "non-amp" ~ segment_value < threshold,
          TRUE ~ FALSE
        )
      ) %>%
      dplyr::filter(is_altered) %>%
      dplyr::select(Start, End)
    
    # Handle case where no altered segments found
    if (nrow(classified_segments) == 0) {
      results[[sample]] <- data.frame(percent = 0, start = 0, end = 0)
      next
    }
    
    # Check telomere/centromere alignment
    arm_type <- ifelse(grepl("p$", arm), "p", "q")
    passes_alignment <- check_telcent_alignment(classified_segments, coord, arm_type, telcent)
    
    if (!passes_alignment) {
      results[[sample]] <- data.frame(percent = 0, start = 0, end = 0)
      next
    }
    
    # Join segments using improved logic
    joined_result <- join_segments_improved(classified_segments, coord, arm_type, telcent)
    
    # Convert to appropriate format
    if (integer) {
      results[[sample]] <- data.frame(
        percent = joined_result$length,
        start = joined_result$start,
        end = joined_result$end
      )
    } else {
      # Normalize to percentages
      start_norm <- (joined_result$start - coord[1]) / arm_length
      end_norm <- (joined_result$end - coord[1]) / arm_length
      percent_norm <- joined_result$length / arm_length
      
      # Handle coordinate system based on arm type and telcent
      if (should_flip_coordinates(arm_type, telcent)) {
        temp_start <- 1 - end_norm
        end_norm <- 1 - start_norm
        start_norm <- temp_start
      }
      
      results[[sample]] <- data.frame(
        percent = max(0, min(1, percent_norm)),
        start = max(0, min(1, start_norm)),
        end = max(0, min(1, end_norm))
      )
    }
  }
  
  # Combine results
  result_df <- dplyr::bind_rows(results, .id = "Sample")
  
  bagel_log(sprintf("Completed processing for arm %s: %d samples processed", 
                   arm, nrow(result_df)), "INFO")
  
  return(result_df)
}

#' Check Telomere/Centromere Alignment
#' 
#' @param segments Classified segments
#' @param coord Arm coordinates
#' @param arm_type p or q arm
#' @param telcent tel or cent
#' @return Logical indicating if alignment passes
#' @keywords internal
check_telcent_alignment <- function(segments, coord, arm_type, telcent) {
  
  if (telcent == "tel") {
    # For telomere alignment, check if segments start from telomere end
    if (arm_type == "p") {
      return(segments$Start[1] == coord[1])  # p-arm telomere is at start
    } else {
      return(segments$End[nrow(segments)] == coord[2])  # q-arm telomere is at end
    }
  } else if (telcent == "cent") {
    # For centromere alignment, check if segments reach centromere
    if (arm_type == "p") {
      return(segments$End[nrow(segments)] == coord[2])  # p-arm centromere is at end
    } else {
      return(segments$Start[1] == coord[1])  # q-arm centromere is at start
    }
  }
  
  # If neither tel nor cent specified, pass alignment check
  return(TRUE)
}

#' Improved Segment Joining Logic
#'
#' @param segments Segments to join
#' @param coord Arm coordinates
#' @param arm_type p or q arm
#' @param telcent tel or cent
#' @return List with length, start, end
#' @keywords internal
join_segments_improved <- function(segments, coord, arm_type, telcent) {
  
  if (nrow(segments) == 1) {
    return(list(
      length = segments$End[1] - segments$Start[1] + 1,
      start = segments$Start[1],
      end = segments$End[1]
    ))
  }
  
  # Sort segments
  segments <- segments[order(segments$Start), ]
  
  # Calculate gaps between segments
  gaps <- segments$Start[-1] - segments$End[-nrow(segments)] - 1
  total_gap_length <- sum(gaps[gaps > 0])
  
  # Calculate total altered segment length
  total_altered_length <- sum(segments$End - segments$Start + 1)
  
  # If gaps are smaller than altered segments, join everything
  if (total_gap_length <= total_altered_length) {
    return(list(
      length = segments$End[nrow(segments)] - segments$Start[1] + 1,
      start = segments$Start[1],
      end = segments$End[nrow(segments)]
    ))
  }
  
  # Otherwise, need to trim segments to make gaps <= altered length
  current_segments <- segments
  
  while (nrow(current_segments) > 1) {
    current_gaps <- current_segments$Start[-1] - current_segments$End[-nrow(current_segments)] - 1
    current_gap_length <- sum(current_gaps[current_gaps > 0])
    current_altered_length <- sum(current_segments$End - current_segments$Start + 1)
    
    if (current_gap_length <= current_altered_length) {
      break
    }
    
    # Remove segments from the appropriate end based on arm type and telcent
    if (should_trim_from_end(arm_type, telcent)) {
      current_segments <- current_segments[-nrow(current_segments), ]
    } else {
      current_segments <- current_segments[-1, ]
    }
  }
  
  if (nrow(current_segments) == 0) {
    return(list(length = 0, start = 0, end = 0))
  }
  
  return(list(
    length = current_segments$End[nrow(current_segments)] - current_segments$Start[1] + 1,
    start = current_segments$Start[1],
    end = current_segments$End[nrow(current_segments)]
  ))
}

#' Determine if coordinates should be flipped
#'
#' @param arm_type p or q
#' @param telcent tel or cent
#' @return Logical
#' @keywords internal
should_flip_coordinates <- function(arm_type, telcent) {
  # Flip coordinates when NOT (p-arm + tel) OR (q-arm + cent)
  return(!((arm_type == "p" && telcent == "tel") || (arm_type == "q" && telcent == "cent")))
}

#' Determine if trimming should happen from end
#'
#' @param arm_type p or q
#' @param telcent tel or cent 
#' @return Logical
#' @keywords internal
should_trim_from_end <- function(arm_type, telcent) {
  # Trim from end for (p-arm + tel) OR (q-arm + cent)
  return((arm_type == "p" && telcent == "tel") || (arm_type == "q" && telcent == "cent"))
}