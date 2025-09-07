#' Improved calculateCopyNumber Function
#'
#' Calculate Copy Number Alterations from Segmentation Data with proper GISTIC-style analysis
#'
#' @param segments A dataframe containing segmentation data with required columns
#' @param breakpoints A list containing tel_bound and cent_bound dataframes
#' @param amp_thres Numeric threshold for amplifications (default=0.25)
#' @param del_thres Numeric threshold for deletions (default=-0.25)
#' @param runCentromere Logical for centromere analysis (default=FALSE)
#'
#' @return List with mean_matrix, label_matrix, and heatmap
#' @export
calculateCopyNumber <- function(segments,
                                breakpoints,
                                amp_thres = 0.25,
                                del_thres = -0.25,
                                runCentromere = FALSE) {
  
  # Input validation
  required_seg_cols <- c("Sample", "Arm", "Start", "End", "Log2Ratios")
  if (!all(required_seg_cols %in% names(segments))) {
    stop("Segments dataframe missing required columns: ",
         paste(setdiff(required_seg_cols, names(segments)), collapse = ", "))
  }
  
  if (!all(c("tel_bound") %in% names(breakpoints))) {
    stop("breakpoints must contain 'tel_bound' dataframe")
  }
  
  if (runCentromere && !"cent_bound" %in% names(breakpoints)) {
    stop("cent_bound required in breakpoints when runCentromere=TRUE")
  }
  
  # Validate thresholds
  if (amp_thres <= 0 || del_thres >= 0) {
    warning("Unusual threshold values: amp_thres should be positive, del_thres should be negative")
  }
  
  excluded_arms <- c("23p", "23q", "24p", "24q", "13p", "14p", "15p")
  
  # Helper function to process alterations with proper error handling
  process_alterations <- function(breakpoint_data, direction_type, threshold) {
    # Validate breakpoint data
    required_bp_cols <- c("arm", "direction", "smallest_start", "largest_end")
    if (!all(required_bp_cols %in% names(breakpoint_data))) {
      stop("Breakpoint data missing required columns: ",
           paste(setdiff(required_bp_cols, names(breakpoint_data)), collapse = ", "))
    }
    
    bp <- breakpoint_data %>%
      dplyr::filter(direction == direction_type,
                    !arm %in% excluded_arms) %>%
      dplyr::arrange(arm, smallest_start)
    
    if (nrow(bp) == 0) {
      message("No breakpoints found for direction: ", direction_type)
      return(data.frame())
    }
    
    # Merge segments with breakpoints
    seg_merged <- segments %>%
      dplyr::inner_join(bp, by = c("Arm" = "arm"), relationship = "many-to-many")
    
    if (nrow(seg_merged) == 0) {
      message("No segments found overlapping with breakpoints for direction: ", direction_type)
      return(data.frame())
    }
    
    # Filter segments based on breakpoint boundaries
    filtered <- seg_merged %>%
      dplyr::group_by(Sample, Arm) %>%
      dplyr::filter(case_when(
        stringr::str_detect(Arm, "p$") ~ Start <= largest_end,
        stringr::str_detect(Arm, "q$") ~ End >= smallest_start,
        TRUE ~ TRUE
      )) %>%
      dplyr::arrange(Sample, Arm, Start)
    
    # Calculate mean log2 ratios and apply thresholds
    result <- filtered %>%
      dplyr::group_by(Sample, Arm) %>%
      dplyr::summarise(
        Mean_Log2Ratios = mean(Log2Ratios, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::filter(
        (direction_type == "amp" & Mean_Log2Ratios >= threshold) |
          (direction_type == "del" & Mean_Log2Ratios <= threshold)
      ) %>%
      dplyr::mutate(
        direction = direction_type,
        label = "user"  # Can be made more sophisticated
      )
    
    return(result)
  }
  
  # Process telomeric events
  message("Processing telomeric amplifications...")
  tel_amp <- process_alterations(breakpoints$tel_bound, "amp", amp_thres)
  
  message("Processing telomeric deletions...")
  tel_del <- process_alterations(breakpoints$tel_bound, "del", del_thres)
  
  # Process centromeric events if requested
  cent_results <- data.frame()
  if (runCentromere) {
    message("Processing centromeric events...")
    
    cent <- breakpoints$cent_bound %>%
      dplyr::filter(!arm %in% excluded_arms) %>%
      dplyr::arrange(arm, smallest_start)
    
    if (nrow(cent) > 0) {
      cent_amp <- process_alterations(cent, "amp", amp_thres)
      cent_del <- process_alterations(cent, "del", del_thres)
      cent_results <- dplyr::bind_rows(cent_amp, cent_del)
    }
  }
  
  # Combine results and resolve duplicates
  seg_total <- dplyr::bind_rows(tel_amp, tel_del, cent_results)
  
  if (nrow(seg_total) == 0) {
    warning("No significant alterations found with current thresholds")
    return(list(
      seg_total = data.frame(),
      mean_matrix = matrix(nrow = 0, ncol = 0),
      label_matrix = data.frame(),
      heatmap = NULL
    ))
  }
  
  # Resolve conflicts by taking the alteration with highest absolute log2 ratio
  seg_total <- seg_total %>%
    dplyr::group_by(Sample, Arm) %>%
    dplyr::slice_max(order_by = abs(Mean_Log2Ratios), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Create output matrices
  mean_matrix <- seg_total %>%
    tidyr::pivot_wider(
      id_cols = Arm,
      names_from = Sample,
      values_from = Mean_Log2Ratios,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("Arm") %>%
    as.matrix()
  
  label_matrix <- seg_total %>%
    dplyr::select(Sample, Arm, direction) %>%
    tidyr::pivot_wider(
      id_cols = Arm,
      names_from = Sample,
      values_from = direction,
      values_fill = "neutral"
    )
  
  # Create heatmap with proper error handling
  heatmap <- NULL
  if (nrow(mean_matrix) > 0 && ncol(mean_matrix) > 0) {
    tryCatch({
      max_val <- max(abs(mean_matrix), na.rm = TRUE)
      if (max_val > 0) {
        heatmap <- pheatmap::pheatmap(
          mean_matrix,
          color = colorRampPalette(c("blue", "white", "red"))(256),
          breaks = seq(-max_val, max_val, length.out = 257),
          main = "Copy Number Alterations",
          silent = TRUE
        )
      }
    }, error = function(e) {
      warning("Could not create heatmap: ", e$message)
    })
  }
  
  return(list(
    seg_total = seg_total,
    mean_matrix = mean_matrix,
    label_matrix = label_matrix,
    heatmap = heatmap
  ))
}