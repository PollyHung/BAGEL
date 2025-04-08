#' Calculate Copy Number Alterations from Segmentation Data
#'
#' This function identifies telomere and centromere-bound copy number alterations (amplifications and deletions)
#' based on segmentation data and predefined breakpoints. It returns summarized results and visualizations.
#'
#' @param segments A dataframe containing segmentation data with columns: Sample, Chromosome, Start, End,
#'        Num_Probes, Log2Ratios, Arm
#' @param breakpoints A list containing two dataframes:
#'        - tel_bound: Telomeric breakpoints with columns: arm, direction, smallest_start, largest_end
#'        - cent_bound: Centromeric breakpoints (required if runCentromere=TRUE)
#' @param thres Numeric threshold for determining amplifications/deletions (default=0.25)
#' @param runCentromere Logical indicating whether to analyze centromere-bound events (default=FALSE)
#'
#' @return A list containing:
#' \itemize{
#'   \item mean_matrix - Matrix of mean log2 ratios (rows=arms, columns=samples)
#'   \item label_matrix - Matrix of alteration labels (rows=arms, columns=samples)
#'   \item heatmap - pheatmap object visualization of alterations
#' }
#'
#' @examples
#' \dontrun{
#' data(example_segments)
#' data(example_breakpoints)
#' results <- calculateCopyNumber(segments = example_segments,
#'                               breakpoints = example_breakpoints,
#'                               thres = 0.3,
#'                               runCentromere = TRUE)
#' }
#'
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import pheatmap pheatmap
#' @export
calculateCopyNumber <- function(segments,
                                breakpoints,
                                amp_thres = 0.25,
                                del_thres = 0.25,
                                runCentromere = FALSE) {

  # Validate input structure
  if (!all(c("tel_bound") %in% names(breakpoints))) {
    stop("breakpoints must contain 'tel_bound' dataframe")
  }
  if (runCentromere && !"cent_bound" %in% names(breakpoints)) {
    stop("cent_bound required in breakpoints when runCentromere=TRUE")
  }

  excluded_arms <- c("23p", "23q", "24p", "24q", "13p", "14p", "15p")

  # Helper function to process alterations
  process_alterations <- function(breakpoint_data, direction_type) {
    bp <- breakpoint_data %>%
      dplyr::filter(direction == direction_type,
                    !arm %in% excluded_arms) %>%
      dplyr::arrange(arm, smallest_start)

    if (nrow(bp) == 0) return(data.frame())

    seg_merged <- segments %>%
      dplyr::inner_join(bp, by = c("Arm" = "arm"))

    filtered <- seg_merged %>%
      dplyr::group_by(Sample, Arm) %>%
      dplyr::filter(case_when(
        stringr::str_detect(Arm, "p$") ~ Start <= largest_end,
        stringr::str_detect(Arm, "q$") ~ End >= smallest_start,
        TRUE ~ TRUE
      )) %>%
      dplyr::arrange(Sample, Arm, Start)

    filtered %>%
      dplyr::group_by(Sample, Arm, label) %>%
      dplyr::summarise(
        Mean_Log2Ratios = mean(Log2Ratios, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::filter(if (direction_type == "amp") {
        Mean_Log2Ratios > amp_thres
      } else {
        Mean_Log2Ratios < del_thres
      })
  }

  # Process telomeric events
  tel_amp <- process_alterations(breakpoints$tel_bound, "amp")
  tel_del <- process_alterations(breakpoints$tel_bound, "del")

  # Process centromeric events if requested
  cent_results <- data.frame()
  if (runCentromere) {
    cent <- breakpoints$cent_bound %>%
      dplyr::filter(!arm %in% excluded_arms,
                    label %in% c("user", "user_ref")) %>%
      dplyr::arrange(arm, smallest_start)

    if (nrow(cent) > 0) {
      cent_merged <- segments %>%
        dplyr::inner_join(cent, by = c("Arm" = "arm"))

      cent_filtered <- cent_merged %>%
        dplyr::group_by(Sample, Arm) %>%
        dplyr::filter(case_when(
          stringr::str_detect(Arm, "p$") ~ Start >= smallest_start,
          stringr::str_detect(Arm, "q$") ~ End <= largest_end,
          TRUE ~ TRUE
        )) %>%
        dplyr::arrange(Sample, Arm, Start)

      cent_results <- cent_filtered %>%
        dplyr::group_by(Sample, Arm, label, direction) %>%
        dplyr::summarise(
          Mean_Log2Ratios = mean(Log2Ratios, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::filter((direction == "del" & Mean_Log2Ratios < del_thres) |
                        (direction == "amp" & Mean_Log2Ratios > amp_thres))
    }
  }

  # Combine results and resolve duplicates
  seg_total <- dplyr::bind_rows(tel_amp, tel_del, cent_results) %>%
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
    tidyr::pivot_wider(
      id_cols = Arm,
      names_from = Sample,
      values_from = label
    )

  # Create heatmap
  heatmap <- pheatmap::pheatmap(
    mean_matrix,
    color = colorRampPalette(c("blue", "white", "red"))(256),
    breaks = seq(-max(abs(mean_matrix)), max(abs(mean_matrix)), length.out = 257),
    main = "Copy Number Alterations"
  )

  return(list(
    seg_total = seg_total,
    mean_matrix = mean_matrix,
    label_matrix = label_matrix,
    heatmap = heatmap
  ))
}
