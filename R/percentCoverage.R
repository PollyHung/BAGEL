#' Calculate Percentage of Altered Segments per Chromosomal Arm
#'
#' This function computes the percentage of amplified or deleted segments in each chromosomal arm
#' based on log2 ratios and coverage metrics. Supports both fixed thresholds and noise-based thresholds.
#'
#' @param segments A dataframe containing segmentation data. Must include columns:
#'   Sample, Arm, Start, End, Log2Ratios
#' @param gaps A dataframe containing coverage information from segmentGaps function. Must include columns:
#'   Sample, Arm, gap_length, sequenced_seg_length, total_seg_length, arm_length
#' @param thres Numeric threshold for amplification/deletion calling. If NULL, uses noise-based thresholds.
#'   (default = 0.1)
#' @param noise Optional dataframe containing noise-based thresholds (not currently implemented)
#'
#' @return A dataframe containing:
#' \itemize{
#'   \item Original segmentation data with coverage metrics
#'   \item aneuploid column indicating alteration status ("amplification", "deletion", or "no_change")
#'   \item Calculated coverage metrics: seg_length, net_seg_coverage, seg_coverage
#' }
#'
#' @examples
#' \dontrun{
#' # Using hard threshold
#' data(segments)
#' data(gaps)
#' results <- percentAlt(segments, gaps, thres = 0.15)
#'
#' # Using noise threshold (requires implementation)
#' data(noise_data)
#' results <- percentAlt(segments, gaps, thres = NULL, noise = noise_data)
#' }
#'
#' @importFrom dplyr inner_join mutate case_when
#' @export
percentCoverage <- function(segments,
                            gaps,
                            amp_thres = 0.1,
                            del_thres = -0.1,
                            noise = NULL) {
  # Input validation
  required_seg_cols <- c("Sample", "Arm", "Start", "End", "Log2Ratios")
  required_gap_cols <- c("Sample", "Arm", "gap_length", "sequenced_seg_length",
                         "total_seg_length", "arm_length")

  if (!all(required_seg_cols %in% names(segments))) {
    stop("Segments dataframe missing required columns: ",
         paste(setdiff(required_seg_cols, names(segments)), collapse = ", "))
  }

  if (!all(required_gap_cols %in% names(gaps))) {
    stop("Gaps dataframe missing required columns: ",
         paste(setdiff(required_gap_cols, names(gaps)), collapse = ", "))
  }

  # Merge data and calculate coverage metrics
  coverage <- segments %>%
    dplyr::inner_join(gaps, by = c("Sample", "Arm")) %>%
    dplyr::mutate(
      seg_length = End - Start,
      net_seg_coverage = seg_length / sequenced_seg_length,
      seg_coverage = seg_length / total_seg_length
    )

  # Calculate alteration fractions
  message("Calculating Fraction of Segments being Amp or Del")

  if (!is.null(thres)) {
    message("Using Hard Threshold: ", thres)

    coverage <- coverage %>%
      dplyr::mutate(
        aneuploid = dplyr::case_when(
          Log2Ratios <= del_thres ~ "deletion",
          Log2Ratios >= amp_thres ~ "amplification",
          TRUE ~ "no_change"
        )
      )
  } else {
    message("Using Noise-Based Threshold")
    if (is.null(noise)) stop("Noise data required when thres = NULL")
    # Implement noise-based threshold logic here
    stop("Noise-based threshold implementation not yet available")
  }

  return(coverage)
}
