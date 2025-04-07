#' Compare User Provided Breakpoint to Reference Used
#'
#' This function compares two dataframes of breakpoint data, merging them to calculate
#' overlap metrics based on the overlap between reference and cut breakpoints.
#'
#' @param reference A dataframe containing reference breakpoint data, including columns such as
#'                  'arm', 'telcent', 'direction', and other relevant metrics.
#' @param breakpoint A dataframe containing cut breakpoint data with a similar structure to 'reference'.
#' @param filename A string representing the directory path where the summary results will be saved as a CSV file.
#'
#' @return A dataframe containing merged data with calculated overlap metrics and deviation scores.
#' @details
#' The function creates unique identifiers for each row based on arm, telcent, and direction.
#' It computes overlap metrics between the reference and cut breakpoints and calculates
#' a combined overlap score that reflects the degree of overlap.
#' A warning is issued if any breakpoints do not overlap.
#'
#' @examples
#' merged_data <- compareCuts(reference_data, cuts_data, "path/to/output")
#' @export
#' @import dplyr

compareCuts <- function(reference,
                        breakpoint,
                        filename) {

  # Create unique row identifiers
  rownames(reference) <- reference$id
  rownames(breakpoint) <- breakpoint$id

  # Intersecting dataframes based on rownames
  reference.cross <- reference[intersect(rownames(reference), rownames(breakpoint)), ]
  breakpoint.cross <- breakpoint[intersect(rownames(reference), rownames(breakpoint)), ]

  # Rename columns for clarity
  reference.cross <- reference.cross %>% dplyr::rename(arm_ref = arm, telcent_ref = telcent, direction_ref = direction,
                                                       smallest_start_ref = smallest_start, largest_end_ref = largest_end)
  breakpoint.cross <- breakpoint.cross %>% dplyr::rename(arm_cut = arm, telcent_cut = telcent, direction_cut = direction,
                                                         smallest_start_cut = smallest_start, largest_end_cut = largest_end)

  # Merge the two dataframes
  merged_df <- inner_join(reference.cross, breakpoint.cross, by = "id")

  # Calculate overlap metrics and midpoint distance
  merged_df <- merged_df %>%
    mutate(
      # Overlap length (same as before)
      overlap_start = pmax(smallest_start_ref, smallest_start_cut),
      overlap_end = pmin(largest_end_ref, largest_end_cut),
      overlap_length = ifelse(overlap_start > overlap_end, 0, overlap_end - overlap_start),

      # Bidirectional overlap scores
      overlap_pct_ref = overlap_length / (largest_end_ref - smallest_start_ref),
      overlap_pct_cut = overlap_length / (largest_end_cut - smallest_start_cut),

      # Combined overlap score (harmonic mean to penalize asymmetry)
      combined_overlap = (overlap_pct_ref + overlap_pct_cut) / 2)  # Simple average

  # Give Out Warning
  non_overlapping_count <- sum(R.utils::isZero(merged_df$combined_overlap))
  if (non_overlapping_count > 0) {
    warning(sprintf("Found %d breakpoints that do not overlap between the user-provided breakpoint file and the reference breakpoint file.",
                    non_overlapping_count))
  }

  # write summary
  write.csv(merged_df, paste0(filename, ".csv"))
  return(merged_df)
}






