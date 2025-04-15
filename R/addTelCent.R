#' Add Telomere or Centromere Labels to Segments
#'
#' This function filters segments based on the specified direction (AMP or DEL)
#' and assigns telomere or centromere labels based on the arm and direction.
#' The resulting dataframe is arranged according to the specified criteria.
#'
#' @param segments A dataframe containing segment data. It must include columns
#'   for `Status`, `Start`, `End`, and `Sample`.
#' @param direction A character string specifying the direction. It must be either
#'   "AMP" or "DEL".
#' @param telcent A character string that specifies whether to label segments
#'   as "TEL" (telomere) or "CENT" (centromere).
#' @param coord A numeric vector of length 2 representing the coordinates of the
#'   telomere and centromere boundaries (for example, `c(Start, End)`).
#' @param arm A character string indicating the arm of the chromosome (e.g., "1p" or "1q").
#'
#' @return A dataframe that includes an additional column `telcent` indicating
#'   the assigned label for each segment.
#'
#' @examples
#' # Example usage
#' result <- addTelCent(segments = segments, direction = "AMP", telcent = "TEL", coord = c(100, 200), arm = "1q")
#'
#' @export
addTelCent <- function(segments,
                       direction,
                       telcent,
                       Arm,
                       coord,
                       fuzzy = FALSE) {  # Added coord parameter
  segments <- segments %>% dplyr::filter(Status == direction)

  if (telcent == "TEL") {
    if (grepl("p", Arm)) {
      segments$telcent <- ifelse(segments$Start <= coord$Start, "TEL", "INTER")
      segments <- segments %>% dplyr::arrange(Sample, Start)
    }
    if (grepl("q", Arm)) {
      segments$telcent <- ifelse(segments$End >= coord$End, "TEL", "INTER")
      segments <- segments %>% dplyr::arrange(Sample, -Start)
    }
  } else if (telcent == "CENT") {
    if (grepl("p", Arm)) {
      segments$telcent <- ifelse(segments$End >= coord$End, "CENT", "INTER")
      segments <- segments %>% dplyr::arrange(Sample, -Start)
    }
    if (grepl("q", Arm)) {
      segments$telcent <- ifelse(segments$Start <= coord$Start, "CENT", "INTER")
      segments <- segments %>% dplyr::arrange(Sample, Start)
    }
  }
  return(segments)
}
