#' Add Telomere or Centromere Labels to Segments
#'
#' This function filters segments based on the specified direction (AMP or DEL)
#' and assigns telomere or centromere labels based on the arm and direction.
#' The resulting dataframe is arranged according to the specified criteria.
#'
#' @param segdf A dataframe containing segment data. It should include columns
#'   for `Status`, `Start`, `End`, and `Sample`.
#' @param direction A character string specifying the direction. It must be either
#'   "AMP" or "DEL".
#' @param telcent A character string that specifies whether to label segments
#'   as "TEL" (telomere) or "CENT" (centromere).
#'
#' @return A dataframe that includes an additional column `telcent` indicating
#'   the assigned label for each segment.
#'
#' @examples
#' # Example usage
#' result <- addTelCent(segments, direction = "AMP", telcent = "TEL")
#'
#' @export
addTelCent <- function(segdf = segments, direction, telcent) {
  ## Filter segments based on the specified direction (AMP or DEL)
  segdf <- segdf %>% dplyr::filter(Status == direction)

  ## Assign telomere or centromere labels based on the arm and direction
  ## Arrange dataframe accordingly
  if (telcent == "TEL") {
    if (grepl("p", arm)) {
      segdf$telcent <- ifelse(segdf$Start <= coord$Start, "TEL", "INTER")
      segdf <- segdf %>% dplyr::arrange(Sample, Start)
    }
    if (grepl("q", arm)) {
      segdf$telcent <- ifelse(segdf$End >= coord$End, "TEL", "INTER")
      segdf <- segdf %>% dplyr::arrange(Sample, -Start)
    }
  } else if (telcent == "CENT") {
    if (grepl("p", arm)) {
      segdf$telcent <- ifelse(segdf$End >= coord$End, "CENT", "INTER")
      segdf <- segdf %>% dplyr::arrange(Sample, -Start)
    }
    if (grepl("q", arm)) {
      segdf$telcent <- ifelse(segdf$Start <= coord$Start, "CENT", "INTER")
      segdf <- segdf %>% dplyr::arrange(Sample, Start)
    }
  }

  return(segdf)
}
