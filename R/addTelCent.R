


addTelCent <- function(segdf = segments, direction, telcent) {
  ## Filter segments based on the specified direction (AMP or DEL)
  segdf <- segdf %>% dplyr::filter(Status == direction)

  ## Assign telomere or centromere labels based on the arm and direction and arrange dataframe accordingly
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
    }
    if (grepl("q", arm)) {
      segdf$telcent <- ifelse(segdf$Start <= coord$Start, "CENT", "INTER")
    }
  }

  return(segdf)
}
