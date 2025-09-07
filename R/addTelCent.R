#' Improved addTelCent Function
#'
#' Add Telomere or Centromere Labels to Segments with proper validation
#'
#' @param segments Dataframe with segment data
#' @param direction Character string: "AMP" or "DEL"
#' @param telcent Character string: "TEL" or "CENT"
#' @param Arm Character string indicating chromosome arm
#' @param coord List or dataframe with Start and End coordinates
#' @param fuzzy Logical for fuzzy matching (default=FALSE)
#'
#' @return Dataframe with telcent column added
#' @export
addTelCent <- function(segments,
                       direction,
                       telcent,
                       Arm,
                       coord,
                       fuzzy = FALSE) {
  
  # Input validation
  required_cols <- c("Status", "Start", "End", "Sample")
  if (!all(required_cols %in% names(segments))) {
    stop("Segments dataframe missing required columns: ",
         paste(setdiff(required_cols, names(segments)), collapse = ", "))
  }
  
  if (!direction %in% c("AMP", "DEL")) {
    stop("direction must be 'AMP' or 'DEL'")
  }
  
  if (!telcent %in% c("TEL", "CENT")) {
    stop("telcent must be 'TEL' or 'CENT'")
  }
  
  # Handle different coord formats
  if (is.data.frame(coord)) {
    coord_start <- coord$Start[1]
    coord_end <- coord$End[1]
  } else if (is.list(coord) && all(c("Start", "End") %in% names(coord))) {
    coord_start <- coord$Start
    coord_end <- coord$End
  } else if (is.numeric(coord) && length(coord) == 2) {
    coord_start <- coord[1]
    coord_end <- coord[2]
  } else {
    stop("coord must be a data.frame with Start/End columns, a list with Start/End elements, or a numeric vector of length 2")
  }
  
  # Filter segments by direction
  filtered_segments <- segments %>% 
    dplyr::filter(Status == direction)
  
  if (nrow(filtered_segments) == 0) {
    warning("No segments found with Status = ", direction)
    return(segments %>% dplyr::mutate(telcent = NA_character_))
  }
  
  # Add telcent labels based on arm type and telcent specification
  if (telcent == "TEL") {
    if (grepl("p$", Arm)) {
      # p-arm: telomere is at the start
      filtered_segments$telcent <- ifelse(
        filtered_segments$Start <= coord_start, 
        "TEL", 
        "INTER"
      )
      filtered_segments <- filtered_segments %>% 
        dplyr::arrange(Sample, Start)
    } else if (grepl("q$", Arm)) {
      # q-arm: telomere is at the end
      filtered_segments$telcent <- ifelse(
        filtered_segments$End >= coord_end, 
        "TEL", 
        "INTER"
      )
      filtered_segments <- filtered_segments %>% 
        dplyr::arrange(Sample, desc(Start))
    } else {
      stop("Arm must end with 'p' or 'q'")
    }
  } else if (telcent == "CENT") {
    if (grepl("p$", Arm)) {
      # p-arm: centromere is at the end
      filtered_segments$telcent <- ifelse(
        filtered_segments$End >= coord_end, 
        "CENT", 
        "INTER"
      )
      filtered_segments <- filtered_segments %>% 
        dplyr::arrange(Sample, desc(Start))
    } else if (grepl("q$", Arm)) {
      # q-arm: centromere is at the start
      filtered_segments$telcent <- ifelse(
        filtered_segments$Start <= coord_start, 
        "CENT", 
        "INTER"
      )
      filtered_segments <- filtered_segments %>% 
        dplyr::arrange(Sample, Start)
    } else {
      stop("Arm must end with 'p' or 'q'")
    }
  }
  
  return(filtered_segments)
}