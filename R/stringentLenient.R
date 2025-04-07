#' Collapse Breakpoints by Arm and Telocentric Position
#'
#' This function collapses a dataframe of breakpoints based on stringent or lenient criteria.
#'
#' @param bpdf A dataframe containing breakpoint data with columns 'arm', 'telcent', 'direction', 'start', and 'end'.
#' @param stringent A logical value indicating whether to use the stringent method (TRUE) or the lenient method (FALSE).
#'
#' @return A dataframe with collapsed breakpoints according to the specified method.
#' @details
#' The stringent method selects the furthest start from the telomere end if the arm is p-arm,
#' and the closest start from the centromere if the arm is q-arm. The lenient method includes
#' all identified breakpoints.
#'
#' @examples
#' bpdf <- data.frame(arm = c("11q", "14", "16q", "17q"),
#'                    telcent = c("tel", "tel", "tel", "tel"),
#'                    direction = c("del", "del", "del", "amp"),
#'                    start = c(65684282, 36767763, 66461199, 33762114),
#'                    end = c(65686531, 36789882, 66529191, 33775856))
#'
#' stringentLenient(bpdf, stringent = TRUE)
#' stringentLenient(bpdf, stringent = FALSE)
#' @export
#' @import dplyr


stringentLenient <- function(bpdf, stringent = FALSE){
  if(isTRUE(stringent)){
    message("----------Collapsing Breakpoints by Arm and Telcent (stringent)-------------")
    bpdf <- bpdf %>%
      group_by(arm, telcent, direction, id) %>%
      summarise(smallest_start = ifelse(any(grepl("p", arm)), max(start), min(start)),
                largest_end = end[which(start == ifelse(any(grepl("p", arm)), max(start), min(start)))][1],
                .groups = 'drop') %>%
      as.data.frame()
    print(bpdf)
  } else {
    message("----------Collapsing Breakpoints by Arm and Telcent (lenient)-------------")
    bpdf <- bpdf %>%
      group_by(arm, telcent, direction, id) %>%
      summarize(smallest_start = min(start), largest_end = max(end), .groups = 'drop') %>%
      as.data.frame()
    print(bpdf)
  }
  return(bpdf)
}


