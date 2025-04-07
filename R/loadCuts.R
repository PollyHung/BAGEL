#' Load and Prepare Breakpoint Data
#'
#' This function reads a breakpoint file and selects specific columns for further analysis.
#' It also converts column names to lowercase for consistency.
#'
#' @param cuts.path A string representing the path to the breakpoint file (in .txt or .delimited format).
#'
#' @return A dataframe containing the selected columns: 'chr', 'arm', 'cytoband', 'start', 'end', 'direction', and 'telcent'.
#' @details
#' The function reads the data from the specified path, ensuring that only the columns of interest are retained.
#' The column names are converted to lowercase to facilitate consistent access in subsequent analyses.
#'
#' @examples
#' cuts_data <- loadCuts("path/to/cuts_file.txt")
#' @export
#' @import dplyr
#' @import readr

loadCuts <- function(cuts.path) {
  # Select interested columns
  bpdf <- read.delim(cuts.path)
  colnames(bpdf) <- tolower(colnames(bpdf))
  bpdf <- bpdf %>% dplyr::select(chr, arm, cytoband, start, end, direction, telcent)
  bpdf$arm <- ifelse(!grepl("p|q", bpdf$arm),
                     paste0(bpdf$arm, gsub("[0-9]|\\.", "", bpdf$cytoband)),
                     bpdf$arm)
  bpdf$id <- paste0(bpdf$arm, "_", bpdf$telcent, "_", bpdf$direction)

  return(bpdf)
}



