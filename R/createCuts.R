#' Create Cuts from Genomic Segments
#'
#' This function processes genomic segments to create necessary cuts based on specified parameters.
#' It loads arm-level coordinates for the specified genome and filters segments accordingly.
#' Additionally, it generates background files for each lineage.
#'
#' @param segments A data frame containing genomic segments. Defaults to `segs`.
#' @param genome A string specifying the genome version (e.g., "hg38", "hg19"). Default is "hg19".
#' @param cutoff A numeric value used as a threshold for filtering segments. Default is 0.25.
#' @param result_dir A string specifying the directory to save output files. Default is "example/pooledOV/".
#' @param cytoband Optional; a path to a custom cytoband CSV file. If provided, overrides the default genome coordinates.
#'
#' @return NULL; the function performs processing and saves output files to the specified result directory.
#'
#' @examples
#' # Example usage
#' createCuts(segments = my_segments, genome = "hg38", cutoff = 0.3, result_dir = "output/")
#'
#' @export
createCuts <- function(segments = segs,
                       genome = "hg19",
                       cutoff = 0.25,
                       result_dir,
                       cytoband = NULL){

  ## Load in the correct arm level coordinate file based on the specified genome
  if (genome == "hg38") {
    coords <- BAGEL::cytoband.hg38
  } else if (genome == "hg19") {
    coords <- BAGEL::cytoband.hg19
  } else if (!is.null(cytoband)) {
    coords <- read.csv(cytoband)
    colnames(coords) <- c("Chromosome", "Arm", "Start", "End")
  } else {
    stop("No valid arm-level coordinates provided, exiting")
  }

  ## Define support functions
  filter_big_small <- function(df) {
    df <- df %>%
      dplyr::filter(percent >= 0.001) %>%
      dplyr::filter(percent <= 0.999) %>%
      dplyr::arrange(percent)
    return(df)
  }

  ## Clean up the segments
  segments <- na.omit(segments)
  segments <- segments %>% dplyr::filter(Arm != "centromere") ## remove segments that are defined in centromere region

  # Call preprocessSeg with ALL parameters
  chr_arms <- unique(segments$Arm)
  parallel::mclapply(chr_arms, function(i) {
    preprocessSeg(
      Segments = segments,   # Explicitly pass
      Coords = coords,       # from createCuts
      Cutoff = cutoff,       # parameters
      Result_dir = result_dir,
      arm = i
    )
  }, mc.cores = parallel::detectCores() - 4)

  background(rd = result_dir)
}






