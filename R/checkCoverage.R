#' Check Coverage of Target Sequencing Markers or SNP Microarray Markers
#'
#' This function checks the coverage of genomic ranges against specified cuts and
#' visualizes the results using an ideogram. It assesses whether the genomic regions
#' have sufficient coverage based on the provided cut information and generates a
#' visual representation of the coverage across chromosomes.
#'
#' @param GRange A GRanges object containing genomic ranges to be checked for coverage.
#' @param cuts A data frame containing cut information, which includes columns for
#'             chromosome arms and corresponding start and end positions.
#' @param result_folder A string representing the directory where the coverage plot will be saved.
#' @param stringent A logical indicating whether to apply stringent criteria (default is FALSE).
#' @param plot_name A string representing the name of the plot file (without extension).
#'
#' @return A plot visualizing the coverage of genomic ranges against the specified cuts.
#' @examples
#' checkCoverage(GRange, cuts, result_folder = "path/to/results", plot_name = "coverage_plot")
#' @export
#' @import dplyr
#' @import RIdeogram


checkCoverage <- function(GRange, cuts, stringent = FALSE,
                          plot_name) {

  message("----------------------Plotting Coverages-----------------------")
  message("!!!!Please check whether your Targeted Sequencing panel includes enough genes to
          support the estimation of arm-level aneuploidy before you proceed!!!")

  # Load in Karyotype
  data(human_karyotype, package = "RIdeogram")

  # Load in the Markers
  grange <- as.data.frame(GRange) %>% dplyr::select(seqnames, start, end)
  colnames(grange) <- c("Chr", "Start", "End")
  grange$Type <- "marker"
  grange$Shape <- "circle"
  grange$color <- "FE9801"
  grange <- grange %>% dplyr::select("Type", "Shape", "Chr", "Start", "End", "color")
  grange$Chr <- gsub("chr", "", grange$Chr) %>% as.numeric()
  grange <- na.omit(grange)

  # Load in the Cuts
  biscut <- cuts
  biscut$Chr <- gsub("p|q", "", biscut$arm)
  biscut <- biscut %>% dplyr::select(Chr, smallest_start, largest_end)
  colnames(biscut) <- c("Chr", "Start", "End")

  # Create GRanges object for cuts and identify gaps
  biscut <- biscut %>% distinct() %>% as.data.frame() %>% toGRanges()
  gap <- gaps(biscut) %>% as.data.frame() %>% dplyr::select(seqnames, start, end)

  biscut <- as.data.frame(biscut) %>% dplyr::select(seqnames, start, end)
  biscut$value <- 100

  gap <- as.data.frame(gap) %>% dplyr::select(seqnames, start, end)
  gap$value <- 0

  # Combine cut regions and gaps
  biscut <- rbind(biscut, gap) %>% arrange(seqnames, start)
  colnames(biscut) <- c("Chr", "Start", "End", "Value")

  # Plot
  ideogram(karyotype = human_karyotype,
           output = paste0(plot_name, ".svg"),
           overlaid = biscut,
           label = grange,
           label_type = "marker",
           colorset1 = c("#e5f5f9", "#2ca25f"))
  convertSVG(paste0(plot_name, ".svg"), device = "png")

}






