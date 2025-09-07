#' Improved checkCoverage Function
#'
#' Check Coverage with consistent parameters and better error handling
#'
#' @param GRange GRanges object with genomic ranges
#' @param cuts Dataframe with cut information
#' @param stringent Logical for stringent criteria (default=FALSE)
#' @param plot_name Character string for plot file name
#' @param result_folder Character string for output directory (default=".")
#'
#' @return Creates coverage plot files
#' @export
checkCoverage <- function(GRange, 
                          cuts, 
                          stringent = FALSE,
                          plot_name,
                          result_folder = ".") {
  
  # Input validation
  if (!is(GRange, "GRanges")) {
    stop("GRange must be a GRanges object")
  }
  
  if (!is.data.frame(cuts)) {
    stop("cuts must be a data.frame")
  }
  
  if (missing(plot_name)) {
    stop("plot_name is required")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(result_folder)) {
    dir.create(result_folder, recursive = TRUE)
  }
  
  message("----------------------Plotting Coverages-----------------------")
  message("!!!!Please check whether your Targeted Sequencing panel includes enough genes to")
  message("     support the estimation of arm-level aneuploidy before you proceed!!!")
  
  # Load in Karyotype
  data(human_karyotype, package = "RIdeogram")
  
  # Prepare markers data
  grange_df <- as.data.frame(GRange) %>% 
    dplyr::select(seqnames, start, end)
  colnames(grange_df) <- c("Chr", "Start", "End")
  
  grange_df <- grange_df %>%
    dplyr::mutate(
      Type = "marker",
      Shape = "circle",
      color = "FE9801",
      Chr = as.numeric(gsub("chr", "", Chr))
    ) %>%
    dplyr::select(Type, Shape, Chr, Start, End, color) %>%
    na.omit()
  
  # Prepare cuts data
  if ("arm" %in% names(cuts) && "smallest_start" %in% names(cuts) && "largest_end" %in% names(cuts)) {
    biscut <- cuts %>%
      dplyr::mutate(Chr = as.numeric(gsub("[pq]", "", arm))) %>%
      dplyr::select(Chr, smallest_start, largest_end) %>%
      dplyr::rename(Start = smallest_start, End = largest_end) %>%
      distinct()
  } else {
    stop("cuts dataframe must contain 'arm', 'smallest_start', and 'largest_end' columns")
  }
  
  # Create coverage visualization data
  tryCatch({
    # Create GRanges for cuts and identify gaps
    cuts_gr <- biscut %>% 
      as.data.frame() %>% 
      regioneR::toGRanges()
    
    gap_gr <- GenomicRanges::gaps(cuts_gr)
    gap_df <- as.data.frame(gap_gr) %>% 
      dplyr::select(seqnames, start, end)
    
    # Prepare data for ideogram
    biscut_plot <- as.data.frame(cuts_gr) %>% 
      dplyr::select(seqnames, start, end) %>%
      dplyr::mutate(value = 100)
    colnames(biscut_plot) <- c("Chr", "Start", "End", "Value")
    
    gap_plot <- gap_df %>%
      dplyr::mutate(value = 0)
    colnames(gap_plot) <- c("Chr", "Start", "End", "Value")
    
    combined_plot <- rbind(biscut_plot, gap_plot) %>%
      dplyr::arrange(Chr, Start)
    
    # Create ideogram
    output_path <- file.path(result_folder, paste0(plot_name, ".svg"))
    RIdeogram::ideogram(
      karyotype = human_karyotype,
      output = output_path,
      overlaid = combined_plot,
      label = grange_df,
      label_type = "marker",
      colorset1 = c("#e5f5f9", "#2ca25f")
    )
    
    # Convert to PNG
    png_path <- file.path(result_folder, paste0(plot_name, ".png"))
    RIdeogram::convertSVG(output_path, device = "png")
    
    message("Coverage plot saved to: ", png_path)
    
  }, error = function(e) {
    stop("Error creating coverage plot: ", e$message)
  })
}