#' Make Markers from Bed file or from Gencode v47
#'
#' This function estimates arm-level copy number based on various sequencing techniques,
#' including targeted sequencing (TS), whole genome sequencing (WGS), whole exome sequencing (WES),
#' and SNP arrays. It generates a GRanges object containing the regions of interest based on the
#' specified parameters.
#'
#' @param bed_file A string representing the path to a BED file (optional). Used if the targeted
#'                 sequencing panel is not predefined.
#' @param technique A character string specifying the sequencing technique to use. Options are
#'                  "TS", "WGS", "WES", or "SNP".
#' @param panel A character string specifying the targeted gene panel to use (default is
#'              "oncoKB_annot"). It should be one of the predefined panels.
#' @param genome_build A character string specifying the genome build to use (default is "hg38").
#' @param gff.path A string representing the path to a GFF file (optional). Used for WGS and WES
#'                  techniques if a custom GFF is required.
#' @param marker.path A string representing the path to a marker file (optional). Required for
#'                    SNP array technique.
#'
#' @return A GRanges object containing the regions of interest derived from the specified
#'         sequencing technique.
#' @examples
#' markers <- makeMarkers(technique = "WGS", gff.path = "path/to/gff")
#' @export
#' @import dplyr


makeMarkers <- function(bed_file = NULL,
                        technique = "TS",
                        genome_build = "hg38",
                        gff.path = NULL,
                        marker.path = NULL) {

  message("----------Make Markers for Calculating Gene Scores-------------")
  # Validate technique input
  valid_techniques <- c("TS", "WGS", "WES", "SNP")
  if (!technique %in% valid_techniques) {
    stop("ERROR: Unsupported technique. Valid options are: ", paste(valid_techniques, collapse = ", "))
  }

  # Initialize message and regions variable
  regions <- NULL

  # Determine the technique
  if (technique == "TS") {
    message("Estimating Arm-Level Copy Number from Targeted Sequencing.")

    # Check if the panel is predefined
    if (!is.null(bed_file)) {
      message("Warning: Please Check Panel Coverage Before Performing Estimation.")
      TS_genes <- read.delim(bed_file, header = FALSE)
      colnames(TS_genes) <- c("chr", "start", "end")
      if (!all(grepl("chr", TS_genes$chr))) {
        TS_genes$chr <- paste0("chr", TS_genes$chr)
      }
      regions <- TS_genes %>% distinct() %>% as.data.frame() %>% regioneR::toGRanges()
    } else {
      stop("ERROR: Please provide a valid sequence panel name or a path to a BED file.")
    }

  } else if (technique %in% c("WGS", "WES")) {
    message(ifelse(technique == "WGS",
                   "Estimating Arm-Level Copy Number from Whole Genome Sequencing.",
                   "Estimating Arm-Level Copy Number from Whole Exome Sequencing."))

    if (is.null(gff.path)) {
      message("Using the default GFF file.")
      regions <-
      genome(regions) <- genome_build
    } else {
      message("Using Customized Version of Gencode GFF")
      gff <- ape::read.gff(gff.path)
      gff <- gff %>% filter(type == ifelse(technique == "WGS", "gene", "exon")) %>% dplyr::select(seqid, start, end)
      colnames(gff) <- c("chr", "start", "end")
      regions <- gff %>% distinct() %>% as.data.frame() %>% regioneR::toGRanges()
      genome(regions) <- genome_build
    }

  } else if (technique == "SNP") {
    message("Estimating Arm-Level Copy Number from SNP arrays")

    if (is.null(marker.path)) {
      stop("ERROR: Marker file does not exist! Please provide a valid marker file path.")
    } else {
      marker <- read.delim(marker.path, header = FALSE)
      colnames(marker) <- c("snp_id", "chr", "pos")
      marker$chr <- paste0("chr", marker$chr)
      regions <- GenomicRanges::GRanges(seqnames = marker$chr,
                         ranges = IRanges(start = marker$pos, end = marker$pos),
                         snp_id = marker$snp_id)
      genome(regions) <- genome_build
    }
  }

  return(regions)
}

