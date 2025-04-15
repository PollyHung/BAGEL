#' Preprocess Segments Function
#'
#' This function preprocesses genomic segments, applying filtering and labeling based on specified cutoffs.
#' It also identifies segments originating from telomeres or centromeres and joins overlapping segments.
#'
#' @param segments A data frame containing genomic segments with columns including Arm and Segment_Mean.
#' @param genome A string specifying the genome version (e.g., "hg38", "hg19"). Default is "hg38".
#' @param arm A string indicating the chromosome arm to be processed (e.g., "1p", "1q").
#' @param cytoband Optional; a path to a custom cytoband CSV file.
#' @param del_cutoff Numeric; threshold for labeling segments as deletions. Default is -0.25.
#' @param amp_cutoff Numeric; threshold for labeling segments as amplifications. Default is 0.25.
#' @param cutoff_SD Numeric; standard deviation used to adjust cutoffs based on segment means.
#' @param result_dir A string specifying the directory to save output files.
#' @return None; output files are saved to the specified result directory.
#' @export
preprocessSeg <- function(segments,
                          genome = "hg19",
                          arm,
                          cytoband = NULL,
                          del_cutoff = -0.25,
                          amp_cutoff = 0.25,
                          cutoff_SD = NULL, ## Try 0.25 if you want to use standard deviation based cutoff
                          result_dir) {

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

  ## Create the result directory if it does not exist
  if (!dir.exists(file.path(result_dir, "breakpoints"))) {
    message("Creating Breakpoint File Directory `/breakpoints`")
    dir.create(file.path(result_dir, "breakpoints"))
  }

  ## Filter segments to only include those from the specified chromosome arm
  segments <- segments %>% dplyr::filter(Arm == arm)
  coord <- coords %>% dplyr::filter(Arm == arm)

  ## Calculate the length of the chromosome arm
  if (grepl("p", arm)) {
    arm_length <- coord$End - 0
  } else {
    arm_length <- coord$End - coord$Start
  }

  ## Adjust cutoffs based on provided standard deviation if applicable
  if (!is.null(cutoff_SD)) {
    message("You provided cutoff standard deviation, thus overriding hard threshold with customized cutoff based on standard deviation")
    segMean <- mean(segments$Segment_Mean)
    segSD <- sd(segments$Segment_Mean)
    amp_cutoff <- segMean + cutoff_SD * segSD
    del_cutoff <- segMean - cutoff_SD * segSD
  }

  ## Label segments as AMP, DEL, or NEUTRAL based on calculated cutoffs
  segments <- segments %>% dplyr::mutate(Status = case_when(Segment_Mean <= del_cutoff ~ "DEL",
                                                            Segment_Mean >= amp_cutoff ~ "AMP",
                                                            TRUE ~ "NEUTRAL"))
  message(paste(names(table(segments$Status)), table(segments$Status), collapse = ", "))

  ## Apply the TELCENT function for both AMP and DEL segments
  seg_list <- list()
  seg_list[["amp_tel"]] <- addTelCent(direction = "AMP", telcent = "TEL")
  seg_list[["del_tel"]] <- addTelCent(direction = "DEL", telcent = "TEL")
  seg_list[["amp_cent"]] <- addTelCent(direction = "AMP", telcent = "CENT")
  seg_list[["del_cent"]] <- addTelCent(direction = "DEL", telcent = "CENT")

  ## Join overlapping segments using the joinSegs_merge function
  seg_list$amp_tel <- joinSegs_merge(segdf = seg_list$amp_tel, TELCENT = "TEL")
  seg_list$del_tel <- joinSegs_merge(segdf = seg_list$del_tel, TELCENT = "TEL")
  seg_list$amp_cent <- joinSegs_merge(segdf = seg_list$amp_cent, TELCENT = "CENT")
  seg_list$del_cent <- joinSegs_merge(segdf = seg_list$del_cent, TELCENT = "CENT")

  ## Save processed segment data to specified output files
  write.table(seg_list$amp_tel, file.path(result_dir, "breakpoints", paste0(arm, "_amp_tel.txt")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(seg_list$del_tel, file.path(result_dir, "breakpoints", paste0(arm, "_del_tel.txt")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(seg_list$amp_cent, file.path(result_dir, "breakpoints", paste0(arm, "_amp_cent.txt")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(seg_list$del_cent, file.path(result_dir, "breakpoints", paste0(arm, "_del_cent.txt")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}
