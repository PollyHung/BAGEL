#' Preprocess Segments Function
#'
#' This function preprocesses genomic segments by applying filtering and labeling based on specified cutoffs.
#' It identifies segments originating from telomeres or centromeres and joins overlapping segments.
#'
#' @param segments A data frame containing genomic segments with columns including `Arm` and `Segment_Mean`.
#' @param genome A string specifying the genome version (e.g., "hg38", "hg19"). Default is "hg38".
#' @param arm A string indicating the chromosome arm to be processed (e.g., "1p", "1q").
#' @param cytoband Optional; a path to a custom cytoband CSV file.
#' @param del_cutoff Numeric; threshold for labeling segments as deletions. Default is -0.25.
#' @param amp_cutoff Numeric; threshold for labeling segments as amplifications. Default is 0.25.
#' @param cutoff_SD Numeric; standard deviation used to adjust cutoffs based on segment means.
#' @param result_dir A string specifying the directory to save output files.
#'
#' @return NULL; output files are saved to the specified result directory.
#'   The files include segments categorized as amplification and deletion
#'   along with their telomere and centromere assignments.
#'
#' @export
preprocessSeg <- function(arm) {

  ## Create the result directory if it does not exist
  if (!dir.exists(file.path(result_dir, "breakpoints"))) {
    message("Creating Breakpoint File Directory `/breakpoints`")
    dir.create(file.path(result_dir, "breakpoints"))
  }

  ## Filter segments to only include those from the specified chromosome arm
  segdf <- segments %>% dplyr::filter(Arm == arm)
  coord <- coords %>% dplyr::filter(Arm == arm)
  coordinates <- unname(unlist(coord[, 3:4]))

  ## Calculate the length of the chromosome arm
  if (grepl("p", arm)) {
    arm_length <- coord$End - 0
  } else {
    arm_length <- coord$End - coord$Start
  }

  ## Label segments as AMP, DEL, or NEUTRAL based on calculated cutoffs
  segdf <- segdf %>% dplyr::mutate(Status = case_when(Segment_Mean <= -cutoff ~ "DEL",
                                                      Segment_Mean >= cutoff ~ "AMP",
                                                      TRUE ~ "NEUTRAL"))
  message(paste(names(table(segdf$Status)), table(segdf$Status), collapse = ", "))

  ## Apply the TELCENT function for both AMP and DEL segments
  seg_list <- list()
  seg_list[["amp_tel"]] <- addTelCent(direction = "AMP", telcent = "TEL")
  seg_list[["del_tel"]] <- addTelCent(direction = "DEL", telcent = "TEL")
  seg_list[["amp_cent"]] <- addTelCent(direction = "AMP", telcent = "CENT")
  seg_list[["del_cent"]] <- addTelCent(direction = "DEL", telcent = "CENT")

  ## Join overlapping segments using the joinSegs_merge function
  amp_tel <- joinSegs(segdf = seg_list$amp_tel, aneu = "AMP", telcent = "TEL", TELCENT = "TEL")
  del_tel <- joinSegs(segdf = seg_list$del_tel, aneu = "DEL", telcent = "TEL", TELCENT = "TEL")
  amp_cent <- joinSegs(segdf = seg_list$del_tel, aneu = "AMP", telcent = "CENT", TELCENT = "CENT")
  del_cent <- joinSegs(segdf = seg_list$del_tel, aneu = "DEL", telcent = "CENT", TELCENT = "CENT")

  ## Save processed segment data to specified output files
  write.table(amp_tel, file.path(result_dir, "breakpoints", paste0(arm, "_amp_tel.txt")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(del_tel, file.path(result_dir, "breakpoints", paste0(arm, "_del_tel.txt")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(amp_cent, file.path(result_dir, "breakpoints", paste0(arm, "_amp_cent.txt")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(del_cent, file.path(result_dir, "breakpoints", paste0(arm, "_del_cent.txt")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}
