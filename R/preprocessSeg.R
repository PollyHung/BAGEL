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
preprocessSeg <- function(Segments,
                          Coords,
                          Cutoff,
                          Result_dir,
                          arm) {

  message(paste0("Processing Chromosome Arm ", arm))

  ## Create the Result Breakpoint Folder if not already created
  if (!dir.exists(file.path(Result_dir, "breakpoints"))) {
    dir.create(file.path(Result_dir, "breakpoints"))
  }
  message(paste0("Files will be stored to ", file.path(Result_dir, "breakpoints")))

  ## Slice the segments in that chromosome arm and slice up the corresponding coordinates
  segdf <- Segments %>% dplyr::filter(Arm == arm)
  coord <- Coords %>% dplyr::filter(Arm == arm)
  coordinates <- unname(unlist(coord[, 3:4]))

  ## Calculate the chromosome arm lengths
  if (grepl("p", arm)) {
    arm_length <- coord$End - 0
  } else {
    arm_length <- coord$End - coord$Start
  }
  message(paste0("Chromosome Arm ", arm, " is ", arm_length, " base pair long"))

  ## Add Amplification and Deletions
  segdf <- segdf %>% dplyr::mutate(Status = case_when(Segment_Mean <= -Cutoff ~ "DEL",
                                                      Segment_Mean >= Cutoff ~ "AMP",
                                                      TRUE ~ "NEUTRAL"))
  if(table(segdf$Status)["AMP"] > table(segdf$Status)["DEL"]){
    message(paste0("There are more amplifications in chromosome arm ", arm))
  } else {
    message(paste0("There are more deletions in chromosome arm ", arm))
  }

  # Pass coord to addTelCent
  seg_list <- list()
  seg_list[["amp_tel"]] <- addTelCent(segments = segdf, direction = "AMP", telcent = "TEL", Arm = arm, coord = coord)
  seg_list[["del_tel"]] <- addTelCent(segments = segdf, direction = "DEL", telcent = "TEL", Arm = arm, coord = coord)
  seg_list[["amp_cent"]] <- addTelCent(segments = segdf, direction = "AMP", telcent = "CENT", Arm = arm, coord = coord)
  seg_list[["del_cent"]] <- addTelCent(segments = segdf, direction = "DEL", telcent = "CENT", Arm = arm, coord = coord)

  # Pass all parameters to joinSegs
  amp_tel <- joinSegs(segments = seg_list$amp_tel, aneu = "AMP", coord = coordinates, ncutoff = Cutoff, telcent = "TEL", TELCENT = "TEL")
  del_tel <- joinSegs(segments = seg_list$del_tel, aneu = "DEL", coord = coordinates, ncutoff = Cutoff, telcent = "TEL", TELCENT = "TEL")
  amp_cent <- joinSegs(segments = seg_list$amp_cent, aneu = "AMP", coord = coordinates, ncutoff = Cutoff, telcent = "CENT", TELCENT = "CENT")
  del_cent <- joinSegs(segments = seg_list$del_cent, aneu = "DEL", coord = coordinates, ncutoff = Cutoff, telcent = "CENT", TELCENT = "CENT")

  ## Save processed segment data to specified output files
  write.table(amp_tel, file.path(Result_dir, "breakpoints", paste0(arm, "_amp_tel.txt")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(del_tel, file.path(Result_dir, "breakpoints", paste0(arm, "_del_tel.txt")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(amp_cent, file.path(Result_dir, "breakpoints", paste0(arm, "_amp_cent.txt")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(del_cent, file.path(Result_dir, "breakpoints", paste0(arm, "_del_cent.txt")), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}


