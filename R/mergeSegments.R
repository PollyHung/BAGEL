#' Merge Segments Function
#'
#' This function merges segments based on the number of probes. Segments with fewer than a specified
#' number of probes are merged into the previous segment.
#'
#' @param seg.path A string representing the file path to a tab-delimited text file
#'                 containing the segmentation data.
#' @param min_probes An integer specifying the minimum number of probes required
#'                   to consider a segment as valid. Segments with fewer probes
#'                   than this threshold will be merged into the preceding segment.
#'
#' @details
#' The function reads the segmentation data from the specified file path, merges segments based on the probe count,
#' and returns a data frame containing the merged segments with columns 'End', 'Num_Probes', and 'Segment_Mean'.
#'
#' @return A data frame containing the merged segments, which includes the columns:
#'         'End', 'Num_Probes', and 'Segment_Mean'.
#' @export
#' @import dplyr

mergeSegments <- function(seg.path, min_probes = 10) {

  message(paste0("----------Merging Segments with min_probe threshold: ", min_probes, "-------------"))
  ## Read in segments
  seg <- read.delim(seg.path)
  colnames(seg) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
  message(paste0("There are ", table(seg$Num_Probes < min_probes)[["TRUE"]],
                 " segments with probe number below min probe threshold. ",
                 "Taking up about ", round((table(seg$Num_Probes < min_probes)[["TRUE"]])/nrow(seg) * 100, digits = 3),
                 "% of all segments"))

  ## Prepare Empty List for Result
  merged_segments <- list()
  current_segment <- NULL

  for (i in 1:nrow(seg)) {
    if (is.null(current_segment)) {
      current_segment <- seg[i, ]
    } else {
      # Check if the current segment can be merged
      if (current_segment$Num_Probes < min_probes) {
        # Merge with the previous segment
        current_segment$End <- seg[i, 'End']
        current_segment$Num_Probes <- current_segment$Num_Probes + seg[i, 'Num_Probes']
        current_segment$Segment_Mean <- mean(c(current_segment$Segment_Mean, seg[i, 'Segment_Mean']))
      } else {
        # If the current segment has enough markers, save it and start a new segment
        merged_segments <- append(merged_segments, list(current_segment))
        current_segment <- seg[i, ]
      }
    }
  }

  # Add the last segment if applicable
  if (!is.null(current_segment)) {
    merged_segments <- append(merged_segments, list(current_segment))
  }

  # Convert the list back to a data frame
  merged_segments_df <- do.call(rbind, merged_segments)
  colnames(merged_segments_df) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Log2Ratios")
  merged_segments_df <- merged_segments_df %>% dplyr::filter(Chromosome %in% c(1:22))
  message(paste0("After merging, ", nrow(merged_segments_df), " segments remain."))

  return(merged_segments_df)
}


