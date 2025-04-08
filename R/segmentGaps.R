#' Calculate Segment Gaps
#'
#' This function calculates the total gaps between segments for each sample,
#' as well as gaps grouped by chromosome and by arm. Gaps in segmentation files
#' can arise from technical issues, such as low probe density, or biological
#' factors where regions may have normal copy numbers. The variability in gap sizes
#' can impact downstream analyses, particularly when comparing segment means or
#' aberrations. This function helps assess the coverage differences, which may
#' warrant correction to ensure reliable comparisons.
#'
#' @param segments A data frame containing segment information with columns:
#'   - Sample: Identifier for the sample
#'   - Chromosome: Chromosome number
#'   - Arm: Chromosome arm
#'   - Start: Start position of the segment
#'   - End: End position of the segment
#'
#' @return A data frame containing:
#'   - gap_length: Total gaps for each sample
#'   - sequenced_seg_length: Total length of sequenced segments
#'   - total_seg_length: Combined length of gaps and sequenced segments
#'   - sequence_coverage: Proportion of the arm length that is covered
#'   - coverage: Proportion of sequenced length to total segment length
#'
#' @import dplyr
#' @export
segmentGaps <- function(segments,
                        genome = "hg38"){

  # read in cytoband
  if(genome == "hg38"){
    arms <- BAGEL::cytoband.hg38
  } else if (genome == "hg19") {
    arms <- BAGEL::cytoband.hg19
  } else if (!is.null(cytoband)) {
    arms <- cytoband
  } else {
    stop ("Please Provide a Cytoband to Proceed!!!!!")
  }

  ## Group by Chromosome Arms to calculate gap length, sequenced_seg_length, and total_seg_length
  gaps <- segments %>%
    arrange(Start) %>%
    group_by(Sample, Chromosome, Arm) %>%
    mutate(Previous_End = lag(End),                      # Get end of the previous segment
           Gap = ifelse(!is.na(Previous_End), Start - Previous_End, NA)) %>%  # Calculate the gap
    summarise(gap_length = sum(Gap, na.rm = TRUE),
              sequenced_seg_length = sum(End - Start + 1),
              total_seg_length = gap_length + sequenced_seg_length) %>%
    as.data.frame()
  arms$arm_length <- arms$End - arms$Start
  arms <- arms %>% dplyr::select(Chromosome, Arm, arm_length)
  gaps <- merge(gaps, arms, by = c("Chromosome", "Arm"))

  ## Calculate Coverage Proportion
  gaps$sequence_coverage <- gaps$total_seg_length/gaps$arm_length

  ## Calculate coverage
  gaps$coverage <- gaps$sequenced_seg_length/gaps$total_seg_length

  return(gaps)
}







