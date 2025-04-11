#' Annotate Segments with Chromosome Arms
#'
#' This function annotates segments with their corresponding chromosome arms
#' based on the provided arm information. It checks for overlaps between
#' segment start and end positions and the defined arms for each chromosome.
#' If a segment crosses the centromere, it is split into two daughter segments:
#' one for the p-arm and one for the q-arm. These daughter segments retain the
#' same log2 ratios, and the number of probes is split based on the proportion
#' of the original segment that lies in each arm.
#'
#' @param segments A data frame containing segment information with columns:
#'   - Sample: Identifier for the sample
#'   - Chromosome: Chromosome number
#'   - Start: Start position of the segment
#'   - End: End position of the segment
#'   - Num_Probes: Number of probes in the segment
#'   - Segment_Mean: Mean value for the segment
#'   - Arm: An initialized column that will be filled with the corresponding arm
#'
#' @param arms A data frame containing chromosome arm information with columns:
#'   - Chromosome: Chromosome number
#'   - Arm: Arm identifier (e.g., "1p", "1q")
#'   - Start: Start position of the chromosome arm
#'   - End: End position of the chromosome arm
#'
#' @return A data frame with the original segments, now including the
#'   corresponding Arm for each segment. Segments that cross the centromere
#'   are split into separate p-arm and q-arm segments, with probes distributed
#'   according to their length in each arm.
#'
#' @import dplyr
#' @export
annotateArms <- function(segments, genome = "hg38", cytoband = NULL) {

  message("Segmenting the Chromosome Arms by Cytoband, excluding centromeres")

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

  # Preprocess arms into p and q for each chromosome
  arms_p <- arms[grep("p$", arms$Arm), ]
  arms_q <- arms[grep("q$", arms$Arm), ]
  arms_pq <- merge(
    arms_p, arms_q,
    by = "Chromosome",
    suffixes = c("_p", "_q")
  )

  # Initialize result list
  result <- list()

  for (i in 1:nrow(segments)) {
    seg <- segments[i, ]
    chrom <- as.character(seg$Chromosome)
    seg_start <- seg$Start
    seg_end <- seg$End

    # Check if chromosome has p/q arms and segment spans both
    if (chrom %in% arms_pq$Chromosome) {
      pq <- arms_pq[arms_pq$Chromosome == chrom, ]
      p_end <- pq$End_p
      q_start <- pq$Start_q

      if (seg_start <= p_end && seg_end >= q_start) {
        # Split into p/q segments
        p_part <- seg
        p_part$Start <- seg_start
        p_part$End <- p_end
        p_part$Arm <- pq$Arm_p

        q_part <- seg
        q_part$Start <- q_start
        q_part$End <- seg_end
        q_part$Arm <- pq$Arm_q

        # Calculate probe distribution
        p_length <- p_end - seg_start + 1
        q_length <- seg_end - q_start + 1
        total_length <- p_length + q_length
        p_prop <- p_length / total_length

        p_part$Num_Probes <- round(seg$Num_Probes * p_prop)
        q_part$Num_Probes <- seg$Num_Probes - p_part$Num_Probes

        # Add to result
        result[[length(result) + 1]] <- p_part
        result[[length(result) + 1]] <- q_part
        next  # Skip further processing
      }
    }

    # If no split needed, find best arm
    max_overlap <- 0
    best_arm <- NA
    chrom_arms <- arms[arms$Chromosome == chrom, ]

    for (j in 1:nrow(chrom_arms)) {
      arm_start <- chrom_arms$Start[j]
      arm_end <- chrom_arms$End[j]

      overlap_start <- max(seg_start, arm_start)
      overlap_end <- min(seg_end, arm_end)
      overlap_length <- if (overlap_start > overlap_end) 0 else overlap_end - overlap_start + 1

      if (overlap_length > max_overlap) {
        max_overlap <- overlap_length
        best_arm <- chrom_arms$Arm[j]
      }
    }

    seg$Arm <- best_arm
    result[[length(result) + 1]] <- seg
  }

  return(do.call(rbind, result))
}




