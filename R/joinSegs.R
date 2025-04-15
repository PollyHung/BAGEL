#' Join Segments Based on Alteration Type
#'
#' This function processes genomic segments based on their alteration type (e.g., amplification or deletion),
#' joining adjacent segments if the gap between them is less than their combined lengths. It also checks
#' for telomere and centromere alignment based on specific coordinates.
#'
#' @param segments A dataframe containing segment information with columns for `Sample`, `Chromosome`,
#'   `Start`, `End`, `Segment_Mean`, and other relevant attributes.
#' @param aneu A character string indicating the alteration type. Options are "amp", "del", "non-del", and "non-amp".
#' @param coord A numeric vector of length 4 representing the coordinates of the chromosome arms.
#' @param ncutoff A numeric cutoff value used to classify segments based on their `Segment_Mean`.
#' @param telcent A character string indicating whether to check for 'telomere' or 'centromere' alignment.
#' @param integer A logical value indicating whether to return integer lengths or normalized percentages.
#' @param TELCENT An additional parameter related to telomere/centromere classification (specific usage not defined).
#'
#' @return A dataframe summarizing the joined segments, including columns for `Sample`, `percent`,
#'   `start`, and `end`.
#'
#' @examples
#' # Example usage
#' result <- joinSegs(segments, aneu = "amp", coord = c(1, 100, 200, 300),
#'                    ncutoff = 0.5, telcent = "tel", integer = FALSE, TELCENT)
#'
#' @export
joinSegs <- function(segments,
                     aneu,
                     coord,          # No default
                     ncutoff,        # No default
                     telcent,
                     integer = FALSE,
                     TELCENT,
                     Arm = arm) {
  # Determine direction based on alteration type
  direction <- switch(aneu,
                      "amp" = 1,
                      "del" = -1,
                      "non-del" = -2,
                      "non-amp" = 2,
                      0)

  # Calculate arm length based on coordinates
  if(length(coord) == 2) {
    arm_length <- coord[2] - coord[1] + 1
  } else {
    arm_length <- (coord[2] - coord[1] + 1) + (coord[4] - coord[3] + 1)
  }

  # Preprocess segments
  preprocess_seg <- function(segments, Arm, coord) {
    segments %>%
      dplyr::filter(Chromosome == sub("[pq]$", "", Arm)) %>%
      dplyr::mutate(Start = pmax(Start, coord[1]),
                    End = pmin(End, if(length(coord) > 2) coord[4] else coord[2]))
  }

  processed_seg <- preprocess_seg(segments, Arm, coord)

  results <- list()

  # Main processing loop
  for(sample in unique(processed_seg$Sample)) {
    sample_seg <- processed_seg %>% dplyr::filter(Sample == sample)

    # Classify segments
    classified_seg <- sample_seg %>%
      dplyr::mutate(
        status = case_when(
          aneu == "del" & Segment_Mean <= -ncutoff ~ direction,
          aneu == "amp" & Segment_Mean >= ncutoff ~ direction,
          aneu == "non-del" & Segment_Mean > -ncutoff ~ direction,
          aneu == "non-amp" & Segment_Mean < ncutoff ~ direction,
          TRUE ~ 0
        )
      ) %>%
      dplyr::select(Start, End, status)

    # Filter relevant segments
    alt_segments <- classified_seg %>% dplyr::filter(status == direction)

    # Early exit conditions
    if(nrow(alt_segments) == 0) {
      results[[sample]] <- data.frame(percent = 0, start = 0, end = 0)
      next
    }

    # Check telomere/centromere alignment
    arm_type <- ifelse(grepl("p$", Arm), "p", "q")
    if((arm_type == "p" && telcent == "tel" && alt_segments$Start[1] != coord[1]) ||
       (arm_type == "p" && telcent == "cent" && alt_segments$End[nrow(alt_segments)] != coord[2]) ||
       (arm_type == "q" && telcent == "tel" && alt_segments$End[nrow(alt_segments)] != coord[2]) ||
       (arm_type == "q" && telcent == "cent" && alt_segments$Start[1] != coord[1])) {
      results[[sample]] <- data.frame(percent = 0, start = 0, end = 0)
      next
    }

    # Segment joining logic
    join_segments <- function(segments) {
      if(nrow(segments) == 1) {
        return(list(
          length = segments$End[1] - segments$Start[1] + 1,
          start = segments$Start[1],
          end = segments$End[1]
        ))
      }

      # Calculate gap lengths
      calc_gaps <- function(s) {
        sum(s$Start[-1] - s$End[-nrow(s)] - 1)
      }

      current_segments <- segments
      repeat {
        gaps <- calc_gaps(current_segments)
        alt_length <- sum(current_segments$End - current_segments$Start + 1)

        if(gaps <= alt_length) break

        # Trim from appropriate end based on arm and telcent
        if((arm_type == "p" && telcent == "tel") || (arm_type == "q" && telcent == "cent")) {
          current_segments <- current_segments[-nrow(current_segments), ]
        } else {
          current_segments <- current_segments[-1, ]
        }

        if(nrow(current_segments) == 0) {
          return(list(percent = 0, start = 0, end = 0))
        }
      }

      list(
        length = current_segments$End[nrow(current_segments)] - current_segments$Start[1] + 1,
        start = current_segments$Start[1],
        end = current_segments$End[nrow(current_segments)]
      )
    }

    joined <- join_segments(alt_segments)

    # Convert to percentages if needed
    if(!integer) {
      start_norm <- (joined$start - coord[1]) / arm_length
      end_norm <- (joined$end - coord[1]) / arm_length
      percent_norm <- joined$length / arm_length

      # Handle coordinate flipping for certain cases
      if(!((arm_type == "p" && telcent == "tel") || (arm_type == "q" && telcent == "cent"))) {
        start_norm <- 1 - end_norm
        end_norm <- 1 - (joined$start - coord[1]) / arm_length
      }

      results[[sample]] <- data.frame(
        percent = pmax(0, pmin(1, percent_norm)),
        start = pmax(0, pmin(1, start_norm)),
        end = pmax(0, pmin(1, end_norm)))
    } else {
      results[[sample]] <- data.frame(
        percent = joined$length,
        start = joined$start,
        end = joined$end
      )
    }
  }

  # Combine results and return
  result_df <- dplyr::bind_rows(results, .id = "Sample")
  return(result_df)
}
