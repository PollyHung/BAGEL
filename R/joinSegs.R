#' Join Segments Function
#'
#' This function processes a data frame of segments, merging them based on specified criteria.
#' It looks for TEL or CENT segments and merges them accordingly.
#'
#' @param segdf A data frame containing segment information, including Sample, Start, End, and other numeric columns.
#' @param TELCENT A string indicating segments to be considered for merging (e.g., "TEL" or "CENT").
#' @return A data frame with merged segment information per sample, including Percent, Start, and End.
#' @export


joinSegs <- function(segdf, TELCENT) {
  results <- list()

  ## Join Segments If they are TELOMERE or CENTROMERE bounded
  joinSegs_telcent <- function(df){

    if(nrow(df) == 0) {return(df)} ## If there are no rows, return emtpy dataframe

    if(sum(grepl(TELCENT, df$telcent)) == 0){
      message("No Telomere or Centromere bound header segment found. This is inter-sCNA")
      ## Write function to join the inter segments
    }

    if(sum(grepl(TELCENT, df$telcent)) > 0){

    }

  }







  joinSegs_telcent <- function(df) {
    # if (nrow(df) == 0) {
    #   message("Empty Dataframe")
    #   return(df)
    # }

    ## Store the first row for later use
    original_first_row <- df[1, , drop = FALSE]
    if (nrow(df) == 1) {
      message("Only One Row Detected")
      return(original_first_row)
    }

    ## Initialize the current row as the first row
    current_row <- df[1, , drop = FALSE]
    for (i in 2:nrow(df)) {
      next_row <- df[i, , drop = FALSE]

      ## Calculate the gap between current and next segment
      gap <- next_row$Start - current_row$End

      ## Calculate the total length of the merged segments
      sum_length <- current_row$seqLength + next_row$seqLength

      ## If segments overlap or touch, merge them
      if (sum_length > gap) {
        # Merge rows by updating the current row's End and seqLength
        merged_row <- current_row
        merged_row$End <- next_row$End
        merged_row$seqLength <- sum_length

        # Handle numeric columns by summing their values
        numeric_cols <- names(df)[sapply(df, is.numeric)]
        numeric_cols <- setdiff(numeric_cols, c("Start", "End", "seqLength"))

        for (col in numeric_cols) {
          merged_row[[col]] <- current_row[[col]] + next_row[[col]]
        }

        # Update the current row to the merged row
        current_row <- merged_row
      } else {
        ## Return the original first row if no merge occurs
        return(original_first_row)
      }
    }

    ## Return the final merged row
    return(current_row)
  }

  ## Process each unique sample in the data frame
  for (sample in unique(segdf$Sample)) {

    ## Slice the data frame for the current sample and calculate sequenced length
    segdf_slice <- segdf %>%
      dplyr::filter(Sample == sample) %>%
      dplyr::mutate(seqLength = End - Start) %>%
      dplyr::select(Sample, Start, End, seqLength, Num_Probes, Segment_Mean, telcent)

    ## Check for presence of TEL or CENT segments
    if (sum(grepl(TELCENT, segdf_slice$telcent)) > 0) {
      segdf_slice <- joinSegs(segdf_slice)
      results[[sample]] <- data.frame(Sample = sample,
                                      Percent = segdf_slice$seqLength / arm_length,
                                      Start = segdf_slice$Start / arm_length,
                                      End = segdf_slice$End / arm_length)
    } else {
      ## If no TEL or CENT segments, return zeros for Percent, Start, and End
      results[[sample]] <- data.frame(Sample = sample, Percent = 0, Start = 0, End = 0)
    }
  }

  ## Combine results into a single data frame and return
  return(do.call(rbind, results))
}
