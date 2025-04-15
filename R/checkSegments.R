#' Check and Clean Segmentation Data
#'
#' This function checks a segmentation DataFrame for several conditions:
#' - Ensures it has exactly 6 columns.
#' - Cleans the Chromosome column by removing "chr" or "Chr" prefixes and converts it to numeric.
#' - Filters Chromosome values to include only integers from 1 to 22.
#' - Converts the Start, End, Num_Probes, and Segment_Mean columns to numeric.
#' - Checks that the End column is greater than the Start column and reports if not.
#' - Identifies and reports outliers in the Segment_Mean column.
#'
#' @param segments A data.frame containing segmentation data with six columns:
#'   Sample, Chromosome, Start, End, Num_Probes, Segment_Mean.
#'
#' @return NULL. The function prints error messages or outliers as needed.
#'
#' @examples
#' # Example usage:
#' checkSegments(segments)
#'
#' @export

checkSegments <- function(segments) {
  # Check if the DataFrame has 6 columns
  if (ncol(segments) != 6) {
    stop("ERROR: The data does not have 6 columns.")
  }

  colnames(segments) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")

  # Check and clean the second column (Chromosome)
  segments$Chromosome <- gsub("^chr|^Chr", "", segments$Chromosome)  # Remove "chr" or "Chr" prefixes
  segments$Chromosome <- as.integer(segments$Chromosome)  # Convert to integer

  # Check if Chromosome is in range 1 to 22
  segments <- segments[segments$Chromosome %in% 1:22, ]

  # Check if columns 3 to 6 are numeric
  for (i in 3:6) {
    segments[, i] <- as.numeric(segments[, i])  # Convert to numeric
  }

  # Check if column 4 (End) is greater than column 3 (Start)
  if (any(segments$End < segments$Start)) {
    stop("ERROR: Reversed segments, start larger than end. Are you sure?")
  }
  if (any(segments$End == segments$Start)){
    warning(paste0("There are ", sum(segments$End == segments$Start),
                   " segments with length of 0. They will be removed"))
    segments <- segments[-which(segments$End == segments$Start), ]
  }

  # Check for outliers in the Segment_Mean (6th column)
  outliers <- which(segments$Segment_Mean <= -20 | segments$Segment_Mean >= 20)
  if(length(outliers) != 0){
    warning("There are ", length(outliers), " segments with Segment Mean outside of -20 and 20 range.")
    summary(segments[outliers, ]$Segment_Mean)
  }

  ## Change column names
  colnames(segments) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")

  return(segments)
}

