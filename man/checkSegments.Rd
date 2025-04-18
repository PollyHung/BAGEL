% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkSegments.R
\name{checkSegments}
\alias{checkSegments}
\title{Check and Clean Segmentation Data}
\usage{
checkSegments(segments)
}
\arguments{
\item{segments}{A data.frame containing segmentation data with six columns:
Sample, Chromosome, Start, End, Num_Probes, Segment_Mean.}
}
\value{
NULL. The function prints error messages or outliers as needed.
}
\description{
This function checks a segmentation DataFrame for several conditions:
\itemize{
\item Ensures it has exactly 6 columns.
\item Cleans the Chromosome column by removing "chr" or "Chr" prefixes and converts it to numeric.
\item Filters Chromosome values to include only integers from 1 to 22.
\item Converts the Start, End, Num_Probes, and Segment_Mean columns to numeric.
\item Checks that the End column is greater than the Start column and reports if not.
\item Identifies and reports outliers in the Segment_Mean column.
}
}
\examples{
# Example usage:
checkSegments(segments)

}
