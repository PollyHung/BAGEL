% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessSeg.R
\name{preprocessSeg}
\alias{preprocessSeg}
\title{Preprocess Segments Function}
\usage{
preprocessSeg(Segments, Coords, Cutoff, Result_dir, arm)
}
\arguments{
\item{arm}{A string indicating the chromosome arm to be processed (e.g., "1p", "1q").}

\item{segments}{A data frame containing genomic segments with columns including \code{Arm} and \code{Segment_Mean}.}

\item{genome}{A string specifying the genome version (e.g., "hg38", "hg19"). Default is "hg38".}

\item{cytoband}{Optional; a path to a custom cytoband CSV file.}

\item{del_cutoff}{Numeric; threshold for labeling segments as deletions. Default is -0.25.}

\item{amp_cutoff}{Numeric; threshold for labeling segments as amplifications. Default is 0.25.}

\item{cutoff_SD}{Numeric; standard deviation used to adjust cutoffs based on segment means.}

\item{result_dir}{A string specifying the directory to save output files.}
}
\value{
NULL; output files are saved to the specified result directory.
The files include segments categorized as amplification and deletion
along with their telomere and centromere assignments.
}
\description{
This function preprocesses genomic segments by applying filtering and labeling based on specified cutoffs.
It identifies segments originating from telomeres or centromeres and joins overlapping segments.
}
