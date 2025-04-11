


preprocessSeg <- function(seg.path,
                          min_probes = 10,
                          genome = "hg38",
                          cytoband = NULL,
                          mergeSeg = FALSE){

  ## Step 1: Load in the Segments and Merge Segments if the number of probes are less than 4
  if(mergeSeg){
    segs <- mergeSegments(seg.path = seg.path, min_probes = min_probes)
  } else {
    message("Skip Merging, proceed with Annotation")
    segs <- read.delim(seg.path)
    colnames(segs) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
  }

  ## Step 2: Cut Arms
  segs <- annotateArms(segments = segs, genome = genome, cytoband = cytoband)
}
