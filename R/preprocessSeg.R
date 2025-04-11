


preprocessSeg <- function(seg.path,
                          min_probes = 10,
                          genome = "hg38",
                          cytoband = NULL){

  ## Step 1: Load in the Segments and Merge Segments if the number of probes are less than 4
  segs <- mergeSegments(seg.path = seg.path, min_probes = min_probes)

  ## Step 2: Cut Arms
  segs <- annotateArms(segs, genome = genome, cytoband = cytoband)
}
