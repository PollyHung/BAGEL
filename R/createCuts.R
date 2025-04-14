
library(parallel)




createCuts <- function(segments = segs){

  ## Clean up the segments
  segments <- na.omit(segments)
  segments <- segments %>% dplyr::filter(Arm != "centromere") ## remove segments that are defined in centromere region

  ## Preprocess to Create All the Cuts needed
  mclapply(unique(segments$Arm), function(i) {
    preprocessSeg(arm = i, result_dir = "example/pooledOV/breakpoints/")
  }, mc.cores = detectCores() - 4)

  ##




}
