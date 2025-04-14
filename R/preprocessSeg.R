## This is a replication of BISCUT Preprocess function in R, the original code is in python


preprocessSeg <- function(segments = segs,
                          genome = "hg38",
                          cytoband = NULL,
                          del_cutoff = -0.2,
                          amp_cutoff = 0.2){

  ## Load in the correct arm level coordinate file
  if(genome == "hg38"){
    coords <- BAGEL::cytoband.hg38
  } else if (genome == "hg19") {
    coords <- BAGEL::cytoband.hg19
  } else if (!is.null(cytoband)) {
    coords <- read.csv(cytoband)
  } else {
    stop("No valid arm-level coordinates provided, exiting")
  }

  ## For each arm
  segs <- segs %>% dplyr::filter(Arm != "centromere") ## remove any centromere segments

  for(arm in unique(segs$Arm)){
    coord <- coords %>% dplyr::filter(coords$Arm == arm)

    segdf <- segs %>% dplyr::filter(Arm == arm)
    segdf <- segdf %>% dplyr::mutate(Status = case_when(Segment_Mean < -del_cutoff ~ "DEL",
                                                        Segment_Mean > amp_cutoff ~ "AMP",
                                                        TRUE ~ "NEUTRAL"))
    for(dir in c("AMP", "DEL", "NEUTRAL")){
      segdel <- segdf %>% dplyr::filter(Status == dir)
    }
  }




}

aneu <- BAGEL::aneu



