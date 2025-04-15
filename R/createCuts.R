

createCuts <- function(segments = segs,
                       genome = "hg19",
                       cutoff = 0.25,
                       result_dir = "example/pooledOV/",
                       cytoband = NULL){

  ## Load in the correct arm level coordinate file based on the specified genome
  if (genome == "hg38") {
    coords <- BAGEL::cytoband.hg38
  } else if (genome == "hg19") {
    coords <- BAGEL::cytoband.hg19
  } else if (!is.null(cytoband)) {
    coords <- read.csv(cytoband)
    colnames(coords) <- c("Chromosome", "Arm", "Start", "End")
  } else {
    stop("No valid arm-level coordinates provided, exiting")
  }

  ## Define support functions
  filter_big_small <- function(df) {
    df <- df %>%
      dplyr::filter(percent >= 0.001) %>%
      dplyr::filter(percent <= 0.999) %>%
      dplyr::arrange(percent)
    return(df)
  }

  ## Clean up the segments
  segments <- na.omit(segments)
  segments <- segments %>% dplyr::filter(Arm != "centromere") ## remove segments that are defined in centromere region

  ## Preprocess to Create All the Cuts needed
  mclapply(unique(segments$Arm), function(i) {
    preprocessSeg(arm = i)
  }, mc.cores = detectCores() - 4)

  ## Generate Backgrounds for each lineage
  background(result_dir)

  ## read in the background files and process
  tel <- read.delim(file.path(result_dir, "backgrounds/background_telomere.txt")) %>% filter_big_small
  telemp <- tel$percent
  cent <- read.delim(file.path(result_dir, "backgrounds/background_centromere.txt")) %>% filter_big_small
  centemp <- cent$percent

  ## Separate the Amp and Del
  amptel <- tel %>% dplyr::filter(amp_del == "amp"); amptelemp <- amptel$Percent
  deltel <- tel %>% dplyr::filter(amp_del == "del"); deltelemp <- deltel$Percent
  ampcent <- cent %>% dplyr::filter(amp_del == "amp"); ampcentemp <- ampcent$Percent
  delcent <- cent %>% dplyr::filter(amp_del == "del"); delcentemp <- delcent$Percent

  ##

}






