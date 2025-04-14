## This is a replication of BISCUT Preprocess function in R, the original code is in python


preprocessSeg <- function(segments = segs,
                          genome = "hg38",
                          arm = "1p",
                          cytoband = NULL,
                          del_cutoff = -0.25,
                          amp_cutoff = 0.25,
                          cutoff_SD = 0.25){

  ## Load in the correct arm level coordinate file
  if(genome == "hg38"){
    coords <- BAGEL::cytoband.hg38
  } else if (genome == "hg19") {
    coords <- BAGEL::cytoband.hg19
  } else if (!is.null(cytoband)) {
    coords <- read.csv(cytoband)
    colnames(coords) <- c("Chromosome", "Arm", "Start", "End")
  } else {
    stop("No valid arm-level coordinates provided, exiting")
  }

  ## Slice out the corresponding chromosome arms and coordinates
  segments <- segments %>% dplyr::filter(Arm == arm)
  coord <- coords %>% dplyr::filter(Arm == arm)

  ## How long is the chromosome arm?
  arm_length <- coord$End - coord$Start

  ## Define the cutoffs
  if(!is.null(cutoff_SD)){
    message("You provided cutoff standard deviation, thus overiding hard threshold with customized cutoff based on standard deviation")
    segMean <- mean(segments$Segment_Mean)
    segSD <- sd(segments$Segment_Mean)
    amp_cutoff <- segMean + cutoff_SD * segSD
    del_cutoff <- segMean - cutoff_SD * segSD
  }


  ## Add the AMP/DEL and NEUTRAL label
  segments <- segments %>% dplyr::mutate(Status = case_when(Segment_Mean <= del_cutoff ~ "DEL",
                                                            Segment_Mean >= amp_cutoff ~ "AMP",
                                                            TRUE ~ "NEUTRAL"))
  message(paste(names(table(segments$Status)), table(segments$Status), collapse = ", "))

  ## For each DEL and AMP, are they originating from telomere?
  telcent <- function(segdf = segments,
                      direction,
                      telcent){

    ## First filter for correct direction
    segdf <- segdf %>% dplyr::filter(Status == direction)

    ## Then add Telomere
    if(telcent == "TEL"){
      if(grepl("p", arm)){segdf$telcent <- ifelse(segdf$Start < coord$Start, "TEL", "INTER")}
      if(grepl("q", arm)){segdf$telcent <- ifelse(segdf$End > coord$End, "TEL", "INTER")}
    } else if (telcent == "CENT"){
      if(grepl("p", arm)){segdf$telcent <- ifelse(segdf$End > coord$End, "CENT", "INTER")}
      if(grepl("q", arm)){segdf$telcent <- ifelse(segdf$Start < coord$Start, "CENT", "INTER")}
    }

    return(segdf)
  }

  ## Apply the TELCENT Function
  seg_list <- list()
  seg_list[["amp_tel"]] <- telcent(direction = "AMP", telcent = "TEL")
  seg_list[["del_tel"]] <- telcent(direction = "DEL", telcent = "TEL")
  seg_list[["amp_cent"]] <- telcent(direction = "AMP", telcent = "CENT")
  seg_list[["del_cent"]] <- telcent(direction = "DEL", telcent = "CENT")

  ## Define the join segments function
  joinSegs <- function(segdf,
                       telcent){
    ## initiate an empty list
    results <- list()

    ## Now for each sample
    for(sample in unique(segdf$Sample)){
      segdf_slice <- segdf %>% dplyr::filter(Sample == sample)
      if(sum(grepl(telcent, segdf_slice$telcent)) > 0){


      } else {
        results[[sample]] <- data.frame(Sample = sample, Percent = 0, Start = 0, End = 0)
      }
    }

    ## bind the list and return
    return(do.call(rbind, results))
  }




}
