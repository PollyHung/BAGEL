
library(parallel)




createCuts <- function(segments = segs,
                       result_dir){

  ## Define support functions
  filter_big_small <- function(df) {
    df <- df %>%
      dplyr::filter(Percent >= 0.001) %>%
      dplyr::filter(Percent <= 0.999) %>%
      dplyr::arrange(Percent)
    return(df)
  }

  ## Clean up the segments
  segments <- na.omit(segments)
  segments <- segments %>% dplyr::filter(Arm != "centromere") ## remove segments that are defined in centromere region

  ## Preprocess to Create All the Cuts needed
  mclapply(unique(segments$Arm), function(i) {
    preprocessSeg(arm = i, result_dir)
  }, mc.cores = detectCores() - 4)

  ## Generate Backgrounds for each lineage
  background(result_dir)

  ## read in the background files and process
  tel <- read.delim(file.path(result_dir, "backgrounds/background_telomere.txt")) %>% filter_big_small
  telemp <- tel$Percent
  cent <- read.delim(file.path(result_dir, "backgrounds/background_centromere.txt")) %>% filter_big_small
  centemp <- cent$Percent

  ## Separate the Amp and Del
  amptel <- tel %>% dplyr::filter(amp_del == "amp"); amptelemp <- amptel$Percent
  deltel <- tel %>% dplyr::filter(amp_del == "del"); deltelemp <- deltel$Percent
  ampcent <- cent %>% dplyr::filter(amp_del == "amp"); ampcentemp <- ampcent$Percent
  delcent <- cent %>% dplyr::filter(amp_del == "del"); delcentemp <- delcent$Percent

  ##

}






