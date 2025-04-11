library(magrittr)
library(dplyr)
library(stringr)


createCuts <- function(segments){

  ## Load cytoband
  cytoband <- BAGEL::cytoband.hg38

  ## segs
  segs <- na.omit(segs)
  segs$Num_Probes[segs$Num_Probes < 4] <- 4

  ## Calculate How many Ends At centromere
  ends_at_centromere <- data.frame(row.names = unique(segs$Arm),
                                   arms = unique(segs$Arm),
                                   centromere_ends = rep(0, length(unique(segs$Arm))))
  for(i in unique(segs$Arm)){
    arm_slice <- segs %>% dplyr::filter(Arm == i)
    cyto_slice <- cytoband %>% dplyr::filter(Arm == i)
    if(grepl("p", i)){
      ends_at_centromere[i, "centromere_ends"] <- sum(arm_slice$End == cyto_slice$End)/(nrow(arm_slice))
      print(hist(arm_slice$End))
      }
    if(grepl("q", i)){
      ends_at_centromere[i, "centromere_ends"] <- sum(arm_slice$Start == cyto_slice$Start)/(nrow(arm_slice))
      print(hist(arm_slice$Start))
      }
  }


}







