library(magrittr)
library(dplyr)
library(stringr)


createCuts <- function(segments){

  ## Load cytoband
  cytoband <- BAGEL::cytoband.hg38

  ## Check NA in segments, remove segments labelled to be in centromere region
  segments <- na.omit(segments)
  segments <- segments %>% dplyr::filter(Arm != "centromere")
  message(paste0("There are ", length(unique(segments$Sample)), " samples."))

  ## First compute how many segments ends at centromere for each chromosoem arm
  ends_at_cent <- data.frame(row.names = unique(segments$Arm),
                             arms = unique(segments$Arm),
                             ends_at_cent = rep(0, length(unique(segments$Arm))))
  for(i in unique(segments$Arm)){
    arm_slice <- segments %>% dplyr::filter(Arm == i)
    cyto_slice <- cytoband %>% dplyr::filter(Arm == i)
    if(grepl("p", i)){
      ends_at_cent[i, "ends_at_cent"] <- sum(arm_slice$End == cyto_slice$End)/nrow(arm_slice)
    } else if (grepl("q", i)){
      ends_at_cent[i, "ends_at_cent"] <- sum(arm_slice$Start == cyto_slice$Start)/(nrow(arm_slice))
    } else {
      message("Mis-labelled arm, neither p or q. Please check your segments")
    }
  }
  ends_at_cent <- ends_at_cent %>% mutate(arms = factor(arms, levels = arms[order(ends_at_cent)]))
  ggplot(ends_at_cent, aes(y = arms, x = ends_at_cent)) +
    geom_col() +
    xlim(0, 1)

  ## Plot Distribution of Breakpoints
  for(i in unique(segments$Arm)){
    arm_slice <- segments %>% dplyr::filter(Arm == i)
  }

}







