#' Get BISCUT-based Chromosome Coordinates
#' 
#' Converts BISCUT chromosome coordinates to BAGEL cytoband format for consistency
#' between chromosome arm definitions used in BISCUT and BAGEL analysis.
#'
#' @return Data frame with BISCUT-based chromosome arm coordinates in cytoband format
#' @import dplyr
#' @importFrom BISCUT get_chromosome_coordinates
#' @export
get_biscut_cytoband <- function() {
  
  # Check if BISCUT is available
  if (!requireNamespace("BISCUT", quietly = TRUE)) {
    stop("BISCUT package is required but not available. Please install it.")
  }
  
  # Get BISCUT chromosome coordinates
  biscut_coords <- BISCUT::get_chromosome_coordinates()
  
  # Convert to cytoband-like format
  cytoband_data <- data.frame()
  
  for (i in 1:nrow(biscut_coords)) {
    chr_num <- biscut_coords$chromosome_info[i]
    
    # Skip sex chromosomes (chr 23)
    if (is.na(chr_num) || chr_num > 22) next
    
    p_start <- biscut_coords$p_start[i]
    p_end <- biscut_coords$p_end[i] 
    q_start <- biscut_coords$q_start[i]
    q_end <- biscut_coords$q_end[i]
    centromere <- biscut_coords$centromere[i]
    
    # Create entries for p arm, centromere, and q arm
    chr_data <- data.frame(
      Chromosome = rep(chr_num, 3),
      Arm = c(paste0(chr_num, "p"), "centromere", paste0(chr_num, "q")),
      Start = c(p_start, p_end, q_start),
      End = c(p_end, q_start, q_end),
      stringsAsFactors = FALSE
    )
    
    cytoband_data <- rbind(cytoband_data, chr_data)
  }
  
  # Sort by chromosome and position
  cytoband_data <- cytoband_data %>%
    arrange(Chromosome, Start)
  
  return(cytoband_data)
}

#' Update Consensus Arm Definitions with BISCUT Coordinates
#' 
#' Updates the consensus arm definitions to use BISCUT chromosome coordinates
#' instead of the default cytoband coordinates for consistency.
#'
#' @param consensus_arm_definitions Data frame with consensus arm definitions
#' @return Updated data frame with BISCUT-based coordinates
#' @import dplyr
#' @importFrom BISCUT get_chromosome_coordinates
#' @export
update_arm_definitions_with_biscut <- function(consensus_arm_definitions) {
  
  # Get BISCUT coordinates
  biscut_coords <- BISCUT::get_chromosome_coordinates()
  
  # Create lookup table for chromosome boundaries
  chr_boundaries <- data.frame()
  
  for (i in 1:nrow(biscut_coords)) {
    chr_num <- biscut_coords$chromosome_info[i]
    
    # Skip sex chromosomes
    if (is.na(chr_num) || chr_num > 22) next
    
    # Add p and q arm boundaries
    boundaries <- data.frame(
      chr_num = rep(chr_num, 2),
      arm_type = c("p", "q"),
      biscut_start = c(biscut_coords$p_start[i], biscut_coords$q_start[i]),
      biscut_end = c(biscut_coords$p_end[i], biscut_coords$q_end[i]),
      stringsAsFactors = FALSE
    )
    
    chr_boundaries <- rbind(chr_boundaries, boundaries)
  }
  
  # Update consensus arm definitions
  updated_definitions <- consensus_arm_definitions %>%
    left_join(chr_boundaries, by = c("chr_num", "arm_type")) %>%
    mutate(
      # Use BISCUT coordinates for full chromosome arm boundaries
      arm_start = case_when(
        arm_type == "p" ~ biscut_start,
        arm_type == "q" ~ biscut_start,
        TRUE ~ arm_start
      ),
      arm_end = case_when(
        arm_type == "p" ~ biscut_end, 
        arm_type == "q" ~ biscut_end,
        TRUE ~ arm_end
      )
    ) %>%
    select(-biscut_start, -biscut_end)
  
  return(updated_definitions)
}