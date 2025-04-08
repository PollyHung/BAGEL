## Technique can be one of TS, WGS, or WES

ploidyThreshold <- function(technique = TS){
  if (! technique %in% c("TS", "WGS", "WES", "SNP")) {
    stop("Technique unsupported, please skip this step.")
  }


  ## If Targeted Sequencing
  if (technique == "TS"){

  }


  ## If it is Whole Genome Sequencing, use GC content
  if (technique == "WGS"){
    # Calculate noise per bin using rolling windows or GC-adjusted baselines
    noise <- calculate_gc_corrected_noise(cna, noise)

    # Set stricter thresholds for high-confidence calls
    amp_thresh <- mean(noise) + 0.5 * sd(noise)
    del_thresh <- -amp_thresh
  }


  ## If it is Whole Exome Sequencing, use exon depth
  if (technique == "WES"){
    # Normalize noise by exon depth
    noise <- noise %>%
      mutate(log2ratio = log2ratio / coverage_depth)

    # Set thresholds
    amp_thresh <- mean(noise) + 1.5 * sd(noise)
    del_thresh <- -amp_thresh
  }


  ## If it is SNP array, use BAF (B-allele frequency) variance alongside log2 ratios
  if(technique == "SNP"){
    # Calculate joint noise from log2 and BAF
    noise <- noise %>%
      mutate(noise = sqrt(sd_log2^2 + sd_baf^2))

    # Set thresholds
    amp_thresh <- mean(noise$log2) + 2 * sd(noise$log2)
    del_thresh <- -amp_thresh
  }
}
