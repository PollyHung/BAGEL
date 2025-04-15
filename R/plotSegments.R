#' Plot Segments of Chromosomes
#'
#' This function generates plots of genomic segments for specified chromosomes or arms,
#' highlighting segments based on their mean values. Plots are saved in a specified
#' directory.
#'
#' @param segments A data frame containing genomic segment information, including
#'   columns for chromosome, arm, start and end positions, and segment mean values.
#' @param genome A character string indicating the genome version to use
#'   ("hg19" or "hg38"). Default is "hg19".
#' @param chrs An optional integer specifying a specific chromosome to plot.
#' @param arms An optional character string specifying a specific chromosome arm to plot.
#' @param all_chrs A logical value indicating whether to plot all chromosomes (1-22).
#' @param result_dir A character string specifying the directory where plots will be saved.
#' @param cutoff A numeric value used to filter segments based on their mean values.
#'   Default is 0.20.
#' @param w A numeric value specifying the width of the output plots in inches.
#'   Default is 6.
#' @param h A numeric value specifying the height of the output plots in inches.
#'   Default is 6.
#'
#' @return The function generates and saves plots of genomic segments to the specified
#'   directory. No value is returned.
#'
#' @import dplyr
#' @import ggplot2
#' @import BAGEL
#' @export
plotSegments <- function(segments,
                         genome = "hg19",
                         chrs = NULL,
                         arms = NULL,
                         all_chrs = FALSE,
                         result_dir,
                         cutoff = 0.20,
                         w = 6,
                         h = 6){
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

  plot_dir = file.path(result_dir, "plots")
  plot_dir1 = file.path(plot_dir, "01.SegsOverCutoff")
  if(!dir.exists(plot_dir)){ dir.create(plot_dir) }
  if(!dir.exists(plot_dir1)){ dir.create(plot_dir1) }

  ## Plot Individual Plots
  if(!is.null(arms)){
    data <- segments %>%
      filter(Arm == arms) %>%
      filter(abs(Segment_Mean) >= cutoff)
    coord <- coords %>% dplyr::filter(Arm == arms)
    start_bp <- coord$Start
    end_bp <- coord$End
    p <- ggplot() + geom_rect(aes(xmin = start_bp, xmax = end_bp, ymin = -1, ymax = 1), fill = "white") +
      geom_rect(data = data, aes(xmin = Start, xmax = End,
                                 ymin = as.numeric(factor(Sample)),
                                 ymax = as.numeric(factor(Sample)) + 0.5,
                                 fill = Segment_Mean > 0), alpha = 0.7) +
      scale_fill_manual(values = c("TRUE" = "#FF8080", "FALSE" = "#8EACCD"),
                        labels = c("TRUE" = "Segment Mean > 0", "FALSE" = "Segment Mean < 0")) +
      geom_vline(xintercept = c(start_bp, end_bp), linewidth = 0.1) +
      theme_minimal(base_line_size = 0) +
      labs(title = paste0("Segments of Chromosome Arm ", arms),
           subtitle = paste0("Number of Samples: ", length(unique(data$Sample)))) +
      ylab("Each Line is One Sample") +
      xlab(paste0("Chromosome Coordinates (", genome, ")")) +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_blank(),
            legend.position = "top")
    plot_name <- paste0("arm_", arms, ".png")
    ggsave(filename = file.path(plot_dir1, plot_name), plot = p, width = w, height = h, units = "in", dpi = 600)
  }

  if (!is.null(chrs)){
    data <- segments %>%
      filter(Chromosome == chrs) %>%
      filter(abs(Segment_Mean) >= cutoff)
    coord <- coords %>% dplyr::filter(Chromosome == chrs)
    start_bp <- coord$Start[1]
    end_bp <- coord$End[3]
    cent_start <- coord$Start[2]
    cent_end <- coord$End[2]
    p <- ggplot() + geom_rect(aes(xmin = start_bp, xmax = end_bp, ymin = -1, ymax = 1), fill = "white") +
      geom_rect(data = data, aes(xmin = Start, xmax = End,
                                 ymin = as.numeric(factor(Sample)),
                                 ymax = as.numeric(factor(Sample)) + 0.5,
                                 fill = Segment_Mean > 0), alpha = 0.7) +
      scale_fill_manual(values = c("TRUE" = "#FF8080", "FALSE" = "#8EACCD"),
                        labels = c("TRUE" = "Amplification", "FALSE" = "Deletion")) +
      geom_vline(xintercept = c(start_bp, end_bp, cent_start, cent_end), linewidth = 0.1) +
      theme_minimal(base_line_size = 0) +
      labs(title = paste0("Segments of Chromosome ", chrs),
           subtitle = paste0("Number of Samples: ", length(unique(data$Sample)))) +
      ylab("Each Line is One Sample") +
      xlab(paste0("Chromosome Coordinates (", genome, ")")) +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_blank(),
            legend.position = "top")
    plot_name <- paste0("chr_", chrs, ".png")
    ggsave(filename = file.path(plot_dir1, plot_name), plot = p, width = w, height = h, units = "in", dpi = 600)
  }

  if (isTRUE(all_chrs)){
    for(chrs in c(1:22)){
      data <- segments %>%
        filter(Chromosome == chrs) %>%
        filter(abs(Segment_Mean) >= cutoff)
      coord <- coords %>% dplyr::filter(Chromosome == chrs)
      start_bp <- coord$Start[1]
      end_bp <- coord$End[3]
      cent_start <- coord$Start[2]
      cent_end <- coord$End[2]
      p <- ggplot() + geom_rect(aes(xmin = start_bp, xmax = end_bp, ymin = -1, ymax = 1), fill = "white") +
        geom_rect(data = data, aes(xmin = Start, xmax = End,
                                   ymin = as.numeric(factor(Sample)),
                                   ymax = as.numeric(factor(Sample)) + 0.5,
                                   fill = Segment_Mean > 0), alpha = 0.7) +
        scale_fill_manual(values = c("TRUE" = "#FF8080", "FALSE" = "#8EACCD"),
                          labels = c("TRUE" = "Amplification", "FALSE" = "Deletion")) +
        geom_vline(xintercept = c(start_bp, end_bp, cent_start, cent_end), linewidth = 0.1) +
        theme_minimal(base_line_size = 0) +
        labs(title = paste0("Segments of Chromosome ", chrs),
             subtitle = paste0("Number of Samples: ", length(unique(data$Sample)))) +
        ylab("Each Line is One Sample") +
        xlab(paste0("Chromosome Coordinates (", genome, ")")) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.title = element_blank(),
              legend.position = "top")
      plot_name <- paste0("chr_", chrs, ".png")
      ggsave(filename = file.path(plot_dir1, plot_name), plot = p, width = w, height = h, units = "in", dpi = 600)
    }
  }
}

