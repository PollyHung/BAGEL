#' Plot Chromosome Ideograms with Functional Arm Definitions
#'
#' Creates ideogram plots showing chromosome arms with functional regions highlighted.
#' Individual plots for each chromosome and a combined panel plot are generated.
#'
#' @param bagel_results List, results from calculateCopyNumber_fixed function
#' @param output_dir Character, directory to save plots (default: current directory)
#' @param cytoband_data Data frame, cytoband data for chromosome boundaries (default: cytoband.hg38)
#' @param plot_width Numeric, width of individual plots in inches (default: 3)
#' @param plot_height Numeric, height of individual plots in inches (default: 3)
#' @param save_plots Logical, whether to save plots to files (default: TRUE)
#' @return List of ggplot objects (individual plots and combined panel)
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_chromosome_ideograms <- function(bagel_results,
                                      output_dir = ".",
                                      cytoband_data = NULL,
                                      plot_width = 3,
                                      plot_height = 3,
                                      save_plots = TRUE) {
  
  # Load required libraries
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("RColorBrewer", quietly = TRUE)
  requireNamespace("patchwork", quietly = TRUE)
  
  # Use BISCUT-based cytoband data if not provided for coordinate consistency
  if (is.null(cytoband_data)) {
    tryCatch({
      cytoband_data <- get_biscut_cytoband()
      cat("Using BISCUT chromosome coordinates for consistency\n")
    }, error = function(e) {
      warning("Could not load BISCUT coordinates, falling back to cytoband.hg38: ", e$message)
      cytoband_data <- BAGEL::cytoband.hg38
    })
  }
  
  # Get arm definitions from BAGEL results
  arm_definitions <- bagel_results$arm_definitions
  
  if (nrow(arm_definitions) == 0) {
    stop("No arm definitions found in bagel_results")
  }
  
  # Get unique chromosomes to plot
  chromosomes_to_plot <- sort(unique(arm_definitions$chr_num))
  set.seed(999)
  chr_colors <- sample(BAGEL::bagel_palette, size = 22, replace = F)
  names(chr_colors) <- c(1:22)
  chr_colors <- chr_colors[chromosomes_to_plot]
  
  # Create output directory
  if (save_plots) {
    ideogram_dir <- file.path(output_dir, "newly_defined_arms_ideogram")
    if (!dir.exists(ideogram_dir)) {
      dir.create(ideogram_dir, recursive = TRUE)
      cat("Created ideogram directory:", ideogram_dir, "\n")
    }
  }
  
  # Function to create individual chromosome ideogram
  create_chromosome_ideogram <- function(chr_num) {
    
    # Get full chromosome boundaries from cytoband data
    chr_cytoband <- cytoband_data %>%
      filter(Chromosome == chr_num) %>%
      arrange(Start)
    
    if (nrow(chr_cytoband) == 0) {
      warning(paste("No cytoband data found for chromosome", chr_num))
      return(NULL)
    }
    
    # Get p and q arm boundaries
    p_arm_data <- chr_cytoband %>% filter(Arm == paste0(chr_num, "p"))
    q_arm_data <- chr_cytoband %>% filter(Arm == paste0(chr_num, "q"))
    centromere_data <- chr_cytoband %>% filter(Arm == "centromere")
    
    if (nrow(p_arm_data) == 0 || nrow(q_arm_data) == 0) {
      warning(paste("Incomplete arm data for chromosome", chr_num))
      return(NULL)
    }
    
    # Calculate arm lengths and positions
    p_start <- min(p_arm_data$Start)
    p_end <- max(p_arm_data$End)
    p_length <- p_end - p_start
    
    q_start <- min(q_arm_data$Start)
    q_end <- max(q_arm_data$End)
    q_length <- q_end - q_start
    
    # Centromere position
    if (nrow(centromere_data) > 0) {
      centromere_start <- min(centromere_data$Start)
      centromere_end <- max(centromere_data$End)
      centromere_length <- abs(centromere_end - centromere_start)
    } else {
      # Approximate centromere as gap between p and q arms
      centromere_start <- p_end
      centromere_end <- q_start
      centromere_length <- 0
    }
    
    # Scale arms for plotting (total height = 10 units)
    total_length <- p_length + q_length + centromere_length
    p_height <- (p_length / total_length) * 10
    q_height <- (q_length / total_length) * 10
    centromere_height <- (centromere_length / total_length) * 10
    
    # Create base plot data
    plot_data <- data.frame(
      arm = c("p", "q"),
      ymin = c(q_height + centromere_height, 0),  # Gap of 0.5 units for centromere
      ymax = c(10, q_height),
      start_pos = c(p_start, q_start),
      end_pos = c(p_end, q_end),
      length = c(p_length, q_length),
      stringsAsFactors = FALSE
    )
    
    # Get functional arm definitions for this chromosome
    chr_arm_defs <- arm_definitions %>% filter(chr_num == !!chr_num)
    
    # Create plot
    chr_color <- chr_colors[as.character(chr_num)]
    
    p <- ggplot() +
      # Background chromosome arms (semi-transparent)
      geom_rect(data = plot_data, 
                aes(xmin = -0.3, xmax = 0.3, ymin = ymin, ymax = ymax),
                fill = chr_color, alpha = 0.5, color = "black", size = 0.5) +
      
      # Coordinate labels - right side for p arm, left side for q arm
      geom_text(data = plot_data[plot_data$arm == "p", ],
                aes(x = 0.4, y = ymax, label = format(start_pos, big.mark = ",")),
                size = 3, hjust = 0) +
      geom_text(data = plot_data[plot_data$arm == "p", ],
                aes(x = 0.4, y = ymin, label = format(end_pos, big.mark = ",")),
                size = 3, hjust = 0) +
      geom_text(data = plot_data[plot_data$arm == "q", ],
                aes(x = -0.4, y = ymax, label = format(start_pos, big.mark = ",")),
                size = 3, hjust = 1) +
      geom_text(data = plot_data[plot_data$arm == "q", ],
                aes(x = -0.4, y = ymin, label = format(end_pos, big.mark = ",")),
                size = 3, hjust = 1) 
    
    # OLD CODE - COMMENTED OUT (replaced with individual subplot approach)
    # # Add functional arm overlays
    # if (nrow(chr_arm_defs) > 0) {
    #   for (i in 1:nrow(chr_arm_defs)) {
    #     arm_def <- chr_arm_defs[i, ]
    #     arm_type <- arm_def$arm_type
    #     
    #     # Get the corresponding background arm data
    #     bg_arm <- plot_data[plot_data$arm == arm_type, ]
    #     
    #     if (nrow(bg_arm) > 0) {
    #       # Calculate functional region position within the arm
    #       func_start <- arm_def$arm_start
    #       func_end <- arm_def$arm_end
    #       
    #       # Calculate proportional position within the arm
    #       arm_total_length <- bg_arm$length
    #       arm_start_genomic <- bg_arm$start_pos
    #       arm_end_genomic <- bg_arm$end_pos
    #       
    #       # Functional region as proportion of total arm
    #       if (arm_type == "p") {
    #         # For p arm, functional region starts from genomic start
    #         func_prop_start <- 1-(func_start - arm_start_genomic) / arm_total_length
    #         func_prop_end <- 1-(func_end - arm_start_genomic) / arm_total_length
    #       } else {
    #         # For q arm, functional region starts from genomic start
    #         func_prop_start <- 1-(func_start - arm_start_genomic) / arm_total_length
    #         func_prop_end <- 1-(func_end - arm_start_genomic) / arm_total_length
    #       }
    #       
    #       # Calculate plot coordinates
    #       arm_height <- bg_arm$ymax - bg_arm$ymin
    #       func_ymin <- bg_arm$ymin + (func_prop_start * arm_height)
    #       func_ymax <- bg_arm$ymin + (func_prop_end * arm_height)
    #       
    #       # Add functional region overlay
    #       p <- p + geom_rect(aes(xmin = -0.3, xmax = 0.3, 
    #                              ymin = func_ymin, ymax = func_ymax),
    #                          fill = chr_color, alpha = 1.0, color = "black", size = 0.5) + 
    #         
    #         # Arm labels
    #         geom_text(data = plot_data,
    #                   aes(x = 0, y = (ymin + ymax)/2, label = arm),
    #                   size = 4, fontface = "bold")
    #       
    #       # Add functional start coordinate label in red
    #       func_coord_y <- if (arm_type == "p") func_ymax else func_ymin
    #       func_coord_x <- if (arm_type == "p") 0.4 else -0.4
    #       func_coord_hjust <- if (arm_type == "p") 0 else 1
    #       func_label <- if (arm_type == "p") func_end else func_start
    #       
    #       p <- p + geom_text(aes(x = func_coord_x, y = func_coord_y, 
    #                              label = format(func_label, big.mark = ",")),
    #                          size = 3, hjust = func_coord_hjust, color = "red", fontface = "bold")
    #       
    #     }
    #   }
    # }
    
    # NEW APPROACH: Check if chromosome has multiple functional regions
    chr_arm_defs$label <- paste0(chr_arm_defs$arm, "_",chr_arm_defs$direction)
    if (nrow(chr_arm_defs) > 1) {
      # Create individual plots for each functional region
      individual_plots <- list()
      chr_arm_defs_split <- split(chr_arm_defs, f = chr_arm_defs$label)
      
      individual_plots <- lapply(chr_arm_defs_split, function(arm_def){
        arm_type <- arm_def$arm_type
        bg_arm <- plot_data[plot_data$arm == arm_type, ]
        
        # Calculate functional region position within the arm
        func_start <- arm_def$arm_start
        func_end <- arm_def$arm_end
        
        # Calculate proportional position within the arm
        arm_total_length <- bg_arm$length
        arm_start_genomic <- bg_arm$start_pos
        
        # Functional region as proportion of total arm
        if (arm_type == "p") {
          func_prop_start <- round(1-(func_start - arm_start_genomic) / arm_total_length)
          func_prop_end <- 1-(func_end - arm_start_genomic) / arm_total_length
        } else {
          func_prop_start <- 1-(func_start - arm_start_genomic) / arm_total_length
          func_prop_end <- round(1-(func_end - arm_start_genomic) / arm_total_length)
        }
        
        # Calculate plot coordinates
        arm_height <- bg_arm$ymax - bg_arm$ymin
        func_ymin <- bg_arm$ymin + (func_prop_start * arm_height)
        func_ymax <- bg_arm$ymin + (func_prop_end * arm_height)
        
        # Create individual subplot for this functional region
        subplot <- ggplot() +
          # Background chromosome arms (semi-transparent)
          geom_rect(data = plot_data, 
                    aes(xmin = -0.3, xmax = 0.3, ymin = ymin, ymax = ymax),
                    fill = chr_color, alpha = 0.5, color = "black", size = 0.5) +
          
          # Add this specific functional overlay
          geom_rect(aes(xmin = -0.3, xmax = 0.3, ymin = func_ymin, ymax = func_ymax),
                    fill = chr_color, alpha = 1.0, color = "black", size = 0.5) +
          
          # Coordinate labels
          geom_text(data = plot_data[plot_data$arm == "p", ],
                    aes(x = 0.4, y = ymax, label = format(start_pos, big.mark = ",")),
                    size = 3, hjust = 0) +
          geom_text(data = plot_data[plot_data$arm == "p", ],
                    aes(x = 0.4, y = ymin, label = format(end_pos, big.mark = ",")),
                    size = 3, hjust = 0) +
          geom_text(data = plot_data[plot_data$arm == "q", ],
                    aes(x = -0.4, y = ymax, label = format(start_pos, big.mark = ",")),
                    size = 3, hjust = 1) +
          geom_text(data = plot_data[plot_data$arm == "q", ],
                    aes(x = -0.4, y = ymin, label = format(end_pos, big.mark = ",")),
                    size = 3, hjust = 1) +
          
          # Arm labels
          geom_text(data = plot_data,
                    aes(x = 0, y = (ymin + ymax)/2, label = arm),
                    size = 4, fontface = "bold") +
          
          # Functional coordinate label in red
          geom_text(aes(x = if(arm_type == "p") 0.4 else -0.4, 
                        y = if(arm_type == "p") func_ymax else func_ymin,
                        label = format(if(arm_type == "p") func_end else func_start, big.mark = ",")),
                    hjust = if(arm_type == "p") 0 else 1,
                    size = 3, color = "red", fontface = "bold") +
          
          # Deletion/Amplification caption in lower right
          geom_text(aes(x = 1.8, y = 0.5, 
                        label = paste(toupper(substring(arm_def$direction, 1, 1)), 
                                      substring(arm_def$direction, 2), sep = "")),
                    hjust = 1, vjust = 0, size = 3, 
                    color = if(arm_def$direction == "amp") "red" else "blue",
                    fontface = "bold") +
          
          # Styling with subtitle showing the functional region
          theme_minimal() + 
          labs(title = paste("Chromosome", chr_num)) +
          theme(axis.text = element_blank(), 
                axis.title = element_blank(), 
                plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                panel.grid = element_blank()) +
          xlim(-2, 2) +
          ylim(0, 10)
        
        return(subplot)
      })
      
      
      # Combine individual plots into a panel
      if (length(individual_plots) <= 2) {
        # Horizontal layout for 2 or fewer
        combined_plot <- wrap_plots(individual_plots, nrow = 1)
      } else {
        # Grid layout for more than 2
        combined_plot <- wrap_plots(individual_plots, ncol = 2)
      }
      
      return(combined_plot)
      
    } else if (nrow(chr_arm_defs) == 1) {
      # Single functional region - use simplified approach
      arm_def <- chr_arm_defs[1, ]
      arm_type <- arm_def$arm_type
      bg_arm <- plot_data[plot_data$arm == arm_type, ]
      
      if (nrow(bg_arm) > 0) {
        # Calculate functional region position within the arm
        func_start <- arm_def$arm_start
        func_end <- arm_def$arm_end
        
        # Calculate proportional position within the arm
        arm_total_length <- bg_arm$length
        arm_start_genomic <- bg_arm$start_pos
        
        # Functional region as proportion of total arm
        if (arm_type == "p") {
          func_prop_start <- round(1-(func_start - arm_start_genomic) / arm_total_length)
          func_prop_end <- 1-(func_end - arm_start_genomic) / arm_total_length
        } else {
          func_prop_start <- 1-(func_start - arm_start_genomic) / arm_total_length
          func_prop_end <- round(1-(func_end - arm_start_genomic) / arm_total_length)
        }
        
        # Calculate plot coordinates
        arm_height <- bg_arm$ymax - bg_arm$ymin
        func_ymin <- bg_arm$ymin + (func_prop_start * arm_height)
        func_ymax <- bg_arm$ymin + (func_prop_end * arm_height)
        
        # Add functional region overlay
        p <- p + geom_rect(aes(xmin = -0.3, xmax = 0.3, 
                               ymin = func_ymin, ymax = func_ymax),
                           fill = chr_color, alpha = 1.0, color = "black", size = 0.5) +
          
          # Arm labels
          geom_text(data = plot_data,
                    aes(x = 0, y = (ymin + ymax)/2, label = arm),
                    size = 4, fontface = "bold") +
          
          # Add functional start coordinate label in red
          geom_text(aes(x = if(arm_type == "p") 0.4 else -0.4, 
                        y = if(arm_type == "p") func_ymax else func_ymin,
                        label = format(if(arm_type == "p") func_end else func_start, big.mark = ",")),
                    hjust = if(arm_type == "p") 0 else 1,
                    size = 3, color = "red", fontface = "bold") +
          
          # Deletion/Amplification caption in lower right
          geom_text(aes(x = 1.8, y = 0.5, 
                        label = paste(toupper(substring(arm_def$direction, 1, 1)), 
                                      substring(arm_def$direction, 2), sep = "")),
                    hjust = 1, vjust = 0, size = 3, 
                    color = if(arm_def$direction == "amp") "red" else "blue",
                    fontface = "bold")
        
      }
    } 
    
    # Final plot styling
    p <- p + 
      theme_minimal() + 
      labs(title = paste("Chromosome", chr_num)) +
      theme(axis.text = element_blank(), 
            axis.title = element_blank(), 
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
            panel.grid = element_blank()) +
      xlim(-2, 2) +
      ylim(0, 10)
    
    return(p)
  }
  
  # Generate individual plots
  individual_plots <- list()
  
  for (chr_num in chromosomes_to_plot) {
    cat("Creating ideogram for chromosome", chr_num, "\n")
    plot <- create_chromosome_ideogram(chr_num)
    
    if (!is.null(plot)) {
      individual_plots[[as.character(chr_num)]] <- plot
      
      # Save individual plot with dynamic sizing
      if (save_plots) {
        # Check if this is a combined plot from multiple functional regions
        if (inherits(plot, "patchwork")) {
          # Calculate dimensions based on number of functional regions for this chromosome
          chr_arm_defs_for_chr <- arm_definitions %>% filter(chr_num == !!chr_num)
          n_subplots <- nrow(chr_arm_defs_for_chr)
          
          if (n_subplots <= 2) {
            # Horizontal layout: nrow = 1
            save_width <- n_subplots * plot_width  # 6x3 for 2 plots
            save_height <- plot_height
          } else {
            # Grid layout: ncol = 2
            save_width <- 2 * plot_width  # Always 6 inches width for grid
            save_height <- plot_height * ceiling(n_subplots / 2)  # 3 * ceiling(n/2)
          }
          
          ggsave(filename = file.path(ideogram_dir, paste0("chromosome_", chr_num, "_ideogram.png")),
                 plot = plot, width = save_width, height = save_height, units = "in", dpi = 300)
        } else {
          # Single plot - use default dimensions
          ggsave(filename = file.path(ideogram_dir, paste0("chromosome_", chr_num, "_ideogram.png")),
                 plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 300)
        }
      }
    }
  }
  
  # Create combined panel plot
  n_plots <- length(individual_plots)
  
  if (n_plots > 0) {
    if (n_plots <= 4) {
      # Single row for 4 or fewer chromosomes
      panel_width <- n_plots * plot_width
      panel_height <- plot_height
      combined_plot <- wrap_plots(individual_plots, nrow = 1)
    } else {
      # Multi-row layout for more chromosomes
      ncol <- 4
      nrow <- ceiling(n_plots / ncol)
      panel_width <- ncol * plot_width
      panel_height <- nrow * plot_height
      combined_plot <- wrap_plots(individual_plots, ncol = ncol)
    }
    
    # Save combined plot
    if (save_plots) {
      ggsave(filename = file.path(ideogram_dir, "all_chromosomes_ideogram_panel.png"),
             plot = combined_plot, width = panel_width, height = panel_height, units = "in", dpi = 300)
    }
    
    return(list(
      individual_plots = individual_plots,
      combined_plot = combined_plot,
      chromosome_colors = chr_colors,
      output_directory = if (save_plots) ideogram_dir else NULL
    ))
    
  } else {
    warning("No valid chromosome ideograms were created")
    return(NULL)
  }
}



