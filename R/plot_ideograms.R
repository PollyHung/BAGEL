#' Plot Chromosome Ideograms with Enhanced Functional Arm Support
#'
#' Creates comprehensive ideogram plots showing chromosome arms with functional regions.
#' Supports both telomere-direction (telomere to breakpoint) and centromere-direction
#' (centromere to breakpoint) arm definitions with proper p/q arm directional handling.
#' Individual plots for each chromosome and combined panel plots are generated.
#'
#' @param bagel_results List, results from calculate_copynumber function containing
#'   arm_definitions with telcent column ("telomere" or "centromere")
#' @param output_dir Character, directory to save plots (default: current directory)
#' @param cytoband_data Data frame, cytoband data for chromosome boundaries (default: BISCUT coordinates)
#' @param plot_width Numeric, width of individual plots in inches (default: 3)
#' @param plot_height Numeric, height of individual plots in inches (default: 3)
#' @param save_plots Logical, whether to save plots to files (default: TRUE)
#' @return List containing individual plots, combined plot, colors, and output directory
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_ideograms <- function(bagel_results,
                           output_dir = ".",
                           cytoband_data = NULL,
                           plot_width = 3,
                           plot_height = 3,
                           save_plots = TRUE) {
  
  # Level 1 Helper Functions ===================================================

  # Load and validate cytoband data
  load_cytoband_data <- function() {
    if (!is.null(cytoband_data)) return(cytoband_data)

    tryCatch({
      cytoband <- get_biscut_cytoband()
      cat("Using BISCUT chromosome coordinates for consistency\n")
      return(cytoband)
    }, error = function(e) {
      warning("Could not load BISCUT coordinates, falling back to cytoband.hg38: ", e$message)
      return(BAGEL::cytoband.hg38)
    })
  }

  # Setup plot environment and directories
  setup_plot_environment <- function() {
    # Load required libraries
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("RColorBrewer", quietly = TRUE)
    requireNamespace("patchwork", quietly = TRUE)

    # Create output directory if saving plots
    ideogram_dir <- NULL
    if (save_plots) {
      ideogram_dir <- file.path(output_dir, "arms_ideogram")
      if (!dir.exists(ideogram_dir)) {
        dir.create(ideogram_dir, recursive = TRUE)
        cat("Created ideogram directory:", ideogram_dir, "\n")
      }
    }

    return(ideogram_dir)
  }

  # Calculate arm positioning with enhanced directional support
  calculate_arm_positions <- function(chr_arm_defs, plot_data) {
    positions <- list()

    for (i in 1:nrow(chr_arm_defs)) {
      arm_def <- chr_arm_defs[i, ]
      arm_type <- arm_def$arm_type
      telcent_direction <- arm_def$telcent  # "telomere" or "centromere"

      # Get background arm data
      bg_arm <- plot_data[plot_data$arm == arm_type, ]
      if (nrow(bg_arm) == 0) next

      # Extract coordinates
      func_start <- arm_def$arm_start
      func_end <- arm_def$arm_end
      arm_total_length <- bg_arm$length
      arm_start_genomic <- bg_arm$start_pos
      arm_end_genomic <- bg_arm$end_pos

      # Calculate proportional positions based on direction and telcent
      if (arm_type == "p") {
        # P arm: telomere at chromosome start, centromere at chromosome end
        if (telcent_direction == "telomere") {
          # Telomere to breakpoint: start from telomere (chromosome start)
          func_prop_start <- 1 - (func_start - arm_start_genomic) / arm_total_length
          func_prop_end <- 1 - (func_end - arm_start_genomic) / arm_total_length
        } else {
          # Centromere to breakpoint: start from centromere (chromosome end)
          func_prop_start <- 1 - (func_end - arm_start_genomic) / arm_total_length
          func_prop_end <- 1 - (func_start - arm_start_genomic) / arm_total_length
        }
      } else {
        # Q arm: centromere at chromosome start, telomere at chromosome end
        if (telcent_direction == "telomere") {
          # Telomere to breakpoint: start from telomere (chromosome end)
          func_prop_start <- (func_end - arm_start_genomic) / arm_total_length
          func_prop_end <- (func_start - arm_start_genomic) / arm_total_length
        } else {
          # Centromere to breakpoint: start from centromere (chromosome start)
          func_prop_start <- (func_start - arm_start_genomic) / arm_total_length
          func_prop_end <- (func_end - arm_start_genomic) / arm_total_length
        }
      }

      # Calculate plot coordinates
      arm_height <- bg_arm$ymax - bg_arm$ymin
      func_ymin <- bg_arm$ymin + (min(func_prop_start, func_prop_end) * arm_height)
      func_ymax <- bg_arm$ymin + (max(func_prop_start, func_prop_end) * arm_height)

      # Determine label position and coordinate based on direction
      if (arm_type == "p") {
        label_x <- 0.4
        label_hjust <- 0
        if (telcent_direction == "telomere") {
          label_y <- func_ymin
          label_coord <- func_end  # Breakpoint coordinate
        } else {
          label_y <- func_ymax
          label_coord <- func_start  # Breakpoint coordinate
        }
      } else {
        label_x <- -0.4
        label_hjust <- 1
        if (telcent_direction == "telomere") {
          label_y <- func_ymax
          label_coord <- func_start  # Breakpoint coordinate
        } else {
          label_y <- func_ymin
          label_coord <- func_end  # Breakpoint coordinate
        }
      }

      positions[[i]] <- list(
        func_ymin = func_ymin,
        func_ymax = func_ymax,
        label_x = label_x,
        label_y = label_y,
        label_hjust = label_hjust,
        label_coord = label_coord,
        arm_def = arm_def
      )
    }

    return(positions)
  }

  # Create single chromosome ideogram plot
  create_single_ideogram <- function(chr_num, arm_definitions, cytoband_data, chr_colors) {
    # Get chromosome cytoband data
    chr_cytoband <- cytoband_data %>%
      filter(Chromosome == chr_num) %>%
      arrange(Start)

    if (nrow(chr_cytoband) == 0) {
      warning(paste("No cytoband data found for chromosome", chr_num))
      return(NULL)
    }

    # Extract arm and centromere data
    p_arm_data <- chr_cytoband %>% filter(Arm == paste0(chr_num, "p"))
    q_arm_data <- chr_cytoband %>% filter(Arm == paste0(chr_num, "q"))
    centromere_data <- chr_cytoband %>% filter(Arm == "centromere")

    if (nrow(p_arm_data) == 0 || nrow(q_arm_data) == 0) {
      warning(paste("Incomplete arm data for chromosome", chr_num))
      return(NULL)
    }

    # Calculate arm dimensions
    p_start <- min(p_arm_data$Start)
    p_end <- max(p_arm_data$End)
    p_length <- p_end - p_start
    q_start <- min(q_arm_data$Start)
    q_end <- max(q_arm_data$End)
    q_length <- q_end - q_start

    # Handle centromere positioning
    if (nrow(centromere_data) > 0) {
      centromere_start <- min(centromere_data$Start)
      centromere_end <- max(centromere_data$End)
      centromere_length <- abs(centromere_end - centromere_start)
    } else {
      centromere_start <- p_end
      centromere_end <- q_start
      centromere_length <- 0
    }

    # Scale for plotting (total height = 10 units)
    total_length <- p_length + q_length + centromere_length
    p_height <- (p_length / total_length) * 10
    q_height <- (q_length / total_length) * 10
    centromere_height <- (centromere_length / total_length) * 10

    # Create plot data structure
    plot_data <- data.frame(
      arm = c("p", "q"),
      ymin = c(q_height + centromere_height, 0),
      ymax = c(10, q_height),
      start_pos = c(p_start, q_start),
      end_pos = c(p_end, q_end),
      length = c(p_length, q_length),
      stringsAsFactors = FALSE
    )

    return(list(plot_data = plot_data, chr_cytoband = chr_cytoband))
  }

  # Generate individual functional region plots
  create_functional_plots <- function(chr_num, chr_arm_defs, plot_data, chr_color) {
    # Calculate functional region positions
    positions <- calculate_arm_positions(chr_arm_defs, plot_data)

    # Create individual plots for each functional region
    functional_plots <- list()

    for (i in seq_along(positions)) {
      pos <- positions[[i]]
      arm_def <- pos$arm_def

      # Create base plot
      subplot <- ggplot() +
        # Background chromosome arms
        geom_rect(data = plot_data,
                  aes(xmin = -0.3, xmax = 0.3, ymin = ymin, ymax = ymax),
                  fill = chr_color, alpha = 0.5, color = "black", size = 0.5) +

        # Functional region overlay
        geom_rect(aes(xmin = -0.3, xmax = 0.3,
                      ymin = pos$func_ymin, ymax = pos$func_ymax),
                  fill = chr_color, alpha = 1.0, color = "black", size = 0.5) +

        # Coordinate labels for chromosome boundaries
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

        # Functional breakpoint coordinate (in red)
        geom_text(aes(x = pos$label_x, y = pos$label_y,
                      label = format(pos$label_coord, big.mark = ",")),
                  hjust = pos$label_hjust, size = 3, color = "red", fontface = "bold") +

        # Direction label with telcent annotation
        geom_text(aes(x = 1.8, y = 0.5,
                      label = paste0(toupper(substring(arm_def$direction, 1, 1)),
                                     substring(arm_def$direction, 2),
                                     " (", arm_def$telcent, ")")),
                  hjust = 1, vjust = 0, size = 2.5,
                  color = if(arm_def$direction == "amp") "red" else "blue",
                  fontface = "bold") +

        # Apply consistent theme
        theme_minimal() +
        labs(title = paste("Chromosome", chr_num, "-", arm_def$arm_type, "arm")) +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
              panel.grid = element_blank()) +
        xlim(-2, 2) + ylim(0, 10)

      # Create unique identifier for this functional region
      region_id <- paste(arm_def$arm_type, arm_def$direction, arm_def$telcent, sep = "_")
      functional_plots[[region_id]] <- subplot
    }

    return(functional_plots)
  }

  # Main Function Execution ====================================================

  # Initialize environment and load data
  cytoband_data <- load_cytoband_data()
  ideogram_dir <- setup_plot_environment()

  # Validate input data
  arm_definitions <- bagel_results$arm_definitions
  if (is.null(arm_definitions) || nrow(arm_definitions) == 0) {
    stop("No arm definitions found in bagel_results")
  }

  # Ensure telcent column exists with proper defaults
  if (!"telcent" %in% names(arm_definitions)) {
    warning("telcent column not found in arm_definitions, defaulting to 'telomere'")
    arm_definitions$telcent <- "telomere"
  }

  # Setup chromosome colors and plotting parameters
  chromosomes_to_plot <- sort(unique(arm_definitions$chr_num))
  set.seed(999)
  chr_colors <- sample(BAGEL::bagel_palette, size = 22, replace = FALSE)
  names(chr_colors) <- as.character(1:22)
  chr_colors <- chr_colors[as.character(chromosomes_to_plot)]
  
  # Generate chromosome ideograms with enhanced directional support
  create_chromosome_ideogram <- function(chr_num) {
    # Get chromosome-specific data
    ideogram_data <- create_single_ideogram(chr_num, arm_definitions, cytoband_data, chr_colors)
    if (is.null(ideogram_data)) return(NULL)

    plot_data <- ideogram_data$plot_data
    chr_color <- chr_colors[as.character(chr_num)]

    # Get functional arm definitions for this chromosome
    chr_arm_defs <- arm_definitions %>% filter(chr_num == !!chr_num)

    # Create base plot with chromosome outline
    base_plot <- ggplot() +
      # Background chromosome arms (semi-transparent)
      geom_rect(data = plot_data,
                aes(xmin = -0.3, xmax = 0.3, ymin = ymin, ymax = ymax),
                fill = chr_color, alpha = 0.5, color = "black", size = 0.5) +

      # Chromosome boundary coordinate labels
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

      # Arm labels (p and q)
      geom_text(data = plot_data,
                aes(x = 0, y = (ymin + ymax)/2, label = arm),
                size = 4, fontface = "bold")

    # Handle multiple functional regions per chromosome
    if (nrow(chr_arm_defs) > 1) {
      # Multiple functional regions - create separate plots
      functional_plots <- create_functional_plots(chr_num, chr_arm_defs, plot_data, chr_color)

      # Combine plots into panel layout
      n_plots <- length(functional_plots)
      if (n_plots <= 2) {
        combined_plot <- wrap_plots(functional_plots, nrow = 1)
      } else {
        combined_plot <- wrap_plots(functional_plots, ncol = 2)
      }

      return(combined_plot)

    } else if (nrow(chr_arm_defs) == 1) {
      # Single functional region - create integrated plot
      positions <- calculate_arm_positions(chr_arm_defs, plot_data)

      if (length(positions) > 0) {
        pos <- positions[[1]]
        arm_def <- pos$arm_def

        # Add functional region overlay to base plot
        final_plot <- base_plot +
          geom_rect(aes(xmin = -0.3, xmax = 0.3,
                        ymin = pos$func_ymin, ymax = pos$func_ymax),
                    fill = chr_color, alpha = 1.0, color = "black", size = 0.5) +

          # Functional breakpoint coordinate (in red)
          geom_text(aes(x = pos$label_x, y = pos$label_y,
                        label = format(pos$label_coord, big.mark = ",")),
                    hjust = pos$label_hjust, size = 3, color = "red", fontface = "bold") +

          # Direction label with telcent annotation
          geom_text(aes(x = 1.8, y = 0.5,
                        label = paste0(toupper(substring(arm_def$direction, 1, 1)),
                                       substring(arm_def$direction, 2),
                                       " (", arm_def$telcent, ")")),
                    hjust = 1, vjust = 0, size = 2.5,
                    color = if(arm_def$direction == "amp") "red" else "blue",
                    fontface = "bold")

        # Apply final styling
        final_plot <- final_plot +
          theme_minimal() +
          labs(title = paste("Chromosome", chr_num)) +
          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                panel.grid = element_blank()) +
          xlim(-2, 2) + ylim(0, 10)

        return(final_plot)
      }
    }

    # Return base plot if no functional regions
    final_plot <- base_plot +
      theme_minimal() +
      labs(title = paste("Chromosome", chr_num)) +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            panel.grid = element_blank()) +
      xlim(-2, 2) + ylim(0, 10)

    return(final_plot)
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



