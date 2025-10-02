#' Plot Chromosome Ideograms with Enhanced Functional Arm Support
#'
#' Creates comprehensive ideogram plots showing chromosome arms with functional regions.
#' Supports both telomere-direction (telomere to breakpoint) and centromere-direction
#' (centromere to breakpoint) arm definitions with proper p/q arm directional handling.
#' Individual plots for each chromosome and combined panel plots are generated.
#'
#' @param arm_definitions arm_definitions with telcent column ("telomere" or "centromere")
#' @param output_dir Character, directory to save plots (default: current directory)
#' @param cytoband_data Data frame, cytoband data for chromosome boundaries (default: BISCUT coordinates)
#' @param plot_width Numeric, width of individual plots in inches (default: 3)
#' @param plot_height Numeric, height of individual plots in inches (default: 3)
#' @param save_plots Logical, whether to save plots to files (default: TRUE)
#' @return List containing individual plots, combined plot, colors, and output directory
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggchicklet geom_rrect
#' @export
plot_ideograms <- function(arm_definitions,
                           output_dir = ".",
                           dir_name = "arms_ideogram", 
                           plot_width = 3,
                           plot_height = 3,
                           save_plots = TRUE) {
  
  
  # Validate input data
  if (is.null(arm_definitions) || nrow(arm_definitions) == 0) {
    stop("No arm definitions found in arm_definitions")
  }
  
  # Create output directory
  ideogram_dir <- NULL
  if (save_plots) {
    ideogram_dir <- file.path(output_dir, dir_name)
    if (!dir.exists(ideogram_dir)) {
      dir.create(ideogram_dir, recursive = TRUE)
      cat("Created ideogram directory:", ideogram_dir, "\n")
    }
  }
  
  # Setup chromosome colors
  chromosomes_to_plot <- sort(unique(arm_definitions$chr))
  set.seed(999)
  chr_colors <- sample(BAGEL::bagel_palette, size = 23, replace = FALSE)
  names(chr_colors) <- paste0("chr", 1:23)
  chr_colors <- chr_colors[as.character(chromosomes_to_plot)]
  
  # Create nested list structure for plots
  plots_by_chromosome <- list()
  combined_by_chromosome <- list()

  for (chr_num in chromosomes_to_plot) {
    cat("Creating ideogram for chromosome", chr_num, "\n")

    # Extract arm data (all events for this chromosome)
    plot_data <- arm_definitions %>%
      dplyr::filter(chr == chr_num)

    # Initialize event plots list for this chromosome
    event_plots_list <- list()
    
    # Make base-plot using direct coordinates ---------------------------------
    background <- plot_data %>%
      dplyr::select(p_start:q_end) %>%
      distinct()

    # Extract genomic coordinates
    p_start <- background$p_start
    p_end <- background$p_end
    p_length <- p_end - p_start
    q_start <- background$q_start
    q_end <- background$q_end
    q_length <- q_end - q_start

    # Centromere coordinates
    centromere_start <- p_end
    centromere_end <- q_start
    centromere_length <- centromere_end - centromere_start

    # get colour
    chr_color <- chr_colors[as.character(chr_num)]

    # Create base plot using direct genomic coordinates
    # P arm: centromere_length (centromere) to (centromere_length + p_length) (telomere) - ABOVE
    # Centromere: 0 to centromere_length - MIDDLE
    # Q arm: -(q_length) (telomere) to 0 (centromere) - BELOW
    base_plot <- ggplot() +
      # Centromere (distinct color, full opacity) - MIDDLE
      ggchicklet::geom_rrect(aes(xmin = -0.3, xmax = 0.3,
                    ymin = 0, ymax = centromere_length),
                fill = "white", alpha = 1.0, color = "white", linewidth = 0.5,
                r = unit(0.3, 'npc')) +

      # P arm background (50% transparency) - ABOVE centromere
      ggchicklet::geom_rrect(aes(xmin = -0.3, xmax = 0.3,
                    ymin = centromere_length, ymax = centromere_length + p_length),
                fill = chr_color, alpha = 0.5, color = "black", linewidth = 0.5,
                r = unit(0.3, 'npc')) +

      # Q arm background (50% transparency) - BELOW centromere
      ggchicklet::geom_rrect(aes(xmin = -0.3, xmax = 0.3, ymin = -q_length, ymax = 0),
                fill = chr_color, alpha = 0.5, color = "black", linewidth = 0.5,
                r = unit(0.3, 'npc')) +

      # P arm boundary labels
      geom_text(aes(x = 0.4, y = centromere_length + p_length,
                    label = format(p_start, big.mark = ",")),
                size = 3, hjust = 0) +  # P telomere (top)
      geom_text(aes(x = 0.4, y = centromere_length,
                    label = format(p_end, big.mark = ",")),
                size = 3, hjust = 0) +  # P-centromere boundary

      # # Centromere boundary labels
      # geom_text(aes(x = -0.4, y = centromere_length,
      #               label = format(centromere_start, big.mark = ",")),
      #           size = 3, hjust = 1, color = "white") +  # Centromere start
      # geom_text(aes(x = -0.4, y = 0,
      #               label = format(centromere_end, big.mark = ",")),
      #           size = 3, hjust = 1, color = "white") +  # Centromere end

      # Q arm boundary labels
      geom_text(aes(x = 0.4, y = 0,
                    label = format(q_start, big.mark = ",")),
                size = 3, hjust = 0) +  # Q-centromere boundary
      geom_text(aes(x = 0.4, y = -q_length,
                    label = format(q_end, big.mark = ",")),
                size = 3, hjust = 0)  # Q telomere (bottom)

    # Store arm label positions for adding on top later
    arm_labels_data <- data.frame(
      x = c(0, 0, 0),
      y = c(centromere_length + p_length / 2, centromere_length / 2, -q_length / 2),
      label = c("p", "cen", "q"),
      size = c(4, 3, 4),
      color = c("black", "black", "black")
    )
    
    # Add functional regions using direct coordinates -------------------------
    if (nrow(plot_data) > 0) {
      for (i in 1:nrow(plot_data)) {
        arm_def <- plot_data[i, ]
        arm_type <- arm_def$arm

        # Get functional region coordinates
        func_start <- arm_def$functional_start
        func_end <- arm_def$functional_end

        # Calculate plot coordinates using direct genomic positions
        if (arm_type == "p") {
          # P arm: centromere_length (bottom) to centromere_length + p_length (top)
          # Convert genomic coords to plot coords
          func_y_start <- centromere_length + (p_end - func_start)
          func_y_end <- centromere_length + (p_end - func_end)

          # Ensure correct order (ymin < ymax)
          func_ymin <- min(func_y_start, func_y_end)
          func_ymax <- max(func_y_start, func_y_end)

        } else {
          # Q arm: -q_length (bottom) to 0 (top at centromere)
          # Convert genomic coords to plot coords
          func_y_start <- -(func_start - q_start)
          func_y_end <- -(func_end - q_start)

          # Ensure correct order (ymin < ymax)
          func_ymin <- min(func_y_start, func_y_end)
          func_ymax <- max(func_y_start, func_y_end)
        }

        # Create data frame for event-specific annotations
        label_coord <- if (arm_def$telcent == "tel") func_end else func_start
        label_x <- if (arm_type == "p") 0.4 else -0.4
        label_hjust <- if (arm_type == "p") 0 else 1
        label_y <- (func_ymin + func_ymax) / 2

        event_data <- data.frame(
          func_ymin = func_ymin,
          func_ymax = func_ymax,
          label_x = label_x,
          label_y = label_y,
          label_coord = label_coord,
          label_hjust = label_hjust,
          direction = arm_def$direction,
          telcent = arm_def$telcent,
          negpos = arm_def$negpos,
          annotation_color = if(arm_def$direction == "amp") "red" else "blue"
        )

        # Create event-specific plot
        event_plot <- base_plot +
          ggchicklet::geom_rrect(data = event_data,
                    aes(xmin = -0.3, xmax = 0.3,
                        ymin = func_ymin, ymax = func_ymax),
                    fill = chr_color, alpha = 1.0, color = "black", linewidth = 0.5,
                    r = unit(0.3, 'npc')) +

          # Add breakpoint coordinate label
          geom_text(data = event_data,
                    aes(x = label_x, y = label_y,
                        label = format(label_coord, big.mark = ",")),
                    hjust = event_data$label_hjust,
                    size = 3, color = "red", fontface = "bold") +

          # Direction and telcent annotation
          geom_text(data = event_data,
                    aes(x = 1.8, y = func_ymin,
                        label = paste0(direction, "\n", telcent, "\n", negpos)),
                    hjust = 1, vjust = 0, size = 3,
                    color = event_data$annotation_color,
                    fontface = "bold") +

          # Add arm labels on top (so they're not covered by event overlay)
          geom_text(data = arm_labels_data,
                    aes(x = x, y = y, label = label),
                    size = arm_labels_data$size,
                    color = arm_labels_data$color,
                    fontface = "bold")

        # Apply theme and styling
        event_plot <- event_plot +
          theme_minimal() +
          labs(title = paste("Chromosome", chr_num, "-", arm_def$peak_id)) +
          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
                panel.grid = element_blank()) +
          coord_cartesian(xlim = c(-2, 2))  # Let y-axis auto-scale

        # Create unique identifier
        event_id <- paste(arm_def$chr_arm, arm_def$direction, arm_def$telcent, arm_def$negpos, sep = "_")

        # Store plot in event_plots_list
        event_plots_list[[event_id]] <- event_plot
      }
    }

    # Store event plots for this chromosome in nested structure
    chr_key <- as.character(chr_num)
    plots_by_chromosome[[chr_key]] <- event_plots_list

    # Create combined grid for this chromosome (n rows × 4 cols) ---------------
    if (length(event_plots_list) > 0) {
      ncol <- 4
      nrow <- ceiling(length(event_plots_list) / ncol)

      chr_combined_plot <- patchwork::wrap_plots(event_plots_list, ncol = ncol)

      # Store combined plot
      combined_by_chromosome[[chr_key]] <- chr_combined_plot

      # Save per-chromosome grid
      if (save_plots) {
        chr_grid_filename <- file.path(ideogram_dir,
                                        paste0(chr_num, "_all_events_grid.png"))
        ggsave(filename = chr_grid_filename,
               plot = chr_combined_plot,
               width = ncol * plot_width,
               height = nrow * plot_height,
               units = "in", dpi = 300)
        cat("Saved chromosome", chr_num, "grid with", length(event_plots_list),
            "events:", chr_grid_filename, "\n")
      }
    }
  }
  
  # Create optional mega-grid of all chromosome grids ------------------------
  all_chromosomes_combined <- NULL
  if (length(combined_by_chromosome) > 0) {
    # Flatten all event plots from all chromosomes
    all_event_plots <- unlist(plots_by_chromosome, recursive = FALSE)
    n_plots <- length(all_event_plots)

    if (n_plots > 0) {
      ncol <- round(sqrt(n_plots))
      nrow <- ceiling(n_plots / ncol)
      all_chromosomes_combined <- patchwork::wrap_plots(all_event_plots, ncol = ncol)

      # # Save mega combined plot
      # if (save_plots) {
      #   mega_filename <- file.path(ideogram_dir, "all_chromosomes_all_events_combined.png")
      #   ggsave(filename = mega_filename,
      #          plot = all_chromosomes_combined,
      #          width = ncol * plot_width,
      #          height = nrow * plot_height,
      #          units = "in", dpi = 300)
      #   cat("Saved mega-grid with all", n_plots, "events:", mega_filename, "\n")
      # }
    }
  }

  # Return nested list structure ----------------------------------------------
  return(list(
    # Level 1: Chromosome → Level 2: Individual event plots
    plots_by_chromosome = plots_by_chromosome,

    # Combined grids: one per chromosome (n rows × 4 cols)
    combined_by_chromosome = combined_by_chromosome,

    # Optional: mega-grid of all chromosomes
    all_chromosomes_combined = all_chromosomes_combined,

    # Metadata
    chromosome_colors = chr_colors,
    output_directory = if (save_plots) ideogram_dir else NULL,
    plot_metadata = list(
      chromosomes_plotted = names(plots_by_chromosome),
      events_per_chromosome = sapply(plots_by_chromosome, length),
      total_events = sum(sapply(plots_by_chromosome, length)),
      total_chromosomes = length(chromosomes_to_plot),
      functional_regions = nrow(arm_definitions),
      telcent_support = "telcent" %in% names(arm_definitions)
    )
  ))
}



