#' Summarise Arm Definitions with Multiple Visualizations
#'
#' Creates multiple visualization plots summarizing defined chromosome arms:
#' 1. Nested pie chart by selection type, direction, and chromosomal location
#' 2. Scatter plot of percentage_length vs percent_sample_support with ellipses
#' 3. Combined plot: ridge plot (percentage_length) + violin/boxplot (percent_sample_support)
#'
#' @param arm_definitions Data frame containing arm definition data with required columns:
#'   peak_id, type, direction, negpos, type_of_selection, percentage_length, percent_sample_support
#' @param output_dir Character, path to output directory (default: NULL, no saving)
#' @param dir_name Character, subdirectory name within output_dir (default: "arms_summary")
#' @param save_plots Logical, whether to save plots to disk (default: FALSE)
#'
#' @return A list containing:
#'   \item{nested_pie_plot}{ggplot2 object - nested pie chart}
#'   \item{scatter_plot}{ggplot2 object - scatter plot with 95% confidence ellipses}
#'   \item{ridge_plot}{patchwork object - ridge plot + violin/boxplot showing distributions by selection type}
#' @export
#'
#' @examples
#' \dontrun{
#' plots <- summarise_arm(arm_definitions)
#' print(plots$nested_pie_plot)
#' print(plots$scatter_plot)
#' print(plots$ridge_plot)
#'
#' # Save plots to disk
#' plots <- summarise_arm(
#'   arm_definitions = arm_definitions,
#'   output_dir = "/path/to/output",
#'   dir_name = "arms_ideogram",
#'   save_plots = TRUE
#' )
#' }
summarise_arm <- function(arm_definitions,
                          output_dir = NULL,
                          dir_name = "arms_summary",
                          save_plots = FALSE) {
  
  # Validate required columns for nested pie plot
  required_cols_pie <- c("peak_id", "type", "direction", "negpos", "type_of_selection")
  missing_cols_pie <- setdiff(required_cols_pie, names(arm_definitions))
  if (length(missing_cols_pie) > 0) {
    stop("Missing required columns for pie plot: ", paste(missing_cols_pie, collapse = ", "))
  }
  
  # Validate required columns for scatter plot
  required_cols_scatter <- c("percentage_length", "percent_sample_support")
  missing_cols_scatter <- setdiff(required_cols_scatter, names(arm_definitions))
  if (length(missing_cols_scatter) > 0) {
    warning("Missing required columns for scatter plot: ", paste(missing_cols_scatter, collapse = ", "),
            "\nScatter plot will not be generated.")
    scatter_available <- FALSE
  } else {
    scatter_available <- TRUE
  }
  
  # Calculate total number of arms
  n_arms <- nrow(arm_definitions)
  
  # Prepare data for plotting
  plot_data <- arm_definitions %>%
    dplyr::select(peak_id, type, direction, negpos, type_of_selection) %>%
    dplyr::mutate(
      # Create combined categories
      level1 = ifelse(negpos == "neg", "Negative Selection", "Positive Selection"),
      level2 = paste0(negpos, "-", direction),
      level2_label = dplyr::case_when(
        negpos == "neg" & direction == "del" ~ "Essential-like",
        negpos == "neg" & direction == "amp" ~ "Toxic-like",
        negpos == "pos" & direction == "del" ~ "TS-like",
        negpos == "pos" & direction == "amp" ~ "Onco-like",
        TRUE ~ "Other"
      ),
      level3 = paste0(level2, "-", type),
      level3_label = type
    )
  
  # Level 1: Negative vs Positive selection
  level1_counts <- plot_data %>%
    dplyr::group_by(level1) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      fraction = n / sum(n),
      ymax = cumsum(fraction),
      ymin = dplyr::lag(ymax, default = 0),
      ymid = (ymin + ymax) / 2
    )
  
  # Level 2: Within each selection type, count amp/del events
  # This should be calculated WITHIN each level1 group
  level2_counts <- plot_data %>%
    dplyr::group_by(level1) %>%
    dplyr::mutate(
      level1_n = dplyr::n(),
      level1_fraction = level1_n / n_arms
    ) %>%
    dplyr::group_by(level1, level2, level2_label, level1_fraction) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(level1) %>%
    dplyr::mutate(
      # Fraction within this level1 group
      fraction_within_level1 = n / sum(n),
      # Adjust to match level1 slice size
      adjusted_fraction = fraction_within_level1 * unique(level1_fraction),
      ymax = cumsum(adjusted_fraction),
      ymin = dplyr::lag(ymax, default = 0),
      ymid = (ymin + ymax) / 2
    ) %>%
    dplyr::ungroup() %>%
    # Shift each level1 group to start at correct position
    dplyr::left_join(
      level1_counts %>% dplyr::select(level1, level1_ymin = ymin),
      by = "level1"
    ) %>%
    dplyr::mutate(
      ymin = ymin + level1_ymin,
      ymax = ymax + level1_ymin,
      ymid = (ymin + ymax) / 2
    )
  
  # Level 3: Within each level2 category, count by type (p_tel, q_tel, etc.)
  level3_counts <- plot_data %>%
    dplyr::group_by(level1, level2, level2_label, level3, level3_label) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(level1, level2) %>%
    dplyr::mutate(
      # Fraction within this level2 group
      fraction_within_level2 = n / sum(n)
    ) %>%
    dplyr::ungroup() %>%
    # Get level2 slice size from level2_counts
    dplyr::left_join(
      level2_counts %>% dplyr::select(level1, level2, level2_ymin = ymin, level2_ymax = ymax),
      by = c("level1", "level2")
    ) %>%
    dplyr::mutate(
      level2_size = level2_ymax - level2_ymin
    ) %>%
    dplyr::group_by(level1, level2) %>%
    dplyr::mutate(
      # Scale fraction to level2 slice size
      adjusted_fraction = fraction_within_level2 * unique(level2_size),
      ymax = cumsum(adjusted_fraction),
      ymin = dplyr::lag(ymax, default = 0)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      # Shift to level2 starting position
      ymin = ymin + level2_ymin,
      ymax = ymax + level2_ymin,
      ymid = (ymin + ymax) / 2
    )
  
  # Define color palettes
  level1_colors <- c("Negative Selection" = "#E41A1C", "Positive Selection" = "#377EB8")
  
  level2_colors <- c(
    "neg-del" = "#A6CEE3",  # Light blue (Essential-like)
    "neg-amp" = "#FB9A99",  # Light red (Toxic-like)
    "pos-del" = "#FDBF6F",  # Light orange (TS-like)
    "pos-amp" = "#B2DF8A"   # Light green (Onco-like)
  )
  
  # Colors for level 3 (type: p_tel, q_tel, p_cent, q_cent)
  level3_colors <- c(
    "p_tel" = "#8DD3C7",
    "q_tel" = "#FFFFB3",
    "p_cent" = "#BEBADA",
    "q_cent" = "#FB8072"
  )
  
  # ============================================================================
  # PLOT 1: Nested Pie Chart
  # ============================================================================
  
  nested_pie_plot <- ggplot2::ggplot() +
    # Level 3 (outermost ring) - by type
    ggplot2::geom_rect(
      data = level3_counts,
      ggplot2::aes(xmin = 3.5, xmax = 4.5, ymin = ymin, ymax = ymax, fill = level3_label),
      color = "white", linewidth = 0.3
    ) +
    # Level 2 (middle ring) - by direction within selection
    ggplot2::geom_rect(
      data = level2_counts,
      ggplot2::aes(xmin = 2.5, xmax = 3.5, ymin = ymin, ymax = ymax, fill = level2),
      color = "white", linewidth = 0.5
    ) +
    # Level 1 (inner ring) - selection type
    ggplot2::geom_rect(
      data = level1_counts,
      ggplot2::aes(xmin = 1.5, xmax = 2.5, ymin = ymin, ymax = ymax, fill = level1),
      color = "white", linewidth = 0.5
    ) +
    # Center circle with count
    ggplot2::annotate(
      "text", x = 0, y = 0.5,
      label = paste0("N Arms\n= ", n_arms),
      size = 3, fontface = "bold", color = "gray30"
    ) +
    # Add labels for level 3 (location types)
    ggplot2::geom_text(
      data = level3_counts,
      ggplot2::aes(x = 4, y = ymid, label = paste0(n)),
      size = 2.5, color = "gray20"
    ) +
    # Add labels for level 2 (selection types)
    ggplot2::geom_text(
      data = level2_counts,
      ggplot2::aes(x = 3, y = ymid, label = paste0(level2_label, "\n(n=", n, ")")),
      size = 2.5, fontface = "bold", color = "gray20"
    ) +
    # Add labels for level 1
    ggplot2::geom_text(
      data = level1_counts,
      ggplot2::aes(x = 2, y = ymid, label = paste0(n)),
      size = 3, fontface = "bold", color = "white"
    ) +
    # Scales and styling
    ggplot2::scale_fill_manual(
      values = c(level1_colors, level2_colors, level3_colors),
      breaks = c(names(level1_colors), names(level2_colors), names(level3_colors)),
      labels = c(
        "Negative Selection", "Positive Selection",
        "Essential-like (neg-del)", "Toxic-like (neg-amp)",
        "TS-like (pos-del)", "Onco-like (pos-amp)",
        "p arm telomere", "q arm telomere", "p arm centromere", "q arm centromere"
      )
    ) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(c(0, 4.5)) +
    ggplot2::labs(
      title = "Summary of Defined Chromosome Arms",
      subtitle = paste0("Total Arms: ", n_arms, " | Nested by Selection Type, Direction, and Location")
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5, color = "gray30"),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 9),
      legend.key.size = unit(0.4, "cm"), 
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        title = "Categories",
        ncol = 1
      )
    )
  
  # ============================================================================
  # PLOT 2: Scatter Plot (percentage_length vs percent_sample_support)
  # ============================================================================
  
  scatter_plot <- NULL
  
  if (scatter_available) {
    # Prepare data with selection type labels
    scatter_data <- arm_definitions %>%
      dplyr::mutate(
        selection_label = dplyr::case_when(
          negpos == "neg" & direction == "del" ~ "Essential-like",
          negpos == "neg" & direction == "amp" ~ "Toxic-like",
          negpos == "pos" & direction == "del" ~ "TS-like",
          negpos == "pos" & direction == "amp" ~ "Onco-like",
          TRUE ~ "Other"
        )
      )
    
    # Create scatter plot
    scatter_plot <- ggplot2::ggplot(
      scatter_data,
      ggplot2::aes(x = percentage_length, y = percent_sample_support, color = selection_label, fill = selection_label)
    ) +
      ggplot2::stat_ellipse(
        geom = "polygon",
        alpha = 0.15,
        level = 0.8,
        show.legend = FALSE
      ) +
      ggplot2::geom_point(size = 2.5, alpha = 1) +
      ggplot2::scale_color_manual(
        values = c(
          "Essential-like" = "#A6CEE3",
          "Toxic-like" = "#FB9A99",
          "TS-like" = "#FDBF6F",
          "Onco-like" = "#B2DF8A",
          "Other" = "gray50"
        )
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "Essential-like" = "#A6CEE3",
          "Toxic-like" = "#FB9A99",
          "TS-like" = "#FDBF6F",
          "Onco-like" = "#B2DF8A",
          "Other" = "gray50"
        ),
        guide = "none"
      ) +
      ggplot2::labs(
        title = "Arm Event Characteristics",
        subtitle = paste0("Functional Length vs Sample Support (n = ", n_arms, " arms)"),
        x = "Percentage of Arm Length",
        y = "Percent Sample Support",
        color = "Selection Type"
      ) +
      ggplot2::scale_x_continuous(
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0", "25", "50", "75", "100"),
        limits = c(-0.25, 
                   max(round(max(scatter_data$percentage_length), digits = 2) + 0.25, 1))
      ) +
      ggplot2::scale_y_continuous(
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0", "25", "50", "75", "100"),
        limits = c(-0.25, 
                   max(round(max(scatter_data$percent_sample_support), digits = 2) + 0.25, 1))
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5, color = "gray30"),
        axis.title = ggplot2::element_text(size = 10, face = "bold"),
        axis.text = ggplot2::element_text(size = 9),
        legend.title = ggplot2::element_text(size = 10, face = "bold"),
        legend.text = ggplot2::element_text(size = 9),
        legend.key.size = ggplot2::unit(0.4, "cm"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "gray30", linewidth = 0.5),
        plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
        panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
        legend.background = ggplot2::element_rect(fill = "transparent", color = NA)
      )
  }
  
  # ============================================================================
  # PLOT 3: Ridge Plots (percentage_length and percent_sample_support)
  # ============================================================================
  
  ridge_plot <- NULL
  
  if (scatter_available) {
    # Check if ggridges is available
    if (!requireNamespace("ggridges", quietly = TRUE)) {
      warning("ggridges package required for ridge plots. Install with: install.packages('ggridges')")
    } else {
      # Reorder selection_label factor for consistent ordering
      scatter_data <- scatter_data %>%
        dplyr::mutate(
          selection_label = factor(selection_label,
                                   levels = c("Onco-like", "TS-like", "Toxic-like", "Essential-like"))
        )
      
      # Calculate medians for percentage_length by selection type
      medians_percentage <- scatter_data %>%
        dplyr::group_by(selection_label) %>%
        dplyr::summarise(median_val = median(percentage_length, na.rm = TRUE), .groups = "drop")
      
      # Perform Kruskal-Wallis test for percentage_length
      kw_p_percentage <- tryCatch({
        kw_test <- kruskal.test(percentage_length ~ selection_label, data = scatter_data)
        kw_test$p.value
      }, error = function(e) {
        if (grepl("all observations are in the same group", e$message)) {
          cat("⚠️ Skipping Kruskal-Wallis test for percentage_length: all observations in same group\n")
          return(NULL)
        } else {
          stop(e)
        }
      })
      
      # Perform pairwise comparisons if Kruskal-Wallis is significant
      pairwise_caption_percentage <- ""
      if (!is.null(kw_p_percentage) && kw_p_percentage < 0.05) {
        pairwise_results <- pairwise.wilcox.test(
          scatter_data$percentage_length,
          scatter_data$selection_label,
          p.adjust.method = "BH"
        )
        
        # Identify which groups differ significantly from others
        p_matrix <- pairwise_results$p.value
        sig_groups <- c()
        
        for (i in rownames(p_matrix)) {
          # Count how many comparisons are significant for this group
          # Check row comparisons
          sig_count <- sum(p_matrix[i, ] < 0.05, na.rm = TRUE)
          
          # Check column comparisons (only if this group appears in column names)
          if (i %in% colnames(p_matrix)) {
            sig_count <- sig_count + sum(p_matrix[, i] < 0.05, na.rm = TRUE)
          }
          
          if (sig_count > 0) {
            sig_groups <- c(sig_groups, i)
          }
        }
        
        # Also check if any column-only groups are significant
        for (j in colnames(p_matrix)) {
          if (!(j %in% sig_groups) && sum(p_matrix[, j] < 0.05, na.rm = TRUE) > 0) {
            sig_groups <- c(sig_groups, j)
          }
        }
        
        if (length(sig_groups) > 0) {
          pairwise_caption_percentage <- paste0(
            "\nSignificant differences: ", paste(sig_groups, collapse = ", ")
          )
        }
      }
      
      caption_percentage <- if (!is.null(kw_p_percentage)) {
        paste0(
          "Kruskal-Wallis p = ", format.pval(kw_p_percentage, digits = 3),
          pairwise_caption_percentage
        )
      } else {
        "Statistical test not applicable (single group)"
      }
      
      base_cols <- c(
        "Essential-like" = "#A6CEE3",
        "Toxic-like"     = "#FB9A99",
        "TS-like"        = "#FDBF6F",
        "Onco-like"      = "#B2DF8A"
      )
      
      scatter_data <- scatter_data %>%
        mutate(
          selection_label = factor(selection_label,
                                   levels = c("Onco-like", "TS-like", "Toxic-like", "Essential-like")),
          telcent         = factor(telcent, levels = c("tel", "cent")),
          fill_key        = interaction(selection_label, telcent, sep = ".", drop = TRUE)
        )
      
      # Ridge plot for percentage_length
      ridge_percentage_length_sum <- ggplot2::ggplot(scatter_data, 
        ggplot2::aes(x = percentage_length,
                     y = selection_label,
                     fill = selection_label)) +
        ggridges::geom_density_ridges(alpha = 0.7, scale = 0.9) +
        ggplot2::geom_segment(
          data = medians_percentage,
          ggplot2::aes(x = median_val, xend = median_val, y = as.numeric(selection_label),
                       yend = as.numeric(selection_label) + 0.9),
          color = "black", linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE
        ) +
        ggplot2::scale_fill_manual(
          values = c("Essential-like" = "#A6CEE3",
                     "Toxic-like"     = "#FB9A99",
                     "TS-like"        = "#FDBF6F",
                     "Onco-like"      = "#B2DF8A")) +
        # ggplot2::scale_x_continuous(
        #   breaks = c(0, 0.25, 0.5, 0.75, 1),
        #   labels = c("0", "25", "50", "75", "100"),
        #   limits = c(0, 1)
        # ) +
        ggplot2::labs(
          title = "Percentage of Arm Length",
          x = "Percentage of Arm Length",
          y = NULL,
          caption = caption_percentage
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
          plot.caption = ggplot2::element_text(size = 8, hjust = 0, color = "gray40"),
          axis.title = ggplot2::element_text(size = 10, face = "bold"),
          axis.text = ggplot2::element_text(size = 9),
          legend.position = "none",
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(color = "gray30", linewidth = 0.5),
          plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
          panel.background = ggplot2::element_rect(fill = "transparent", color = NA)
        ) 
      
      
      # Ridge plot for percentage_length
      ridge_percentage_length <- ggplot2::ggplot(scatter_data) +
        # ggplot2::aes(x = percentage_length, 
        #              y = selection_label, 
        #              fill = selection_label, 
        #              group = interaction(selection_label, telcent))) +
        # ggridges::geom_density_ridges(alpha = 0.7, scale = 0.9) +
        ggridges::stat_density_ridges(aes(x = percentage_length,
                                          y = selection_label,
                                          fill = fill_key,
                                          group = interaction(selection_label, telcent)),
                                      quantiles       = 0.5,
                                      quantile_lines  = TRUE,
                                      geom            = "density_ridges",
                                      scale           = 0.9,
                                      rel_min_height  = 0.001) + 
        # ggplot2::geom_segment(
        #   data = medians_percentage,
        #   ggplot2::aes(x = median_val, xend = median_val, y = as.numeric(selection_label),
        #                yend = as.numeric(selection_label) + 0.9),
        #   color = "black", linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE
        # ) +
        ggplot2::scale_fill_manual(
          values = c(setNames(alpha(darken(base_cols, amount = 0.35), 0.7), paste0(names(base_cols), ".tel")),
                     setNames(alpha(base_cols, 0.7), paste0(names(base_cols), ".cent")))) +
        # ggplot2::scale_x_continuous(
        #   breaks = c(0, 0.25, 0.5, 0.75, 1),
        #   labels = c("0", "25", "50", "75", "100"),
        #   limits = c(0, 1)
        # ) +
        ggplot2::labs(
          title = "Percentage of Arm Length",
          x = "Percentage of Arm Length",
          y = NULL,
          subtitle = "LightShade=CENT;DarkShade=TEL"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
          plot.caption = ggplot2::element_text(size = 8, hjust = 0, color = "gray40"),
          axis.title = ggplot2::element_text(size = 10, face = "bold"),
          axis.text = ggplot2::element_text(size = 9),
          legend.position = "none",
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(color = "gray30", linewidth = 0.5),
          plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
          panel.background = ggplot2::element_rect(fill = "transparent", color = NA)
        ) 
      
      
      # Perform Kruskal-Wallis test for percent_sample_support
      kw_p_support <- tryCatch({
        kw_test <- kruskal.test(percent_sample_support ~ selection_label, data = scatter_data)
        kw_test$p.value
      }, error = function(e) {
        if (grepl("all observations are in the same group", e$message)) {
          cat("⚠️ Skipping Kruskal-Wallis test for percent_sample_support: all observations in same group\n")
          return(NULL)
        } else {
          stop(e)
        }
      })
      
      # Perform pairwise comparisons if Kruskal-Wallis is significant
      pairwise_caption_support <- ""
      if (!is.null(kw_p_support) && kw_p_support < 0.05) {
        pairwise_results <- pairwise.wilcox.test(
          scatter_data$percent_sample_support,
          scatter_data$selection_label,
          p.adjust.method = "BH"
        )
        
        # Identify which groups differ significantly from others
        p_matrix <- pairwise_results$p.value
        sig_groups <- c()
        
        for (i in rownames(p_matrix)) {
          # Count how many comparisons are significant for this group
          # Check row comparisons
          sig_count <- sum(p_matrix[i, ] < 0.05, na.rm = TRUE)
          
          # Check column comparisons (only if this group appears in column names)
          if (i %in% colnames(p_matrix)) {
            sig_count <- sig_count + sum(p_matrix[, i] < 0.05, na.rm = TRUE)
          }
          
          if (sig_count > 0) {
            sig_groups <- c(sig_groups, i)
          }
        }
        
        # Also check if any column-only groups are significant
        for (j in colnames(p_matrix)) {
          if (!(j %in% sig_groups) && sum(p_matrix[, j] < 0.05, na.rm = TRUE) > 0) {
            sig_groups <- c(sig_groups, j)
          }
        }
        
        if (length(sig_groups) > 0) {
          pairwise_caption_support <- paste0(
            "\nSignificant differences: ", paste(sig_groups, collapse = ", ")
          )
        }
      }
      
      caption_support <- if (!is.null(kw_p_support)) {
        paste0(
          "Kruskal-Wallis p = ", format.pval(kw_p_support, digits = 3),
          pairwise_caption_support
        )
      } else {
        "Statistical test not applicable (single group)"
      }
      
      # Violin plot with boxplot for percent_sample_support
      violin_boxplot_support <- ggplot2::ggplot(
        scatter_data,
        ggplot2::aes(x = selection_label, y = percent_sample_support, fill = selection_label)
      ) +
        ggplot2::geom_violin(alpha = 0.7, trim = FALSE) +
        ggplot2::geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.size = 1) +
        ggplot2::scale_fill_manual(
          values = c(
            "Essential-like" = "#A6CEE3",
            "Toxic-like" = "#FB9A99",
            "TS-like" = "#FDBF6F",
            "Onco-like" = "#B2DF8A"
          )
        ) +
        ggplot2::scale_y_continuous(
          breaks = c(0, 0.25, 0.5, 0.75, 1),
          labels = c("0", "25", "50", "75", "100"),
          limits = c(0, 1)
        ) +
        ggplot2::labs(
          title = "Percent Sample Support",
          x = NULL,
          y = "Percent Sample Support",
          caption = caption_support
        ) +
        ggplot2::coord_flip() +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
          plot.caption = ggplot2::element_text(size = 8, hjust = 0, color = "gray40"),
          axis.title = ggplot2::element_text(size = 10, face = "bold"),
          axis.text = ggplot2::element_text(size = 9),
          legend.position = "none",
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(color = "gray30", linewidth = 0.5),
          plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
          panel.background = ggplot2::element_rect(fill = "transparent", color = NA)
        )
      
      # Combine using patchwork
      if (!requireNamespace("patchwork", quietly = TRUE)) {
        warning("patchwork package required for combined plots. Install with: install.packages('patchwork')")
        ridge_plot <- list(
          percentage_length = ridge_percentage_length,
          sample_support = violin_boxplot_support
        )
      } else {
        ridge_plot <- ridge_percentage_length_sum + ridge_percentage_length + violin_boxplot_support + 
          patchwork::plot_annotation(
            title = "Distribution of Arm Event Characteristics by Selection Type",
            theme = ggplot2::theme(
              plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5)
            )
          )
      }
    }
  }
  
  # ============================================================================
  # Save plots if requested
  # ============================================================================
  
  if (save_plots) {
    if (is.null(output_dir)) {
      stop("output_dir must be provided when save_plots = TRUE")
    }
    
    # Create output directory
    plot_dir <- file.path(output_dir, dir_name)
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Save nested pie plot
    ggplot2::ggsave(
      filename = file.path(plot_dir, "nested_pie_plot.png"),
      plot = nested_pie_plot,
      width = 5.1,
      height = 2.5,
      units = "in",
      dpi = 600,
      bg = "white"
    )
    cat("✅ Saved nested_pie_plot.png\n")
    
    # Save scatter plot (if available)
    if (!is.null(scatter_plot)) {
      ggplot2::ggsave(
        filename = file.path(plot_dir, "scatter_plot.png"),
        plot = scatter_plot,
        width = 4.2,
        height = 3,
        units = "in",
        dpi = 600,
        bg = "white"
      )
      cat("✅ Saved scatter_plot.png\n")
    }
    
    # Save ridge plot (if available)
    if (!is.null(ridge_plot)) {
      ggplot2::ggsave(
        filename = file.path(plot_dir, "ridge_plot.png"),
        plot = ridge_plot,
        width = 8.5,
        height = 3.7,
        units = "in",
        dpi = 600,
        bg = "white"
      )
      cat("✅ Saved ridge_plot.png\n")
    }
    
    cat("📁 Plots saved to:", plot_dir, "\n")
  }
  
  # ============================================================================
  # Return list of plots
  # ============================================================================
  
  return(list(
    nested_pie_plot = nested_pie_plot,
    scatter_plot = scatter_plot,
    ridge_plot = ridge_plot
  ))
}
