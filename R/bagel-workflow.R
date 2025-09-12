#' BAGEL v2.0 Workflow Functions
#' 
#' Complete workflow functions for BAGELv2 analysis including data loading,
#' matrix creation, and report generation.
#' 
#' @name bagel-workflow

library(dplyr)
library(readr)
library(tidyr)
library(stringr)

#' Load Pan-Cancer Breakpoint Data
#'
#' @return List containing all breakpoint datasets
#' @export
load_breakpoint_data <- function() {
    data_dir <- system.file("data", package = "BAGEL")
    
    # Load all breakpoint datasets
    load(file.path(data_dir, "consensus_arm_definitions.rda"))
    load(file.path(data_dir, "cancer_specific_breakpoints.rda"))
    load(file.path(data_dir, "arm_breakpoints.rda"))
    load(file.path(data_dir, "consensus_breakpoints.rda"))
    
    return(list(
        consensus_arm_definitions = consensus_arm_definitions,
        cancer_specific_breakpoints = cancer_specific_breakpoints,
        arm_breakpoints = arm_breakpoints,
        consensus_breakpoints = consensus_breakpoints
    ))
}

#' Get Available Cancer Types
#'
#' @param data_dir Character, path to data directory
#' @return Character vector of available cancer types
#' @export
get_available_cancer_types <- function(data_dir = NULL) {
    breakpoint_data <- load_breakpoint_data()
    available_types <- names(breakpoint_data$cancer_specific_breakpoints)
    
    # Also check for segmentation files if data_dir provided
    if (!is.null(data_dir) && dir.exists(data_dir)) {
        seg_dirs <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
        available_types <- intersect(available_types, seg_dirs)
    }
    
    return(available_types)
}

#' Load and Clean Segmentation Data for a Cancer Type
#'
#' @param cancer_type Character, cancer type name
#' @param data_dir Character, path to data directory
#' @return Cleaned segmentation data frame
#' @export
load_segmentation_data <- function(cancer_type, data_dir) {
    seg_file <- file.path(data_dir, cancer_type, "segmentation.seg")
    
    if (!file.exists(seg_file)) {
        stop(paste("Segmentation file not found for", cancer_type, "at", seg_file))
    }
    
    # Load data
    segments <- read_tsv(seg_file, show_col_types = FALSE)
    
    # Convert column name for consistency
    if ("Segment_Mean" %in% names(segments) && !"Log2Ratios" %in% names(segments)) {
        segments$Log2Ratios <- segments$Segment_Mean
    }
    
    # Clean data
    original_rows <- nrow(segments)
    segments_clean <- segments %>%
        filter(
            !is.na(Start) & !is.na(End) & !is.na(Log2Ratios),
            Start < End,
            Chromosome %in% 1:22
        ) %>%
        mutate(
            Chromosome = as.numeric(Chromosome),
            Segment_Length = End - Start + 1
        )
    
    cleaned_rows <- nrow(segments_clean)
    bagel_log(paste("Loaded", cancer_type, ":", original_rows, "->", cleaned_rows, "segments (",
        original_rows - cleaned_rows, "removed)"))
    
    return(segments_clean)
}

#' Get Arm Definitions for Cancer Type
#'
#' @param cancer_type Character, cancer type name
#' @param breakpoint_data List, loaded breakpoint data
#' @return Arm definitions data frame
#' @export
get_arm_definitions <- function(cancer_type, breakpoint_data) {
    
    # Get BISCUT chromosome coordinates for consistent boundaries
    tryCatch({
        biscut_coords <- BISCUT::get_chromosome_coordinates()
        bagel_log("Using BISCUT coordinates for chromosome arm boundaries")
    }, error = function(e) {
        warning("Could not load BISCUT coordinates: ", e$message)
        biscut_coords <- NULL
    })
    
    # Try cancer-specific first
    if (cancer_type %in% names(breakpoint_data$cancer_specific_breakpoints)) {
        bagel_log(paste("Using", cancer_type, "specific breakpoints"))
        cancer_breakpoints <- breakpoint_data$cancer_specific_breakpoints[[cancer_type]]
        
        # Convert to arm definitions
        if (nrow(cancer_breakpoints) > 0 && "peak_start" %in% names(cancer_breakpoints)) {
            
            if (!is.null(biscut_coords)) {
                # Use BISCUT coordinates for arm boundaries
                arm_definitions <- cancer_breakpoints %>%
                    mutate(
                        chr_num = as.numeric(chr),
                        arm_type = ifelse(grepl("p$", arm), "p", "q")
                    ) %>%
                    left_join(
                        biscut_coords %>%
                            select(chromosome_info, p_start, p_end, q_start, q_end) %>%
                            rename(chr_num = chromosome_info),
                        by = "chr_num"
                    ) %>%
                    mutate(
                        arm_start = case_when(
                            arm_type == "p" ~ p_start,
                            arm_type == "q" ~ peak_start  # Use functional breakpoint for q arm start
                        ),
                        arm_end = case_when(
                            arm_type == "p" ~ peak_end,   # Use functional breakpoint for p arm end
                            arm_type == "q" ~ q_end       # Use BISCUT end for q arm
                        )
                    ) %>%
                    filter(!is.na(arm_start) & !is.na(arm_end) & arm_start < arm_end) %>%
                    select(arm, chr_num, arm_type, arm_start, arm_end, direction) %>%
                    arrange(chr_num, arm_type)
                
            } else {
                # Fallback to original hardcoded coordinates
                arm_definitions <- cancer_breakpoints %>%
                    mutate(
                        chr_num = as.numeric(chr),
                        arm_type = ifelse(grepl("p$", arm), "p", "q"),
                        arm_start = case_when(
                            arm_type == "p" ~ 1,
                            arm_type == "q" ~ peak_start
                        ),
                        arm_end = case_when(
                            arm_type == "p" ~ peak_end,
                            arm_type == "q" ~ case_when(
                                chr_num == 1 ~ 249250621, chr_num == 2 ~ 243199373, chr_num == 3 ~ 198022430,
                                chr_num == 4 ~ 191154276, chr_num == 5 ~ 180915260, chr_num == 6 ~ 171115067,
                                chr_num == 7 ~ 159138663, chr_num == 8 ~ 146364022, chr_num == 9 ~ 141213431,
                                chr_num == 10 ~ 135534747, chr_num == 11 ~ 135006516, chr_num == 12 ~ 133851895,
                                chr_num == 13 ~ 115169878, chr_num == 14 ~ 107349540, chr_num == 15 ~ 102531392,
                                chr_num == 16 ~ 90354753, chr_num == 17 ~ 81195210, chr_num == 18 ~ 78077248,
                                chr_num == 19 ~ 59128983, chr_num == 20 ~ 63025520, chr_num == 21 ~ 48129895,
                                chr_num == 22 ~ 51304566, TRUE ~ NA_real_
                            )
                        )
                    ) %>%
                    filter(!is.na(arm_start) & !is.na(arm_end) & arm_start < arm_end) %>%
                    select(arm, chr_num, arm_type, arm_start, arm_end, direction) %>%
                    arrange(chr_num, arm_type)
            }
            
            return(arm_definitions)
        }
    }
    
    # Fallback to consensus (which now uses BISCUT coordinates)
    bagel_log(paste("Using consensus breakpoints for", cancer_type))
    return(breakpoint_data$consensus_arm_definitions)
}

#' Create Chromosome Arm Copy Number Matrix
#'
#' @param arm_summaries Data frame with arm-level summaries
#' @param output_dir Character, output directory
#' @return List with matrix and summary statistics
#' @export
create_arm_matrix <- function(arm_summaries, output_dir = NULL) {
    
    # Handle duplicates by taking the mean for each Arm-Sample combination
    arm_summary_clean <- arm_summaries %>%
        select(Sample, Arm, mean_log2ratio) %>%
        group_by(Sample, Arm) %>%
        summarise(mean_log2ratio = mean(mean_log2ratio, na.rm = TRUE), .groups = "drop")
    
    # Create arm-by-sample matrix
    arm_matrix_df <- arm_summary_clean %>%
        pivot_wider(names_from = Sample, values_from = mean_log2ratio, names_sort = TRUE)
    
    # Convert to matrix - handle potential NA values
    arm_names <- arm_matrix_df$Arm
    arm_matrix <- as.matrix(arm_matrix_df[, -1])
    
    # Check if conversion was successful
    if (any(sapply(arm_matrix_df[, -1], function(x) is.list(x)))) {
        stop("Matrix conversion failed: list columns detected after duplicate handling")
    }
    
    # Ensure all values are numeric
    if (!is.numeric(arm_matrix)) {
        warning("Non-numeric values detected in arm_matrix, attempting conversion")
        arm_matrix <- apply(arm_matrix, 2, as.numeric)
    }
    
    rownames(arm_matrix) <- arm_names
    
    # =============================================================================
    # GISTIC2-Style Copy Number Calculation Fix
    # 
    # ISSUE: Original BAGEL used absolute copy number conversion: 2 * (2^log2ratio)
    # This produces values centered around 2 (diploid) and always positive, which:
    # - Does not follow GISTIC2 conventions
    # - Makes survival analysis difficult (no negative values for deletions)
    # - Incompatible with standard copy number analysis tools
    #
    # FIX: Use GISTIC2-style calculation where:
    # - Log2 ratios are used directly (centered at 0)
    # - 0 = diploid (normal)
    # - Positive values = amplifications  
    # - Negative values = deletions
    # - Matches GISTIC2 parameters: amp_threshold=0.25, del_threshold=0.25
    # =============================================================================
    
    # OLD CODE - COMMENTED OUT:
    # Original BAGEL approach (absolute copy numbers, always positive)
    # copy_number_matrix <- 2 * (2^arm_matrix)
    
    # NEW CODE - GISTIC2-Style:
    # Use log2 ratios directly (GISTIC2 convention)
    copy_number_matrix <- arm_matrix  # Keep log2 ratios as-is, centered at 0
    
    # =============================================================================
    # GISTIC2-Style Summary Statistics Fix
    # 
    # ISSUE: Original summary used absolute copy number statistics (Mean_CN, etc.)
    # which don't make sense for log2-centered values
    #
    # FIX: Update to GISTIC2-style statistics:
    # - Use log2 ratio terminology
    # - Add alteration frequency calculations
    # - Apply GISTIC2 thresholds from parameter file (amp=0.25, del=0.25)
    # =============================================================================
    
    # OLD CODE - COMMENTED OUT:
    # Original absolute copy number statistics
    # arm_stats <- data.frame(
    #     Arm = rownames(copy_number_matrix),
    #     Mean_CN = round(rowMeans(copy_number_matrix, na.rm = TRUE), 3),
    #     Median_CN = round(apply(copy_number_matrix, 1, median, na.rm = TRUE), 3),
    #     Min_CN = round(apply(copy_number_matrix, 1, min, na.rm = TRUE), 3),
    #     Max_CN = round(apply(copy_number_matrix, 1, max, na.rm = TRUE), 3),
    #     SD_CN = round(apply(copy_number_matrix, 1, sd, na.rm = TRUE), 3)
    # )
    
    # NEW CODE - GISTIC2-Style:
    # GISTIC2 thresholds from parameter file
    gistic_amp_threshold <- 0.25  # From GISTIC2 parameter file
    gistic_del_threshold <- -0.25 # From GISTIC2 parameter file (negative for deletions)
    
    # Calculate GISTIC2-style summary statistics
    arm_stats <- data.frame(
        Arm = rownames(copy_number_matrix),
        Mean_Log2Ratio = round(rowMeans(copy_number_matrix, na.rm = TRUE), 4),
        Median_Log2Ratio = round(apply(copy_number_matrix, 1, median, na.rm = TRUE), 4),
        Min_Log2Ratio = round(apply(copy_number_matrix, 1, min, na.rm = TRUE), 4),
        Max_Log2Ratio = round(apply(copy_number_matrix, 1, max, na.rm = TRUE), 4),
        SD_Log2Ratio = round(apply(copy_number_matrix, 1, sd, na.rm = TRUE), 4),
        Pct_Amplified = round(rowSums(copy_number_matrix >= gistic_amp_threshold, na.rm = TRUE) / ncol(copy_number_matrix) * 100, 1),
        Pct_Deleted = round(rowSums(copy_number_matrix <= gistic_del_threshold, na.rm = TRUE) / ncol(copy_number_matrix) * 100, 1),
        Pct_Altered = round(rowSums(abs(copy_number_matrix) >= gistic_amp_threshold, na.rm = TRUE) / ncol(copy_number_matrix) * 100, 1)
    )
    
    # Save matrices if output_dir provided
    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        # =============================================================================
        # GISTIC2-Style Output Files Fix
        # 
        # ISSUE: Original file names and content structure don't reflect GISTIC2 convention
        #
        # FIX: Update file names and data structure:
        # - arm_copynumber_matrix.csv now contains log2 ratios (GISTIC2 style)
        # - Add discrete calls matrix (-1, 0, +1) for compatibility
        # - Update column names to reflect log2 ratio nature
        # =============================================================================
        
        # Save log2 ratio matrix (this is now the "copy number" matrix in GISTIC2 style)
        write.csv(arm_matrix, file.path(output_dir, "arm_log2ratio_matrix.csv"))
        write.csv(copy_number_matrix, file.path(output_dir, "arm_copynumber_matrix.csv"))
        write.csv(arm_stats, file.path(output_dir, "arm_copynumber_summary.csv"), row.names = FALSE)
        
        # NEW: Create discrete calls matrix (-1, 0, +1) following GISTIC2 convention
        calls_matrix <- matrix(0, nrow = nrow(copy_number_matrix), ncol = ncol(copy_number_matrix))
        rownames(calls_matrix) <- rownames(copy_number_matrix)
        colnames(calls_matrix) <- colnames(copy_number_matrix)
        calls_matrix[copy_number_matrix >= gistic_amp_threshold] <- 1   # Amplification
        calls_matrix[copy_number_matrix <= gistic_del_threshold] <- -1  # Deletion
        write.csv(calls_matrix, file.path(output_dir, "arm_calls_matrix_gistic_style.csv"))
        
        # Long format - updated to reflect GISTIC2 style
        long_format <- as.data.frame(copy_number_matrix) %>%
            mutate(Arm = rownames(.)) %>%
            pivot_longer(-Arm, names_to = "Sample", values_to = "Log2Ratio") %>%  # Changed from Copy_Number to Log2Ratio
            arrange(Arm, Sample)
        
        # Add discrete calls to long format
        calls_long <- as.data.frame(calls_matrix) %>%
            mutate(Arm = rownames(.)) %>%
            pivot_longer(-Arm, names_to = "Sample", values_to = "Call") %>%
            arrange(Arm, Sample)
        
        # Combine continuous and discrete data
        long_format <- long_format %>%
            left_join(calls_long, by = c("Arm", "Sample")) %>%
            mutate(
                Call_Label = case_when(
                    Call == 1 ~ "Amplification",
                    Call == -1 ~ "Deletion", 
                    Call == 0 ~ "Normal"
                )
            )
        
        write.csv(long_format, file.path(output_dir, "arm_copynumber_long.csv"), row.names = FALSE)
    } else {
        # Long format for return - GISTIC2 style
        calls_matrix <- matrix(0, nrow = nrow(copy_number_matrix), ncol = ncol(copy_number_matrix))
        rownames(calls_matrix) <- rownames(copy_number_matrix)
        colnames(calls_matrix) <- colnames(copy_number_matrix)
        calls_matrix[copy_number_matrix >= gistic_amp_threshold] <- 1
        calls_matrix[copy_number_matrix <= gistic_del_threshold] <- -1
        
        long_format <- as.data.frame(copy_number_matrix) %>%
            mutate(Arm = rownames(.)) %>%
            pivot_longer(-Arm, names_to = "Sample", values_to = "Log2Ratio") %>%  # Changed from Copy_Number
            arrange(Arm, Sample)
        
        calls_long <- as.data.frame(calls_matrix) %>%
            mutate(Arm = rownames(.)) %>%
            pivot_longer(-Arm, names_to = "Sample", values_to = "Call") %>%
            arrange(Arm, Sample)
        
        long_format <- long_format %>%
            left_join(calls_long, by = c("Arm", "Sample")) %>%
            mutate(
                Call_Label = case_when(
                    Call == 1 ~ "Amplification",
                    Call == -1 ~ "Deletion",
                    Call == 0 ~ "Normal"
                )
            )
    }
    
    # =============================================================================
    # GISTIC2-Style Return Object Fix
    # 
    # FIX: Update return object to include GISTIC2-style components:
    # - cn_matrix now contains log2 ratios (GISTIC2 style) instead of absolute copy numbers
    # - Add calls_matrix for discrete classifications
    # - Add GISTIC2 thresholds for reference
    # =============================================================================
    
    return(list(
        log2_matrix = arm_matrix,                    # Raw log2 ratios
        cn_matrix = copy_number_matrix,             # GISTIC2-style copy numbers (same as log2_matrix)
        calls_matrix = calls_matrix,                # NEW: Discrete calls (-1, 0, +1)
        summary_stats = arm_stats,                  # Updated GISTIC2-style statistics
        long_format = long_format,                  # Updated long format with calls
        gistic_thresholds = list(                   # NEW: GISTIC2 parameters used
            amp_threshold = gistic_amp_threshold,
            del_threshold = gistic_del_threshold
        )
    ))
}

#' Generate Analysis Report
#'
#' @param cancer_type Character, cancer type name
#' @param results List, analysis results
#' @param output_dir Character, output directory
#' @export
generate_analysis_report <- function(cancer_type, results, output_dir) {
    
    report_file <- file.path(output_dir, "BAGEL_V2_ANALYSIS_REPORT.md")
    
    # Generate report content
    report_content <- c(
        paste("# BAGEL v2.0 Analysis Report -", stringr::str_to_title(gsub("_", " ", cancer_type))),
        "",
        "## Analysis Summary",
        paste("- **Analysis Date**:", Sys.Date()),
        paste("- **Cancer Type**:", cancer_type),
        paste("- **BAGEL Version**: 2.0.0"),
        paste("- **Samples Analyzed**:", length(unique(results$segments$Sample))),
        paste("- **Chromosome Arms**:", length(unique(results$arm_summaries$Arm))),
        "",
        "## Key Findings",
        ""
    )
    
    # Add significant arms if available
    if (!is.null(results$significant_arms) && nrow(results$significant_arms) > 0) {
        report_content <- c(
            report_content,
            paste("### Significant Chromosome Arms (q < 0.25)"),
            ""
        )
        
        for (i in 1:min(10, nrow(results$significant_arms))) {
            arm_info <- results$significant_arms[i, ]
            report_content <- c(
                report_content,
                paste("- **", arm_info$Arm, "**: z-score =", round(arm_info$mean_z_score, 2),
                     ", samples =", arm_info$n_samples)
            )
        }
        report_content <- c(report_content, "")
    }
    
    # Add alteration frequencies
    if (!is.null(results$matrices)) {
      alteration_summary <- matrices$long_format %>%
        count(Arm, Call_Label) %>% 
        pivot_wider(names_from = Call_Label, values_from = n, values_fill = 0) %>% 
        dplyr::mutate(Total_Samples = Normal + Amplification + Deletion,
                      Amp_Freq = round(Amplification / Total_Samples * 100, 1),
                      Del_Freq = round(Deletion / Total_Samples * 100, 1),
                      Alt_Freq = Amp_Freq + Del_Freq) %>%
        arrange(desc(Alt_Freq))
        
        report_content <- c(
            report_content,
            "### Top Altered Arms",
            ""
        )
        
        for (i in 1:min(10, nrow(alteration_summary))) {
            arm_info <- alteration_summary[i, ]
            report_content <- c(
                report_content,
                paste("- **", arm_info$Arm, "**: ",
                     arm_info$Amp_Freq, "% amplified, ",
                     arm_info$Del_Freq, "% deleted")
            )
        }
    }
    
    # Add file listing
    report_content <- c(
        report_content,
        "",
        "## Output Files",
        "",
        "### Matrices",
        "- `arm_copynumber_matrix.csv` - Copy number matrix (arms × samples)",
        "- `arm_log2ratio_matrix.csv` - Log2 ratio matrix (as used in analysis)",
        "- `arm_copynumber_long.csv` - Long format for statistical analysis",
        "- `arm_copynumber_summary.csv` - Summary statistics per arm",
        "",
        "### Analysis Results", 
        "- `BAGEL_results_*/arm_level_summaries.txt` - Detailed arm summaries",
        "- `BAGEL_results_*/arm_definitions.txt` - Chromosome arm definitions used in analysis",
        "- `BAGEL_results_*/gistic_results.txt` - GISTIC statistical results",
        "- `BAGEL_results_*/significant_arms.txt` - Significant chromosome arms",
        "- `bagel_v2_complete_results.RData` - Complete R results object",
        "",
        "---",
        paste("**Analysis completed**:", Sys.time()),
        "**BAGEL Version**: 2.0.0"
    )
    
    # Write report
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(report_content, report_file)
    bagel_log(paste("Generated analysis report:", report_file))
}

#' Run Complete BAGEL v2.0 Analysis
#'
#' @param cancer_type Character, cancer type name
#' @param data_dir Character, path to data directory
#' @param output_dir Character, output directory
#' @param amp_threshold Numeric, amplification threshold (default log2(2.5/2))
#' @param del_threshold Numeric, deletion threshold (default log2(1.5/2))
#' @param stringent_threshold Numeric, stringent threshold (default 0.9)
#' @param use_gistic Logical, whether to use GISTIC analysis (default TRUE)
#' @param save_results Logical, whether to save results (default TRUE)
#' @return List with analysis results or TRUE/FALSE if used as workflow
#' @export
# =============================================================================
# GISTIC2-Style Workflow Parameters Fix
# 
# ISSUE: Workflow function used inconsistent thresholds with main function
#
# FIX: Update to use same GISTIC2 parameters as calculateCopyNumber_fixed
# =============================================================================

run_bagel_workflow <- function(cancer_type, 
                               data_dir,
                               output_dir = NULL,
                               # OLD CODE - COMMENTED OUT:
                               # amp_threshold = log2(2.5/2),
                               # del_threshold = log2(1.5/2),
                               
                               # NEW CODE - GISTIC2 Standard:
                               amp_threshold = 0.25,     # From GISTIC2 parameter file
                               del_threshold = -0.25,    # From GISTIC2 parameter file
                               stringent_threshold = 0.9,
                               use_gistic = TRUE,
                               save_results = TRUE) {
    
    bagel_log(paste("=== BAGEL v2.0 Analysis:", cancer_type, "==="))
    
    # Setup output directory
    if (is.null(output_dir)) {
        output_dir <- file.path("BAGEL_results", cancer_type, "bagel_v2_analysis")
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Setup logging
    log_file <- file.path(output_dir, "analysis.log")
    setup_bagel_logging(log_level = "INFO", log_file = log_file)
    
    tryCatch({
        
        # Step 1: Load breakpoint data
        bagel_log("Step 1: Loading breakpoint data...")
        breakpoint_data <- load_breakpoint_data()
        
        # Step 2: Load segmentation data  
        bagel_log("Step 2: Loading segmentation data...")
        segments <- load_segmentation_data(cancer_type, data_dir)
        
        if (nrow(segments) == 0) {
            stop("No valid segments found for ", cancer_type)
        }
        
        # Step 3: Get arm definitions
        bagel_log("Step 3: Getting arm definitions...")
        arm_definitions <- get_arm_definitions(cancer_type, breakpoint_data)
        
        bagel_log(paste("Arms available:", nrow(arm_definitions)))
        
        # Step 4: Run complete BAGEL analysis
        bagel_log("Step 4: Running BAGEL v2.0 analysis...")
        
        bagel_results <- calculateCopyNumber_fixed(
            segments = segments,
            breakpoints = arm_definitions,
            amp_threshold = amp_threshold,
            del_threshold = del_threshold,
            stringent_threshold = stringent_threshold,
            output_dir = output_dir,
            cancer_type = cancer_type,
            use_gistic = use_gistic,
            save_results = save_results
        )
        
        # Step 5: Create chromosome arm matrices
        bagel_log("Step 5: Creating chromosome arm matrices...")
        
        matrices <- create_arm_matrix(bagel_results$arm_summaries, output_dir)
        bagel_results$matrices <- matrices
        
        # Step 6: Save complete results
        if (save_results) {
            bagel_log("Step 6: Saving complete results...")
            save(bagel_results, file = file.path(output_dir, "bagel_v2_complete_results.RData"))
        }
        
        # Step 7: Generate analysis report
        bagel_log("Step 7: Generating analysis report...")
        generate_analysis_report(cancer_type, bagel_results, output_dir)
        
        bagel_log("=== ANALYSIS COMPLETED SUCCESSFULLY ===")
        bagel_log(paste("Cancer Type:", cancer_type))
        bagel_log(paste("Samples:", length(unique(bagel_results$segments$Sample))))
        bagel_log(paste("Arms:", length(unique(bagel_results$arm_summaries$Arm))))
        bagel_log(paste("Output directory:", output_dir))
        
        return(bagel_results)
        
    }, error = function(e) {
        
        error_msg <- paste("ERROR in", cancer_type, "analysis:", e$message)
        bagel_log(error_msg, "ERROR")
        
        # Create error summary
        error_summary <- c(
            paste("BAGEL v2.0 Analysis FAILED -", cancer_type),
            paste("Error occurred:", Sys.time()),
            paste("Error message:", e$message),
            "",
            "Check analysis.log for detailed error information."
        )
        
        writeLines(error_summary, file.path(output_dir, "ANALYSIS_ERROR.txt"))
        
        stop(e)
    })
}