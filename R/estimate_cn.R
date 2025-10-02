#' Calculate Copy Number Analysis with Functional Aneuploidy Events
#'
#' @description
#' Performs comprehensive copy number analysis using biologically meaningful functional
#' aneuploidy regions defined by \code{\link{define_arm}}. Each unique functional event
#' (identified by peak_id) is analyzed independently to generate arm-level copy number
#' summaries, GISTIC2-style statistical significance, and comprehensive matrix outputs.
#'
#' This function implements a complete pipeline from segment-to-event assignment through
#' statistical testing and matrix generation, following BAGEL coding standards with
#' maximum 2-level function nesting.
#'
#' @section Workflow Overview:
#' The analysis follows a six-stage pipeline:
#' \enumerate{
#'   \item \strong{Input Validation}: Validates segments and arm_definitions structure,
#'         ensuring peak_id uniqueness and required columns
#'   \item \strong{Segment-to-Event Assignment}: Assigns each copy number segment to
#'         the functional aneuploidy event with maximum overlap
#'   \item \strong{Event-Level Summarization}: Calculates weighted mean log2 ratios
#'         for each peak_id in each sample
#'   \item \strong{GISTIC2 Statistical Testing}: Performs permutation-based significance
#'         testing with FDR correction (Benjamini-Hochberg)
#'   \item \strong{Matrix Generation}: Creates peak_id × sample matrices in multiple
#'         formats (log2, copy number, discrete calls)
#'   \item \strong{Results Export}: Saves comprehensive outputs including statistics,
#'         summaries, and visualizations
#' }
#'
#' @section Input Data Structure:
#' \strong{Segments} (copy number segmentation data):
#' \describe{
#'   \item{Sample}{Character. Sample identifier}
#'   \item{Chromosome}{Numeric. Chromosome number (1-22)}
#'   \item{Start}{Numeric. Segment start position (bp)}
#'   \item{End}{Numeric. Segment end position (bp)}
#'   \item{Log2Ratios or Segment_Mean}{Numeric. Copy number in log2 scale}
#'   \item{NumMarkers}{Numeric. Number of probes/markers (optional, defaults to 1)}
#' }
#'
#' \strong{Arm Definitions} (from \code{define_arm()}):
#' \describe{
#'   \item{peak_id}{Character. \strong{UNIQUE} identifier for each functional event
#'         (e.g., "10q_amp_tel_neg", "11p_del_tel_neg")}
#'   \item{chr}{Character or Numeric. Chromosome number (1-22)}
#'   \item{arm_start}{Numeric. Functional region start position (bp)}
#'   \item{arm_end}{Numeric. Functional region end position (bp)}
#'   \item{direction}{Character. Expected alteration direction ("amp" or "del")}
#'   \item{telcent}{Character. Boundary type ("tel" for telomere, "cent" for centromere)}
#'   \item{type_of_selection}{Character. Selection classification
#'         (e.g., "onco-like", "TS-like", "toxic-like", "essential-like")}
#'   \item{functional_length}{Numeric. Length of functional region (bp, optional)}
#'   \item{percentage_length}{Numeric. Percentage of arm covered (optional)}
#'   \item{combined_sig, n_events}{Numeric. BISCUT statistical metrics (optional)}
#' }
#'
#' @section Algorithm Details:
#' \subsection{Segment Assignment}{
#'   Each segment is assigned to the functional event with \strong{maximum overlap}:
#'   \itemize{
#'     \item Overlap = min(segment_end, arm_end) - max(segment_start, arm_start) + 1
#'     \item Segments overlapping multiple events are assigned to the event with
#'           greatest overlap length
#'     \item Only segments on autosomes (chr 1-22) are processed
#'   }
#' }
#'
#' \subsection{Event-Level Summarization}{
#'   For each Peak_ID in each Sample:
#'   \itemize{
#'     \item \strong{mean_log2ratio}: Weighted mean by segment length
#'     \item \strong{median_log2ratio}: Unweighted median
#'     \item \strong{n_segments}: Number of segments assigned
#'     \item \strong{total_length}: Sum of segment lengths
#'   }
#' }
#'
#' \subsection{GISTIC2 Statistical Analysis}{
#'   Implements full GISTIC2 methodology:
#'   \enumerate{
#'     \item \strong{Arm Medians}: Calculate median log2 ratio per event across samples
#'     \item \strong{Threshold Filtering}: Separate amp (≥0.25) and del (≤-0.25) events
#'     \item \strong{Frequency Calculation}: Proportion of samples with significant alterations
#'     \item \strong{Permutation Testing}: Generate null distribution by randomly
#'           permuting samples (default: 1000 permutations)
#'     \item \strong{P-value Calculation}: Empirical two-sided test comparing observed
#'           G-scores to null distribution
#'     \item \strong{FDR Correction}: Benjamini-Hochberg procedure for multiple testing
#'     \item \strong{Significance Classification}: q<0.01 (highly significant),
#'           q<0.05 (significant), q<0.10 (suggestive), q<0.25 (nominally significant)
#'   }
#' }
#'
#' \subsection{Matrix Formats}{
#'   Three matrix formats are generated (all with Peak_ID as row names):
#'   \itemize{
#'     \item \strong{log2_matrix}: Raw log2 ratios (continuous values)
#'     \item \strong{cn_matrix}: Copy number values (same as log2, GISTIC2 convention)
#'     \item \strong{calls_matrix}: Discrete calls (-1=Deletion, 0=Neutral, +1=Amplification)
#'   }
#' }
#'
#' @section Key Features:
#' \itemize{
#'   \item \strong{Unique Event Identification}: Uses peak_id to distinguish functional
#'         events (e.g., "1q_amp_tel_pos" vs "1q_del_cent_pos" are separate rows)
#'   \item \strong{GISTIC2 Compliance}: Standard thresholds (amp=0.25, del=-0.25),
#'         proper statistical testing, and FDR correction
#'   \item \strong{Comprehensive Output}: Statistical summaries, matrix formats,
#'         long-format data, and metadata
#'   \item \strong{Flexible Thresholds}: Customizable amp/del thresholds and
#'         stringency parameters
#'   \item \strong{Memory Efficient}: Chunked processing for large datasets
#'   \item \strong{Selection Type Integration}: Preserves evolutionary selection
#'         classifications from BISCUT analysis
#' }
#'
#' @param segments Data frame containing segmentation data. \strong{Required columns}:
#'   \itemize{
#'     \item Sample: Sample identifier (character)
#'     \item Chromosome: Chromosome number 1-22 (numeric)
#'     \item Start: Segment start position in bp (numeric)
#'     \item End: Segment end position in bp (numeric)
#'     \item Log2Ratios or Segment_Mean: Copy number in log2 scale (numeric)
#'     \item NumMarkers: Number of probes (numeric, optional - defaults to 1)
#'   }
#' @param breakpoints Data frame with functional arm definitions from \code{\link{define_arm}}.
#'   \strong{Required columns}: peak_id, chr, arm_start, arm_end, direction, telcent,
#'   type_of_selection. See "Input Data Structure" section for details.
#' @param amp_threshold Numeric. Log2 ratio threshold for amplifications (default: 0.25).
#'   GISTIC2 standard: 0.25 corresponds to ~1.19x copy number ratio.
#' @param del_threshold Numeric. Log2 ratio threshold for deletions (default: -0.25).
#'   GISTIC2 standard: -0.25 corresponds to ~0.84x copy number ratio.
#' @param stringent_threshold Numeric. Multiplier for high-confidence calls (default: 0.9).
#'   Strong alterations defined as |log2ratio| ≥ threshold × stringent_threshold.
#' @param output_dir Character. Directory path for saving results (default: ".").
#'   Creates subdirectories for different output types.
#' @param cancer_type Character. Cancer type identifier for logging and metadata (default: "Unknown").
#'   Used in output file naming and analysis reports.
#' @param use_gistic Logical. Perform GISTIC2-style statistical analysis (default: TRUE).
#'   When FALSE, skips permutation testing and significance calculations.
#' @param chunk_size Numeric. Number of samples per processing chunk (default: 100).
#'   For memory-efficient processing of large datasets.
#' @param save_results Logical. Save results to files (default: TRUE).
#'   When FALSE, returns results in memory only without file output.
#' @param create_matrices Logical. Generate matrix outputs (default: TRUE).
#'   Creates log2_matrix, cn_matrix, calls_matrix, and summary statistics.
#' @param matrix_formats Character. Matrix format selection (default: "all").
#'   Options: "all" (all formats), "basic" (log2 and calls), "statistical" (with stats).
#'
#' @return Named list with comprehensive analysis results:
#' \describe{
#'   \item{segments}{Data frame. Segments with assigned_peak_id indicating which
#'         functional event each segment belongs to}
#'   \item{arm_summaries}{Data frame. Event-level summaries with columns:
#'         \itemize{
#'           \item Peak_ID: Unique functional event identifier
#'           \item Sample: Sample identifier
#'           \item mean_log2ratio: Weighted mean copy number
#'           \item median_log2ratio: Median copy number
#'           \item n_segments: Number of segments assigned to this event
#'           \item total_length: Total genomic length covered (bp)
#'           \item alteration_call: Classification (Amplification/Deletion/Neutral)
#'           \item alteration_strength: Strength (Strong/Moderate/Weak)
#'           \item matches_expected_direction: Logical, TRUE if observed matches expected
#'         }}
#'   \item{arm_definitions}{Data frame. Validated arm definitions used in analysis}
#'   \item{gistic_results}{Data frame. GISTIC2 statistical results (when use_gistic=TRUE):
#'         \itemize{
#'           \item Peak_ID: Functional event identifier
#'           \item median_alteration: Median log2 ratio across samples
#'           \item frequency_amplified: Percentage of samples amplified
#'           \item frequency_deleted: Percentage of samples deleted
#'           \item g_score: GISTIC G-score (mean alteration)
#'           \item p_value: Empirical p-value from permutation testing
#'           \item q_value: FDR-corrected q-value (Benjamini-Hochberg)
#'           \item is_significant_p05, is_significant_q25, is_significant_q10: Logical flags
#'         }}
#'   \item{significant_arms}{Data frame. Subset of gistic_results with p<0.05 or q<0.25,
#'         includes significance_level and effect_size classifications}
#'   \item{matrices}{List (when create_matrices=TRUE):
#'         \itemize{
#'           \item log2_matrix: Matrix (Peak_ID × Sample) with log2 ratios
#'           \item cn_matrix: Matrix (Peak_ID × Sample) with copy numbers (GISTIC2 format)
#'           \item calls_matrix: Matrix (Peak_ID × Sample) with discrete calls (-1/0/+1)
#'         }}
#'   \item{matrix_summary}{Data frame (when create_matrices=TRUE). Per-event statistics:
#'         Mean/Median/SD of log2 ratios, counts and percentages of Amp/Del/Neutral calls}
#'   \item{matrix_long}{Data frame (when create_matrices=TRUE). Long-format data with
#'         columns: Peak_ID, Sample, Log2Ratio, Call, Call_Label, Pct_Altered}
#'   \item{parameters}{List. Analysis parameters and settings used}
#'   \item{metadata}{List. Summary information:
#'         \itemize{
#'           \item n_samples: Number of samples analyzed
#'           \item n_arms: Number of functional events (unique peak_ids)
#'           \item has_enhanced_arms: Logical, TRUE for updated arm_definitions structure
#'           \item n_significant_arms: Number of significantly altered events
#'           \item matrices_created: Logical, whether matrices were generated
#'           \item matrix_dimensions: Character, dimensions of output matrices
#'         }}
#' }
#'
#' @section Output Files:
#' When save_results=TRUE, the following files are created in output_dir:
#' \itemize{
#'   \item \strong{enhanced_arm_summaries.tsv}: Complete event-level summaries
#'   \item \strong{gistic_statistics.tsv}: GISTIC2 statistical results
#'   \item \strong{significant_arms.tsv}: Significantly altered events only
#'   \item \strong{analysis_parameters.tsv}: Parameters used
#'   \item \strong{arm_log2ratio_matrix.csv}: Log2 ratio matrix
#'   \item \strong{arm_copynumber_matrix.csv}: Copy number matrix (GISTIC2 format)
#'   \item \strong{arm_calls_matrix_gistic_style.csv}: Discrete calls matrix
#'   \item \strong{arm_copynumber_summary.csv}: Per-event frequency statistics
#'   \item \strong{arm_copynumber_long.csv}: Long-format data for visualization
#'   \item \strong{final_results.RData}: Complete R object with all results
#' }
#'
#' @examples
#' \dontrun{
#' # Standard workflow with define_arm() output
#' arm_defs <- define_arm(
#'   custom_biscut_file = "/path/to/all_BISCUT_results.txt",
#'   cancer_type = "BRCA"
#' )
#'
#' results <- calculate_copynumber(
#'   segments = segments_df,
#'   breakpoints = arm_defs,
#'   cancer_type = "BRCA",
#'   output_dir = "./analysis/brca"
#' )
#'
#' # Access results
#' head(results$arm_summaries)  # Event-level summaries
#' results$matrices$log2_matrix[1:5, 1:5]  # First 5 events × 5 samples
#' significant_events <- results$significant_arms  # q < 0.25
#'
#' # High-stringency analysis
#' results_strict <- calculate_copynumber(
#'   segments = segments_df,
#'   breakpoints = arm_defs,
#'   amp_threshold = 0.3,
#'   del_threshold = -0.3,
#'   stringent_threshold = 0.95
#' )
#'
#' # Quick analysis without saving files
#' results_quick <- calculate_copynumber(
#'   segments = segments_df,
#'   breakpoints = arm_defs,
#'   save_results = FALSE,
#'   create_matrices = FALSE
#' )
#' }
#'
#' @seealso \code{\link{define_arm}}, \code{\link{create_matrix}}
#'
#' @import dplyr
#' @importFrom readr write_tsv
#' @importFrom stats median sd
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom rlang sym
#' @export
calculate_copynumber <- function(segments,
                                 breakpoints,
                                 amp_threshold = 0.25,      # GISTIC2 standard
                                 del_threshold = -0.25,     # GISTIC2 standard
                                 stringent_threshold = 0.9,
                                 output_dir = ".",
                                 cancer_type = "Unknown",
                                 use_gistic = TRUE,
                                 chunk_size = 100,
                                 save_results = TRUE,
                                 create_matrices = TRUE,    # NEW: Matrix generation
                                 matrix_formats = "all") {   # NEW: Matrix format control
  
  cat("Cancer Type:", cancer_type, "\n")
  cat("Samples:", length(unique(segments$Sample)), "\n")
  
  # ============================================================================
  # LEVEL 1 NESTED FUNCTIONS: Core analysis components
  # ============================================================================
  
  # Level 1: Validate and prepare input data
  prepare_input_data <- function(segs, arms, params) {
    
    # Level 2: Comprehensive segment validation
    validate_segments_enhanced <- function(segment_data) {
      required_cols <- c("Sample", "Chromosome", "Start", "End")
      
      # Check for Log2Ratios or Segment_Mean
      if (!"Log2Ratios" %in% names(segment_data)) {
        if ("Segment_Mean" %in% names(segment_data)) {
          segment_data$Log2Ratios <- segment_data$Segment_Mean
        } else {
          stop("Missing Log2Ratios or Segment_Mean column in segments")
        }
      }
      
      # Add NumMarkers if missing (for compatibility)
      if (!"NumMarkers" %in% names(segment_data)) {
        segment_data$NumMarkers <- 1  # Default value
      }
      
      # Basic validation
      missing_cols <- setdiff(required_cols, names(segment_data))
      if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
      }
      
      # Clean and filter data
      clean_data <- segment_data %>%
        dplyr::filter(
          !is.na(Start), !is.na(End), !is.na(Log2Ratios),
          Start < End,
          Chromosome %in% 1:22
        ) %>%
        dplyr::mutate(
          Chromosome = as.numeric(Chromosome),
          Segment_Length = End - Start + 1
        ) %>%
        dplyr::arrange(Sample, Chromosome, Start)
      
      cat("Segments validation: ", nrow(segment_data), "->", nrow(clean_data), "segments\n")
      return(clean_data)
    }
    
    # Level 2: Enhanced arm definitions validation
    # =================================================================
    # POLLY HUNG UPDATE - 2025-10-02
    # CHANGED: Updated to work with new arm_definitions structure from define_arm()
    # NEW STRUCTURE: Each row = unique functional aneuploidy event with peak_id
    # KEY COLUMNS: peak_id, chr, arm_start, arm_end, direction, telcent, type_of_selection
    # =================================================================
    validate_arms_enhanced <- function(arm_data) {
      # NEW: Required columns for updated arm_definitions format
      required_cols <- c("peak_id", "chr", "arm_start", "arm_end",
                         "direction", "telcent", "type_of_selection")

      # Check for required columns
      missing_cols <- setdiff(required_cols, names(arm_data))
      if (length(missing_cols) > 0) {
        stop("Missing required columns in arm_definitions: ",
             paste(missing_cols, collapse = ", "),
             "\nExpected columns from define_arm(): ",
             paste(required_cols, collapse = ", "))
      }

      # NEW: Check for peak_id uniqueness
      if (any(duplicated(arm_data$peak_id))) {
        duplicated_ids <- arm_data$peak_id[duplicated(arm_data$peak_id)]
        stop("Duplicate peak_id found in arm_definitions: ",
             paste(unique(duplicated_ids), collapse = ", "))
      }

      cat("Using updated arm definitions with functional regions\n")

      # NEW: Validate updated format
      validated_arms <- arm_data %>%
        dplyr::mutate(
          # Ensure chr is numeric (remove "chr" prefix if present)
          chr_num = as.numeric(gsub("chr", "", as.character(chr)))
        ) %>%
        dplyr::filter(
          chr_num %in% 1:22,
          !is.na(arm_start), !is.na(arm_end),
          arm_start < arm_end,
          direction %in% c("amp", "del"),
          telcent %in% c("tel", "cent"),
          !is.na(peak_id)
        ) %>%
        dplyr::mutate(is_enhanced = TRUE) %>%
        dplyr::arrange(chr_num, arm_start)

      cat("Arm definitions validated:", nrow(validated_arms), "functional events\n")
      cat("Sample peak_ids:", paste(head(validated_arms$peak_id, 3), collapse=", "), "...\n")

      return(validated_arms)
    }

    # # OLD CODE - COMMENTED OUT FOR REFERENCE
    # validate_arms_enhanced <- function(arm_data) {
    #   required_base_cols <- c("arm", "chr", "arm_type", "arm_start", "arm_end", "direction")
    #   enhanced_cols <- c("telcent", "type_of_selection", "functional_length_mb",
    #                      "combined_sig", "n_events", "n_right", "n_left")
    #
    #   # Check if we have enhanced format
    #   has_enhanced <- all(enhanced_cols %in% names(arm_data))
    #
    #   if (has_enhanced) {
    #     cat("Using enhanced arm definitions with functional regions\n")
    #
    #     # Validate enhanced format
    #     validated_arms <- arm_data %>%
    #       dplyr::filter(
    #         chr %in% 1:22,
    #         arm_type %in% c("p", "q"),
    #         !is.na(arm_start), !is.na(arm_end),
    #         arm_start < arm_end,
    #         direction %in% c("amp", "del"),
    #         telcent %in% c("tel", "cent"),
    #         !is.na(functional_length_mb)
    #       ) %>%
    #       dplyr::mutate(
    #         chr_num = as.numeric(chr),
    #         is_enhanced = TRUE
    #       ) %>%
    #       dplyr::arrange(chr_num, arm_type, arm_start)
    #
    #   } else {
    #     cat("Using basic arm definitions\n")
    #
    #     # Handle basic format
    #     validated_arms <- arm_data %>%
    #       dplyr::mutate(
    #         chr_num = if ("chr_num" %in% names(.)) as.numeric(chr_num) else as.numeric(chr),
    #         direction = if ("direction" %in% names(.)) direction else "amp",
    #         is_enhanced = FALSE
    #       ) %>%
    #       dplyr::filter(
    #         chr_num %in% 1:22,
    #         arm_type %in% c("p", "q"),
    #         !is.na(arm_start), !is.na(arm_end),
    #         arm_start < arm_end
    #       ) %>%
    #       dplyr::arrange(chr_num, arm_type, arm_start)
    #   }
    #
    #   cat("Arm definitions validated:", nrow(validated_arms), "regions\n")
    #   return(validated_arms)
    # }
    
    # Main Level 1 logic
    validated_segments <- validate_segments_enhanced(segs)
    validated_arms <- validate_arms_enhanced(arms)
    
    # Parameter validation
    if (params$amp_threshold <= 0 || params$del_threshold >= 0) {
      warning("Non-standard thresholds: amp should be positive, del should be negative")
    }
    
    return(list(
      segments = validated_segments,
      arm_definitions = validated_arms,
      n_samples = length(unique(validated_segments$Sample)),
      has_enhanced_arms = all(validated_arms$is_enhanced)
    ))
  }
  
  # Level 1: Perform arm-level analysis
  analyze_arm_levels <- function(validated_data, analysis_params) {
    
    # Level 2: Annotate segments with arm assignments
    # =================================================================
    # POLLY HUNG UPDATE - 2025-10-02
    # CHANGED: Assign segments to peak_id instead of arm
    # WHY: Each peak_id represents a unique functional aneuploidy event
    # =================================================================
    annotate_segments_with_arms <- function(segs, arms) {
      annotated_segments <- segs %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          assigned_peak_id = {  # NEW: Changed from assigned_arm
            # Find overlapping functional regions for this segment
            chr_arms <- arms[arms$chr_num == Chromosome, ]
            overlaps <- which(
              chr_arms$arm_start <= End &
                chr_arms$arm_end >= Start
            )

            if (length(overlaps) > 0) {
              # Take the functional region with maximum overlap
              overlap_lengths <- pmin(chr_arms$arm_end[overlaps], End) -
                pmax(chr_arms$arm_start[overlaps], Start) + 1
              best_overlap <- which.max(overlap_lengths)
              chr_arms$peak_id[overlaps[best_overlap]]  # NEW: Changed from arm
            } else {
              NA_character_
            }
          }
        ) %>%
        dplyr::ungroup() %>%
        dplyr::filter(!is.na(assigned_peak_id))  # NEW: Changed from assigned_arm

      return(annotated_segments)
    }

    # # OLD CODE - COMMENTED OUT FOR REFERENCE
    # annotate_segments_with_arms <- function(segs, arms) {
    #   annotated_segments <- segs %>%
    #     dplyr::rowwise() %>%
    #     dplyr::mutate(
    #       assigned_arm = {
    #         # Find overlapping arms for this segment
    #         chr_arms <- arms[arms$chr_num == Chromosome, ]
    #         overlaps <- which(
    #           chr_arms$arm_start <= End &
    #             chr_arms$arm_end >= Start
    #         )
    #
    #         if (length(overlaps) > 0) {
    #           # Take the arm with maximum overlap
    #           overlap_lengths <- pmin(chr_arms$arm_end[overlaps], End) -
    #             pmax(chr_arms$arm_start[overlaps], Start) + 1
    #           best_overlap <- which.max(overlap_lengths)
    #           chr_arms$arm[overlaps[best_overlap]]
    #         } else {
    #           NA_character_
    #         }
    #       }
    #     ) %>%
    #     dplyr::ungroup() %>%
    #     dplyr::filter(!is.na(assigned_arm))
    #
    #   return(annotated_segments)
    # }
    
    # Level 2: Calculate arm-level summaries
    # =================================================================
    # POLLY HUNG UPDATE - 2025-10-02
    # CHANGED: Group by peak_id directly, no need to create functional_arm_id
    # WHY: peak_id already contains all necessary information
    # =================================================================
    calculate_arm_summaries <- function(annotated_segs, arms, params) {

      # NEW: Basic summaries grouped by peak_id
      basic_summaries <- annotated_segs %>%
        dplyr::group_by(Sample, assigned_peak_id) %>%  # NEW: Changed from assigned_arm
        dplyr::summarise(
          mean_log2ratio = weighted.mean(Log2Ratios, Segment_Length, na.rm = TRUE),
          median_log2ratio = median(Log2Ratios, na.rm = TRUE),
          n_segments = n(),
          total_length = sum(Segment_Length),
          .groups = "drop"
        ) %>%
        dplyr::rename(Peak_ID = assigned_peak_id) %>%  # NEW: Use Peak_ID
        # NEW: Add functional event information by joining on peak_id
        dplyr::left_join(
          arms %>% dplyr::select(peak_id, chr_num, arm, direction, telcent,
                                 type_of_selection,
                                 dplyr::any_of(c("functional_length", "percentage_length",
                                                 "combined_sig", "n_events"))),
          by = c("Peak_ID" = "peak_id")
        )

      # Add alteration classifications
      enhanced_summaries <- basic_summaries %>%
        dplyr::mutate(
          alteration_call = case_when(
            mean_log2ratio >= params$amp_threshold ~ "Amplification",
            mean_log2ratio <= params$del_threshold ~ "Deletion",
            TRUE ~ "Neutral"
          ),
          alteration_strength = case_when(
            abs(mean_log2ratio) >= max(abs(params$amp_threshold), abs(params$del_threshold)) * params$stringent_threshold ~ "Strong",
            abs(mean_log2ratio) >= max(abs(params$amp_threshold), abs(params$del_threshold)) ~ "Moderate",
            TRUE ~ "Weak"
          ),
          matches_expected_direction = case_when(
            is.na(direction) ~ TRUE,  # No expectation
            direction == "amp" & mean_log2ratio >= params$amp_threshold ~ TRUE,
            direction == "del" & mean_log2ratio <= params$del_threshold ~ TRUE,
            TRUE ~ FALSE
          )
        )

      return(enhanced_summaries)
    }

    # # OLD CODE - COMMENTED OUT FOR REFERENCE
    # calculate_arm_summaries <- function(annotated_segs, arms, params) {
    #
    #   # Basic arm summaries
    #   basic_summaries <- annotated_segs %>%
    #     dplyr::group_by(Sample, assigned_arm) %>%
    #     dplyr::summarise(
    #       mean_log2ratio = weighted.mean(Log2Ratios, Segment_Length, na.rm = TRUE),
    #       median_log2ratio = median(Log2Ratios, na.rm = TRUE),
    #       n_segments = n(),
    #       total_length = sum(Segment_Length),
    #       .groups = "drop"
    #     ) %>%
    #     dplyr::rename(Arm = assigned_arm) %>%
    #     # Add arm information
    #     dplyr::left_join(
    #       arms %>% dplyr::select(arm, chr_num, arm_type, direction,
    #                              dplyr::any_of(c("telcent", "type_of_selection",
    #                                              "functional_length_mb", "combined_sig"))),
    #       by = c("Arm" = "arm")
    #     ) %>%
    #     # =================================================================
    #   # POLLY HUNG FIX - 2025-09-26
    #   # ISSUE: Basic arm summaries lack functional arm identifiers
    #   # SOLUTION: Create functional_arm_id for matrix row names
    #   # WHY: User requires functional arms like "1q_del_tel" not just "1q"
    #   # =================================================================
    #   dplyr::mutate(
    #     functional_arm_id = if ("telcent" %in% names(.)) {
    #       case_when(
    #         !is.na(direction) & !is.na(telcent) ~ paste0(Arm, "_", direction, "_", substr(telcent, 1, 3)),
    #         !is.na(direction) ~ paste0(Arm, "_", direction),
    #         TRUE ~ Arm
    #       )
    #     } else {
    #       case_when(
    #         !is.na(direction) ~ paste0(Arm, "_", direction),
    #         TRUE ~ Arm
    #       )
    #     }
    #   )
    #
    #   # Add alteration classifications
    #   # =================================================================
    #   # POLLY HUNG MODIFICATION - 2025-09-26
    #   # CHANGED: basic_summaries now includes functional_arm_id from join above
    #   # =================================================================
    #   enhanced_summaries <- basic_summaries %>%
    #     dplyr::mutate(
    #       alteration_call = case_when(
    #         mean_log2ratio >= params$amp_threshold ~ "Amplification",
    #         mean_log2ratio <= params$del_threshold ~ "Deletion",
    #         TRUE ~ "Neutral"
    #       ),
    #       alteration_strength = case_when(
    #         abs(mean_log2ratio) >= max(abs(params$amp_threshold), abs(params$del_threshold)) * params$stringent_threshold ~ "Strong",
    #         abs(mean_log2ratio) >= max(abs(params$amp_threshold), abs(params$del_threshold)) ~ "Moderate",
    #         TRUE ~ "Weak"
    #       ),
    #       matches_expected_direction = case_when(
    #         is.na(direction) ~ TRUE,  # No expectation
    #         direction == "amp" & mean_log2ratio >= params$amp_threshold ~ TRUE,
    #         direction == "del" & mean_log2ratio <= params$del_threshold ~ TRUE,
    #         TRUE ~ FALSE
    #       )
    #     )
    #
    #   return(enhanced_summaries)
    # }
    
    # Main Level 1 logic
    cat("Annotating segments with arm assignments...\n")
    annotated_segments <- annotate_segments_with_arms(
      validated_data$segments,
      validated_data$arm_definitions
    )
    
    cat("Calculating arm-level summaries...\n")
    arm_summaries <- calculate_arm_summaries(
      annotated_segments,
      validated_data$arm_definitions,
      analysis_params
    )
    
    return(list(
      segments = annotated_segments,
      arm_summaries = arm_summaries,
      n_annotated_segments = nrow(annotated_segments)
    ))
  }
  
  # Level 1: Perform statistical analysis
  perform_statistical_analysis <- function(analysis_results, params) {
    
    # Level 2: GISTIC2-style significance testing
    # =================================================================
    # POLLY HUNG MAJOR ENHANCEMENT - 2025-09-26
    # ISSUE: Oversimplified z-score calculation lacked GISTIC2 rigor
    # SOLUTION: Implement full GISTIC2 statistical methodology
    # WHY: Need proper permutation testing, FDR correction, G-scores
    # BASED ON: GISTIC2 source code analysis (arm_medians.m, gistic_broad_analysis.m)
    # =================================================================
    # =================================================================
    # POLLY HUNG UPDATE - 2025-10-02
    # CHANGED: Use Peak_ID column instead of functional_arm_id
    # =================================================================
    calculate_gistic_statistics <- function(arm_data, analysis_params) {

      # NEW: Handle duplicates by taking mean for each Peak_ID-Sample combination
      clean_arm_data <- arm_data %>%
        dplyr::select(Sample, Peak_ID, mean_log2ratio) %>%  # NEW: Changed from functional_arm_id
        dplyr::group_by(Sample, Peak_ID) %>%  # NEW: Changed from functional_arm_id
        dplyr::summarise(mean_log2ratio = mean(mean_log2ratio, na.rm = TRUE), .groups = "drop")

      # NEW: Create peak_id × sample matrix for statistical testing
      arm_matrix <- clean_arm_data %>%
        tidyr::pivot_wider(names_from = Sample, values_from = mean_log2ratio, values_fill = NA) %>%
        tibble::column_to_rownames("Peak_ID") %>%  # NEW: Changed from functional_arm_id
        as.matrix()
      
      # Ensure matrix is numeric
      if (!is.numeric(arm_matrix)) {
        arm_matrix <- apply(arm_matrix, 2, as.numeric)
      }
      
      n_arms <- nrow(arm_matrix)
      n_samples <- ncol(arm_matrix)
      
      cat("GISTIC2-style analysis on", n_arms, "functional arms ×", n_samples, "samples\n")
      
      # =================================================================
      # STEP 1: Calculate GISTIC2-style arm medians and frequencies
      # Based on gistic_broad_analysis.m lines 54-75
      # =================================================================
      
      # Calculate median alteration per arm (equivalent to GISTIC2 arm_medians)
      arm_medians <- apply(arm_matrix, 1, median, na.rm = TRUE)
      
      # Apply GISTIC2 significance thresholds (separate for amp/del)
      arm_medians_amp <- pmax(arm_medians, 0)  # Keep only positive values
      arm_medians_del <- pmax(-arm_medians, 0)  # Keep only negative values (as positive)
      
      # Threshold filtering (GISTIC2 approach)
      arm_medians_amp[arm_medians_amp < analysis_params$amp_threshold] <- 0
      arm_medians_del[arm_medians_del < abs(analysis_params$del_threshold)] <- 0
      
      # Exclude opposing alterations (GISTIC2 lines 66-67)
      arm_medians_amp[arm_medians_del > 0] <- NA
      arm_medians_del[arm_medians_amp > 0] <- NA
      
      # Calculate GISTIC2-style frequencies (lines 70-71)
      # Frequency = altered samples / (total - opposite_alteration_samples)
      n_amp_samples <- rowSums(arm_matrix >= analysis_params$amp_threshold, na.rm = TRUE)
      n_del_samples <- rowSums(arm_matrix <= analysis_params$del_threshold, na.rm = TRUE)
      
      freq_amp <- n_amp_samples / (n_samples - n_del_samples)
      freq_del <- n_del_samples / (n_samples - n_amp_samples)
      
      # Handle division by zero
      freq_amp[is.nan(freq_amp) | is.infinite(freq_amp)] <- 0
      freq_del[is.nan(freq_del) | is.infinite(freq_del)] <- 0
      
      # =================================================================
      # STEP 2: Permutation-based significance testing
      # Based on zigg_score_permutations.m
      # =================================================================
      
      n_permutations <- min(1000, n_samples * 10)  # Computational efficiency
      cat("Running", n_permutations, "permutations for significance testing...\n")
      
      # Calculate observed G-scores (mean across samples)
      observed_g_scores <- rowMeans(arm_matrix, na.rm = TRUE)
      
      # Generate null distribution through permutation
      null_g_scores <- matrix(NA, n_arms, n_permutations)
      
      for (perm in 1:n_permutations) {
        # Permute each sample independently (GISTIC2 approach)
        permuted_matrix <- arm_matrix
        for (s in 1:n_samples) {
          if (sum(!is.na(arm_matrix[, s])) > 1) {
            valid_indices <- which(!is.na(arm_matrix[, s]))
            permuted_matrix[valid_indices, s] <- sample(arm_matrix[valid_indices, s])
          }
        }
        null_g_scores[, perm] <- rowMeans(permuted_matrix, na.rm = TRUE)
      }
      
      # Calculate empirical p-values
      p_values <- numeric(n_arms)
      for (i in 1:n_arms) {
        if (!is.na(observed_g_scores[i])) {
          # Two-sided test
          null_abs <- abs(null_g_scores[i, ])
          obs_abs <- abs(observed_g_scores[i])
          p_values[i] <- (sum(null_abs >= obs_abs, na.rm = TRUE) + 1) / (n_permutations + 1)
        } else {
          p_values[i] <- 1.0
        }
      }
      
      # =================================================================
      # STEP 3: FDR correction (Benjamini-Hochberg)
      # Based on make_final_Qs.m approach
      # =================================================================
      
      # Calculate q-values using Benjamini-Hochberg FDR correction
      valid_p <- !is.na(p_values) & p_values > 0
      q_values <- rep(1.0, n_arms)
      
      if (sum(valid_p) > 0) {
        p_sorted_indices <- order(p_values[valid_p])
        p_sorted <- p_values[valid_p][p_sorted_indices]
        m <- length(p_sorted)
        
        # Benjamini-Hochberg correction
        q_sorted <- numeric(m)
        q_sorted[m] <- p_sorted[m]
        for (i in (m-1):1) {
          q_sorted[i] <- min(q_sorted[i+1], p_sorted[i] * m / i)
        }
        
        # Map back to original order
        q_values[valid_p][p_sorted_indices] <- q_sorted
      }
      
      # =================================================================
      # STEP 4: Compile comprehensive statistics
      # =================================================================
      # NEW: Changed column name from Functional_Arm to Peak_ID
      arm_stats <- data.frame(
        Peak_ID = rownames(arm_matrix),  # NEW: Changed from Functional_Arm
        
        # GISTIC2-style medians and frequencies
        median_alteration = arm_medians,
        median_amplification = ifelse(is.na(arm_medians_amp), 0, arm_medians_amp),
        median_deletion = ifelse(is.na(arm_medians_del), 0, arm_medians_del),
        frequency_amplified = round(freq_amp * 100, 1),
        frequency_deleted = round(freq_del * 100, 1),
        
        # Traditional statistics for compatibility
        mean_alteration = rowMeans(arm_matrix, na.rm = TRUE),
        sd_alteration = apply(arm_matrix, 1, sd, na.rm = TRUE),
        n_samples = rowSums(!is.na(arm_matrix)),
        n_amplified = n_amp_samples,
        n_deleted = n_del_samples,
        pct_amplified = round(n_amp_samples / n_samples * 100, 1),
        pct_deleted = round(n_del_samples / n_samples * 100, 1),
        pct_altered = round((n_amp_samples + n_del_samples) / n_samples * 100, 1),
        
        # GISTIC2-style significance scores
        g_score = observed_g_scores,
        p_value = p_values,
        q_value = q_values,
        log10_p = -log10(pmax(p_values, 1e-10)),  # Avoid log(0)
        log10_q = -log10(pmax(q_values, 1e-10)),
        
        # Significance classifications
        is_significant_p05 = p_values < 0.05,
        is_significant_q25 = q_values < 0.25,
        is_significant_q10 = q_values < 0.10,
        
        stringsAsFactors = FALSE
      ) %>%
        dplyr::arrange(p_value, desc(abs(g_score)))
      
      cat("GISTIC2 analysis complete:\n")
      cat("- Significant arms (p < 0.05):", sum(arm_stats$is_significant_p05), "\n")
      cat("- Significant arms (q < 0.25):", sum(arm_stats$is_significant_q25), "\n")
      cat("- Significant arms (q < 0.10):", sum(arm_stats$is_significant_q10), "\n")
      
      return(arm_stats)
    }
    
    # Level 2: Identify significant arms using GISTIC2 criteria
    # =================================================================
    # POLLY HUNG ENHANCEMENT - 2025-09-26
    # CHANGED: Use GISTIC2 p-values and q-values instead of z-scores
    # WHY: Proper statistical significance based on permutation testing
    # =================================================================
    identify_significant_arms <- function(gistic_stats, 
                                          p_threshold = 0.05, 
                                          q_threshold = 0.25) {
      
      # Apply multiple significance criteria (GISTIC2 standard)
      significant <- gistic_stats %>%
        dplyr::filter(
          p_value < p_threshold | q_value < q_threshold
        ) %>%
        dplyr::mutate(
          # GISTIC2-style significance classification
          significance_level = case_when(
            q_value < 0.01 ~ "q < 0.01 (highly significant)",
            q_value < 0.05 ~ "q < 0.05 (significant)", 
            q_value < 0.10 ~ "q < 0.10 (suggestive)",
            q_value < 0.25 ~ "q < 0.25 (nominally significant)",
            p_value < 0.05 ~ "p < 0.05 (uncorrected)",
            TRUE ~ "not significant"
          ),
          
          # Alteration direction classification
          alteration_type = case_when(
            median_amplification > 0 & median_deletion == 0 ~ "Amplification",
            median_deletion > 0 & median_amplification == 0 ~ "Deletion",
            median_amplification > 0 & median_deletion > 0 ~ "Mixed",
            TRUE ~ "Neutral"
          ),
          
          # Effect size categories
          effect_size = case_when(
            abs(g_score) >= 0.5 ~ "Large",
            abs(g_score) >= 0.25 ~ "Medium", 
            abs(g_score) >= 0.1 ~ "Small",
            TRUE ~ "Minimal"
          )
        ) %>%
        dplyr::arrange(q_value, p_value, desc(abs(g_score)))
      
      return(significant)
    }
    
    # Main Level 1 logic
    gistic_results <- NULL
    significant_arms <- NULL
    
    if (params$use_gistic && nrow(analysis_results$arm_summaries) > 0) {
      cat("Performing GISTIC-style statistical analysis...\n")
      
      gistic_results <- calculate_gistic_statistics(
        analysis_results$arm_summaries,
        params
      )
      
      significant_arms <- identify_significant_arms(gistic_results)
      
      cat("Statistical analysis completed:", nrow(significant_arms), "significant arms\n")
    }
    
    return(list(
      gistic_results = gistic_results,
      significant_arms = significant_arms
    ))
  }
  
  # Level 1: Create comprehensive matrix outputs (INTEGRATED FROM create_matrix.R)
  create_matrices_output <- function(arm_summaries, analysis_params) {
    
    # Use analysis_params directly (no alias needed)
    
    # Level 2: Prepare matrix data (adapted from create_matrix.R)
    # =================================================================
    # POLLY HUNG UPDATE - 2025-10-02
    # CHANGED: Require Peak_ID column, remove fallback logic
    # =================================================================
    prepare_matrix_data_integrated <- function(arm_data, params) {
      # NEW: Require Peak_ID column
      if (!"Peak_ID" %in% names(arm_data)) {
        stop("Peak_ID column not found in arm summaries. ",
             "Expected from updated arm_definitions structure.")
      }

      arm_column <- "Peak_ID"
      cat("Using Peak_ID identifiers for matrix rows\n")

      # Check required columns
      required_cols <- c("Sample", "Peak_ID", "mean_log2ratio")
      missing_cols <- setdiff(required_cols, names(arm_data))
      if (length(missing_cols) > 0) {
        stop("Missing required columns for matrix creation: ", paste(missing_cols, collapse = ", "))
      }

      # NEW: Handle duplicates by taking mean for each Peak_ID-Sample combination
      clean_data <- arm_data %>%
        dplyr::select(Sample, Peak_ID, mean_log2ratio) %>%
        dplyr::group_by(Sample, Peak_ID) %>%
        dplyr::summarise(mean_log2ratio = mean(mean_log2ratio, na.rm = TRUE), .groups = "drop") %>%
        dplyr::filter(!is.na(mean_log2ratio))

      # NEW: Create peak_id-by-sample matrix
      arm_matrix_df <- clean_data %>%
        tidyr::pivot_wider(names_from = Sample, values_from = mean_log2ratio,
                           names_sort = TRUE, values_fill = NA)

      # NEW: Convert to matrix using peak_id as row names
      peak_ids <- arm_matrix_df$Peak_ID  # Contains peak_ids like "10q_amp_tel_neg"
      arm_matrix <- as.matrix(arm_matrix_df[, -1])

      # Ensure matrix is numeric
      if (!is.numeric(arm_matrix)) {
        arm_matrix <- apply(arm_matrix, 2, as.numeric)
      }

      rownames(arm_matrix) <- peak_ids

      cat("Matrix created:", nrow(arm_matrix), "events ×", ncol(arm_matrix), "samples\n")
      return(list(
        log2_matrix = arm_matrix,
        n_arms = nrow(arm_matrix),
        n_samples = ncol(arm_matrix)
      ))
    }

    # # OLD CODE - COMMENTED OUT FOR REFERENCE
    # prepare_matrix_data_integrated <- function(arm_data, params) {
    #   # Check for functional_arm_id column (preferred) or Arm column (fallback)
    #   if ("functional_arm_id" %in% names(arm_data)) {
    #     arm_column <- "functional_arm_id"
    #     cat("Using functional arm identifiers for matrix rows\n")
    #   } else {
    #     arm_column <- "Arm"
    #     cat("Warning: functional_arm_id not found, using basic Arm identifiers\n")
    #   }
    #
    #   # Check required columns
    #   required_cols <- c("Sample", arm_column, "mean_log2ratio")
    #   missing_cols <- setdiff(required_cols, names(arm_data))
    #   if (length(missing_cols) > 0) {
    #     stop("Missing required columns for matrix creation: ", paste(missing_cols, collapse = ", "))
    #   }
    #
    #   # Handle duplicates by taking mean for each functional_arm-Sample combination
    #   clean_data <- arm_data %>%
    #     dplyr::select(Sample, !!rlang::sym(arm_column), mean_log2ratio) %>%
    #     dplyr::rename(Arm_ID = !!rlang::sym(arm_column)) %>%  # Standardize column name
    #     dplyr::group_by(Sample, Arm_ID) %>%
    #     dplyr::summarise(mean_log2ratio = mean(mean_log2ratio, na.rm = TRUE), .groups = "drop") %>%
    #     dplyr::filter(!is.na(mean_log2ratio))
    #
    #   # Create functional-arm-by-sample matrix
    #   arm_matrix_df <- clean_data %>%
    #     tidyr::pivot_wider(names_from = Sample, values_from = mean_log2ratio,
    #                        names_sort = TRUE, values_fill = NA)
    #
    #   # Convert to matrix using functional arm identifiers
    #   arm_names <- arm_matrix_df$Arm_ID  # Contains functional IDs like "1q_del_tel"
    #   arm_matrix <- as.matrix(arm_matrix_df[, -1])
    #
    #   # Ensure matrix is numeric
    #   if (!is.numeric(arm_matrix)) {
    #     arm_matrix <- apply(arm_matrix, 2, as.numeric)
    #   }
    #
    #   rownames(arm_matrix) <- arm_names
    #
    #   cat("Matrix created:", nrow(arm_matrix), "arms ×", ncol(arm_matrix), "samples\n")
    #   return(list(
    #     log2_matrix = arm_matrix,
    #     n_arms = nrow(arm_matrix),
    #     n_samples = ncol(arm_matrix)
    #   ))
    # }
    
    # Level 2: Generate matrix formats (adapted from create_matrix.R)
    generate_matrices_integrated <- function(matrix_data, params) {
      # GISTIC2-style copy number matrix (same as log2)
      cn_matrix <- matrix_data$log2_matrix
      
      # Generate discrete calls matrix
      calls_matrix <- matrix(0, nrow = nrow(cn_matrix), ncol = ncol(cn_matrix))
      rownames(calls_matrix) <- rownames(cn_matrix)
      colnames(calls_matrix) <- colnames(cn_matrix)
      
      # Apply thresholds
      calls_matrix[cn_matrix >= analysis_params$amp_threshold] <- 1   # Amplification
      calls_matrix[cn_matrix <= analysis_params$del_threshold] <- -1  # Deletion
      # Neutral remains 0
      
      cat("Matrix formats generated: log2, copy number, and discrete calls\n")
      
      return(list(
        log2_matrix = matrix_data$log2_matrix,
        cn_matrix = cn_matrix,
        calls_matrix = calls_matrix
      ))
    }
    
    # Level 2: Calculate matrix statistics (adapted from create_matrix.R)
    # =================================================================
    # POLLY HUNG UPDATE - 2025-10-02
    # CHANGED: Use Peak_ID column name instead of Functional_Arm
    # =================================================================
    calculate_matrix_stats_integrated <- function(matrices, params) {
      cn_matrix <- matrices$cn_matrix
      calls_matrix <- matrices$calls_matrix
      n_samples <- ncol(cn_matrix)

      # NEW: Generate comprehensive event statistics using Peak_ID
      arm_stats <- data.frame(
        Peak_ID = rownames(cn_matrix),  # NEW: Changed from Functional_Arm
        Mean_Log2Ratio = round(rowMeans(cn_matrix, na.rm = TRUE), 4),
        Median_Log2Ratio = round(apply(cn_matrix, 1, median, na.rm = TRUE), 4),
        SD_Log2Ratio = round(apply(cn_matrix, 1, sd, na.rm = TRUE), 4),
        N_Samples = n_samples,
        N_Amplified = rowSums(calls_matrix == 1, na.rm = TRUE),
        N_Deleted = rowSums(calls_matrix == -1, na.rm = TRUE),
        N_Neutral = rowSums(calls_matrix == 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::mutate(
          Pct_Amplified = round(N_Amplified / N_Samples * 100, 1),
          Pct_Deleted = round(N_Deleted / N_Samples * 100, 1),
          Pct_Altered = round((N_Amplified + N_Deleted) / N_Samples * 100, 1),
          Pct_Neutral = round(N_Neutral / N_Samples * 100, 1)
        ) %>%
        dplyr::arrange(desc(Pct_Altered), Peak_ID)  # NEW: Changed from Functional_Arm

      cat("Matrix statistics calculated:", nrow(arm_stats), "events with frequency data\n")
      return(arm_stats)
    }
    
    # Level 2: Create long format data (adapted from create_matrix.R)
    # =================================================================
    # POLLY HUNG UPDATE - 2025-10-02
    # CHANGED: Use Peak_ID column name instead of Functional_Arm
    # =================================================================
    create_long_format_integrated <- function(matrices, summary_stats, params) {
      cn_matrix <- matrices$cn_matrix
      calls_matrix <- matrices$calls_matrix

      # NEW: Create continuous data in long format with Peak_ID
      continuous_long <- as.data.frame(cn_matrix) %>%
        dplyr::mutate(Peak_ID = rownames(.)) %>%  # NEW: Changed from Functional_Arm
        tidyr::pivot_longer(-Peak_ID, names_to = "Sample", values_to = "Log2Ratio") %>%
        dplyr::arrange(Peak_ID, Sample)

      # NEW: Create discrete calls in long format with Peak_ID
      calls_long <- as.data.frame(calls_matrix) %>%
        dplyr::mutate(Peak_ID = rownames(.)) %>%  # NEW: Changed from Functional_Arm
        tidyr::pivot_longer(-Peak_ID, names_to = "Sample", values_to = "Call") %>%
        dplyr::arrange(Peak_ID, Sample) %>%
        dplyr::mutate(
          Call_Label = case_when(
            Call == 1 ~ "Amplification",
            Call == -1 ~ "Deletion",
            Call == 0 ~ "Normal",
            TRUE ~ "Unknown"
          )
        )

      # NEW: Combine continuous and discrete data using Peak_ID
      combined_long <- continuous_long %>%
        dplyr::left_join(calls_long %>% dplyr::select(Peak_ID, Sample, Call, Call_Label),
                         by = c("Peak_ID", "Sample")) %>%
        dplyr::left_join(summary_stats %>% dplyr::select(Peak_ID, Pct_Altered, Pct_Amplified, Pct_Deleted),
                         by = "Peak_ID") %>%
        dplyr::arrange(Peak_ID, Sample)

      cat("Long format created:", nrow(combined_long), "event-sample observations\n")
      return(combined_long)
    }
    
    # Level 2: Save matrix files (adapted from create_matrix.R)
    save_matrices_integrated <- function(matrices, stats, long_data, output_path, params) {
      if (!is.null(output_path) && analysis_params$save_results) {
        # Save all matrix formats to same directory as TSV files
        write.csv(matrices$log2_matrix,
                  file.path(output_path, "arm_log2ratio_matrix.csv"))
        write.csv(matrices$cn_matrix,
                  file.path(output_path, "arm_copynumber_matrix.csv"))
        write.csv(matrices$calls_matrix,
                  file.path(output_path, "arm_calls_matrix_gistic_style.csv"))
        write.csv(stats,
                  file.path(output_path, "arm_copynumber_summary.csv"), row.names = FALSE)
        write.csv(long_data,
                  file.path(output_path, "arm_copynumber_long.csv"), row.names = FALSE)
        
        cat("Matrix files saved to:", output_path, "\n")
        return(output_path)
      } else {
        cat("Matrix file output skipped\n")
        return(NULL)
      }
    }
    
    # Main Level 1 logic for matrix creation
    if (!analysis_params$create_matrices) {
      cat("Matrix creation skipped (create_matrices=FALSE)\n")
      return(NULL)
    }
    
    if (is.null(arm_summaries) || nrow(arm_summaries) == 0) {
      cat("No arm summaries available for matrix creation\n")
      return(NULL)
    }
    
    cat("=== Creating Comprehensive Matrix Outputs ===\n")
    
    # Step 1: Prepare matrix data
    matrix_data <- prepare_matrix_data_integrated(arm_summaries, analysis_params)
    
    # Step 2: Generate matrix formats
    matrices <- generate_matrices_integrated(matrix_data, analysis_params)
    
    # Step 3: Calculate statistics
    summary_stats <- calculate_matrix_stats_integrated(matrices, analysis_params)
    
    # Step 4: Create long format
    long_format <- create_long_format_integrated(matrices, summary_stats, analysis_params)
    
    # Step 5: Save matrix files
    matrix_output_dir <- if (!is.null(analysis_params$output_dir) && analysis_params$save_results) {
      analysis_params$output_dir
    } else {
      NULL
    }
    
    save_matrices_integrated(matrices, summary_stats, long_format, matrix_output_dir, analysis_params)
    
    cat("Matrix creation complete:\n")
    cat("- Formats:", matrix_data$n_arms, "arms ×", matrix_data$n_samples, "samples\n")
    cat("- Outputs: log2, copy number, discrete calls, summary, long format\n")
    
    return(list(
      matrices = list(
        log2_matrix = matrices$log2_matrix,
        cn_matrix = matrices$cn_matrix,
        calls_matrix = matrices$calls_matrix
      ),
      summary_stats = summary_stats,
      long_format = long_format,
      gistic_thresholds = list(
        amp_threshold = analysis_params$amp_threshold,
        del_threshold = analysis_params$del_threshold
      ),
      metadata = list(
        n_arms = matrix_data$n_arms,
        n_samples = matrix_data$n_samples,
        creation_time = Sys.time()
      )
    ))
  }
  
  # Level 1: Save comprehensive results
  save_comprehensive_results <- function(all_results, params) {
    
    # Level 2: Create output structure
    create_output_files <- function(results, output_path, cancer_name) {
      
      # Create cancer-specific directory
      results_dir <- file.path(output_path)
      dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Save arm-level summaries
      if (!is.null(results$arm_summaries)) {
        readr::write_tsv(results$arm_summaries,
                         file.path(results_dir, "enhanced_arm_summaries.tsv"))
      }
      
      # Save statistical results
      if (!is.null(results$gistic_results)) {
        readr::write_tsv(results$gistic_results,
                         file.path(results_dir, "gistic_statistics.tsv"))
      }
      
      if (!is.null(results$significant_arms)) {
        readr::write_tsv(results$significant_arms,
                         file.path(results_dir, "significant_arms.tsv"))
      }
      
      # Save parameters
      params_df <- data.frame(
        parameter = names(results$parameters),
        value = sapply(results$parameters, as.character),
        stringsAsFactors = FALSE
      )
      readr::write_tsv(params_df, file.path(results_dir, "analysis_parameters.tsv"))
      
      return(results_dir)
    }
    
    # Main Level 1 logic
    if (params$save_results) {
      output_path <- create_output_files(all_results, params$output_dir, params$cancer_type)
      cat("Results saved to:", output_path, "\n")
    }
  }
  
  # ============================================================================
  # MAIN ANALYSIS WORKFLOW
  # ============================================================================
  
  # Prepare analysis parameters
  analysis_params <- list(
    amp_threshold = amp_threshold,
    del_threshold = del_threshold,
    stringent_threshold = stringent_threshold,
    output_dir = output_dir,
    cancer_type = cancer_type,
    use_gistic = use_gistic,
    chunk_size = chunk_size,
    save_results = save_results,
    create_matrices = create_matrices,
    matrix_formats = matrix_formats,
    analysis_date = Sys.time()
  )
  
  # Step 1: Prepare and validate input data
  prepared_data <- prepare_input_data(segments, 
                                      breakpoints, 
                                      analysis_params)
  
  # Step 2: Perform arm-level analysis
  analysis_results <- analyze_arm_levels(prepared_data, 
                                         analysis_params)
  
  # Step 3: Perform statistical analysis
  statistical_results <- perform_statistical_analysis(analysis_results, 
                                                      analysis_params)
  
  # Step 4: Compile enhanced final results
  final_results <- list(
    segments = analysis_results$segments,
    arm_summaries = analysis_results$arm_summaries,
    arm_definitions = prepared_data$arm_definitions,
    gistic_results = statistical_results$gistic_results,
    significant_arms = statistical_results$significant_arms,
    parameters = analysis_params,
    metadata = list(
      n_samples = prepared_data$n_samples,
      n_arms = nrow(prepared_data$arm_definitions),
      has_enhanced_arms = prepared_data$has_enhanced_arms,
      n_significant_arms = if(!is.null(statistical_results$significant_arms)) nrow(statistical_results$significant_arms) else 0
    )
  )
  
  # Step 5: Save results
  save_comprehensive_results(final_results, analysis_params)
  
  # Step 6: Create matrix outputs (NEW INTEGRATION)
  matrix_results <- create_matrices_output(analysis_results$arm_summaries, 
                                           analysis_params)
  
  # Step 7: Compile Again 
  final_results <- list(
    segments = analysis_results$segments,
    arm_summaries = analysis_results$arm_summaries,
    arm_definitions = prepared_data$arm_definitions,
    gistic_results = statistical_results$gistic_results,
    significant_arms = statistical_results$significant_arms,
    
    # NEW MATRIX OUTPUTS (when create_matrices=TRUE)
    matrices = if (!is.null(matrix_results)) matrix_results$matrices else NULL,
    matrix_summary = if (!is.null(matrix_results)) matrix_results$summary_stats else NULL,
    matrix_long = if (!is.null(matrix_results)) matrix_results$long_format else NULL,
    
    parameters = analysis_params,
    metadata = list(
      n_samples = prepared_data$n_samples,
      n_arms = nrow(prepared_data$arm_definitions),
      has_enhanced_arms = prepared_data$has_enhanced_arms,
      n_significant_arms = if(!is.null(statistical_results$significant_arms)) nrow(statistical_results$significant_arms) else 0,
      matrices_created = !is.null(matrix_results),
      matrix_dimensions = if (!is.null(matrix_results)) paste(matrix_results$metadata$n_arms, "×", matrix_results$metadata$n_samples) else "None"
    )
  )
  
  # Step 8: Save Again 
  save_comprehensive_results(final_results, 
                             analysis_params)
  
  # Step 9: Save Object 
  save(final_results, 
       file = file.path(output_dir, "final_results.RData"))
  
  
  return(final_results)
}