#' Display Data Dictionary for BAGEL Analysis Results
#'
#' Provides comprehensive documentation of column definitions and statistical measures
#' used in BAGEL (Breakpoint-based Arm-level Genomic Evaluation and Linkage) analysis results.
#' This function serves as an interactive reference guide for interpreting output data.
#'
#' @param define Character string specifying the dictionary type to display.
#'   Currently supported:
#'   - "BISCUT": BISCUT peak identification and BAGEL analysis results (default)
#'   - "ARM": Arm-level copy number analysis results
#'   - "GISTIC": GISTIC-style statistical analysis results
#'
#' @details
#' The BISCUT dictionary covers six main groups of variables:
#' \itemize{
#'   \item Peak Identification & Location: Genomic coordinates and peak boundaries
#'   \item Gene Annotations: Gene symbols, RefSeq IDs, and cytogenetic locations
#'   \item Statistical Significance: KS test statistics and multiple testing corrections
#'   \item Selection Analysis: Positive/negative selection classification
#'   \item Analysis Parameters: Search limits and algorithm settings
#' }
#'
#' Statistical measures include Kolmogorov-Smirnov test statistics for measuring
#' deviation from background distribution, with Benjamini-Yekutieli correction
#' for multiple testing.
#'
#' @return Invisible NULL. Function is called for its side effect of printing
#'   formatted dictionary information to the console.
#'
#' @examples
#' # Display BISCUT analysis dictionary
#' dictionary()
#' dictionary("BISCUT")
#'
#' # Display arm-level analysis dictionary
#' dictionary("ARM")
#'
#' @author Polly Hung
#' @export
definition <- function(define = "BISCUT"){
  
  # Input validation
  if (!is.character(define) || length(define) != 1) {
    stop("Parameter 'define' must be a single character string")
  }

  define <- toupper(trimws(define))

  if (define == "BISCUT") {
    # Pre-build dictionary content for better performance
    dictionary_content <- paste(
      "",
      "----- Group 1: Peak Identification & Location -------",
      "",
      "peak_id              # Unique identifier for each significant peak",
      "                     # Format: [arm]_[direction]_[telcent]_[negpos]_[number]",
      "                     # e.g., '8p_del_tel_p_1' = chromosome 8p, deletion, telomere-bounded, positive selection, peak #1",
      "",
      "Chr                  # Chromosome number (8, 19)",
      "arm                  # Chromosome arm (8p, 19p)",
      "direction            # Type of alteration: 'del' (deletion) or 'amp' (amplification)",
      "telcent              # Peak boundary type: 'tel' (telomere-bounded) or 'cent' (centromere-bounded)",
      "",
      "peak_interval        # Genomic coordinates of the peak region",
      "                     # Format: 'Chr:Start-End' (e.g., '8:26259620-27069416')",
      "Peak.Start           # Peak start position in base pairs (26259620)",
      "Peak.End             # Peak end position in base pairs (27069416)",
      "Peak.Start.1         # Peak start as fraction of arm length (0.6059678)",
      "Peak.End.1           # Peak end as fraction of arm length (0.6251048)",
      "Peak.Position        # Peak center in base pairs (26579830)",
      "Peak.Position.1      # Peak center as fraction of arm length (0.6135349)",
      "",
      "----- Group 2: Gene Annotations -----",
      "",
      "Gene                 # Gene symbol within the peak (BNIP3L, PNMA2, etc.)",
      "RefSeqName           # RefSeq transcript identifier (NM_004331, NM_007257, etc.)",
      "Start.1              # Gene start as fraction of arm length (0.6055164)",
      "End.1                # Gene end as fraction of arm length (0.6062283)",
      "Cytoband             # Cytogenetic band location (p21.2, p13.12)",
      "",
      "----- Group 3: Statistical Significance -----",
      "",
      "ks_stat              # Kolmogorov-Smirnov test statistic (0.4576758)",
      "                     # Measures effect size - how much the observed distribution",
      "                     # deviates from background",
      "",
      "ks_p                 # Raw KS test p-value (1.071357e-51)",
      "log10_ks_p           # -log10 of raw p-value (50.97007)",
      "",
      "ksby                 # Benjamini-Yekutieli corrected p-value (5.031175e-49)",
      "log10_ksby           # -log10 of corrected p-value (48.29833)",
      "",
      "combined_sig         # Combined significance score (22.10498)",
      "                     # Formula: ks_stat × log10_ksby",
      "                     # Higher values = more significant peaks",
      "",
      "----- Group 4: Selection Analysis -----",
      "",
      "code                 # Selection direction: 'p' (positive) or 'n' (negative)",
      "negpos               # Same as code - selection type indicator",
      "",
      "type_of_selection    # Biological interpretation:",
      "                     # 'TS-like' = Tumor suppressor-like (deletions under positive selection)",
      "                     # 'onco-like' = Oncogene-like (amplifications under positive selection)",
      "                     # 'essential-like' = Essential gene-like (deletions under negative selection)",
      "                     # 'toxic-like' = Toxic amplification-like (amplifications under negative selection)",
      "",
      "----- Group 5: Analysis Parameters -----",
      "",
      "search_lowlim        # Lower search limit for peak detection (0)",
      "search_highlim       # Upper search limit for peak detection (1)",
      "",
      sep = "\n"
    )

    cat(dictionary_content)

  } else if (define == "ARM") {

    arm_content <- paste(
      "",
      "----- ARM-LEVEL COPY NUMBER ANALYSIS DICTIONARY -----",
      "",
      "Sample               # Sample identifier from input data",
      "Chromosome           # Chromosome number (1-22, X, Y)",
      "Arm                  # Chromosome arm (p or q)",
      "Chr_Arm              # Combined chromosome and arm (e.g., '1p', '2q')",
      "",
      "Copy_Number          # Arm-level copy number estimate",
      "Log2_Ratio           # Log2 ratio of copy number relative to diploid",
      "Call                 # Discrete copy number call (-2, -1, 0, 1, 2)",
      "                     # -2: Deep deletion, -1: Shallow deletion",
      "                     # 0: Neutral, 1: Gain, 2: Amplification",
      "",
      "Segment_Count        # Number of segments contributing to arm estimate",
      "Arm_Coverage         # Fraction of arm covered by segments",
      "Breakpoint_Used      # Whether BISCUT breakpoint was used (TRUE/FALSE)",
      "",
      sep = "\n"
    )

    cat(arm_content)

  } else if (define == "GISTIC") {

    gistic_content <- paste(
      "",
      "----- GISTIC-STYLE STATISTICAL ANALYSIS DICTIONARY -----",
      "",
      "Chromosome           # Chromosome number",
      "Arm                  # Chromosome arm (p or q)",
      "Start                # Arm start position (bp)",
      "End                  # Arm end position (bp)",
      "",
      "Frequency_Altered    # Fraction of samples with alteration",
      "Frequency_Amp        # Fraction of samples with amplification",
      "Frequency_Del        # Fraction of samples with deletion",
      "",
      "Q_Value              # False discovery rate corrected p-value",
      "G_Score              # GISTIC G-score (frequency × amplitude)",
      "Peak_Boundaries      # Peak region boundaries from BISCUT",
      "",
      "Residual_Q           # Residual q-value after peak removal",
      "Wide_Peak_Start      # Wide peak start coordinate",
      "Wide_Peak_End        # Wide peak end coordinate",
      "",
      sep = "\n"
    )

    cat(gistic_content)

  } else {
    available_types <- c("BISCUT", "ARM", "GISTIC")
    stop(sprintf("Dictionary type '%s' not recognized. Available types: %s",
                 define, paste(available_types, collapse = ", ")))
  }

  invisible(NULL)
}
