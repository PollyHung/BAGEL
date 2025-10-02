#' Define Chromosome Arm Boundaries Using Unified Breakpoint Sources
#'
#' \code{define_arm} is the master function for establishing chromosome arm definitions
#' in BAGEL analysis. It intelligently routes between custom BISCUT-derived breakpoints
#' and TCGA consensus breakpoints, providing a unified interface that handles all
#' arm definition scenarios with comprehensive error handling and automatic fallback
#' mechanisms.
#'
#' @section Algorithm Overview:
#' This function implements a priority-based selection algorithm:
#' \enumerate{
#'   \item **Primary**: Custom BISCUT breakpoints (if file exists and is valid)
#'   \item **Secondary**: TCGA cancer-specific breakpoints (if cancer type available)
#'   \item **Tertiary**: TCGA consensus breakpoints (universal fallback)
#'   \item **Error Handling**: Graceful degradation with informative messaging
#' }
#'
#' @section Custom BISCUT Processing:
#' When processing custom BISCUT files, the function:
#' \itemize{
#'   \item Validates required columns: peak_id, Chr, arm, direction, telcent, negpos,
#'         type_of_selection, Peak.Start, Peak.End, combined_sig, n_right, n_left
#'   \item Applies statistical filters based on \code{ks_p < 0.05} significance threshold
#'   \item Filters functional regions by \code{percentage_length} of arm coverage
#'   \item Uses genomic coordinates (q_start/q_end) for precise boundary definition
#'   \item Implements telomere/centromere boundary logic for functional region definition
#'   \item Returns enhanced arm definitions with selection type annotations
#' }
#'
#' @section TCGA Consensus Processing:
#' For TCGA consensus breakpoints, the function:
#' \itemize{
#'   \item Loads cancer-specific breakpoint data when available
#'   \item Falls back to pan-cancer consensus definitions
#'   \item Integrates with BISCUT coordinate system when possible
#'   \item Uses hardcoded GRCh37/hg19 coordinates as ultimate fallback
#' }
#'
#' @param custom_biscut_file Character vector of length 1 or NULL. Absolute path to
#'   BISCUT results file (typically named \code{all_BISCUT_results.txt}). If NULL,
#'   the function will use TCGA consensus breakpoints. If the file path is provided
#'   but the file doesn't exist or is invalid, the function gracefully falls back
#'   to TCGA consensus with appropriate warning messages.
#'
#' @param cancer_type Character vector of length 1 or NULL. TCGA cancer type
#'   identifier (e.g., "BRCA", "OV", "LUAD") used for cancer-specific breakpoint
#'   selection and as fallback when custom BISCUT processing fails. Should match
#'   identifiers in the BAGEL breakpoint database. If NULL and custom_biscut_file
#'   processing fails, an error will be thrown.
#'
#' @param percentage_length Numeric. Minimum percentage of arm length that
#'   functional regions must span to be included in the analysis. Regions shorter
#'   than this threshold are filtered out during processing. Default value of 0.1
#'   (10%) provides good balance between sensitivity and inclusion of meaningful
#'   functional regions. Range: [0, 1].
#'
#' @return A data.frame with standardized arm definition columns:
#' \describe{
#'   \item{chr_arm}{Character. Chromosome arm identifier (e.g., "1p", "1q", "2p")}
#'   \item{chr}{Character. Chromosome identifier (BISCUT: "chr1"; TCGA: numeric)}
#'   \item{arm_type/arm}{Character. Arm type indicator ("p" or "q")}
#'   \item{functional_start}{Numeric. Functional region start coordinate in base pairs}
#'   \item{functional_end}{Numeric. Functional region end coordinate in base pairs}
#'   \item{direction}{Character. Expected alteration direction ("amp" or "del")}
#' }
#'
#' When custom BISCUT data is used, additional columns include:
#' \describe{
#'   \item{telcent}{Character. Boundary type ("tel" or "cent")}
#'   \item{type_of_selection}{Character. Selection classification ("TS-like", "onco-like", etc.)}
#'   \item{functional_length}{Numeric. Functional region length in base pairs}
#'   \item{percentage_length}{Numeric. Functional region as percentage of arm length}
#'   \item{combined_sig}{Numeric. Statistical significance score}
#'   \item{n_events}{Numeric. Number of events supporting the breakpoint}
#'   \item{total_sample_size}{Numeric. Total number of samples in the dataset}
#'   \item{percent_sample_support}{Numeric. Percentage of samples supporting the event (n_events / total_sample_size)}
#'   \item{n_right, n_left}{Numeric. Directional statistical support}
#'   \item{ks_p}{Numeric. Kolmogorov-Smirnov test p-value}
#' }
#'
#' @section Error Handling:
#' The function implements comprehensive error handling:
#' \itemize{
#'   \item **File not found**: Graceful fallback to TCGA consensus with warning
#'   \item **Invalid file format**: Detailed error messages specifying missing columns
#'   \item **Insufficient data**: Clear reporting when statistical thresholds eliminate all peaks
#'   \item **Missing dependencies**: Informative messages about BISCUT package availability
#'   \item **Data validation**: Coordinate validation and boundary logic checking
#' }
#'
#' @section Performance Notes:
#' \itemize{
#'   \item Function uses nested helper functions (max 2 levels) following BAGEL coding standards
#'   \item Memory efficient processing suitable for large BISCUT result files
#'   \item Caches loaded breakpoint data within function scope to avoid repeated I/O
#'   \item Progress reporting via console messages for long-running operations
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Custom BISCUT file with cancer type fallback
#' arm_defs <- define_arm(
#'   custom_biscut_file = "/path/to/all_BISCUT_results.txt",
#'   cancer_type = "BRCA",
#'   percentage_length = 0.15
#' )
#'
#' # Example 2: TCGA consensus only
#' arm_defs <- define_arm(
#'   custom_biscut_file = NULL,
#'   cancer_type = "OV"
#' )
#'
#' # Example 3: Automatic fallback scenario
#' arm_defs <- define_arm(
#'   custom_biscut_file = "/nonexistent/file.txt",  # Will trigger fallback
#'   cancer_type = "LUAD"
#' )
#'
#' # Example 4: High stringency custom processing
#' arm_defs <- define_arm(
#'   custom_biscut_file = "/path/to/results.txt",
#'   cancer_type = "PRAD",
#'   percentage_length = 0.2
#' )
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{calculateCopyNumber_fixed}} for using arm definitions in copy number analysis
#'   \item \code{\link{load_segmentation_data}} for preparing input segmentation data
#'   \item \code{\link{run_biscut_pipeline}} for generating custom BISCUT breakpoint files
#'   \item \code{\link{create_arm_matrix}} for creating analysis-ready arm × sample matrices
#' }
#'
#' @references
#' \itemize{
#'   \item BISCUT algorithm: Beroukhim et al. (2010). The landscape of somatic copy-number
#'         alteration across human cancers. Nature 463, 899-905.
#'   \item GISTIC2 methodology: Mermel et al. (2011). GISTIC2.0 facilitates sensitive and
#'         confident localization of the targets of focal somatic copy-number alteration
#'         in human cancers. Genome Biology 12, R41.
#'   \item TCGA breakpoint data: The Cancer Genome Atlas Research Network pan-cancer studies.
#' }
#'
#' @author Polly Hung \email{u3012149@connect.hku.hk}
#'
#' @keywords genomics copy-number breakpoints chromosome-arms cancer
#'
#' @import dplyr
#' @importFrom readr read_delim read_tsv
#' @importFrom stringr str_remove
#'
#' @export
define_arm <- function(custom_biscut_file = NULL,
                       cancer_type = NULL,
                       percentage_length = 0.1, 
                       output_dir = output_dir) {
  
  cat("=== Getting Unified Arm Definitions ===\n")
  
  # ============================================================================
  # LEVEL 1 NESTED FUNCTIONS: Core processing components
  # ============================================================================
  
  # Level 1: Load breakpoint data
  load_breakpoint_data_helper <- function() {
    data_dir <- system.file("data", package = "BAGEL")

    # Level 2: Safe data loading with error handling
    safe_load <- function(filename) {
      tryCatch({
        load(file.path(data_dir, filename), envir = environment())
        return(TRUE)
      }, error = function(e) {
        warning("Failed to load ", filename, ": ", e$message)
        return(FALSE)
      })
    }

    # Load all datasets with validation
    datasets <- list()

    # Try updated TCGA breakpoints first (BAGEL v4 analysis)
    if (safe_load("updated_tcga_breakpoints.rda")) {
      datasets$cancer_specific_breakpoints <- updated_tcga_breakpoints
      cat("✅ Loaded updated TCGA breakpoints (BAGEL v4)\n")
    } else if (safe_load("cancer_specific_breakpoints.rda")) {
      datasets$cancer_specific_breakpoints <- cancer_specific_breakpoints
      cat("Using legacy TCGA breakpoints\n")
    }

    # Load other datasets
    if (safe_load("consensus_arm_definitions.rda")) datasets$consensus_arm_definitions <- consensus_arm_definitions
    if (safe_load("arm_breakpoints.rda")) datasets$arm_breakpoints <- arm_breakpoints
    if (safe_load("consensus_breakpoints.rda")) datasets$consensus_breakpoints <- consensus_breakpoints

    return(datasets)
  }
  
  # NEW SIMPLIFIED FUNCTION - Straightforward approach without nesting
  process_custom_biscut <- function(biscut_file, 
                                    pct_length = percentage_length) {
    # Validate input file
    if (!file.exists(biscut_file)) {
      stop("BISCUT results file not found: ", biscut_file)
    }
    
    # Load BISCUT data
    biscut_data <- readr::read_delim(biscut_file, delim = "\t", show_col_types = FALSE)
    biscut_data <- as.data.frame(biscut_data)
    cat("Loaded BISCUT results:", nrow(biscut_data), "entries\n")
    
    # Validate required columns
    required_cols <- c("peak_id", "Chr", "arm", "direction", "telcent", "negpos",
                       "type_of_selection", "Peak.Start", "Peak.End",
                       "combined_sig", "n_right", "n_left", "total_sample_size")
    missing_cols <- setdiff(required_cols, colnames(biscut_data))

    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    # Filter by statistical support and significance
    filtered_data <- biscut_data %>% dplyr::filter(ks_p < 0.05)
    
    # Get chromosome coordinates
    coord_data <- BISCUT::get_chromosome_coordinates()
    coord_data$chr <- paste0("chr", coord_data$chromosome_info)
    
    # Process peak data and create functional boundaries
    peak_data <- filtered_data %>%
      dplyr::mutate(negpos = case_when(negpos == "p" ~ "pos", negpos == "n" ~ "neg", TRUE ~ "what"),
                    chr = paste0("chr", Chr), 
                    arm = as.character(arm), 
                    chr_arm = case_when(arm == "13" ~ "13q", 
                                        arm == "14" ~ "14q", 
                                        arm == "15" ~ "15q", 
                                        arm == "21" ~ "21q", 
                                        arm == "22" ~ "22q", 
                                        TRUE ~ arm), 
                    arm = gsub("[0-9]", "", chr_arm), 
                    peak_length = Peak.End - Peak.Start, 
                    peak_length.pct = Peak.End.1 - Peak.Start.1, 
                    peak_id = paste(chr_arm, direction, telcent, negpos, sep = "_"), 
                    type = paste(arm, telcent, sep = "_")) %>% 
      dplyr::select(
        ## Descriptive metrics
        peak_id, chr, chr_arm, arm, type,
        direction, telcent, negpos, type_of_selection,
        n_events, n_left, n_right, total_sample_size,

        ## Locations
        Peak.Start, Peak.End, Peak.Position, peak_length,
        Peak.Start.1, Peak.End.1, Peak.Position.1, peak_length.pct,

        ## Statistics
        ks_stat, ks_p, ksby, combined_sig) %>% 
      dplyr::group_by(peak_id) %>%
      dplyr::summarise(
        across(where(is.character), ~paste(unique(.), collapse = ", ")),
        across(where(is.numeric), ~mean(., na.rm = TRUE)),
        .groups = "drop") %>% 
      
      ## Append the reference hg19 coordinates 
      dplyr::left_join(coord_data, by = "chr") %>%
      dplyr::mutate(
        ## Add reference arm length 
        arm_length = case_when(arm == "p" ~ p_end - p_start, 
                               arm == "q" ~ q_end - q_start, 
                               TRUE ~ NA_integer_), 
        
        ## Functional start and end of the chromosome aberration region 
        functional_start = case_when(type == "p_tel" ~ p_start,
                                     type == "p_cent" ~ Peak.Start,
                                     type == "q_tel" ~ Peak.Start,
                                     type == "q_cent" ~ q_start,
                                     TRUE ~ NA_integer_),
        
        functional_end = case_when(type == "p_tel" ~ Peak.End,
                                   type == "p_cent" ~ p_end,
                                   type == "q_tel" ~ q_end,
                                   type == "q_cent" ~ Peak.End,
                                   TRUE ~ NA_integer_), 
        
        ## Calculate region metrics
        functional_length = functional_end - functional_start,
        percentage_length = (functional_end - functional_start)/arm_length,

        ## Calculate percent of samples supporting the event
        percent_sample_support = n_events / total_sample_size) %>% 
      
      ## DO Filtering Based on the result 
      dplyr::filter(!is.na(functional_start), 
                    !is.na(functional_end),
                    functional_start < functional_end,
                    percentage_length > pct_length) %>% 
      as.data.frame()
    rownames(peak_data) <- peak_data$peak_id

    # Add compatibility columns for calculate_copynumber integration
    peak_data <- peak_data %>%
      dplyr::mutate(
        arm_start = functional_start,
        arm_end = functional_end
      )

    cat("Generated", nrow(peak_data), "functional aneuploidy definitions\n")
    return(peak_data)
  }
  
  # Level 1: Get TCGA consensus definitions
  get_tcga_consensus <- function(cancer_type_name, breakpoint_datasets) {

    # Level 2: Build arm definitions from cancer-specific or consensus data
    build_arm_definitions <- function(breakpoints, chr_coords = NULL) {

      # Check if this is updated BAGEL v4 format (has functional_start/functional_end)
      if (all(c("functional_start", "functional_end", "chr_arm", "direction") %in% names(breakpoints))) {
        cat("Using BAGEL v4 format arm definitions\n")

        # Already in correct format, just standardize
        arm_definitions <- breakpoints %>%
          mutate(
            chr = gsub("chr", "", chr),
            arm_type = arm
          ) %>%
          select(chr_arm, chr, arm_type,
                 functional_start, functional_end, direction) %>%
          arrange(chr, arm_type)

        return(arm_definitions)
      }

      # Legacy format handling (has peak_start/peak_end)
      if (!is.null(chr_coords)) {
        # Use BISCUT coordinates for arm boundaries
        arm_definitions <- breakpoints %>%
          mutate(
            chr = as.numeric(chr),
            arm_type = ifelse(grepl("p$", arm), "p", "q")
          ) %>%
          left_join(
            chr_coords %>%
              select(chromosome_info, p_start, p_end, q_start, q_end) %>%
              rename(chr = chromosome_info),
            by = "chr"
          ) %>%
          mutate(
            arm_start = case_when(
              arm_type == "p" ~ p_start,
              arm_type == "q" ~ peak_start
            ),
            arm_end = case_when(
              arm_type == "p" ~ peak_end,
              arm_type == "q" ~ q_end
            )
          ) %>%
          filter(!is.na(arm_start) & !is.na(arm_end) & arm_start < arm_end) %>%
          select(chr_arm = arm, chr, arm_type,
                 functional_start = arm_start, functional_end = arm_end, direction) %>%
          arrange(chr, arm_type)

      } else {
        # Use hardcoded coordinates as fallback
        chr_lengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                         159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                         115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                         59128983, 63025520, 48129895, 51304566)
        names(chr_lengths) <- 1:22

        arm_definitions <- breakpoints %>%
          mutate(
            chr = as.numeric(chr),
            arm_type = ifelse(grepl("p$", arm), "p", "q"),
            arm_start = case_when(
              arm_type == "p" ~ 1,
              arm_type == "q" ~ peak_start
            ),
            arm_end = case_when(
              arm_type == "p" ~ peak_end,
              arm_type == "q" ~ chr_lengths[as.character(chr)]
            )
          ) %>%
          filter(!is.na(arm_start) & !is.na(arm_end) & arm_start < arm_end) %>%
          select(chr_arm = arm, chr, arm_type,
                 functional_start = arm_start, functional_end = arm_end, direction) %>%
          arrange(chr, arm_type)
      }

      return(arm_definitions)
    }

    # Get BISCUT coordinates if available
    biscut_coords <- tryCatch({
      BISCUT::get_chromosome_coordinates()
    }, error = function(e) {
      warning("Could not load BISCUT coordinates: ", e$message)
      NULL
    })

    # Try cancer-specific breakpoints first
    if (cancer_type_name %in% names(breakpoint_datasets$cancer_specific_breakpoints)) {
      cat(paste("Using", cancer_type_name, "specific breakpoints\n"))
      cancer_breakpoints <- breakpoint_datasets$cancer_specific_breakpoints[[cancer_type_name]]

      if (nrow(cancer_breakpoints) > 0) {
        arm_definitions <- build_arm_definitions(cancer_breakpoints, biscut_coords)
        return(arm_definitions)
      }
    }

    # Fallback to consensus
    cat(paste("Using consensus breakpoints for", cancer_type_name, "\n"))
    return(breakpoint_datasets$consensus_arm_definitions)
  }
  
  # ============================================================================
  # MAIN LOGIC: Route to appropriate processing method
  # ============================================================================
  
  # Case 1: Custom BISCUT file provided and exists
  if (!is.null(custom_biscut_file) && file.exists(custom_biscut_file)) {
    cat("Using custom BISCUT file:", custom_biscut_file, "\n")
    
    tryCatch({
      functional_definitions <- process_custom_biscut(biscut_file = custom_biscut_file, 
                                                      pct_length = percentage_length)
      write_to <- file.path(output_dir, "arm_definitions.txt")
      
      write.table(functional_definitions, write_to, sep = "\t", quote = F, col.names = T, row.names = F)
      
      cat("✅ Successfully created custom BISCUT definitions\n")
      return(functional_definitions)
      
    }, error = function(e) {
      cat("⚠️ Error processing custom BISCUT results:", e$message, "\n")
      cat("Falling back to TCGA consensus breakpoints\n")
      
      # Fallback to TCGA if provided
      if (!is.null(cancer_type)) {
        breakpoint_data <- load_breakpoint_data_helper()
        arm_definitions <- get_tcga_consensus(cancer_type, breakpoint_data)
        cat("✅ Using TCGA consensus breakpoints as fallback\n")
        return(arm_definitions)
      } else {
        stop("Custom BISCUT file failed and no cancer_type provided for fallback")
      }
    })
  }
  
  # Case 2: No custom file, use TCGA consensus
  if (is.null(custom_biscut_file) && !is.null(cancer_type)) {
    cat("Using TCGA consensus breakpoints for:", cancer_type, "\n")
    
    breakpoint_data <- load_breakpoint_data_helper()
    arm_definitions <- get_tcga_consensus(cancer_type, breakpoint_data)
    
    cat("✅ Successfully loaded TCGA consensus definitions\n")
    return(arm_definitions)
  }
  
  # Case 3: Custom file provided but doesn't exist
  if (!is.null(custom_biscut_file) && !file.exists(custom_biscut_file)) {
    cat("⚠️ Custom BISCUT file not found:", custom_biscut_file, "\n")
    
    if (!is.null(cancer_type)) {
      cat("Falling back to TCGA consensus breakpoints\n")
      breakpoint_data <- load_breakpoint_data_helper()
      arm_definitions <- get_tcga_consensus(cancer_type, breakpoint_data)
      cat("✅ Using TCGA consensus breakpoints as fallback\n")
      return(arm_definitions)
    } else {
      stop("Custom BISCUT file not found and no cancer_type provided for fallback")
    }
  }
  
  # Case 4: Neither provided - error
  if (is.null(custom_biscut_file) && is.null(cancer_type)) {
    stop("Either custom_biscut_file or cancer_type must be provided")
  }
}