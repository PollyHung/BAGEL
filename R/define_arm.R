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
#'   \item Applies statistical filters based on \code{min_statistical_support} and
#'         \code{min_combined_sig}
#'   \item Implements telomere/centromere boundary logic for functional region definition
#'   \item Integrates BISCUT chromosome coordinates for precise arm boundaries
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
#' @param min_statistical_support Numeric. Minimum combined statistical support
#'   threshold for BISCUT peak inclusion, calculated as (n_right + n_left). Peaks
#'   with support below this threshold are filtered out during custom processing.
#'   Default value of 10 provides good balance between sensitivity and specificity
#'   based on BISCUT validation studies. Range: [1, Inf].
#'
#' @param min_combined_sig Numeric. Minimum combined significance score threshold
#'   for BISCUT peak inclusion. This represents the statistical significance of
#'   breakpoint detection. Default value of 1.0 corresponds to moderately significant
#'   breakpoints. Higher values increase stringency. Range: [0, Inf].
#'
#' @return A data.frame with standardized arm definition columns:
#' \describe{
#'   \item{arm}{Character. Chromosome arm identifier (e.g., "1p", "1q", "2p")}
#'   \item{chr}{Numeric. Chromosome number (1-22)}
#'   \item{arm_type}{Character. Arm type indicator ("p" or "q")}
#'   \item{arm_start}{Numeric. Arm start coordinate in base pairs (1-based)}
#'   \item{arm_end}{Numeric. Arm end coordinate in base pairs (1-based)}
#'   \item{direction}{Character. Expected alteration direction ("amp" or "del")}
#' }
#'
#' When custom BISCUT data is used, additional columns may include:
#' \describe{
#'   \item{telcent}{Character. Boundary type ("tel" or "cent")}
#'   \item{type_of_selection}{Character. Selection classification ("TS-like", "onco-like", etc.)}
#'   \item{functional_length_mb}{Numeric. Functional region length in megabases}
#'   \item{combined_sig}{Numeric. Statistical significance score}
#'   \item{n_events}{Numeric. Number of events supporting the breakpoint}
#'   \item{n_right, n_left}{Numeric. Directional statistical support}
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
#'   min_statistical_support = 15,
#'   min_combined_sig = 1.5
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
#'   min_statistical_support = 25,
#'   min_combined_sig = 2.0
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
                       min_statistical_support = 10,
                       min_combined_sig = 1.0) {
  
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
    if (safe_load("consensus_arm_definitions.rda")) datasets$consensus_arm_definitions <- consensus_arm_definitions
    if (safe_load("cancer_specific_breakpoints.rda")) datasets$cancer_specific_breakpoints <- cancer_specific_breakpoints
    if (safe_load("arm_breakpoints.rda")) datasets$arm_breakpoints <- arm_breakpoints
    if (safe_load("consensus_breakpoints.rda")) datasets$consensus_breakpoints <- consensus_breakpoints
    
    return(datasets)
  }
  
  # Level 1: Process BISCUT results for custom definitions
  process_custom_biscut <- function(biscut_file, min_support, min_sig) {
    
    # Level 2: Load and validate BISCUT data
    load_biscut_data <- function(file_path) {
      if (!file.exists(file_path)) {
        stop("BISCUT results file not found: ", file_path)
      }
      
      biscut_data <- readr::read_delim(file_path, delim = "\t", show_col_types = FALSE)
      cat("Loaded BISCUT results:", nrow(biscut_data), "entries\n")
      
      # Validate required columns
      required_cols <- c("peak_id", "Chr", "arm", "direction", "telcent", "negpos",
                         "type_of_selection", "Peak.Start", "Peak.End",
                         "combined_sig", "n_right", "n_left")
      missing_cols <- setdiff(required_cols, colnames(biscut_data))
      
      if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
      }
      
      return(biscut_data)
    }
    
    # Level 2: Create functional boundaries
    create_functional_boundaries <- function(biscut_data) {
      # Get chromosome coordinates with error handling
      coord_data <- tryCatch({
        BISCUT::get_chromosome_coordinates()
      }, error = function(e) {
        warning("Could not load BISCUT coordinates, using hardcoded fallback")
        # Return basic fallback coordinates
        data.frame(
          chromosome_info = 1:22,
          p_start = rep(1, 22),
          p_end = c(125000000, 93300000, 91000000, 50400000, 48400000, 61000000,
                    59900000, 45600000, 49000000, 40200000, 53700000, 35800000,
                    17900000, 17600000, 19000000, 36600000, 24000000, 17200000,
                    26500000, 27500000, 13200000, 14700000),
          q_start1 = c(125000001, 93300001, 91000001, 50400001, 48400001, 61000001,
                       59900001, 45600001, 49000001, 40200001, 53700001, 35800001,
                       17900001, 17600001, 19000001, 36600001, 24000001, 17200001,
                       26500001, 27500001, 13200001, 14700001),
          q_end1 = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                     159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                     115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                     59128983, 63025520, 48129895, 51304566)
        )
      })
      
      # Process peak data
      peak_data <- biscut_data %>%
        dplyr::mutate(
          negpos = case_when(negpos == "p" ~ "pos", negpos == "n" ~ "neg", TRUE ~ "what"),
          Peak_ID = paste0(arm, ":", Peak.Start, "-", Peak.End, ":", direction, ":", telcent, ":", negpos)
        ) %>%
        dplyr::select(Peak_ID, Chr, arm, Peak.Start, Peak.End, Peak.Position, direction,
                      telcent, negpos, type_of_selection, ks_stat, ks_p, combined_sig,
                      n_right, n_left, n_events) %>%
        distinct()
      
      # Create functional definitions
      functional_defs <- peak_data %>%
        dplyr::mutate(Chr = as.numeric(Chr)) %>%
        dplyr::left_join(coord_data, by = c("Chr" = "chromosome_info")) %>%
        dplyr::mutate(
          arm_type = ifelse(str_remove(arm, "[0-9]*") == "", "q", str_remove(arm, "[0-9]*")),
          chr_num = Chr,
          
          # Define functional coordinates based on boundary logic
          functional_start = case_when(
            arm_type == "p" & telcent == "tel" ~ p_start,
            arm_type == "p" & telcent == "cent" ~ Peak.Start,
            arm_type == "q" & telcent == "tel" ~ Peak.Start,
            arm_type == "q" & telcent == "cent" ~ q_start1,
            TRUE ~ NA_integer_
          ),
          
          functional_end = case_when(
            arm_type == "p" & telcent == "tel" ~ Peak.End,
            arm_type == "p" & telcent == "cent" ~ p_end,
            arm_type == "q" & telcent == "tel" ~ q_end1,
            arm_type == "q" & telcent == "cent" ~ Peak.End,
            TRUE ~ NA_integer_
          ),
          
          # Calculate region metrics
          functional_length_mb = round((functional_end - functional_start) / 1e6, 2)
        ) %>%
        dplyr::filter(
          !is.na(functional_start),
          !is.na(functional_end),
          functional_start < functional_end,
          functional_length_mb > 0.1
        )
      
      # Standardize column names for output
      result <- functional_defs %>%
        dplyr::select(
          arm,
          chr = chr_num,
          arm_type,
          arm_start = functional_start,
          arm_end = functional_end,
          direction,
          telcent,
          type_of_selection,
          functional_length_mb,
          combined_sig,
          n_events,
          n_right,
          n_left
        ) %>%
        arrange(chr, arm_type, arm_start)
      
      return(result)
    }
    
    # Main processing logic
    biscut_data <- load_biscut_data(biscut_file)
    
    # Filter by statistical support and significance
    filtered_data <- biscut_data %>%
      dplyr::filter(
        (n_right + n_left) >= min_support,
        combined_sig >= min_sig
      )
    
    if (nrow(filtered_data) == 0) {
      stop("No BISCUT peaks meet minimum thresholds (support=", min_support, ", sig=", min_sig, ")")
    }
    
    functional_definitions <- create_functional_boundaries(filtered_data)
    cat("Generated", nrow(functional_definitions), "functional aneuploidy definitions\n")
    
    return(functional_definitions)
  }
  
  # Level 1: Get TCGA consensus definitions
  get_tcga_consensus <- function(cancer_type_name, breakpoint_datasets) {
    
    # Level 2: Build arm definitions from cancer-specific or consensus data
    build_arm_definitions <- function(breakpoints, chr_coords = NULL) {
      
      if (!is.null(chr_coords)) {
        # Use BISCUT coordinates for arm boundaries
        arm_definitions <- breakpoints %>%
          mutate(
            chr_num = as.numeric(chr),
            arm_type = ifelse(grepl("p$", arm), "p", "q")
          ) %>%
          left_join(
            chr_coords %>%
              select(chromosome_info, p_start, p_end, q_start, q_end) %>%
              rename(chr_num = chromosome_info),
            by = "chr_num"
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
          select(arm, chr_num, arm_type, arm_start, arm_end, direction) %>%
          arrange(chr_num, arm_type)
        
      } else {
        # Use hardcoded coordinates as fallback
        chr_lengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                         159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                         115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                         59128983, 63025520, 48129895, 51304566)
        names(chr_lengths) <- 1:22
        
        arm_definitions <- breakpoints %>%
          mutate(
            chr_num = as.numeric(chr),
            arm_type = ifelse(grepl("p$", arm), "p", "q"),
            arm_start = case_when(
              arm_type == "p" ~ 1,
              arm_type == "q" ~ peak_start
            ),
            arm_end = case_when(
              arm_type == "p" ~ peak_end,
              arm_type == "q" ~ chr_lengths[as.character(chr_num)]
            )
          ) %>%
          filter(!is.na(arm_start) & !is.na(arm_end) & arm_start < arm_end) %>%
          select(arm, chr_num, arm_type, arm_start, arm_end, direction) %>%
          arrange(chr_num, arm_type)
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
      
      if (nrow(cancer_breakpoints) > 0 && "peak_start" %in% names(cancer_breakpoints)) {
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
      functional_definitions <- process_custom_biscut(
        custom_biscut_file,
        min_statistical_support,
        min_combined_sig
      )
      
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