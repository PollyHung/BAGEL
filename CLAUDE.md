# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is the **BAGEL R Package** (v2.0.0) - a comprehensive cancer genomics analysis package for breakpoint-based arm-level copy number analysis. BAGEL uses biologically meaningful chromosome arm breakpoints (identified through BISCUT) instead of fixed percentage thresholds for more accurate aneuploidy detection.

### Core Architecture

The package follows a three-tier analysis architecture:

1. **Arm Definition Layer** (`R/arm_definition.R`): Unified module handling all breakpoint data sources
2. **Copy Number Analysis Layer** (`R/calculateCopyNumber-fixed.R`): GISTIC2-style statistical analysis
3. **Workflow Integration Layer** (`R/bagel-workflow.R`): Complete pipeline orchestration

### Key Modules

- **`arm_definition.R`**: **NEWEST** - Consolidated arm definition functions with unified wrapper `get_unified_arm_definitions()`
- **`calculateCopyNumber-fixed.R`**: Enhanced copy number analysis with proper GISTIC2 implementation
- **`bagel-workflow.R`**: Complete workflow functions for data loading and matrix creation
- **`process_biscut_results.R`**: BISCUT integration for custom breakpoint definitions
- **`plot-ideograms.R`**: Chromosome visualization with arm-level alterations
- **`utils.R`**: Logging, validation, and utility functions

## Development Commands

### Package Development (ARM64 macOS)

```bash
# IMPORTANT: Always use absolute R path for ARM64 compatibility
cd /Users/polly_hung/Desktop/BAGEL/packages/BAGEL

# Build package
/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/R CMD build .

# Install package
/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/R CMD INSTALL BAGEL_2.0.0.tar.gz

# Run tests
/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/Rscript -e "devtools::test()"

# Check package
/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/R CMD check BAGEL_2.0.0.tar.gz

# Update documentation (after @export changes)
/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/Rscript -e "roxygen2::roxygenise()"
```

### Running Main Analysis

```bash
# Run complete analysis pipeline
cd /Users/polly_hung/Desktop/BAGEL/packages/BAGEL
/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/Rscript runMe.R

# Run with timeout for long analyses
timeout 30m /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/Rscript runMe.R
```

## Core Data Flow

### Input Requirements

**Segmentation Data** (`segmentation.seg`):
```
Sample    Chromosome    Start       End       NumMarkers     Segment_Mean
TCGA-01   1            1           1000000     4				0.23
TCGA-01   1            1000001     2000000        8				  -0.45
```

**BISCUT Results** (optional `all_BISCUT_results.txt`):
- Peak coordinates with statistical significance
- Selection type classifications (TS-like, onco-like, etc.)
- Telomere/centromere boundary information

### Main Analysis Workflow

1. **Arm Definition Resolution**:
   ```r
   # NEW: Single unified function handles all scenarios
   arm_definitions <- get_unified_arm_definitions(
     custom_biscut_file = custom_biscut_file,    # NULL if not available
     cancer_type = cancer_type,                  # TCGA cancer type for fallback
     min_statistical_support = 10,
     min_combined_sig = 1.0
   )
   ```

2. **Copy Number Analysis**:
   ```r
   bagel_results <- calculateCopyNumber_fixed(
     segments = segments,
     breakpoints = arm_definitions,
     amp_threshold = 0.25,              # GISTIC2 standard
     del_threshold = -0.25,             # GISTIC2 standard
     stringent_threshold = 0.9,
     use_gistic = TRUE
   )
   ```

3. **Matrix Creation & Export**:
   ```r
   matrices <- create_arm_matrix(bagel_results$arm_summaries, output_dir)
   ```

### Output Structure

```
output_dir/
├── arm_copynumber_matrix.csv          # Log2 ratios (GISTIC2 style)
├── arm_calls_matrix_gistic_style.csv  # Discrete calls (-1, 0, +1)
├── arm_copynumber_summary.csv         # Statistical summaries
├── arm_copynumber_long.csv           # Long format with calls
├── BAGEL_V2_ANALYSIS_REPORT.md       # Comprehensive report
└── newly_defined_arms_ideogram/       # Chromosome visualizations
```

## Key Design Patterns

### Unified Arm Definition Pattern

**Problem Solved**: Previous code had complex conditional logic scattered across multiple files for handling custom BISCUT files vs TCGA consensus breakpoints.

**Solution**: Single entry point with automatic fallback:
```r
# Handles all scenarios in one call:
# - Custom BISCUT file exists → use custom definitions
# - Custom file missing/fails → fallback to TCGA consensus
# - No custom file provided → use TCGA consensus
arm_definitions <- get_unified_arm_definitions(custom_biscut_file, cancer_type)
```

### GISTIC2-Style Analysis Pattern

**Fixed Issue**: Original BAGEL used absolute copy numbers `2 * (2^log2ratio)` which were incompatible with GISTIC2 standards.

**Current Implementation**:
- Uses log2 ratios directly (centered at 0)
- Standard thresholds: amp_threshold=0.25, del_threshold=-0.25
- Generates both continuous values and discrete calls (-1, 0, +1)

### Error Handling Pattern

All major functions include comprehensive error handling:
```r
tryCatch({
  # Main operation
  result <- main_function()
}, error = function(e) {
  cat("⚠️ Error occurred:", e$message, "\n")
  # Graceful fallback or informative error
})
```

## Dependencies

### Required R Packages
- **Core**: `dplyr`, `tidyr`, `ggplot2`, `readr`, `stringr`, `magrittr`
- **Genomics**: `GenomicRanges`, `IRanges`, `S4Vectors` (Bioconductor)
- **Optional**: `BISCUT` (for custom breakpoint analysis)

### External Data
- **Built-in breakpoint data**: Pan-cancer TCGA consensus breakpoints in `data/`
- **Custom breakpoint files**: BISCUT-generated `all_BISCUT_results.txt`

## Testing

**Current Status**: Limited test framework (`tests/testthat/test-bagel-v2.R`)

**Testing Strategy**:
```r
# Basic function testing
devtools::test()

# Integration testing via runMe.R
# Use small test datasets for validation
```

## Code Style Guidelines

### Function Documentation
All exported functions use roxygen2 documentation:
```r
#' Function Title
#'
#' Detailed description of what the function does.
#'
#' @param param1 Type, description
#' @param param2 Type, description
#' @return Description of return value
#' @export
function_name <- function(param1, param2) {
  # Implementation
}
```

### Error Messages
Use informative error messages with context:
```r
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file,
       "\nExpected location: ", dirname(input_file))
}
```

### Progress Reporting
Use consistent progress reporting:
```r
cat("=== Step 1: Loading Data ===\n")
# operations
cat("✅ Data loaded successfully\n")
```

## Critical Development Notes

- **ARM64 Compatibility**: Always use absolute R path `/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/Rscript`
- **BISCUT Integration**: BISCUT package is optional but recommended for custom breakpoint analysis
- **Memory Management**: Large datasets require chunked processing (see `chunk_size` parameters)
- **GISTIC2 Compliance**: All copy number outputs follow GISTIC2 conventions for downstream compatibility

## Common Development Tasks

### Adding New Functions
1. Write function in appropriate `R/` file
2. Add `@export` tag for public functions
3. Update `NAMESPACE`: `roxygen2::roxygenise()`
4. Add tests in `tests/testthat/`
5. Update documentation

### Modifying Arm Definitions
- **Primary location**: `R/arm_definition.R`
- **Entry point**: `get_unified_arm_definitions()`
- Test with both custom BISCUT and TCGA consensus scenarios

### Updating Analysis Pipeline
- **Main pipeline**: `runMe.R`
- **Workflow functions**: `R/bagel-workflow.R`
- **Core analysis**: `R/calculateCopyNumber-fixed.R`

### Troubleshooting
- **R issues**: Check package installation and library paths
- **BISCUT issues**: Verify BISCUT package availability and input file formats
- **Memory issues**: Adjust chunk_size parameters for large datasets
- **Output issues**: Check file permissions and directory creation

### R Function Organization Style:
- Maximum 2 levels of function nesting allowed
- Prefer comprehensive, self-contained functions over scattered helper functions
- Keep related logic grouped within the main function scope
- Use nested functions for intermediate steps rather than separate top-level functions
- Naming of a function should be kept as a_b.R, only one _ underscore allowed. Keep it to 2 words. Preferably verb_noun.R combination
- The function name should be consistent as the script name. So if the function is called define_arm then the script should be correspondingly define_arm.R

Example Pattern 
```
main_function <- function(params) {
  # Nested helper functions (level 1)
  helper_step1 <- function(x) {
    # Optional deeper nesting (level 2) if absolutely needed
    inner_calc <- function(y) { return(y * 2) }
    return(inner_calc(x) + 1)
  }
  
  helper_step2 <- function(x) {
    return(x^2)
  }
  
  # Main logic using helpers
  result1 <- helper_step1(params$a)
  result2 <- helper_step2(result1)
  return(result2)
}
```