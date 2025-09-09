# BAGEL v2.0: Breakpoint-based Arm-level Genomic Evaluation and Linkage

[![R build status](https://img.shields.io/badge/build-passing-brightgreen.svg)](https://github.com/your-username/BAGEL)
[![License: GPL-3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

BAGEL (Breakpoint-based Arm-level Genomic Evaluation and Linkage) is an R package for comprehensive copy number analysis using biologically meaningful chromosome arm breakpoints. Unlike traditional methods that rely on fixed percentage thresholds (e.g., 98% coverage), BAGEL uses breakpoints identified through the BISCUT algorithm to define chromosome arm boundaries for more accurate aneuploidy detection.

### Key Features

- **Breakpoint-based arm definitions**: Uses BISCUT-derived chromosome arm breakpoints instead of arbitrary thresholds
- **Pan-cancer breakpoint database**: Built-in breakpoint data from 29 TCGA cancer types
- **GISTIC-style statistical analysis**: Proper background modeling and significance testing
- **Memory-efficient processing**: Chunked processing for large datasets
- **Comprehensive error handling**: Robust input validation and informative error messages
- **Flexible input formats**: Supports multiple breakpoint data formats

## Installation

```r
# Install from GitHub
devtools::install_github("your-username/BAGEL")

# Or install from local source
devtools::install_local("path/to/BAGEL")
```

## Quick Start

### Complete BAGEL v2.0 Analysis Pipeline

Follow this step-by-step guide for a complete BAGEL analysis using both BISCUT breakpoint detection and BAGEL copy number analysis:

#### Step 1: Setup and Data Preparation

```r
# Load required libraries
library(BAGEL)
library(BISCUT)
library(dplyr)
library(readr)

# Set up your analysis parameters
data_dir <- "/path/to/your/data"
cancer_type <- "ovarian_serous_cystadenocarcinoma"  # Use your cancer type
output_dir <- file.path(data_dir, cancer_type, "bagel_v2_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Enable logging for tracking analysis progress
setup_bagel_logging(log_level = "INFO", log_file = file.path(output_dir, "analysis.log"))
```

#### Step 2: Run BISCUT Breakpoint Detection

```r
# Run BISCUT pipeline to identify cancer-specific breakpoints
# This requires a segmentation.seg file in your cancer_type directory
run_biscut_pipeline(
  cancer_folder = cancer_type, 
  results_dir = data_dir, 
  cores = 4
)
```

#### Step 3: Load Data and Breakpoints

```r
# Load segmentation data
segments <- load_segmentation_data(cancer_type, data_dir)

# Load custom BISCUT breakpoints (if available)
custom_biscut_file <- list.files(
  path = file.path(data_dir, cancer_type), 
  pattern = "all_BISCUT_results.txt", 
  recursive = TRUE, 
  full.names = TRUE
)

if (length(custom_biscut_file) > 0 && file.exists(custom_biscut_file)) {
  # Use cancer-specific BISCUT breakpoints
  arm_definitions <- create_custom_arm_definitions(custom_biscut_file)
} else {
  # Fall back to TCGA consensus breakpoints
  arm_definitions <- get_arm_definitions(cancer_type)
}
```

#### Step 4: Run BAGEL Copy Number Analysis

```r
# Run the main BAGEL analysis
bagel_results <- calculateCopyNumber_fixed(
  segments = segments,
  breakpoints = arm_definitions,
  amp_threshold = log2(2.5/2),       # Amplification threshold
  del_threshold = log2(1.5/2),       # Deletion threshold
  stringent_threshold = 0.9,         # Stringent threshold
  output_dir = output_dir,
  cancer_type = cancer_type,
  use_gistic = TRUE,                 # Enable GISTIC statistical analysis
  save_results = TRUE                # Save intermediate results
)
```

#### Step 5: Generate Analysis Matrices and Reports

```r
# Create chromosome arm copy number matrices
matrices <- create_arm_matrix(bagel_results$arm_summaries, output_dir)
bagel_results$matrices <- matrices

# Generate comprehensive analysis report
generate_analysis_report(cancer_type, bagel_results, output_dir)

# Save complete results
save(bagel_results, file = file.path(output_dir, "bagel_v2_complete_results.RData"))
```

#### Step 6: Create Chromosome Ideograms

```r
# Generate chromosome ideogram visualizations
ideogram_results <- plot_chromosome_ideograms(
  bagel_results = bagel_results,
  output_dir = output_dir,
  save_plots = TRUE
)
```

#### Step 7: Verify Analysis Completion

```r
# Summarize all output files and verify analysis completion
summarize_bagel_outputs(output_dir = output_dir)

# Quick completeness check
if (check_bagel_completeness(output_dir)) {
  cat("✅ Analysis completed successfully!\n")
} else {
  cat("⚠️ Analysis may be incomplete - check logs\n")
}
```

### Alternative: Simplified Analysis

For a quick analysis using pre-defined breakpoints:

```r
library(BAGEL)

# Load your segmentation data
segments <- read.table("your_segments.txt", header = TRUE)

# Run complete BAGEL workflow with pan-cancer consensus breakpoints
results <- bagel_workflow(
  segments = segments,
  cancer_type = "BRCA",
  output_dir = "./BAGEL_results"
)

# View significant arms
print(results$significant_arms)
```

### Using Cancer-Specific Breakpoints

```r
# Get cancer-specific breakpoint data
brca_breakpoints <- get_breakpoint_data("breast_invasive_carcinoma")

# Run analysis with cancer-specific breakpoints
results <- calculateCopyNumber_fixed(
  segments = segments,
  breakpoints = brca_breakpoints,
  cancer_type = "BRCA",
  use_gistic = TRUE
)
```

### Large Dataset Processing

```r
# For large datasets, use chunked processing
results <- calculateCopyNumber_fixed(
  segments = large_segments,
  breakpoints = "consensus",
  chunk_size = 50,  # Process 50 samples at a time
  cancer_type = "PANCAN"
)
```

## Data Format Requirements

### Segmentation Data

Your segmentation data should be a tab-delimited file named `segmentation.seg` with the following columns:

```r
segments <- data.frame(
  Sample = c("Sample1", "Sample1", "Sample2"),
  Chromosome = c(1, 1, 2),
  Start = c(1000, 50000, 1000),
  End = c(49999, 100000, 60000),
  Log2Ratios = c(0.1, -0.3, 0.5)  # or Segment_Mean
)
```

**Required columns:**
- `Sample`: Sample identifier (character)
- `Chromosome`: Chromosome number (1-22, numeric)
- `Start`: Segment start position in base pairs (numeric)
- `End`: Segment end position in base pairs (numeric)  
- `Log2Ratios` or `Segment_Mean`: Copy number values in log2 scale (numeric)

### Directory Structure

For the complete pipeline, organize your data as follows:

```
your_data_directory/
├── cancer_type_name/
│   ├── segmentation.seg          # Required: your segmentation data
│   └── bagel_v2_analysis/        # Created: output directory
│       ├── analysis.log
│       ├── arm_copynumber_matrix.csv
│       ├── BAGEL_V2_ANALYSIS_REPORT.md
│       ├── BAGEL_results_cancer_type/
│       └── newly_defined_arms_ideogram/
```

### Data Quality Requirements

- **Minimum samples**: At least 10 samples recommended for statistical analysis
- **Chromosome coverage**: Include autosomes 1-22 (sex chromosomes excluded)
- **Segment quality**: Remove low-quality segments with extreme log2 ratios (|log2ratio| > 3)
- **Sample coverage**: Ensure adequate genome coverage per sample (>80% recommended)

## Built-in Breakpoint Data

BAGEL v2.0 includes comprehensive breakpoint data from TCGA:

### Available Datasets

- `consensus_arm_definitions`: Pan-cancer consensus arm boundaries
- `cancer_specific_breakpoints`: Cancer-type specific breakpoints  
- `all_breakpoints`: Complete BISCUT results from all cancer types
- `consensus_breakpoints`: Consensus breakpoints across cancer types

### Accessing Breakpoint Data

```r
# Get pan-cancer consensus breakpoints
consensus <- get_breakpoint_data("consensus", "arm_definitions")

# Get cancer-specific breakpoints
brca <- get_breakpoint_data("breast_invasive_carcinoma")

# List available cancer types
names(cancer_specific_breakpoints)
```

### Available Cancer Types

The package includes breakpoint data for 29 TCGA cancer types:

- `adrenocortical_cancer`
- `bladder_urothelial_carcinoma`
- `brain_lower_grade_glioma`
- `breast_invasive_carcinoma`
- `cervical_endocervical_cancer`
- `cholangiocarcinoma`
- `colon_adenocarcinoma`
- `diffuse_large_b_cell_lymphoma`
- `esophageal_carcinoma`
- `glioblastoma_multiforme`
- `head_neck_squamous_cell_carcinoma`
- `kidney_clear_cell_carcinoma`
- `kidney_papillary_cell_carcinoma`
- `liver_hepatocellular_carcinoma`
- `lung_adenocarcinoma`
- `lung_squamous_cell_carcinoma`
- `mesothelioma`
- `ovarian_serous_cystadenocarcinoma`
- `pancreatic_adenocarcinoma`
- `pheochromocytoma_paraganglioma`
- `prostate_adenocarcinoma`
- `rectum_adenocarcinoma`
- `sarcoma`
- `skin_cutaneous_melanoma`
- `stomach_adenocarcinoma`
- `testicular_germ_cell_tumor`
- `thyroid_carcinoma`
- `uterine_carcinosarcoma`
- `uterine_corpus_endometrioid_carcinoma`

## Advanced Usage

### Custom Thresholds

```r
results <- calculateCopyNumber_fixed(
  segments = segments,
  breakpoints = "consensus",
  amp_threshold = log2(3/2),     # Higher amplification threshold
  del_threshold = log2(1/2),     # More stringent deletion threshold
  stringent_threshold = 0.8      # Stringent threshold multiplier
)
```

### GISTIC Analysis Options

```r
# Run GISTIC analysis separately
gistic_results <- gistic_analysis(
  segments = segments,
  arm_definitions = consensus_arm_definitions,
  q_threshold = 0.1  # More stringent q-value threshold
)
```

### Memory Management

```r
# Enable logging for monitoring
setup_bagel_logging(log_level = "INFO", log_file = "bagel.log")

# Process very large datasets in chunks
results <- process_segments_chunked(
  segments = huge_segments,
  arm_definitions = consensus_arm_definitions,
  chunk_size = 25  # Smaller chunks for memory constraints
)
```

## Expected Output Files

Following the complete BAGEL v2.0 pipeline, you can expect the following output files in your analysis directory:

### Core Analysis Files

- **`arm_copynumber_matrix.csv`**: Copy number calls matrix with samples as rows and chromosome arms as columns
- **`arm_log2ratio_matrix.csv`**: Log2 ratio values matrix for all samples and arms
- **`arm_copynumber_long.csv`**: Long format copy number data for further analysis
- **`arm_copynumber_summary.csv`**: Summary statistics for each chromosome arm
- **`bagel_v2_complete_results.RData`**: Complete R data object with all analysis results
- **`BAGEL_V2_ANALYSIS_REPORT.md`**: Comprehensive analysis report in Markdown format

### BAGEL Results Subdirectory

In the `BAGEL_results_[cancer_type]/` subdirectory:

- **`analysis_parameters.txt`**: Parameters used in the analysis
- **`arm_level_summaries.txt`**: Detailed arm-level copy number summaries
- **`arm_definitions.txt`**: Chromosome arm definitions used
- **`gistic_results.txt`**: GISTIC statistical analysis results
- **`significant_arms.txt`**: Statistically significant chromosome arms (q < 0.25)
- **`stringent_calls.txt`**: High-confidence alteration calls

### Visualization Files

In the `newly_defined_arms_ideogram/` subdirectory:

- **Individual chromosome PDFs**: `chromosome_1_ideogram.pdf`, `chromosome_2_ideogram.pdf`, etc.
- **Combined panel plot**: `all_chromosomes_ideogram_panel.pdf`
- **PNG versions**: High-resolution PNG versions of all plots

### BISCUT Output Files

If BISCUT analysis was run:

- **`breakpoint_files_[cancer_type]/`**: Directory containing breakpoint files for each sample
- **`results_[cancer_type]/`**: Directory containing BISCUT analysis results
- **`all_BISCUT_results.txt`**: Comprehensive BISCUT breakpoint results file

## Results Interpretation

### Understanding the Output

The BAGEL analysis provides multiple layers of information:

1. **Raw Copy Number Values**: Log2 ratios representing copy number relative to diploid
2. **Binary Calls**: Classified as amplification (+1), deletion (-1), or normal (0)
3. **Statistical Significance**: GISTIC q-values for each chromosome arm
4. **Alteration Frequencies**: Percentage of samples with alterations in each arm

### Key Metrics

- **Amplification threshold**: Default `log2(2.5/2)` ≈ 0.32
- **Deletion threshold**: Default `log2(1.5/2)` ≈ -0.42
- **Significance threshold**: q-value < 0.25 (adjustable)
- **Stringent threshold**: 90% of arm length must be altered

### Visualization Guide

The chromosome ideograms show:
- **Red regions**: Amplifications
- **Blue regions**: Deletions  
- **White regions**: Normal copy number
- **Intensity**: Alteration frequency across samples

## Troubleshooting

### Common Issues and Solutions

#### 1. BISCUT Installation Issues
```r
# If BISCUT is not available, install from appropriate source
# Check BISCUT documentation for installation instructions
if (!requireNamespace("BISCUT", quietly = TRUE)) {
  cat("BISCUT package required but not found\n")
  cat("Please install BISCUT before running the pipeline\n")
}
```

#### 2. Segmentation File Not Found
```bash
Error: Segmentation file not found: /path/to/cancer_type/segmentation.seg
```
**Solution**: Ensure your segmentation file is:
- Named exactly `segmentation.seg`
- Located in the correct cancer type directory
- Tab-delimited with proper column headers

#### 3. Memory Issues with Large Datasets
```r
# Use chunked processing for large datasets
results <- calculateCopyNumber_fixed(
  segments = segments,
  breakpoints = arm_definitions,
  chunk_size = 25,  # Reduce chunk size
  # ... other parameters
)
```

#### 4. No Significant Arms Found
This is normal for some datasets. Check:
- Sample size (need ≥10 samples for robust statistics)
- Data quality (extreme outliers can affect results)
- Threshold settings (try more lenient q-value thresholds)

#### 5. GISTIC Analysis Fails
```r
# Run without GISTIC if needed
results <- calculateCopyNumber_fixed(
  segments = segments,
  breakpoints = arm_definitions,
  use_gistic = FALSE  # Disable GISTIC
)
```

### Getting Help

- **Error logs**: Check `analysis.log` in your output directory
- **Function documentation**: Use `?function_name` for detailed help
- **Package vignettes**: Browse with `browseVignettes("BAGEL")`
- **Output verification**: Use `summarize_bagel_outputs()` to check completion

### Performance Tips

1. **Use appropriate cores**: Set `cores = 4` or based on your system
2. **Monitor memory**: Use chunked processing for >1000 samples
3. **Check data quality**: Remove extreme outliers before analysis
4. **Enable logging**: Always use `setup_bagel_logging()` for debugging

## Comparison with BAGEL v1.0

### Major Improvements in v2.0

1. **GISTIC-style Statistical Analysis**: Proper background modeling and significance testing
2. **Pan-cancer Breakpoint Database**: Built-in breakpoints from 29 cancer types
3. **Improved Error Handling**: Comprehensive input validation and informative errors
4. **Memory Efficiency**: Chunked processing for large datasets
5. **Fixed Logic Errors**: Resolved coordinate system issues and parameter inconsistencies
6. **Better Documentation**: Comprehensive function documentation and examples

### Migration from v1.0

Most v1.0 functions are still supported through legacy format conversion:

```r
# v1.0 style (still works)
old_breakpoints <- list(tel_bound = your_breakpoint_df)
results <- calculateCopyNumber_fixed(segments, old_breakpoints)

# v2.0 style (recommended)
results <- bagel_workflow(segments, "consensus")
```

## Citation

If you use BAGEL in your research, please cite:

```
[Your Citation Here - Replace with actual publication details]
```

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

## Support

- **Issues**: Report bugs and request features on [GitHub Issues](https://github.com/your-username/BAGEL/issues)
- **Documentation**: See function documentation with `?function_name`
- **Vignettes**: Browse vignettes with `browseVignettes("BAGEL")`

## License

BAGEL is licensed under the GPL-3 License. See [LICENSE](LICENSE) for details.

## References

1. BISCUT Algorithm: [Nature paper reference]
2. GISTIC Method: Beroukhim et al., Nature Genetics, 2007
3. TCGA Data: [TCGA reference]