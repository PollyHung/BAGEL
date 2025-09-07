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

### Basic Analysis

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

Your segmentation data should be a data frame with the following columns:

```r
segments <- data.frame(
  Sample = c("Sample1", "Sample1", "Sample2"),
  Chromosome = c(1, 1, 2),
  Start = c(1000, 50000, 1000),
  End = c(49999, 100000, 60000),
  Log2Ratios = c(0.1, -0.3, 0.5)  # or Segment_Mean
)
```

Required columns:
- `Sample`: Sample identifier
- `Chromosome`: Chromosome number (1-22)
- `Start`: Segment start position (bp)
- `End`: Segment end position (bp)  
- `Log2Ratios` or `Segment_Mean`: Copy number values in log2 scale

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

## Output Description

The analysis returns a comprehensive results object containing:

- `segments`: Processed segments with arm annotations
- `arm_summaries`: Arm-level copy number summaries
- `stringent_summaries`: High-confidence alteration calls
- `gistic_results`: Statistical analysis results (if enabled)
- `significant_arms`: Statistically significant chromosome arms
- `parameters`: Analysis parameters used

### Results Files

When `save_results = TRUE`, the following files are generated:

- `arm_level_summaries.txt`: Arm-level copy number summaries
- `stringent_calls.txt`: High-confidence alteration calls
- `gistic_results.txt`: GISTIC statistical analysis results
- `significant_arms.txt`: Significantly altered chromosome arms
- `analysis_parameters.txt`: Parameters used in the analysis

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