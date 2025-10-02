# Updating TCGA Fallback Breakpoints

This directory contains scripts for updating the TCGA fallback breakpoint data used in `define_arm.R`.

## Overview

The `define_arm()` function uses a fallback hierarchy:
1. **Primary**: Custom BISCUT breakpoints (if provided)
2. **Secondary**: TCGA cancer-specific breakpoints (updated BAGEL v4)
3. **Tertiary**: TCGA consensus breakpoints (legacy)

## Update Process

### Step 1: Run the consolidation script

```bash
cd /Users/polly_hung/Desktop/BAGEL/packages/BAGEL/data-raw
/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/Rscript update_tcga_fallback.R
```

This script:
- Reads all 29 TCGA arm_definitions.txt files from `~/Desktop/BAGEL/data/biscut+bagel/tcga_*/`
- Consolidates them into a single R list object
- Saves to `packages/BAGEL/data/updated_tcga_breakpoints.rda`
- Creates metadata CSV with arm counts per cancer type

### Step 2: Rebuild the package

```bash
cd /Users/polly_hung/Desktop/BAGEL/packages/BAGEL
/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/R CMD build .
/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/R CMD INSTALL BAGEL_2.0.0.tar.gz
```

### Step 3: Verify the update

```r
library(BAGEL)

# Test fallback with a TCGA cancer type
arm_defs <- define_arm(
  custom_biscut_file = NULL,
  cancer_type = "breast_invasive_carcinoma"
)

# Should see: "✅ Loaded updated TCGA breakpoints (BAGEL v4)"
# Should see: "Using BAGEL v4 format arm definitions"
```

## Data Format

### Input files (arm_definitions.txt)
```
chr    chr_arm  arm  functional_start  functional_end  direction  type_of_selection  ...
chr1   1p       p    1                 121535434       del        TS-like           ...
chr1   1q       q    142535434         249250621       amp        onco-like         ...
```

### Output object (updated_tcga_breakpoints)
```r
list(
  breast_invasive_carcinoma = <data.frame with arm definitions>,
  lung_adenocarcinoma = <data.frame>,
  ...
)
```

## Cancer Types Included

29 TCGA cancer types:
- adrenocortical_cancer
- bladder_urothelial_carcinoma
- brain_lower_grade_glioma
- breast_invasive_carcinoma
- cervical_endocervical_cancer
- cholangiocarcinoma
- colon_adenocarcinoma
- diffuse_large_b_cell_lymphoma
- esophageal_carcinoma
- glioblastoma_multiforme
- head_neck_squamous_cell_carcinoma
- kidney_clear_cell_carcinoma
- kidney_papillary_cell_carcinoma
- liver_hepatocellular_carcinoma
- lung_adenocarcinoma
- lung_squamous_cell_carcinoma
- mesothelioma
- ovarian_serous_cystadenocarcinoma
- pancreatic_adenocarcinoma
- pheochromocytoma_paraganglioma
- prostate_adenocarcinoma
- rectum_adenocarcinoma
- sarcoma
- skin_cutaneous_melanoma
- stomach_adenocarcinoma
- testicular_germ_cell_tumor
- thyroid_carcinoma
- uterine_carcinosarcoma
- uterine_corpus_endometrioid_carcinoma

## Modified Code

### define_arm.R changes

1. **load_breakpoint_data_helper()**: Now tries to load `updated_tcga_breakpoints.rda` first
2. **build_arm_definitions()**: Added BAGEL v4 format detection and handling
3. **get_tcga_consensus()**: Updated to handle new data structure

### Key improvements

- **Automatic format detection**: Detects BAGEL v4 format vs legacy format
- **Graceful fallback**: Falls back to legacy breakpoints if updated ones unavailable
- **Consistent output**: Returns standardized arm_definitions regardless of source
- **Informative logging**: Clear console messages about which data source is used

## Notes

- The updated breakpoints are based on BAGEL v4 analysis with direct coordinate mapping
- Each cancer type has its own custom arm definitions derived from BISCUT analysis
- The fallback still supports legacy TCGA consensus breakpoints for backward compatibility
- Cancer type names are standardized (lowercase with underscores, no "tcga_" prefix)