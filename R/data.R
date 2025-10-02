#' @title cytoband.hg19
#'
#' @description Cytoband version hg19
#'
#' @format ## cytoband.hg19
#' \describe{
#'   \item{chr}{chromosome}
#'   \item{arm}{chromosome arms}
#'   \item{start}{arm start}
#'   \item{end}{arm end}
#' }
#' @source <https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/>
"cytoband.hg19"

#' @title cytoband.hg38
#'
#' @description Cytoband version hg38, with centromeres
#'
#' @format ## cytoband.hg38
#' \describe{
#'   \item{chr}{chromosome}
#'   \item{arm}{chromosome arms}
#'   \item{start}{arm start}
#'   \item{end}{arm end}
#' }
#' @source <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/>
"cytoband.hg38"

#' @title cytoband.T2T
#'
#' @description Cytoband version T2T, with centromeres
#'
#' @format ## cytoband.T2T
#' \describe{
#'   \item{chr}{chromosome}
#'   \item{arm}{chromosome arms}
#'   \item{start}{arm start}
#'   \item{end}{arm end}
#' }
#' @source <https://github.com/marbl/CHM13?tab=readme-ov-file>
"cytoband.T2T"

#' @title aneu
#'
#' @description Aneuploidy dataframe
#'
#' @format ## aneu
#' \describe{
#'   \item{direction}{1, -1, -2, 2}
#' }
#' @source NA
"aneu"

#' Pan-Cancer Breakpoint Data from BISCUT Analysis
#'
#' Complete BISCUT breakpoint detection results across 29 cancer types from TCGA.
#' This dataset contains detailed peak information for chromosome arm breakpoints
#' identified using the BISCUT algorithm.
#'
#' @format A data frame with variables:
#' \describe{
#'   \item{peak_id}{Unique identifier for each breakpoint peak}
#'   \item{arm}{Chromosome arm (e.g., "1p", "1q")}
#'   \item{direction}{Peak direction: "amp" (amplification) or "del" (deletion)}
#'   \item{Chr}{Chromosome number}
#'   \item{Peak.Start}{Peak start position (bp)}
#'   \item{Peak.End}{Peak end position (bp)}
#'   \item{n_events}{Number of events contributing to the peak}
#'   \item{ks_stat}{Kolmogorov-Smirnov statistic}
#'   \item{ks_p}{K-S test p-value}
#'   \item{log10_ks_p}{-log10(p-value)}
#'   \item{cancer_type}{Cancer type from TCGA}
#'   \item{...}{Additional BISCUT analysis columns}
#' }
#' @source TCGA pan-cancer segmentation data processed with BISCUT algorithm
"all_breakpoints"

#' Arm-Level Breakpoint Summaries by Cancer Type
#'
#' Summarized breakpoint information at the chromosome arm level for each cancer type.
#' This dataset provides a more condensed view of breakpoint patterns across cancer types.
#'
#' @format A data frame with variables:
#' \describe{
#'   \item{arm}{Chromosome arm (e.g., "1p", "1q")}
#'   \item{direction}{Breakpoint direction: "amp" or "del"}
#'   \item{cancer_type}{Cancer type}
#'   \item{chr}{Chromosome number}
#'   \item{n_peaks}{Number of peaks detected for this arm}
#'   \item{peak_start}{Minimum peak start position}
#'   \item{peak_end}{Maximum peak end position}
#' }
#' @source Derived from all_breakpoints data
"arm_breakpoints"

#' Pan-Cancer Consensus Breakpoints
#'
#' Consensus breakpoint definitions based on breakpoints found in multiple cancer types.
#' These represent the most robust and reproducible chromosome arm breakpoints across
#' different cancer contexts.
#'
#' @format A data frame with variables:
#' \describe{
#'   \item{arm}{Chromosome arm (e.g., "1p", "1q")}
#'   \item{direction}{Consensus direction: "amp" or "del"}
#'   \item{chr}{Chromosome number}
#'   \item{n_cancer_types}{Number of cancer types supporting this breakpoint}
#'   \item{median_start}{Median start position across cancer types}
#'   \item{median_end}{Median end position across cancer types}
#'   \item{cancer_types}{Semicolon-separated list of contributing cancer types}
#' }
#' @source Consensus analysis of arm_breakpoints across cancer types
"consensus_breakpoints"

#' Cancer-Specific Breakpoint Data
#'
#' A named list containing breakpoint data for each individual cancer type.
#' Each element is a data frame with breakpoint information specific to that cancer type.
#'
#' @format A named list where each element is a data frame with variables:
#' \describe{
#'   \item{arm}{Chromosome arm}
#'   \item{direction}{Breakpoint direction}
#'   \item{chr}{Chromosome number}
#'   \item{n_peaks}{Number of peaks}
#'   \item{peak_start}{Peak start position}
#'   \item{peak_end}{Peak end position}
#'   \item{cancer_type}{Cancer type (same for all rows in each list element)}
#' }
#' @source Derived from arm_breakpoints, split by cancer type
"cancer_specific_breakpoints"

#' Pan-Cancer Consensus Arm Definitions
#'
#' Standardized chromosome arm boundary definitions based on pan-cancer consensus
#' breakpoints. This dataset provides the arm boundaries used for copy number analysis
#' in BAGEL, replacing traditional cytogenetic band definitions with biologically
#' meaningful breakpoint-derived boundaries.
#'
#' @format A data frame with variables:
#' \describe{
#'   \item{arm}{Chromosome arm identifier (e.g., "1p", "1q")}
#'   \item{chr_num}{Chromosome number (1-22)}
#'   \item{arm_type}{Arm type: "p" (short arm) or "q" (long arm)}
#'   \item{arm_start}{Start position of arm (1 for p-arms, breakpoint for q-arms)}
#'   \item{arm_end}{End position of arm (breakpoint for p-arms, chromosome end for q-arms)}
#'   \item{direction}{Primary alteration direction for this arm}
#' }
#' @details 
#' The arm boundaries are derived from BISCUT breakpoint analysis across 29 cancer types.
#' For p-arms, the boundary extends from chromosome start (position 1) to the breakpoint.
#' For q-arms, the boundary extends from the breakpoint to the chromosome end.
#' Chromosome end positions are based on hg19 reference genome.
#' 
#' @source Pan-cancer TCGA analysis with BISCUT algorithm
#' @seealso \code{\link{get_breakpoint_data}} for accessing this data programmatically
"consensus_arm_definitions"

#' BAGEL Color Palette
#'
#' A comprehensive color palette based on RColorBrewer palettes, designed for 
#' BAGEL visualization outputs. This 29-color palette combines three well-established
#' ColorBrewer palettes to provide distinct, colorblind-friendly colors suitable 
#' for chromosome ideograms, copy number plots, and other genomic visualizations.
#'
#' @format A character vector of 29 hexadecimal color codes:
#' \describe{
#'   \item{Colors 1-9}{Set1 palette: High-contrast qualitative colors (red, blue, green, purple, orange, yellow, brown, pink, grey)}
#'   \item{Colors 10-17}{Set2 palette: Pastel qualitative colors (light green, orange, blue, pink, lime, yellow, brown, grey)}
#'   \item{Colors 18-29}{Set3 palette: Extended qualitative colors (12 additional distinct hues)}
#' }
#' 
#' @details 
#' The palette combines three RColorBrewer qualitative palettes:
#' \itemize{
#'   \item \strong{Set1 (9 colors)}: High-contrast, saturated colors for primary distinctions
#'   \item \strong{Set2 (8 colors)}: Softer, pastel versions for secondary categories  
#'   \item \strong{Set3 (12 colors)}: Extended range for fine-grained categorization
#' }
#' 
#' Benefits of this approach:
#' \itemize{
#'   \item Scientifically validated color schemes from ColorBrewer
#'   \item Optimized for colorblind accessibility (Set1 and Set2)
#'   \item Print-friendly with good contrast ratios
#'   \item Widely recognized in scientific visualizations
#'   \item Sufficient colors for human chromosome visualization (22 + sex chromosomes)
#' }
#' 
#' @examples
#' # Use full palette
#' plot(1:29, col = bagel_palette, pch = 19, cex = 2)
#' 
#' # Use specific ColorBrewer sets
#' set1_colors <- bagel_palette[1:9]     # High contrast
#' set2_colors <- bagel_palette[10:17]   # Pastel
#' set3_colors <- bagel_palette[18:29]   # Extended
#' brewer.pal(12, "Paired"), 
#' brewer.pal(8, "Dark2")
#' 
#' # Display palette
#' library(RColorBrewer)
#' display.brewer.all(type = "qual")
#' 
#' @source Combination of RColorBrewer Set1, Set2, and Set3 qualitative palettes
#' @references Brewer, C. A. (2003). ColorBrewer 2.0. http://www.colorbrewer2.org
#' @seealso \code{\link{plot_chromosome_ideograms}} for usage in ideogram visualization
#' @importFrom RColorBrewer brewer.pal
"bagel_palette"

#' Updated TCGA Breakpoint Definitions (BAGEL v4 Analysis)
#'
#' Cancer-specific chromosome arm boundary definitions derived from BAGEL v4 analysis
#' of 29 TCGA cancer types. This dataset represents the most current breakpoint data,
#' replacing legacy TCGA consensus breakpoints with empirically-derived functional
#' arm definitions from comprehensive BISCUT + BAGEL analysis.
#'
#' @format A named list where each element corresponds to a cancer type.
#'   Each cancer type contains a data frame with standardized arm definition columns:
#' \describe{
#'   \item{peak_id}{Unique identifier for the breakpoint event (e.g., "1p_amp_tel_neg")}
#'   \item{chr}{Chromosome identifier with "chr" prefix (e.g., "chr1", "chr2")}
#'   \item{chr_arm}{Chromosome arm identifier (e.g., "1p", "1q", "2p")}
#'   \item{arm}{Arm type: "p" (short arm) or "q" (long arm)}
#'   \item{type}{Boundary type combining arm and location (e.g., "p_tel", "q_cent")}
#'   \item{direction}{Alteration direction: "amp" (amplification) or "del" (deletion)}
#'   \item{telcent}{Boundary location: "tel" (telomere) or "cent" (centromere)}
#'   \item{negpos}{Selection direction: "pos" (positive) or "neg" (negative)}
#'   \item{type_of_selection}{Functional classification: "TS-like", "onco-like", "toxic-like", "essential-like"}
#'   \item{n_events}{Number of events supporting this breakpoint}
#'   \item{n_left}{Left-side statistical support}
#'   \item{n_right}{Right-side statistical support}
#'   \item{Peak.Start}{Breakpoint peak start position (bp, hg19)}
#'   \item{Peak.End}{Breakpoint peak end position (bp, hg19)}
#'   \item{Peak.Position}{Peak center position (bp, hg19)}
#'   \item{peak_length}{Physical length of peak region (bp)}
#'   \item{Peak.Start.1}{Normalized peak start (0-1 scale relative to arm length)}
#'   \item{Peak.End.1}{Normalized peak end (0-1 scale relative to arm length)}
#'   \item{Peak.Position.1}{Normalized peak center position}
#'   \item{peak_length.pct}{Peak length as percentage of arm length}
#'   \item{ks_stat}{Kolmogorov-Smirnov test statistic}
#'   \item{ks_p}{K-S test p-value}
#'   \item{ksby}{K-S test adjusted by sample size}
#'   \item{combined_sig}{Combined significance score}
#'   \item{chromosome_info}{Chromosome number (1-22)}
#'   \item{size}{Total chromosome size (bp, hg19)}
#'   \item{p_start}{P arm start position (bp)}
#'   \item{p_end}{P arm end position (bp)}
#'   \item{q_start}{Q arm start position (bp)}
#'   \item{q_end}{Q arm end position (bp)}
#'   \item{centromere}{Centromere position (bp)}
#'   \item{offset}{Genomic coordinate offset for visualization}
#'   \item{middle_of_arm}{Midpoint of arm (bp)}
#'   \item{arm_length}{Total arm length (bp)}
#'   \item{functional_start}{Functional region start coordinate (bp, hg19)}
#'   \item{functional_end}{Functional region end coordinate (bp, hg19)}
#'   \item{functional_length}{Length of functional region (bp)}
#'   \item{percentage_length}{Functional region as percentage of arm length}
#'   \item{arm_start}{Arm boundary start (for compatibility)}
#'   \item{arm_end}{Arm boundary end (for compatibility)}
#' }
#'
#' @details
#' \strong{Data Generation Process:}
#' \enumerate{
#'   \item BISCUT breakpoint detection on TCGA segmentation data (29 cancer types)
#'   \item BAGEL v4 functional arm definition with direct coordinate mapping
#'   \item Statistical filtering (ks_p < 0.05, percentage_length > 0.1)
#'   \item Quality control and validation
#'   \item Consolidation into cancer-specific datasets
#' }
#'
#' \strong{Key Differences from Legacy Data:}
#' \itemize{
#'   \item Uses \code{functional_start} and \code{functional_end} (not peak_start/peak_end)
#'   \item Includes detailed selection type classifications
#'   \item Contains complete statistical support metrics
#'   \item Based on empirical BISCUT analysis (not consensus averaging)
#'   \item Each cancer type has unique, biologically-relevant breakpoints
#' }
#'
#' \strong{Cancer Types Included (29):}
#' \itemize{
#'   \item adrenocortical_cancer, bladder_urothelial_carcinoma
#'   \item brain_lower_grade_glioma, breast_invasive_carcinoma
#'   \item cervical_endocervical_cancer, cholangiocarcinoma
#'   \item colon_adenocarcinoma, diffuse_large_b_cell_lymphoma
#'   \item esophageal_carcinoma, glioblastoma_multiforme
#'   \item head_neck_squamous_cell_carcinoma, kidney_clear_cell_carcinoma
#'   \item kidney_papillary_cell_carcinoma, liver_hepatocellular_carcinoma
#'   \item lung_adenocarcinoma, lung_squamous_cell_carcinoma
#'   \item mesothelioma, ovarian_serous_cystadenocarcinoma
#'   \item pancreatic_adenocarcinoma, pheochromocytoma_paraganglioma
#'   \item prostate_adenocarcinoma, rectum_adenocarcinoma
#'   \item sarcoma, skin_cutaneous_melanoma
#'   \item stomach_adenocarcinoma, testicular_germ_cell_tumor
#'   \item thyroid_carcinoma, uterine_carcinosarcoma
#'   \item uterine_corpus_endometrioid_carcinoma
#' }
#'
#' @section Usage in define_arm():
#' This dataset is automatically loaded as the primary TCGA fallback when:
#' \itemize{
#'   \item Custom BISCUT file is not provided or fails to load
#'   \item \code{cancer_type} parameter matches one of the 29 cancer types
#'   \item The function \code{define_arm()} is called with \code{custom_biscut_file = NULL}
#' }
#'
#' The fallback hierarchy is:
#' \enumerate{
#'   \item Custom BISCUT file (if provided and valid)
#'   \item \strong{updated_tcga_breakpoints} (BAGEL v4, cancer-specific)
#'   \item consensus_arm_definitions (legacy pan-cancer consensus)
#' }
#'
#' @examples
#' # Access specific cancer type
#' breast_arms <- updated_tcga_breakpoints$breast_invasive_carcinoma
#'
#' # Get all cancer types available
#' available_cancers <- names(updated_tcga_breakpoints)
#'
#' # Check number of arms defined per cancer
#' sapply(updated_tcga_breakpoints, nrow)
#'
#' # Use in define_arm() function
#' arm_defs <- define_arm(
#'   custom_biscut_file = NULL,
#'   cancer_type = "lung_adenocarcinoma"
#' )
#'
#' @source BAGEL v4 analysis of TCGA pan-cancer segmentation data
#' @references
#' \itemize{
#'   \item BISCUT: Beroukhim et al. (2010) Nature 463, 899-905
#'   \item BAGEL v4: Direct coordinate mapping with functional arm definitions
#'   \item TCGA: The Cancer Genome Atlas Research Network
#' }
#' @seealso
#' \itemize{
#'   \item \code{\link{define_arm}} for using this data in arm definition
#'   \item \code{\link{cancer_specific_breakpoints}} for legacy breakpoint data
#'   \item \code{\link{consensus_arm_definitions}} for pan-cancer consensus
#' }
"updated_tcga_breakpoints"
