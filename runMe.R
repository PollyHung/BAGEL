# devtools::install_github("PollyHung/BAGEL")
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(GenomicRanges)
library(ape)
library(readxl)
library(regioneR)
library(RIdeogram)
library(BAGEL)
library(data.table)


## Step 1: Define File Paths ---------------------------------------------------
base_dir <- "/Users/polly_hung/Desktop/BAGEL/example/qmh/"
setwd(base_dir)
cuts.path <- "all_BISCUT_results.txt"
ref.path <- "/Users/polly_hung/Desktop/BAGEL/docs/breakpoints/ovarian.txt"
seg.path <- "QMH_TS_hg38_403.seg"
bed.file <- "CPOS_812_HPV_EBV_targets_hg38_sorted_merged_mod.bed"


## Step 2: Process the Breakpoints ---------------------------------------------
breakpoints <- processCuts(cuts.path = "all_BISCUT_results.txt",
                           cutoff = 0.98,
                           model = "human",
                           genome = "hg38",
                           ref.path = ref.path,
                           ref2.type = "pancancer",
                           stringent = FALSE)

markers <- makeMarkers(bed_file = bed.file, technique = "TS")
checkCoverage(GRange = markers, cuts = breakpoints$amp, plot_name = "amp")
checkCoverage(GRange = markers, cuts = breakpoints$del, plot_name = "del")


## Step 3: Segmentation to Copy Number -----------------------------------------
segments <- mergeSegments(seg.path = seg.path, min_probes = 4)
seg2 <- annotateArms(segments = segments, genome = "hg38")
gaps <- segmentGaps(segments = seg2, genome = "hg38")
coverage <- percentCoverage(segments = seg2, gaps = gaps, thres = 0.1)
raw_copyNumber <- calculateCopyNumber(segments = seg2,
                                      breakpoints = breakpoints,
                                      amp_thres = 0.1,
                                      del_thres = -0.1,
                                      runCentromere = TRUE)




