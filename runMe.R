# devtools::install_github("PollyHung/BAGEL")
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(GenomicRanges)
library(ggbio)
library(ape)
library(readxl)
library(regioneR)
library(RIdeogram)


## Step 1: Define File Paths ---------------------------------------------------
base_dir <- "/Users/polly_hung/Desktop/BAGEL/example/qmh/"
setwd(base_dir)
cuts.path <- "all_BISCUT_results.txt"
ref.path <- "/Users/polly_hung/Desktop/BAGEL/docs/breakpoints/ovarian.txt"
seg.path <- "QMH_TS_hg38_403.seg"
bed.file <- "CPOS_812_HPV_EBV_targets_hg38_sorted_merged_mod.bed"



## Step 2: Load in Libraries and Functions -------------------------------------
setwd(base_dir)
functions <- list.files(path = file.path(base_dir, "R/"), pattern = ".R", full.names = TRUE)
lapply(functions, function(x){source(x)})


## Step 2: Preprocess ----------------------------------------------------------
segments <- mergeSegments(seg.path = seg.path, min_probes = 4)

breakpoints <- processCuts(cuts.path = "all_BISCUT_results.txt",
                           cutoff = 0.98,
                           model = "human",
                           genome = "hg38",
                           ref.path = ref.path,
                           ref2.type = "pancancer",
                           stringent = FALSE)

markers <- makeMarkers(bed_file = bed.file, technique = "TS")

checkCoverage(GRange = markers, cuts = breakpoints$amp, )
checkCoverage(GRange = markers, cuts = breakpoints$del, result_folder = "/Users/polly_hung/Desktop/BESAA/results/")

