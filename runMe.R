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
# library(BAGEL)
library(data.table)


## Step 1: Define File Paths ---------------------------------------------------
base_dir <- "/Users/polly_hung/Desktop/BAGEL/example/pooledOV/"
setwd(base_dir)
seg.path <- "pooled_2787.seg"


## Step 2: Preprocess the Segments ---------------------------------------------
segs <- preprocessSeg(seg.path = seg.path,
                      min_probes = 4,
                      genome = "hg38")

## Step 3: Create Cuts ---------------------------------------------------------



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
raw_copyNumber <- calculateCopyNumber(segments = seg2,
                                      breakpoints = breakpoints,
                                      amp_thres = 0.25,
                                      del_thres = -0.25,
                                      runCentromere = TRUE)


## Step 4: Adjust Copy Number based on ...--------------------------------------
gaps <- segmentGaps(segments = seg2, genome = "hg38")
coverage <- percentCoverage(segments = seg2,
                            gaps = gaps,
                            amp_thres = 0.1,
                            del_thres = 0.1)

