# Create Pan-Cancer Breakpoint Data for BAGEL 2.0
# This script consolidates all BISCUT results into built-in package data

library(dplyr)
library(readr)

# List all cancer types and their BISCUT results
cancer_types <- c(
  "adrenocortical_cancer", "bladder_urothelial_carcinoma", "brain_lower_grade_glioma",
  "breast_invasive_carcinoma", "cervical_endocervical_cancer", "cholangiocarcinoma",
  "colon_adenocarcinoma", "diffuse_large_b_cell_lymphoma", "esophageal_carcinoma",
  "glioblastoma_multiforme", "head_neck_squamous_cell_carcinoma", "kidney_clear_cell_carcinoma",
  "kidney_papillary_cell_carcinoma", "liver_hepatocellular_carcinoma", "lung_adenocarcinoma",
  "lung_squamous_cell_carcinoma", "mesothelioma", "ovarian_serous_cystadenocarcinoma",
  "pancreatic_adenocarcinoma", "pheochromocytoma_paraganglioma", "prostate_adenocarcinoma",
  "rectum_adenocarcinoma", "sarcoma", "skin_cutaneous_melanoma", "stomach_adenocarcinoma",
  "testicular_germ_cell_tumor", "thyroid_carcinoma", "uterine_carcinosarcoma",
  "uterine_corpus_endometrioid_carcinoma"
)

# Function to read and process a single cancer type's breakpoints
read_cancer_breakpoints <- function(cancer_type) {
  file_path <- file.path("/Users/polly_hung/Desktop/BAGEL/results", 
                        cancer_type, 
                        paste0("results_", cancer_type),
                        "all_BISCUT_results.txt")
  
  if (!file.exists(file_path)) {
    warning(sprintf("File not found: %s", file_path))
    return(NULL)
  }
  
  # Read the BISCUT results
  results <- read_delim(file_path, delim = "\t", show_col_types = FALSE)
  
  # Add cancer type annotation
  results$cancer_type <- cancer_type
  
  return(results)
}

# Read all cancer types
cat("Reading BISCUT results for all cancer types...\n")
all_breakpoints <- do.call(rbind, lapply(cancer_types, read_cancer_breakpoints))

# Create arm-level breakpoint summary
arm_breakpoints <- all_breakpoints %>%
  group_by(arm, direction, cancer_type) %>%
  summarise(
    chr = unique(Chr)[1],
    n_peaks = n(),
    peak_start = min(Peak.Start),
    peak_end = max(Peak.End),
    .groups = "drop"
  )

# Create pan-cancer consensus breakpoints (used across multiple cancer types)
consensus_breakpoints <- arm_breakpoints %>%
  group_by(arm, direction) %>%
  summarise(
    chr = unique(chr)[1],
    n_cancer_types = n(),
    median_start = median(peak_start),
    median_end = median(peak_end),
    cancer_types = paste(cancer_type, collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(n_cancer_types >= 3) # Require at least 3 cancer types for consensus

# Create cancer-specific breakpoint data
cancer_specific_breakpoints <- split(arm_breakpoints, arm_breakpoints$cancer_type)

# Create standardized arm definitions based on breakpoints
create_arm_definitions <- function(breakpoints_df) {
  breakpoints_df %>%
    mutate(
      chr_num = as.numeric(chr),
      arm_type = ifelse(grepl("p$", arm), "p", "q"),
      arm_start = case_when(
        arm_type == "p" ~ 1,
        arm_type == "q" ~ median_start
      ),
      arm_end = case_when(
        arm_type == "p" ~ median_end,
        arm_type == "q" ~ case_when(
          chr_num == 1 ~ 249250621,
          chr_num == 2 ~ 243199373,
          chr_num == 3 ~ 198022430,
          chr_num == 4 ~ 191154276,
          chr_num == 5 ~ 180915260,
          chr_num == 6 ~ 171115067,
          chr_num == 7 ~ 159138663,
          chr_num == 8 ~ 146364022,
          chr_num == 9 ~ 141213431,
          chr_num == 10 ~ 135534747,
          chr_num == 11 ~ 135006516,
          chr_num == 12 ~ 133851895,
          chr_num == 13 ~ 115169878,
          chr_num == 14 ~ 107349540,
          chr_num == 15 ~ 102531392,
          chr_num == 16 ~ 90354753,
          chr_num == 17 ~ 81195210,
          chr_num == 18 ~ 78077248,
          chr_num == 19 ~ 59128983,
          chr_num == 20 ~ 63025520,
          chr_num == 21 ~ 48129895,
          chr_num == 22 ~ 51304566,
          TRUE ~ NA_real_
        )
      )
    ) %>%
    select(arm, chr_num, arm_type, arm_start, arm_end, direction) %>%
    arrange(chr_num, arm_type)
}

# Create consensus arm definitions
consensus_arm_definitions <- create_arm_definitions(consensus_breakpoints)

# Save the data objects that will be included in the package
usethis::use_data(all_breakpoints, overwrite = TRUE)
usethis::use_data(arm_breakpoints, overwrite = TRUE)  
usethis::use_data(consensus_breakpoints, overwrite = TRUE)
usethis::use_data(cancer_specific_breakpoints, overwrite = TRUE)
usethis::use_data(consensus_arm_definitions, overwrite = TRUE)

cat("Created package data objects:\n")
cat("- all_breakpoints: Complete BISCUT results from all cancer types\n")
cat("- arm_breakpoints: Summary by cancer type and arm\n") 
cat("- consensus_breakpoints: Pan-cancer consensus breakpoints\n")
cat("- cancer_specific_breakpoints: List of cancer-specific breakpoints\n")
cat("- consensus_arm_definitions: Standardized arm definitions\n")