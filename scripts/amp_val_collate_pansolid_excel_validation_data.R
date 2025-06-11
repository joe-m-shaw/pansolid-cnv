# Load and Collate Data for PanSolid Gene Amplifications Cohort

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(here)

# Functions -------------------------------------------------------------------------

source(here("functions/pansolid_cnv_excel_functions.R"))

# Files -----------------------------------------------------------------------------

amp_cohort_filepaths <- list.files(path = paste0(config::get("data_folderpath"),
                                                 "validation/",
                                                 "DOC6283_amplifications/",
                                                 "raw/",
                                                 "pansolid_ngs_amplifications/"), 
                               full.names = TRUE,
                               recursive = TRUE,
                               pattern = "Annotated_.+.xlsx")

if(length(amp_cohort_filepaths) == 0) {
  stop("No files in location")
}

# Collate results -------------------------------------------------------------------

pos_cnv_results_collated <- amp_cohort_filepaths |> 
  map(\(amp_cohort_filepaths) read_pos_cnv_results(file = amp_cohort_filepaths,
                                                   sheet = get_amp_sheetname(
                                                     amp_cohort_filepaths))) |> 
  list_rbind() |> 
  # The limit of detection normal control was listed as 24039973 (CNVMix12CopiesSERASEQ
  # at 0%) but can be changed to 24039975 (DNAWTmixSERASEQ at 100%)
  mutate(
    labno = case_when(
    
      labno_suffix_worksheet == "24039973d_WS144291" ~"24039975",
      TRUE ~labno),
    
    patient_name = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"DNAWTmixSERASEQ",
      TRUE ~patient_name),
    
    labno_suffix = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"24039975d",
      TRUE ~labno_suffix),
    
    labno_suffix_worksheet = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"24039975d_WS144291",
      TRUE ~labno_suffix_worksheet)
    )

fold_change_threshold <- 2.8

all_amp_gene_results_collated <- amp_cohort_filepaths |> 
  map(\(amp_cohort_filepaths) read_all_amp_genes_results(file = amp_cohort_filepaths,
                                                         sheet = get_amp_sheetname(
                                                           amp_cohort_filepaths))) |> 
  list_rbind() |> 
  mutate(pansolid_call = case_when(
    max_region_fold_change >= fold_change_threshold ~"amplification",
    max_region_fold_change < fold_change_threshold ~"normal result"),
    labno = case_when(
      
      labno_suffix_worksheet == "24039973d_WS144291" ~"24039975",
      TRUE ~labno),
    
    patient_name = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"DNAWTmixSERASEQ",
      TRUE ~patient_name),
    
    labno_suffix = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"24039975d",
      TRUE ~labno_suffix),
    
    labno_suffix_worksheet = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"24039975d_WS144291",
      TRUE ~labno_suffix_worksheet)
  )

stdev_results_collated <-  amp_cohort_filepaths |> 
  map(\(amp_cohort_filepaths) read_stdev_results(file = amp_cohort_filepaths,
                                                 sheet = get_amp_sheetname(
                                                   amp_cohort_filepaths))) |> 
  list_rbind() |> 
  mutate(
    labno = case_when(
      
      labno_suffix_worksheet == "24039973d_WS144291" ~"24039975",
      TRUE ~labno),
    
    patient_name = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"DNAWTmixSERASEQ",
      TRUE ~patient_name),
    
    labno_suffix = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"24039975d",
      TRUE ~labno_suffix),
    
    labno_suffix_worksheet = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"24039975d_WS144291",
      TRUE ~labno_suffix_worksheet)
  )

percent_138_collated <- amp_cohort_filepaths |> 
  map(\(amp_cohort_filepaths) 
      read_percent_138_results(file = amp_cohort_filepaths,
                               sheet = get_amp_sheetname(amp_cohort_filepaths))) |> 
  list_rbind() |> 
  mutate(
    labno = case_when(
      
      labno_suffix_worksheet == "24039973d_WS144291" ~"24039975",
      TRUE ~labno),
    
    patient_name = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"DNAWTmixSERASEQ",
      TRUE ~patient_name),
    
    labno_suffix = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"24039975d",
      TRUE ~labno_suffix),
    
    labno_suffix_worksheet = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"24039975d_WS144291",
      TRUE ~labno_suffix_worksheet)
  )

# Check all results collated --------------------------------------------------------

if(
  
  unique(unique(pos_cnv_results_collated$filepath) == amp_cohort_filepaths) != TRUE |
  
  unique(unique(all_amp_gene_results_collated$filepath) == amp_cohort_filepaths) != TRUE |
  
  unique(stdev_results_collated$filepath == amp_cohort_filepaths) != TRUE |
  
  unique(percent_138_collated$filepath == amp_cohort_filepaths) != TRUE
  
) {
  
  stop("Not all filepaths have been read in") 
  
} else {
  
  message("All PanSolid files successfully read in")
  
} 

# Save collated results -------------------------------------------------------------

processed_validation_data_folder <- paste0(config::get("data_folderpath"), 
                                          "validation/DOC6283_amplifications/",
                                          "processed/")

write.csv(x = pos_cnv_results_collated,
          file = paste0(processed_validation_data_folder, 
                        "amp_validation_pos_cnv_results_collated.csv"),
          row.names = FALSE)

write.csv(x = all_amp_gene_results_collated,
          file = paste0(processed_validation_data_folder, 
                        "amp_validation_all_amp_gene_results_collated.csv"),
          row.names = FALSE)

write.csv(x = stdev_results_collated,
          file = paste0(processed_validation_data_folder, 
                        "amp_validation_stdev_results_collated.csv"),
          row.names = FALSE)

write.csv(x = percent_138_collated,
          file = paste0(processed_validation_data_folder, 
                        "amp_validation_percent_138_collated.csv"),
          row.names = FALSE)
