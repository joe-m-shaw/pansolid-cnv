# Load and Collate Data for PanSolid Gene Amplifications Cohort

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(here)

# Functions -------------------------------------------------------------------------

source(here("functions/pansolid_excel_functions.R"))

source(here("scripts/set_shared_drive_filepath.R"))

# Files -----------------------------------------------------------------------------

amp_cohort_data_location <- "S:/central shared/Genetics/NGS/Bioinformatics/1_Pan-solid-Cancer/CNV/Amplifications_final_with_MYC/"

amp_cohort_filepaths <- list.files(path = amp_cohort_data_location, 
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
  list_rbind()

fold_change_threshold <- 2.8

all_amp_gene_results_collated <- amp_cohort_filepaths |> 
  map(\(amp_cohort_filepaths) read_all_amp_genes_results(file = amp_cohort_filepaths,
                                                         sheet = get_amp_sheetname(
                                                           amp_cohort_filepaths))) |> 
  list_rbind() |> 
  mutate(pansolid_call = case_when(
    max_region_fold_change >= fold_change_threshold ~"amplification",
    max_region_fold_change < fold_change_threshold ~"normal result"))

stdev_results_collated <-  amp_cohort_filepaths |> 
  map(\(amp_cohort_filepaths) read_stdev_results(file = amp_cohort_filepaths,
                                                 sheet = get_amp_sheetname(
                                                   amp_cohort_filepaths))) |> 
  list_rbind() 

percent_138_collated <- amp_cohort_filepaths |> 
  map(\(amp_cohort_filepaths) 
      read_percent_138_results(file = amp_cohort_filepaths,
                               sheet = get_amp_sheetname(amp_cohort_filepaths))) |> 
  list_rbind()


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

processed_validation_data_folder <- paste0(data_folder, 
                                          "validation/processed/")

write.csv(x = pos_cnv_results_collated,
          file = paste0(processed_validation_data_folder, 
                        "validation_pos_cnv_results_collated.csv"),
          row.names = FALSE)

write.csv(x = all_amp_gene_results_collated,
          file = paste0(processed_validation_data_folder, 
                        "validation_all_amp_gene_results_collated.csv"),
          row.names = FALSE)

write.csv(x = stdev_results_collated,
          file = paste0(processed_validation_data_folder, 
                        "validation_stdev_results_collated.csv"),
          row.names = FALSE)

write.csv(x = percent_138_collated,
          file = paste0(processed_validation_data_folder, 
                        "validation_percent_138_collated.csv"),
          row.names = FALSE)
