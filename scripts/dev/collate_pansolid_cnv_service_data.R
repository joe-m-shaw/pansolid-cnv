# Collate PanSolidv2 Copy Number Variant Data

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

# Functions -------------------------------------------------------------------------

source(here::here("functions/cnv_functions.R"))

# S drive filepaths -----------------------------------------------------------------

pansolidv2_worksheets <- read_excel(here::here("data/pansolid_live_service_worksheets.xlsx"))

worksheet_list <- list(pansolidv2_worksheets$worksheet)

s_drive_filepaths <- worksheet_list |> 
  map(\(worksheet_list) get_annotated_filepaths(worksheet_list)) |> 
  flatten()

s_drive_file_df <- tibble(
  filepath = unlist(s_drive_filepaths)) |> 
  mutate(filename = str_extract(string = filepath, 
                                pattern = "Annotated_WS\\d{6}_.+.xlsx"))

# Copy to local drive ---------------------------------------------------------------

file.copy(from = s_drive_filepaths, 
          to = here::here("data/live_service_annotated_files/"))

# Local drive filepaths -------------------------------------------------------------

local_drive_file_df <- tibble(
  filepath = unlist(list.files(here::here("data/live_service_annotated_files/"),
                               full.names = TRUE))) |> 
  mutate(filename = str_extract(string = filepath, 
                                pattern = "Annotated_WS\\d{6}_.+.xlsx"),
         labno = str_extract(string = filename, 
                             pattern = "\\d{8}")) |> 
  filter(!labno %in% c("24023280", "24025207"))

local_filepaths <- list(local_drive_file_df$filepath) |> 
  flatten()

# Collate local data ----------------------------------------------------------------

amp_gene_collated <- local_filepaths |> 
  map(\(local_filepaths) read_all_amp_genes_results(file = local_filepaths, 
                                            sheet = get_amp_sheetname(local_filepaths))) |> 
  list_rbind()

std_dev_collated <- local_filepaths |> 
  map(\(local_filepaths) read_stdev_results(local_filepaths, 
                                            sheet = get_amp_sheetname(local_filepaths))) |> 
  list_rbind()

pos_cnv_collated <- local_filepaths |> 
  map(\(local_filepaths) read_pos_cnv_results(local_filepaths, 
                                              sheet = get_amp_sheetname(local_filepaths))) |> 
  list_rbind()

percent_138_collated <- local_filepaths |> 
  map(\(local_filepaths) read_percent_138_results(local_filepaths, 
                                                  sheet = get_amp_sheetname(local_filepaths))) |> 
  list_rbind()

# Save collated data ----------------------------------------------------------------

write.csv(x = amp_gene_collated, 
          file = here::here("data/live_service_collated_data/live_service_amp_gene_results_collated.csv"),
          row.names = FALSE)

write.csv(std_dev_collated, here::here("data/live_service_collated_data/live_service_std_dev_results_collated.csv"),
          row.names = FALSE)

write.csv(pos_cnv_collated, here::here("data/live_service_collated_data/live_service_pos_cnv_results_collated.csv"),
          row.names = FALSE)

write.csv(percent_138_collated, here::here("data/live_service_collated_data/live_service_percent_138_results_collated.csv"),
          row.names = FALSE)
