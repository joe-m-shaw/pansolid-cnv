# Update PanSolidv2 Copy Number Variant Data

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
                                pattern = "Annotated_WS\\d{6}_.+.xlsx"),
         labno = str_extract(string = filename, 
                             pattern = "\\d{8}"))

# Local drive filepaths -------------------------------------------------------------

local_drive_file_df <- tibble(
  filepath = unlist(list.files(here::here("data/live_service_annotated_files/"),
                               full.names = TRUE))) |> 
  mutate(filename = str_extract(string = filepath, 
                                pattern = "Annotated_WS\\d{6}_.+.xlsx"),
         labno = str_extract(string = filename, 
                             pattern = "\\d{8}"))

# Identify and copy new files -------------------------------------------------------

new_files <- s_drive_file_df |> 
  filter(!filename %in% local_drive_file_df$filename &
           labno != "24023280")

if (nrow(new_files) > 0) {
  
  file.copy(from = new_files$filepath, 
            to = here::here("data/live_service_annotated_files/"))
  
}

# Get new file local filepaths ------------------------------------------------------

local_drive_file_df <- tibble(
  filepath = unlist(list.files(here::here("data/live_service_annotated_files/"),
                               full.names = TRUE))) |> 
  mutate(filename = str_extract(string = filepath, 
                                pattern = "Annotated_WS\\d{6}_.+.xlsx"),
         labno = str_extract(string = filename, 
                             pattern = "\\d{8}"))

new_file_local_paths_df <- local_drive_file_df |> 
  filter(filename %in% new_files$filename & 
           # Remove samples without "Amplifications" tab
           !labno %in% c("24023280", "24025207"))

new_file_local_paths <- list(new_file_local_paths_df$filepath) |> 
  flatten()

# Collate new file data -------------------------------------------------------------

new_amp_gene_collated <- new_file_local_paths |> 
  map(\(new_file_local_paths) read_annotated_file_all_amp(new_file_local_paths)) |> 
  list_rbind()

new_std_dev_collated <- new_file_local_paths |> 
  map(\(new_file_local_paths) read_annotated_file_stdev(new_file_local_paths)) |> 
  list_rbind()

# Load previously collated data -----------------------------------------------------

amp_gene_results <- read_csv(here::here("data/live_service_collated_data/live_service_amp_gene_results_collated.csv"))

std_dev_results <- read_csv(here::here("data/live_service_collated_data/live_service_std_dev_results_collated.csv"))

# Check columns ---------------------------------------------------------------------

# Scientists have been adding comments as extra columns, which makes binding the data
# tricky.

all_amp_cols <- c("worksheet", "labno", "suffix", "patient_name", 
                 "labno_suffix", "labno_suffix_worksheet", "filepath", 
                 "gene", "max_region_fold_change", "min_region_fold_change")

std_dev_cols <- c("worksheet", "labno", "suffix", "patient_name", 
                  "labno_suffix", "labno_suffix_worksheet", "filepath", 
                  "st_dev_signal_adjusted_log2_ratios")

amp_gene_results_cols_defined <- amp_gene_results |> 
  select(all_of(all_amp_cols))

std_dev_results_cols_defined <- std_dev_results |> 
  select(all_of(std_dev_cols))

new_amp_gene_collated_cols_defined <- new_amp_gene_collated |> 
  select(all_of(all_amp_cols))

new_std_dev_collated_cols_defined <- new_std_dev_collated |> 
  select(all_of(std_dev_cols))

stopifnot(colnames(amp_gene_results_cols_defined) == 
            colnames(new_amp_gene_collated_cols_defined))

stopifnot(colnames(std_dev_results_cols_defined) == 
            colnames(new_std_dev_collated_cols_defined))

# Add new data to collated data -----------------------------------------------------

amp_gene_results_updated <- rbind(amp_gene_results_cols_defined,
                                  new_amp_gene_collated_cols_defined)

std_dev_results_updated <- rbind(std_dev_results_cols_defined,
                                 new_std_dev_collated_cols_defined)

# Archive previous collated data ----------------------------------------------------

write.csv(amp_gene_results_cols_defined,
          here::here(str_c("data/live_service_collated_data/archived_collated_data/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_amp_gene_results_collated.csv")),
          row.names = FALSE)

write.csv(std_dev_results_cols_defined,
          here::here(str_c("data/live_service_collated_data/archived_collated_data/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_std_dev_results_collated.csv")),
          row.names = FALSE)

# Save updated collated data --------------------------------------------------------

write.csv(amp_gene_results_updated,
          here::here(str_c("data/live_service_collated_data/",
                           "live_service_amp_gene_results_collated.csv")),
          row.names = FALSE)

write.csv(std_dev_results_updated,
          here::here(str_c("data/live_service_collated_data/",
                           "live_service_std_dev_results_collated.csv")),
          row.names = FALSE)
