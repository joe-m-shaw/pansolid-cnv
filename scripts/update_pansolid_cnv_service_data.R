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

worksheet_labno_regex <- "(WS\\d{6})_(\\d{6,8})_"

panel_regex <-".+WorksheetAnalysedData/WS\\d{6}/(\\w{1,30})/(Ann.+|Genotyped/Ann.+)"

s_drive_file_df <- tibble(
  filepath = unlist(s_drive_filepaths)) |> 
  mutate(filename = str_extract(string = filepath, 
                         pattern = str_replace(string = pansolidv2_excel_regex, 
                                               pattern = "\\^", 
                                               replacement = "")),
         worksheet = str_extract(string = filename, 
                                 pattern = worksheet_labno_regex,
                                 group = 1),
         labno = str_extract(string = filename, 
                             pattern = worksheet_labno_regex,
                             group = 2),
         panel = str_extract(string = filepath,
                             pattern = panel_regex,
                             group = 1))

if (any(grepl(pattern = "/", x = s_drive_file_df$filename))) {
  stop("Error: filenames contain backslashes") 
}

if (anyNA(s_drive_file_df)) {
  stop("Error: there are NA values in the filepath table")
}

write.csv(s_drive_file_df, 
          here::here(str_c("data/live_service_collated_data/",
                           "pansolidv2_sample_worksheet_panel_information.csv")),
          row.names = FALSE)

# Local drive filepaths -------------------------------------------------------------

local_drive_file_df <- tibble(
  filepath = unlist(list.files(here::here("data/live_service_annotated_files/"),
                               full.names = TRUE))) |> 
  mutate(filename = str_extract(string = filepath, 
                                pattern = str_replace(string = pansolidv2_excel_regex, 
                                                      pattern = "\\^", 
                                                      replacement = "")),
         labno = str_extract(string = filename, 
                             pattern = "WS\\d{6}_(\\d{6,8})_",
                             group = 1))

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
                                pattern = str_replace(string = pansolidv2_excel_regex, 
                                                      pattern = "\\^", 
                                                      replacement = "")),
         labno = str_extract(string = filename, 
                             pattern = "WS\\d{6}_(\\d{6,8})_",
                             group = 1))

new_file_local_paths_df <- local_drive_file_df |> 
  filter(filename %in% new_files$filename & 
           # Remove samples without "Amplifications" tab
           !labno %in% c("24023280", "24025207", "24027566", "24033006",
                         "24033959"))

new_file_local_paths <- list(new_file_local_paths_df$filepath) |> 
  flatten()

# Collate new file data -------------------------------------------------------------

new_amp_gene_collated <- new_file_local_paths |> 
  map(\(new_file_local_paths) read_all_amp_genes_results(file = new_file_local_paths,
                                                         sheet = get_amp_sheetname(new_file_local_paths))) |> 
  list_rbind()

new_pos_cnv_collated <- new_file_local_paths |> 
  map(\(new_file_local_paths) read_pos_cnv_results(file = new_file_local_paths,
                                                   sheet = get_amp_sheetname(new_file_local_paths))) |> 
  list_rbind()

new_std_dev_collated <- new_file_local_paths |> 
  map(\(new_file_local_paths) read_stdev_results(file = new_file_local_paths,
                                                        sheet = get_amp_sheetname(new_file_local_paths))) |> 
  list_rbind()

new_percent_138_collated <- new_file_local_paths |> 
  map(\(new_file_local_paths) read_percent_138_results(file = new_file_local_paths,
                                                        sheet = get_amp_sheetname(new_file_local_paths))) |> 
  list_rbind()

# Load previously collated data -----------------------------------------------------

amp_gene_results <- read_csv(here::here("data/live_service_collated_data/live_service_amp_gene_results_collated.csv"))

std_dev_results <- read_csv(here::here("data/live_service_collated_data/live_service_std_dev_results_collated.csv"))

percent_138_results <- read_csv(here::here("data/live_service_collated_data/live_service_percent_138_results_collated.csv"))

pos_cnv_results <- read_csv(here::here("data/live_service_collated_data/live_service_pos_cnv_results_collated.csv"))

# Check columns ---------------------------------------------------------------------

# Scientists have been adding comments as extra columns, which makes binding the data
# tricky.

id_cols <- c("worksheet", "labno", "suffix", "patient_name", 
             "labno_suffix", "labno_suffix_worksheet", "filepath")

all_amp_cols <- c(id_cols, "gene", "max_region_fold_change", "min_region_fold_change")

std_dev_cols <- c(id_cols, "st_dev_signal_adjusted_log2_ratios")

pos_cnv_cols <- c(id_cols, "gene",	"chromosome",	"cnv_co_ordinates",	"cnv_length",
                  "consequence", "fold_change",	"p_value",	
                  "no_targets",	"start",	"end")

percent_138_cols <- c(id_cols, "percent_whole_panel_covered_at_138x")

# Add new data to collated data -----------------------------------------------------

amp_gene_results_updated <- rbind(amp_gene_results |> 
                                    select(all_of(all_amp_cols)),
                                  new_amp_gene_collated |> 
                                    select(all_of(all_amp_cols)))

pos_cnv_results_updated <- rbind(pos_cnv_results |> 
                                   select(all_of(pos_cnv_cols)),
                                 new_pos_cnv_collated |> 
                                   select(all_of(pos_cnv_cols)))

std_dev_results_updated <- rbind(std_dev_results |> 
                                   select(all_of(std_dev_cols)),
                                 new_std_dev_collated |> 
                                   select(all_of(std_dev_cols)))

percent_138_results_updated <- rbind(percent_138_results |> 
                               select(all_of(percent_138_cols)),
                             new_percent_138_collated |> 
                               select(all_of(percent_138_cols)))

# Archive previous collated data ----------------------------------------------------

write.csv(amp_gene_results,
          here::here(str_c("data/live_service_collated_data/archived_collated_data/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_amp_gene_results_collated.csv")),
          row.names = FALSE)

write.csv(std_dev_results,
          here::here(str_c("data/live_service_collated_data/archived_collated_data/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_std_dev_results_collated.csv")),
          row.names = FALSE)

write.csv(pos_cnv_results,
          here::here(str_c("data/live_service_collated_data/archived_collated_data/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_pos_cnv_results_collated.csv")),
          row.names = FALSE)

write.csv(percent_138_results,
          here::here(str_c("data/live_service_collated_data/archived_collated_data/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_percent_138_results_collated.csv")),
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

write.csv(pos_cnv_results_updated,
          here::here(str_c("data/live_service_collated_data/",
                           "live_service_pos_cnv_results_collated.csv")),
          row.names = FALSE)

write.csv(percent_138_results_updated,
          here::here(str_c("data/live_service_collated_data/",
                           "live_service_percent_138_results_collated.csv")),
          row.names = FALSE)
