# Update PanSolidv2 Copy Number Variant Data

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

# Functions and filepaths -----------------------------------------------------------

data_folder <- config::get("data_filepath")

source(here("functions/pansolid_excel_functions.R"))

# S drive filepaths -----------------------------------------------------------------

pansolidv2_worksheets <- read_excel(paste0(data_folder,
                                           "live_service/pansolid_live_service_worksheets.xlsx"))

message("PanSolid worksheet list read")

worksheet_list <- list(pansolidv2_worksheets$worksheet)

s_drive_filepaths <- worksheet_list |> 
  map(\(worksheet_list) get_annotated_filepaths(worksheet = worksheet_list)) |> 
  flatten()

worksheet_labno_regex <- "(WS\\d{6})_(\\d{6,8})(|a|b|c|d)_"

panel_regex <-".+WorksheetAnalysedData/WS\\d{6}/(\\w{1,30})/.+"

pansolidv2_excel_regex <- "^Annotated(_|_v2.+_)WS\\d{6}_.+.xlsx"

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
  warning("Error: there are NA values in the filepath table")
}

message("Sample filepaths compiled")

write_csv(s_drive_file_df, 
          paste0(data_folder, "live_service/collated/",
                 "pansolidv2_sample_worksheet_panel_information.csv"))

# Single folder filepaths -----------------------------------------------------------

single_folder_file_df <- tibble(
  filepath = unlist(list.files(paste0(data_folder, "live_service/raw/"),
                               full.names = TRUE))) |> 
  mutate(filename = str_extract(string = filepath, 
                                pattern = str_replace(string = pansolidv2_excel_regex, 
                                                      pattern = "\\^", 
                                                      replacement = "")),
         labno = str_extract(string = filename, 
                             pattern = "WS\\d{6}_(\\d{6,8})",
                             group = 1))


if(anyNA(single_folder_file_df$filepath)){
  warning("NA values in filepath column")
} else {
  message("Filepath column has no NA values")
}

if(anyNA(single_folder_file_df$filename)){
  warning("NA values in filename column")
} else {
  message("Filename column has no NA values")
}

if(anyNA(single_folder_file_df$labno)){
  warning("NA values in labno column")
} else {
  message("Labno column has no NA values")
}

# Identify and copy new files -------------------------------------------------------

new_files <- s_drive_file_df |> 
  filter(!filename %in% single_folder_file_df$filename)

if (nrow(new_files) > 0) {
  
  file.copy(from = new_files$filepath, 
            to = paste0(data_folder, "live_service/raw/"))
  
  message(paste0("Copying ", nrow(new_files), " new files to raw_data folder"))
  
}

# Get new single folder filepaths ---------------------------------------------------

single_folder_file_df <- tibble(
  filepath = unlist(list.files(paste0(data_folder, "live_service/raw/"),
                               full.names = TRUE))) |> 
  mutate(filename = str_extract(string = filepath, 
                                pattern = str_replace(string = pansolidv2_excel_regex, 
                                                      pattern = "\\^", 
                                                      replacement = "")),
         labno = str_extract(string = filename, 
                             pattern = "WS\\d{6}_(\\d{6,8})",
                             group = 1))

samples_without_amp_tabs <- c("24023280", "24025207", "24027566", "24033006",
                              "24033959", "24038848")

new_file_paths_df <- single_folder_file_df |> 
  filter(filename %in% new_files$filename & 
           # Remove samples without "Amplifications" tab
           !labno %in% samples_without_amp_tabs)

new_file_paths <- list(new_file_paths_df$filepath) |> 
  flatten()

# Collate new file data -------------------------------------------------------------

new_amp_gene_collated <- new_file_paths |> 
  map(\(new_file_paths) 
      read_all_amp_genes_results(file = new_file_paths,
                                 sheet = get_amp_sheetname(new_file_paths))) |> 
  list_rbind()

new_pos_cnv_collated <- new_file_paths |> 
  map(\(new_file_paths) 
      read_pos_cnv_results(file = new_file_paths,
                           sheet = get_amp_sheetname(new_file_paths))) |> 
  list_rbind()

new_std_dev_collated <- new_file_paths |> 
  map(\(new_file_paths) 
      read_stdev_results(file = new_file_paths,
                         sheet = get_amp_sheetname(new_file_paths))) |> 
  list_rbind()

new_percent_138_collated <- new_file_paths |> 
  map(\(new_file_paths) 
      read_percent_138_results(file = new_file_paths,
                               sheet = get_amp_sheetname(new_file_paths))) |> 
  list_rbind()

# Load previously collated data -----------------------------------------------------

amp_gene_results <- read_csv(str_c(data_folder, "live_service/collated/",
                                   "/live_service_amp_gene_results_collated.csv"),
                             col_types = list(
                               worksheet = col_character(),
                               labno = col_character(),
                               suffix = col_character(),
                               patient_name = col_character(),
                               labno_suffix = col_character(),
                               labno_suffix_worksheet = col_character(),
                               filepath = col_character(),
                               gene = col_character(),
                               max_region_fold_change = col_double(),
                               min_region_fold_change = col_double()
                             ))

std_dev_results <- read_csv(str_c(data_folder, "live_service/collated/", 
                                  "/live_service_std_dev_results_collated.csv"),
                            col_types = list(
                              worksheet = col_character(),
                              labno = col_character(),
                              suffix = col_character(),
                              patient_name = col_character(),
                              labno_suffix = col_character(),
                              labno_suffix_worksheet = col_character(),
                              filepath = col_character(),
                              st_dev_signal_adjusted_log2_ratios = col_double()
                            ))

percent_138_results <- read_csv(str_c(data_folder, "live_service/collated/", 
                                      "/live_service_percent_138_results_collated.csv"),
                                col_types = list(
                                  worksheet = col_character(),
                                  labno = col_character(),
                                  suffix = col_character(),
                                  patient_name = col_character(),
                                  labno_suffix = col_character(),
                                  labno_suffix_worksheet = col_character(),
                                  filepath = col_character(),
                                  percent_whole_panel_covered_at_138x = col_double()
                                ))

pos_cnv_results <- read_csv(str_c(data_folder, "live_service/collated/", 
                                  "/live_service_pos_cnv_results_collated.csv"),
                            col_types = list(
                              worksheet = col_character(),
                              labno = col_character(),
                              suffix = col_character(),
                              patient_name = col_character(),
                              labno_suffix = col_character(),
                              labno_suffix_worksheet = col_character(),
                              filepath = col_character(),
                              gene = col_character(),
                              chromosome = col_character(),
                              cnv_co_ordinates = col_character(),
                              cnv_length = col_double(),
                              consequence = col_character(),
                              fold_change = col_double(),
                              p_value = col_double(),
                              start = col_double(),
                              end = col_double()
                            ))

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

if(all(nrow(new_amp_gene_collated) > 0,
       nrow(new_amp_gene_collated) > 0,
       nrow(new_std_dev_collated) > 0,
       nrow(new_percent_138_collated) > 0)) {
  
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
  message("New data added")
  
} else {
  
  amp_gene_results_updated <- amp_gene_results
  
  pos_cnv_results_updated <- pos_cnv_results
  
  std_dev_results_updated <- std_dev_results
  
  percent_138_results_updated <- percent_138_results
  
  message("No new data added")
  
}

# Checks ----------------------------------------------------------------------------

# Each PanSolid worksheet should have 48 samples on
expected_file_number <- (length(unique(pansolidv2_worksheets$worksheet)) * 48)

file_number <- length(list.files(paste0(data_folder, "live_service/raw/"), 
                                 pattern = ".xlsx"))

message(str_c("Check: ", expected_file_number, " files predicted and ",
                file_number, " files found."))

# Archive previous collated data ----------------------------------------------------

write.csv(amp_gene_results,
          paste0(data_folder, "live_service/collated/archive/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_amp_gene_results_collated.csv"),
          row.names = FALSE)

write.csv(std_dev_results,
          paste0(data_folder, "live_service/collated/archive/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_std_dev_results_collated.csv"),
          row.names = FALSE)

write.csv(pos_cnv_results,
          paste0(data_folder, "live_service/collated/archive/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_pos_cnv_results_collated.csv"),
          row.names = FALSE)

write.csv(percent_138_results,
          paste0(data_folder, "live_service/collated/archive/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_percent_138_results_collated.csv"),
          row.names = FALSE)

# Save updated collated data --------------------------------------------------------

write_csv(amp_gene_results_updated,
          paste0(data_folder, "live_service/collated/",
                           "live_service_amp_gene_results_collated.csv"))

write_csv(std_dev_results_updated,
          paste0(data_folder, "live_service/collated/",
                           "live_service_std_dev_results_collated.csv"))

write_csv(pos_cnv_results_updated,
          paste0(data_folder, "live_service/collated/",
                           "live_service_pos_cnv_results_collated.csv"))

write_csv(percent_138_results_updated,
          paste0(data_folder, "live_service/collated/",
                           "live_service_percent_138_results_collated.csv"))
