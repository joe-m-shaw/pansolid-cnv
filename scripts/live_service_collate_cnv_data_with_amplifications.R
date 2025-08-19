# Update PanSolidv2 Copy Number Variant Data for Amplifications Service 
# (April 2024-May 2025)

# Packages ---------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

# Functions and filepaths ------------------------------------------------------

data_folderpath <- config::get("data_folderpath")

source(here("functions/pansolid_cnv_excel_functions.R"))

source(here::here("scripts/connect_to_dna_db.R"))

# Find worksheets --------------------------------------------------------------

all_worksheets <- dna_db_worksheets |> 
  select(pcrid, date, description) |> 
  collect() |> 
  mutate(worksheet = paste0("WS", pcrid))

stopifnot(nrow(all_worksheets) > 0)

ps_ws_info <- all_worksheets |> 
  # Jewish BRCA samples were only tested using the CNV pipeline since
  # WS147693
  # New CNV Excel layout with deletions and LOH started with WS152758
  filter(pcrid >= 147693 & pcrid < 152758) |> 
  filter(grepl(pattern = "pansolid|pan-solid|pan_solid|pan solid", 
               x = description,
               ignore.case = TRUE))

ps_worksheets <- ps_ws_info$worksheet

# S drive filepaths ------------------------------------------------------------

s_drive_filepaths <- ps_worksheets |> 
  map(\(ps_worksheets) get_annotated_filepaths(worksheet = ps_worksheets)) |> 
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
                             group = 1)) |> 
  arrange(worksheet, labno)

if (any(grepl(pattern = "/", x = s_drive_file_df$filename))) {
  stop("Error: filenames contain backslashes") 
}

if (anyNA(s_drive_file_df)) {
  warning("Error: there are NA values in the filepath table")
}

message("Sample filepaths compiled")

# Load collated data -----------------------------------------------------------

amp_service_collated_data_folder <- paste0(config::get("data_folderpath"),
                               "live_service/collated/pansolid_amplifications_live_service/")

stdev_collated <- read_csv(paste0(amp_service_collated_data_folder,
                              "live_service_std_dev_results_collated.csv"))

pos_cnv_collated <- read_csv(paste0(amp_service_collated_data_folder,
                                    "live_service_pos_cnv_results_collated.csv"))

amp_genes_collated <- read_csv(paste0(amp_service_collated_data_folder,
                                      "live_service_amp_gene_results_collated.csv"))

percent_138_collated <- read_csv(paste0(amp_service_collated_data_folder,
                                        "live_service_percent_138_results_collated.csv"))

# Get filenames from collated data ---------------------------------------------

stdev_with_filenames <- stdev_collated |> 
  mutate(filename = str_extract(string = filepath, 
                                pattern = str_replace(string = pansolidv2_excel_regex, 
                                                      pattern = "\\^", 
                                                      replacement = "")))

stopifnot(anyNA(stdev_with_filenames$filename) == FALSE)

# Identify files without CNV tabs ----------------------------------------------

# Some samples have such low coverage that a CNV tab is not included in 
# the output Excel. These files will cause an error message at the data 
# collation stage if they are included.

samples_without_amp_tabs <- c("24023280", "24025207", "24027566", "24033006",
                              "24033959", "24038848")

# Identify files not already collated ------------------------------------------

message("Identifying new files for collation")

# Some samples were reanalysed with the new PanSolid pipeline
# with v2b panels. These results are not compatible with the functions for
# reading old-format PanSolid Excels.
reanalysed_samples <- paste(c("v2PANSOLID_WS152548_25022367",
                              "v2NF1_PS_WS152693_25010089",
                              "v2M1_tLYNCH_PS_WS152171_25015829",
                              "v2b"), collapse = "|")

ps_new_filepath_df <- s_drive_file_df |> 
  filter(!filename %in% stdev_with_filenames$filename &
           !labno %in% samples_without_amp_tabs &

           !grepl(pattern = reanalysed_samples, 
                  x = filename))

if(length(ps_new_filepath_df) > 0) {
  message(paste0(length(ps_new_filepath_df$filepath),
                 " new files identified"))
} else {
  stop("No new files identified")
}

# Prepare raw data folder ------------------------------------------------------

raw_folder_path <- paste0(config::get("data_folderpath"),
                          "live_service/raw/")

if(length(list.files(raw_folder_path)) != 0){
  stop("Raw file folder is not empty")
} else {
  message("Raw data folder is empty")
}

# Copy new files to raw data folder --------------------------------------------

new_file_paths_to_copy <-  ps_new_filepath_df$filepath

file.copy(from = new_file_paths_to_copy,
          to = raw_folder_path)

new_filepaths <- list.files(path = raw_folder_path,
                            full.names = TRUE,
                            pattern = "Annotated.*.xlsx")

message(paste0(length(new_filepaths), " files copied into raw data folder"))

# Collate new file data --------------------------------------------------------

message("Collating new data files")

new_std_dev_collated <- new_filepaths |> 
  map(\(new_filepaths) 
      read_stdev_results(file = new_filepaths,
                         sheet = get_sheetname(filepath = new_filepaths,
                                               sheet_regex = "Amplifications|CNVs_"))) |> 
  list_rbind()

new_percent_138_collated <- new_filepaths |> 
  map(\(new_filepaths) 
      read_percent_138_results(file = new_filepaths,
                               sheet = get_sheetname(filepath = new_filepaths,
                                                     sheet_regex = "Amplifications|CNVs_"))) |> 
  list_rbind()

new_amp_gene_collated <- new_filepaths |> 
  map(\(new_filepaths) 
      read_all_amp_genes_results(file = new_filepaths,
                                 sheet = get_sheetname(filepath = new_filepaths,
                                                       sheet_regex = "Amplifications|CNVs_"))) |> 
  list_rbind()

new_pos_cnv_collated <- new_filepaths |> 
  map(\(new_filepaths) 
      read_pos_cnv_results(file = new_filepaths,
                           sheet = get_sheetname(filepath = new_filepaths,
                                                 sheet_regex = "Amplifications|CNVs_"))) |> 
  list_rbind()

# Check data -------------------------------------------------------------------

stopifnot(anyNA.data.frame(new_std_dev_collated) == FALSE)

stopifnot(anyNA.data.frame(new_percent_138_collated) == FALSE)

stopifnot(anyNA.data.frame(new_amp_gene_collated) == FALSE)

stopifnot(anyNA.data.frame(new_pos_cnv_collated) == FALSE)

stopifnot(nrow(new_std_dev_collated) > 0)

stopifnot(nrow(new_percent_138_collated) > 0)

stopifnot(nrow(new_amp_gene_collated) > 0)

stopifnot(nrow(new_pos_cnv_collated) > 0)

# Check columns ----------------------------------------------------------------

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

# Add new data to collated data ------------------------------------------------

std_dev_results_updated <- rbind(stdev_collated |> 
                                   select(all_of(std_dev_cols)),
                                 new_std_dev_collated |> 
                                   select(all_of(std_dev_cols)))

percent_138_results_updated <- rbind(percent_138_collated |> 
                                       select(all_of(percent_138_cols)),
                                     new_percent_138_collated |> 
                                       select(all_of(percent_138_cols)))

amp_gene_results_updated <- rbind(amp_genes_collated |> 
                                    select(all_of(all_amp_cols)),
                                  new_amp_gene_collated |> 
                                    select(all_of(all_amp_cols)))
  
pos_cnv_results_updated <- rbind(pos_cnv_collated |> 
                                   select(all_of(pos_cnv_cols)),
                                 new_pos_cnv_collated |> 
                                   select(all_of(pos_cnv_cols)))

# Archive previous collated data -----------------------------------------------

write.csv(amp_genes_collated,
          paste0(data_folderpath, "live_service/collated/archive/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_amp_gene_results_collated.csv"),
          row.names = FALSE)

write.csv(stdev_collated,
          paste0(data_folderpath, "live_service/collated/archive/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_std_dev_results_collated.csv"),
          row.names = FALSE)

write.csv(pos_cnv_collated,
          paste0(data_folderpath, "live_service/collated/archive/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_pos_cnv_results_collated.csv"),
          row.names = FALSE)

write.csv(percent_138_collated,
          paste0(data_folderpath, "live_service/collated/archive/",
                           format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           "live_service_percent_138_results_collated.csv"),
          row.names = FALSE)

# Save updated collated data ---------------------------------------------------

write_csv(amp_gene_results_updated,
          paste0(data_folderpath, "live_service/collated/pansolid_amplifications_live_service/",
                           "live_service_amp_gene_results_collated.csv"))

write_csv(std_dev_results_updated,
          paste0(data_folderpath, "live_service/collated/pansolid_amplifications_live_service/",
                           "live_service_std_dev_results_collated.csv"))

write_csv(pos_cnv_results_updated,
          paste0(data_folderpath, "live_service/collated/pansolid_amplifications_live_service/",
                           "live_service_pos_cnv_results_collated.csv"))

write_csv(percent_138_results_updated,
          paste0(data_folderpath, "live_service/collated/pansolid_amplifications_live_service/",
                           "live_service_percent_138_results_collated.csv"))

# Delete raw files -------------------------------------------------------------

message("Deleting new raw files")

file.remove(new_filepaths)

if(length(list.files(raw_folder_path)) == 0){
  message("Raw file folder is empty")
}
