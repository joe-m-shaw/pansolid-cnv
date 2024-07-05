# Deletion cohort samples

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

# Functions and filepaths -----------------------------------------------------------

source(here("scripts/set_shared_drive_filepath.R"))
source(here("functions/dna_database_connection.R"))
source(here("functions/dna_database_functions.R"))
source(here("functions/cnv_functions.R"))

# PanSolid filepaths ----------------------------------------------------------------

pansolidv2_worksheets <- read_excel(paste0(data_folder,
                                           "pansolid_live_service_worksheets.xlsx"))

worksheet_list <- list(pansolidv2_worksheets$worksheet)

s_drive_filepaths <- worksheet_list |> 
  map(\(worksheet_list) get_annotated_filepaths(worksheet = worksheet_list)) |> 
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

# Deletion cohort samples -----------------------------------------------------------

deletion_cohort <- read_excel(path = paste0(data_folder, "deletion_cohort_samples.xlsx"),
                              col_types = c("text", "text"))

sample_extractions <- get_extraction_method(sample_vector = deletion_cohort$labno) |> 
  select(labno, method_name)

worksheet_info <- s_drive_file_df |> 
  filter(labno %in% deletion_cohort$labno) |> 
  select(labno, worksheet, panel)

deletion_samples <- deletion_cohort$labno

sample_info <- sample_tbl |> 
  filter(labno %in% deletion_samples) |> 
  select(labno, firstname, surname, comments) |> 
  collect()

deletion_cohort_mod <- deletion_cohort |> 
  left_join(worksheet_info, by = "labno") |> 
  left_join(sample_info, by = "labno") |> 
  left_join(sample_extractions, by = "labno") |> 
  rename(category = note,
         extraction_method = method_name) |> 
  arrange(worksheet)

write.csv(deletion_cohort_mod, paste0(outputs_folder, "deletion_cohort_info.csv"),
          row.names = FALSE)
