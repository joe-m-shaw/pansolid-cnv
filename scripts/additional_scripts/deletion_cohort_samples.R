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

# Deletion cohort samples -----------------------------------------------------------

deletion_cohort <- read_excel(path = paste0(data_folder, "deletion_cohort_samples.xlsx"),
                              col_types = c("text", "text"))

sample_extractions <- get_extraction_method(sample_vector = deletion_cohort$labno) |> 
  select(labno, method_name) |> 
  distinct()

if(any(duplicated(sample_extractions$labno))) stop()

deletion_samples <- deletion_cohort$labno

sample_info <- sample_tbl |> 
  filter(labno %in% deletion_samples) |> 
  select(labno, firstname, surname, nhsno, comments) |> 
  collect()

if(any(duplicated(sample_info$labno))) stop()

# Panel information -----------------------------------------------------------------

panel_df <- read_csv(file = paste0(collated_data_path, 
                                     "pansolidv2_sample_worksheet_panel_information.csv")) |> 
  mutate(labno = as.character(labno))

worksheet_info <- panel_df |> 
  filter(labno %in% deletion_cohort$labno) |> 
  select(labno, worksheet, panel) 

if(any(duplicated(worksheet_info$labno))) stop()

# WGS patient identifiers -----------------------------------------------------------

wgs_htmls <- list.files(path = paste0(data_folder, "wgs_result_htmls/"),
                        full.names = TRUE,
                        pattern = "*.html")

wgs_pids <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_pid_text(wgs_htmls)) |> 
  list_rbind()

wgs_headers <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_header(wgs_htmls)) |> 
  list_rbind()

wgs_info <- inner_join(wgs_pids, wgs_headers, by = "filepath") 

# Collate information ---------------------------------------------------------------

deletion_cohort_mod <- deletion_cohort |> 
  inner_join(worksheet_info, by = "labno") |> 
  inner_join(sample_info, by = "labno") |> 
  inner_join(sample_extractions, by = "labno") |> 
  left_join(wgs_info |> 
              select(nhs_no_clean, wgs_r_no, wgs_p_no),
            join_by("nhsno" == "nhs_no_clean")) |> 
  rename(category = note,
         extraction_method = method_name) |> 
  arrange(worksheet) |> 
  relocate(category)

if(any(duplicated(deletion_cohort_mod$labno))) stop()

# Export ----------------------------------------------------------------------------

time <- format(Sys.time(), "%Y-%m-%d")

write.csv(deletion_cohort_mod, paste0(outputs_folder, time, "_deletion_cohort_info.csv"),
          row.names = FALSE)
