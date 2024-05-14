# Find BRAIN MATRIX samples

library(tidyverse)
library(readxl)

clin_trial_filepath <- "S:/central shared/Genetics/Mol_Shared/Cancer Team/100kGP Cancer Program Validation & Feedback/Live Clinical Trials/"

# Brain Matrix Excel ----------------------------------------------------------------

brain_matrix_data_1 <- read_excel(path = str_c(clin_trial_filepath, "Brain Matrix Data.xlsx"),
                                  sheet = "Sheet1",
                                  col_types = c("text", "text", "text", "text",
                                                "date", "text", "date", "date",
                                                "date", "date", "text", "text")) |> 
  janitor::clean_names()

brain_matrix_data_2 <- read_excel(path = str_c(clin_trial_filepath, "Brain Matrix Data.xlsx"),
                                  sheet = "Sheet2",
                                  col_types = c("text", "text", "text", "date",
                                                "text", "date", "date", "date", "date")) |> 
  janitor::clean_names()

shared_cols <- intersect(colnames(brain_matrix_data_2), colnames(brain_matrix_data_1)) 

brain_matrix_data_join <- brain_matrix_data_1 |> 
  select(all_of(shared_cols)) |> 
  rbind(brain_matrix_data_2 |> 
          select(all_of(shared_cols))) |> 
  filter(!is.na(brain_matrix_id) & brain_matrix_id != "Brain Matrix ID") |> 
  filter(!duplicated(brain_matrix_id)) |> 
  mutate(gel_id_clean = str_replace_all(string = gel_id, 
                                        pattern = " ", replacement = "")) |> 
  arrange(brain_matrix_id)

# WGS Tracker Excel -----------------------------------------------------------------

wgs_tracker <- read_excel(path = str_c(clin_trial_filepath, "WGS tracker.xlsx")) |> 
  janitor::clean_names() |> 
  mutate(referral_id_clean = str_replace_all(string = referral_id, 
                                             pattern = " ", replacement = ""),
         
         patient_id_clean = str_replace_all(string = patient_id, 
                                            pattern = " ", replacement = ""))


# Find HTML files -------------------------------------------------------------------

file_regex <- "\\d{10}_p\\d{11}_LP\\d{7}-DNA_\\w{1}\\d{2}_LP\\d{7}-DNA_\\w{1}\\d{2}.+(supplementary|supplementary\\s\\(1\\)|supplementary\\s\\(2\\)).html$"

all_htmls <- list.files(path = str_c(clin_trial_filepath, "MTB sheets/"),
                        pattern = file_regex,
                        recursive = TRUE,
                        full.names = TRUE,
                        ignore.case = TRUE)

all_html_df <- data.frame(
  s_drive_filepath = all_htmls) |> 
  mutate(wgs_p_no = str_extract(string = s_drive_filepath,
                                patter = "p\\d{11}"))

# Brain Matrix patient identifiers --------------------------------------------------

brain_matrix_patient_info <- wgs_tracker |> 
  filter(referral_id_clean %in% brain_matrix_data_join$gel_id |
           patient_id_clean %in% brain_matrix_data_join$gel_id) |> 
  filter(!is.na(patient_id_clean))

brain_matrix_p_nos <- list(brain_matrix_patient_info$patient_id_clean)

# Find Brain Matrix html paths ------------------------------------------------------

brain_matrix_filepath_df <- all_html_df |> 
  filter(wgs_p_no %in% brain_matrix_patient_info$patient_id_clean)
  
brain_matrix_filepaths <- brain_matrix_filepath_df$s_drive_filepath

# Copy Brain Matrix htmls -----------------------------------------------------------

file.copy(from = brain_matrix_filepaths, 
          to = here::here("data/brain_matrix_htmls/"))
