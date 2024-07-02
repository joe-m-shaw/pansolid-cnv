# Find WGS Samples with QIASymphony Extractions

# This is a script to identify case which have had/are having whole genome sequencing
# and also have a DNA sample extracted by the QIASymphony method from the same or 
# similar pathology block to the whole genome sequencing sample.

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

# Functions -------------------------------------------------------------------------

source(here::here("functions/dna_database_connection.R"))

source(here::here("functions/dna_database_functions.R"))

# WGS Pathway Tracker ---------------------------------------------------------------

wgs_pathway_tracker <- read_excel(path = here("data/WGS pathway tracker_copy_2024-07-01.xlsx"),
                                  sheet = "Cancer") |> 
  janitor::clean_names()

wgs_pathway_tracker_dna_no_df <- wgs_pathway_tracker |> 
  filter(!is.na(mol_db_number)) |> 
  mutate(labno = str_extract(string = mol_db_number,
                             pattern = "\\d{8}"))

wgs_dna_numbers <- wgs_pathway_tracker_dna_no_df$mol_db_number

# Get NHS numbers -------------------------------------------------------------------

wgs_nhs_numbers_df <- sample_tbl |> 
  filter(labno %in% wgs_dna_numbers) |> 
  select(labno, nhsno) |> 
  collect() |> 
  filter(!is.na(nhsno))

wgs_nhs_numbers <- wgs_nhs_numbers_df$nhsno

# Get DNA numbers -------------------------------------------------------------------

all_samples_from_wgs_patients_df <-  sample_tbl |> 
  filter(nhsno %in% wgs_nhs_numbers) |> 
  select(labno, nhsno, pathno) |> 
  collect()

all_samples_from_wgs_patients <- all_samples_from_wgs_patients_df$labno

all_samples_from_wgs_patients_extractions <- extraction_tbl |> 
  filter(lab_no %in% all_samples_from_wgs_patients) |> 
  select(lab_no, extraction_id, extraction_batch_fk) |> 
  collect()

extraction_batches <- all_samples_from_wgs_patients_extractions$extraction_batch_fk

extraction_types <- extraction_batch_tbl |> 
  filter(extraction_batch_id %in% extraction_batches) |> 
  collect() |> 
  left_join(extraction_method_key, join_by(extraction_method_fk == extraction_method_id ))

extraction_ids_with_method <- all_samples_from_wgs_patients_extractions |> 
  left_join(extraction_types |> 
              select(extraction_batch_id, method_name),
            join_by(extraction_batch_fk == extraction_batch_id))

lab_nos_with_method <- all_samples_from_wgs_patients_df |> 
  left_join(extraction_ids_with_method, join_by(labno == lab_no)) |> 
  filter(method_name %in% c("Fresh tissue", "QIAsymphony_DNA_FFPE" )) |> 
  arrange(nhsno) |> 
  left_join(wgs_pathway_tracker_dna_no_df |> 
              select(labno, ngis_patient_id, ngis_referral_id),
            by = "labno") |> 
  arrange(nhsno)

qiasymphony_lab_nos <- lab_nos_with_method |> 
  filter(method_name == "QIAsymphony_DNA_FFPE")

wgs_samples_with_qiasymphony_dna <- lab_nos_with_method |> 
  filter(nhsno %in% qiasymphony_lab_nos$nhsno &
           !is.na(ngis_patient_id))

# Get HTML files --------------------------------------------------------------------

patient_ids_to_get <- wgs_samples_with_qiasymphony_dna$ngis_patient_id

clin_trial_filepath <- "S:/central shared/Genetics/Mol_Shared/Cancer Team/100kGP Cancer Program Validation & Feedback/Live Clinical Trials/"

file_regex <- "\\d{10}_p\\d{11}_LP\\d{7}-DNA_\\w{1}\\d{2}_LP\\d{7}-DNA_\\w{1}\\d{2}.+(supplementary|supplementary\\s\\(1\\)|supplementary\\s\\(2\\)).html$"

all_htmls <- list.files(path = str_c(clin_trial_filepath, "MTB sheets/"),
                        pattern = file_regex,
                        recursive = TRUE,
                        full.names = TRUE,
                        ignore.case = TRUE)

all_html_df <- data.frame(
  s_drive_filepath = all_htmls) |> 
  mutate(wgs_p_no = str_extract(string = s_drive_filepath,
                                pattern = "p\\d{11}"))

file_path_df_samples_to_get <- all_html_df |> 
  filter(wgs_p_no %in% patient_ids_to_get)

file_paths_to_copy <- file_path_df_samples_to_get$s_drive_filepath

file.copy(from = file_paths_to_copy, to = here::here("data/wgs_result_htmls/"))
