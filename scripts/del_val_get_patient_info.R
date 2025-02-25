# Get Patient Information for Sample Cohort

# Packages ---------------------------------------------------------------------

library(tidyverse)
library(here)

# Source scripts ---------------------------------------------------------------

source(here("scripts/connect_to_dna_db.R"))

source(here("functions/dna_db_functions.R"))

# Get patient information ------------------------------------------------------

del_val_collated_stdev <- read_csv(paste0(config::get("data_folderpath"), 
                                          "validation/DOC6567_deletions/processed/",
                                          "del_val_collated_stdev.csv"),
                                              col_types = list(
                                                "labno" = col_character()))

del_val_ddpcr_collated <- read_csv(paste0(config::get("data_folderpath"), 
                                          "validation/DOC6567_deletions/processed/",
                                             "del_val_ddpcr_collated.csv")) |> 
  filter(grepl(pattern = "\\d{8}", x = sample))

del_val_wgs_html_ids <- read_csv(paste0(config::get("data_folderpath"), 
                                        "validation/DOC6567_deletions/processed/",
                                "del_val_wgs_html_ids.csv"),
                         col_types = list(
                           "labno" = col_character()
                         ))

sample_labnos <- unique(c(del_val_collated_stdev$labno, 
                          del_val_ddpcr_collated$sample, 
                          del_val_wgs_html_ids$labno))

patient_info <- sample_tbl |> 
  filter(labno %in% sample_labnos) |> 
  select(labno, firstname, surname, dob, date_in, tissue, nhsno, pathno,
         comments) |> 
  collect() |> 
  mutate(dob = as_date(substr(dob, start = 1, stop = 10),
                       format = "%Y-%m-%d"),
         date_sample_received = as_date(substr(date_in, start = 1, stop = 10),
                                        format = "%Y-%m-%d"),
         years_at_sample_receipt = round(interval(dob, date_sample_received) / years(1), 
                                         2),
         ncc = parse_ncc(comments),
         tissue = as.numeric(tissue)) |> 
  left_join(tissue_types |> 
              select(tissue_type_id, tissue_type),
            join_by("tissue" == "tissue_type_id")) |> 
  # WGS frozen tissue is coded as "Other" on DNA Database
  mutate(tissue_type = ifelse(tissue_type == "Other", "Frozen tissue",
                              tissue_type))

# Add DNA extraction method ----------------------------------------------------

extraction_methods <- get_extraction_method(sample_vector = sample_labnos) |> 
  filter(!duplicated(labno)) |> 
  mutate(method_name = ifelse(extraction_batch_fk == 81056,
                              "QIAsymphony_DNA_FFPE",
                              method_name))

patient_info_extraction_method <- patient_info |> 
  left_join(extraction_methods, by = "labno") |> 
  mutate(method_name = case_when(
    is.na(method_name) & 
      tissue_type == "DNA" ~"DNA received. No DNA extraction performed",
    method_name == "QIAsymphony_DNA_FFPE" ~"QIAsymphony",
    TRUE ~method_name)) |> 
  rename(extraction_method = method_name)

# Get worksheet details --------------------------------------------------------

pansolid_worksheets <- data.frame(
  "worksheet" = c(unique(del_val_collated_stdev$worksheet))) |> 
  mutate(pcrid = str_extract(string = worksheet, 
                      pattern = "WS(\\d{6})",
                      group = 1))

pansolid_pcrids <- pansolid_worksheets$pcrid

pansolid_worksheet_details <- dna_db_worksheets |> 
  filter(pcrid %in% pansolid_pcrids) |> 
  select(pcrid, description, date) |> 
  collect()

# Export results ---------------------------------------------------------------

del_val_sample_patient_info <- patient_info_extraction_method |> 
  select(labno, firstname, surname, nhsno, dob, date_sample_received,
         years_at_sample_receipt, tissue, tissue_type, pathno, ncc,
         extraction_method, extraction_batch_fk, comments)

write_csv(x = del_val_sample_patient_info,
          file = paste0(config::get("data_folderpath"), 
                        "validation/DOC6567_deletions/processed/", 
                        "del_val_sample_patient_info.csv"))

write_csv(x = pansolid_worksheet_details,
          file = paste0(config::get("data_folderpath"), 
                        "validation/DOC6567_deletions/processed/", 
                        "del_val_pansolid_worksheet_details.csv"))
