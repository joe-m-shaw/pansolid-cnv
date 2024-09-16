# Get Patient Information for Sample Cohort

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(here)

# Source scripts --------------------------------------------------------------------

source(here("scripts/set_shared_drive_filepath.R"))

source(here("scripts/connect_to_dna_db.R"))

source(here("functions/dna_db_functions.R"))

# Get patient information -----------------------------------------------------------

amp_validation_stdev_results_collated <- read_csv(paste0(data_folder, 
                                                     "validation/processed/",
                                                     "amp_validation_stdev_results_collated.csv"),
                                              col_types = list(
                                                "labno" = col_character()))

validation_ddpcr_collated <- read_csv(paste0(data_folder,
                                             "validation/processed/",
                                             "validation_ddpcr_collated.csv")) |> 
  filter(grepl(pattern = "\\d{8}", x = sample))

wgs_html_ids <- read_csv(paste0(data_folder, 
                                "validation/processed/",
                                "wgs_html_ids.csv"),
                         col_types = list(
                           "labno" = col_character()
                         ))

sample_labnos <- unique(c(amp_validation_stdev_results_collated$labno, 
                          validation_ddpcr_collated$sample, 
                          wgs_html_ids$labno))

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
                                         2))

# Add NCC and tissue type -----------------------------------------------------------

patient_info_ncc <- patient_info |> 
  mutate(ncc = parse_ncc(comments),
         tissue = as.numeric(tissue)) |> 
  left_join(tissue_types |> 
              select(tissue_type_id, tissue_type),
            join_by("tissue" == "tissue_type_id")) |> 
  # WGS frozen tissue is coded as "Other" on DNA Database
  mutate(tissue_type = ifelse(tissue_type == "Other", "Frozen tissue",
                              tissue_type))

patient_ncc_manual <- patient_info_ncc |> 
  filter(is.na(ncc)) |> 
  select(labno) |> 
  mutate(ncc_manual = "")

write.csv(x = patient_ncc_manual,
          file = paste0(data_folder, "validation/processed/", 
                        "patient_ncc_manual.csv"),
          row.names = FALSE)

manual_ncc_values <- read_csv(file = paste0(data_folder, "validation/processed/",
                       "patient_ncc_manual_edit.csv"),
                       col_types = list(
                         "labno" = col_character(),
                         "ncc_manual" = col_character()))

patient_info_ncc_edit <- patient_info_ncc |> 
  left_join(manual_ncc_values, by = "labno") |> 
  mutate(neoplastic_cell_content = case_when(
    
    is.na(ncc) & !is.na(ncc_manual) ~ncc_manual,
    !is.na(ncc) & is.na(ncc_manual) ~ncc,
    TRUE ~"Error"))

# Extract cancer type from comments field -------------------------------------------

cancer_type_vector <- c("(?:O|o)varian",
                         "(?:C|c)olorectal",
                         "(?:L|l)ung",
                         "(?:B|b)ladder",
                         "CNS",
                         "(?:E|e)ndometrial",
                         "(?:M|m)elanoma",
                         "Oligodendroglioma",
                         "(?:P|p)rostate",
                         "Glioma",
                         "Pancreatic",
                         "Glioblastoma")

cancer_type_regex <- paste0("(",
                            paste(cancer_type_vector, collapse  = "|"),
                            ")")

patient_info_cancer_type <- patient_info_ncc_edit |> 
  mutate(cancer_comment = tolower(str_extract(pattern = cancer_type_regex,
                        string = comments,
                        group = 1)),
         
         cancer_group = case_when(
           
           cancer_comment %in% c("cns", "glioblastoma", "glioma",
                                 "oligodendroglioma") ~"central nervous system",
           
           TRUE ~cancer_comment))

cancer_type_manual <- patient_info_cancer_type |> 
  filter(is.na(cancer_group)) |> 
  select(labno) |> 
  mutate(cancer_group_manual = "")

write.csv(x = cancer_type_manual,
          file = paste0(data_folder, "validation/processed/", 
                        "cancer_type_manual.csv"),
          row.names = FALSE)

manual_cancer_types <- read_csv(file = paste0(data_folder, "validation/processed/",
                                            "cancer_type_manual_edit.csv"),
                              col_types = list(
                                "labno" = col_character(),
                                "cancer_group_manual" = col_character()))

patient_info_cancer_type_edit <- patient_info_cancer_type |> 
  left_join(manual_cancer_types, by = "labno") |> 
  mutate(cancer_tissue_source = case_when(
    
    is.na(cancer_group) & !is.na(cancer_group_manual) ~cancer_group_manual,
    !is.na(cancer_group) & is.na(cancer_group_manual) ~cancer_group,
    TRUE ~"Error"))

# Add DNA extraction method ---------------------------------------------------------

extraction_methods <- get_extraction_method(sample_vector = sample_labnos) |> 
  # Some samples were on a QIAsymphony run with an instrument failure
  # so they have multiple extraction entries for the same lab number
  filter(!duplicated(labno)) |> 
  # The method for extraction batch 81056 was incorrectly entered
  # as "QIAsymphony_RNA_FFPE", but it should have been "QIAsymphony_DNA_FFPE"
  mutate(method_name = ifelse(extraction_batch_fk == 81056,
                              "QIAsymphony_DNA_FFPE",
                              method_name))

patient_info_extraction_method <- patient_info_cancer_type_edit |> 
  left_join(extraction_methods, by = "labno") |> 
  mutate(method_name = case_when(
    is.na(method_name) & surname == "Seraseq" ~"No DNA extraction performed",
    TRUE ~method_name))

# Get worksheet details -------------------------------------------------------------

pansolid_worksheets <- data.frame(
  "worksheet" = c(unique(amp_validation_stdev_results_collated$worksheet))) |> 
  mutate(pcrid = str_extract(string = worksheet, 
                      pattern = "WS(\\d{6})",
                      group = 1))

pansolid_pcrids <- pansolid_worksheets$pcrid

pansolid_worksheet_details <- dna_db_worksheets |> 
  filter(pcrid %in% pansolid_pcrids) |> 
  select(pcrid, description, date) |> 
  collect()

# Export results --------------------------------------------------------------------

validation_sample_patient_info <- patient_info_extraction_method |> 
  select(-c(tissue))

write.csv(x = validation_sample_patient_info,
          file = paste0(data_folder, "validation/processed/", 
                        "amp_validation_sample_patient_info.csv"),
          row.names = FALSE)

write.csv(x = pansolid_worksheet_details,
          file = paste0(data_folder, "validation/processed/", 
                        "pansolid_worksheet_details.csv"),
          row.names = FALSE)
