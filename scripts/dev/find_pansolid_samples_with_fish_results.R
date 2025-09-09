# Find PanSolid NGS results for cases with FISH results

# Packages and connection -------------------------------------------------

library(tidyverse)
library(readxl)

source(here::here("scripts/connect_to_dna_db.R"))

# Load data ---------------------------------------------------------------

fish_data <- read_excel(paste0(config::get("data_folderpath"),
                               "validation/DOC6791_chromosome_arms/",
                               "raw/fish/",
                               "Glioma FISH results jun24_June25 for NGS CNV.xlsx"),
                        sheet = "Sheet1") |> 
  janitor::clean_names() |> 
  mutate(nhsno = str_replace_all(nhs_number, "\\s", ""))

wgs_tracker <- read_excel(paste0(config::get("data_folderpath"),
                                 "validation/DOC6791_chromosome_arms/",
                                 "2025-07-29_WGS pathway tracker.xlsx"),
                          sheet = "Cancer") |> 
  janitor::clean_names() |> 
  mutate(nhsno = str_replace_all(nhs_number, "\\s", ""))


# Find PanSolid worksheet samples -----------------------------------------

all_worksheets <- dna_db_worksheets |> 
  select(pcrid, date, description) |> 
  collect() |> 
  mutate(worksheet = paste0("WS", pcrid))

stopifnot(nrow(all_worksheets) > 0)

ps_ws_info <- all_worksheets |> 
  filter(pcrid >= 140721) |> 
  filter(grepl(pattern = "pansolid|pan-solid|pan_solid|pan solid", 
               x = description,
               ignore.case = TRUE)) |> 
  mutate(ps_category = case_when(
    grepl(pattern = "jBRCA|j_BRCA|j-BRCA|jew",
          x = description,
          ignore.case = TRUE) ~"PanSolid Jewish BRCA",
    TRUE ~"PanSolid FFPE"
  )) |> 
  filter(ps_category == "PanSolid FFPE")

ps_ws_vector <- unique(ps_ws_info$pcrid)

all_ps_samples <- dna_db_pcr_records |> 
  filter(pcrid %in% ps_ws_vector) |> 
  select(pcrid, sample) |> 
  collect()

ps_labno_vector <- unique(all_ps_samples$sample)

ps_nhsno_tbl <- sample_tbl |> 
  filter(labno %in% ps_labno_vector) |> 
  select(labno, nhsno) |> 
  collect()

ps_nhsno_tbl_clean <- ps_nhsno_tbl |> 
  mutate(nhsno = str_replace_all(nhsno, "\\s", ""))

all_ps_samples_with_nhsno <- all_ps_samples |> 
  left_join(ps_nhsno_tbl_clean, join_by("sample" == "labno"))

# Find PanSolid samples with FISH results ---------------------------------

ps_samples_with_fish_results <- all_ps_samples_with_nhsno |> 
  filter(nhsno %in% fish_data$nhsno) 

ps_fish_result_sample_vector <- ps_samples_with_fish_results$sample

sample_ids <- sample_tbl |> 
  select(labno, nhsno, firstname, surname, dob) |> 
  filter(labno %in% ps_fish_result_sample_vector) |> 
  collect()

ps_samples_with_fish_results_and_names <- ps_samples_with_fish_results |> 
  left_join(sample_ids, join_by("sample" == "labno"))
