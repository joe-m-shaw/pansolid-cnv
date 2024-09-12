# 1p19q codeletion samples from Liverpool

# Packages ---------------------------------------------------------------------------------------

library(here)
library(tidyverse)
library(readxl)
library(janitor)

# Source scripts ---------------------------------------------------------------------------------

source(here("scripts/set_shared_drive_filepath.R"))
source(here("scripts/connect_to_dna_db.R"))
source(here("scripts/load_processed_live_service_data.R"))
source(here("functions/join_dna_submission_sheets.R"))
source(here("functions/dna_db_functions.R"))

# Load FISH results from Liverpool ---------------------------------------------------------------

liverpool_fish <- read_excel(path = paste0(data_folder,
                                           "liverpool_glioma_FISH_results.xlsx"),
                             sheet = "Sheet2") |> 
  clean_names() |> 
  mutate(nhs_no_clean = str_replace_all(nhs_no, pattern = " ", replacement = ""))

# Have any samples already been tested on PanSolid? ----------------------------------------------

labnos_tested_on_pansolid <- live_service_std_dev_results_collated$labno

pansolid_info <- sample_tbl |> 
  filter(labno %in% labnos_tested_on_pansolid) |> 
  select(labno, nhsno, firstname, surname) |> 
  collect()

# Only finds one sample
base::intersect(pansolid_info$nhsno, liverpool_fish$nhs_no_clean)

# What are the labnos for the FISH samples tested on QIAseq Core? --------------------------------

fish_nhs_nos <- liverpool_fish$nhs_no_clean

fish_manchester_info <- sample_tbl |> 
  filter(nhsno %in% fish_nhs_nos) |> 
  select(labno, nhsno) |> 
  collect()

fish_manchester_labnos <- fish_manchester_info$labno

fish_manchester_ngs_results <- results_tbl |> 
  filter(labno %in% fish_manchester_labnos) |> 
  select(labno, surname, test, genotype) |> 
  collect() |> 
  filter(grepl(pattern = "core", x = test, 
               ignore.case = TRUE)) |> 
  filter(!duplicated(labno)) |> 
  left_join(fish_manchester_info, by = "labno") |> 
  left_join(liverpool_fish |> 
              select(nhs_no_clean, genotypesummarycomment), join_by("nhsno" == "nhs_no_clean"))

fish_manchester_extractions <- get_extraction_method(sample_vector = fish_manchester_labnos)

# Get DNA concentrations -------------------------------------------------------------------------

qiaseq_core_submission_sheets <- join_qiaseq_core_submission_sheets() |> 
  filter(!duplicated(labno))

fish_manchester_ngs_results_with_info <- fish_manchester_ngs_results |> 
  left_join(qiaseq_core_submission_sheets |> 
              select(labno, stock_qubit),
            by = "labno") |> 
  left_join(fish_manchester_extractions |> 
              select(labno, method_name),
            by = "labno")
