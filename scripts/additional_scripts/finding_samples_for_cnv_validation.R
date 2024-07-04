# Finding Samples for the CNV Validation

# The aim of this script is to find existing samples with data which could be used for the 
# PanSolid CNV validation.

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

# Functions -------------------------------------------------------------------------

source(here::here("functions/dna_database_connection.R"))

source(here::here("functions/dna_database_functions.R"))

# Identify samples tested on PanSolid -----------------------------------------------

pansolidv2_worksheets <- read_excel(here("data/pansolid_live_service_worksheets.xlsx")) |> 
  mutate(pcrid = str_replace(string = worksheet,
                             pattern = "WS",
                             replacement = ""))

pcrid_list <- pansolidv2_worksheets$pcrid

pansolid_worksheet_samples <- dna_db_pcr_records |> 
  filter(pcrid %in% pcrid_list) |> 
  select(pcrid, sample, name) |> 
  collect() |> 
  rename(labno = sample)

samples_on_ps_worksheets <- pansolid_worksheet_samples$labno

pansolid_sample_info <- sample_tbl |> 
  filter(labno %in% samples_on_ps_worksheets) |> 
  select(labno, firstname, surname, pathno, nhsno, comments) |> 
  collect()

# Finding WGS samples with QIAsymphony extractions ----------------------------------

# We want to identify samples which have had/are having whole genome sequencing
# and also have a DNA sample extracted by the QIASymphony method from the same or 
# similar pathology block to the whole genome sequencing sample.

## Find WGS samples -----------------------------------------------------------------

wgs_pathway_tracker <- read_excel(path = here("data/WGS pathway tracker_copy_2024-07-01.xlsx"),
                                  sheet = "Cancer") |> 
  janitor::clean_names()

wgs_pathway_tracker_dna_no_df <- wgs_pathway_tracker |> 
  filter(!is.na(mol_db_number)) |> 
  mutate(labno = str_extract(string = mol_db_number,
                             pattern = "\\d{8}"))

wgs_labnos <- wgs_pathway_tracker_dna_no_df$mol_db_number

## Find WGS NHS numbers -------------------------------------------------------------

wgs_nhs_numbers_df <- sample_tbl |> 
  filter(labno %in% wgs_labnos) |> 
  select(labno, nhsno) |> 
  collect() |> 
  filter(!is.na(nhsno))

wgs_nhsnos <- wgs_nhs_numbers_df$nhsno

## Find WGS samples tested on PanSolid ----------------------------------------------

all_samples_from_wgs_patients_df <-  sample_tbl |> 
  filter(nhsno %in% wgs_nhsnos) |> 
  select(labno, nhsno, firstname, surname, pathno) |> 
  collect()

wgs_samples_tested_on_pansolid <- pansolid_worksheet_samples |> 
  filter(labno %in% all_samples_from_wgs_patients_df$labno)

## Find WGS samples with QIAsymphony extractions ------------------------------------

all_samples_from_wgs_patients <- all_samples_from_wgs_patients_df$labno

wgs_samples_extraction_info <- get_extraction_method(sample_vector =
                                                       all_samples_from_wgs_patients)

all_samples_from_wgs_patients_df_with_extraction <- all_samples_from_wgs_patients_df |> 
  left_join(wgs_samples_extraction_info |> 
              select(labno, method_name), by = "labno") |> 
  # Remove blood sample extractions
  filter(!method_name %in% c("Chemagen 360", "Saliva Chemagic 360-D", "FFPE RNA")) |> 
  arrange(nhsno)

# I manually went through and identified samples which:
# - had a QiaSymphony extraction
# - were from the same (or similar) pathology block to the WGS samples

# I couldn't think of a non-manual way which didn't involve committing patient 
# identifiable information.

# Results saved here: "data/qiasymphony_extraction_samples_with_wgs_testing.xlsx"

# Find samples with PanSolid and NF2 LOH testing ------------------------------------

# NF2 microsatellite testing can be used to confirm loss of heterozygosity (LOH).

dlims_results <- results_tbl |> 
  filter(genodate > "2024-03-01 00:00:00") |> 
  select(labno, surname, test, exon, genotype, genotype2,
         genocomm) |> 
  collect() |> 
  filter(!labno %in% c("water", "Water", "WATER"))

loh_results <- dlims_results |> 
  filter(grepl(pattern = "LOH", x = test, ignore.case = TRUE))

pansolid_samples_with_loh_results <- loh_results |> 
  filter(labno %in% pansolid_worksheet_samples$labno) |> 
  arrange(surname, labno) |>  
  left_join(pansolid_worksheet_samples |> 
              select(-name), by = "labno") |> 
  rename(pansolid_worksheet = pcrid)

# Find EQA samples tested on PanSolid -----------------------------------------------

pansolid_eqa_samples <- pansolid_sample_info |> 
  filter(grepl(pattern  = "emqn|eqa|genqa",
               x = comments, ignore.case = TRUE)) 

eqa_samples_pansolid <- recent_samples |> 
  filter(labno %in% eqa_samples$labno &
           labno %in% pansolid_sample_list) |> 
  left_join(recent_results |> 
              select(-genodate), by = "labno") |> 
  filter(test == "NGS DNA QIAseq Pansolid V2")

# Find samples tested with PanSolid and MLPA ----------------------------------------

all_tests_on_pansolid_samples <- results_tbl |> 
  filter(labno %in% samples_on_ps_worksheets) |> 
  select(labno, test, genotype) |> 
  collect()

pansolid_sample_mlpa_results <- all_tests_on_pansolid_samples |> 
  filter(grepl(pattern = "P\\d{2,3}|MLPA",
               x = test, ignore.case = TRUE) &
           genotype != "Fail") |> 
  distinct()
