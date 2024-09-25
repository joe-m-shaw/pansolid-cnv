# Get ddPCR CNV results for Leigh

library(tidyverse)
library(here)

source(here("scripts/set_shared_drive_filepath.R"))
source(here("scripts/load_processed_validation_data.R"))
source(here("scripts/connect_to_dna_db.R"))
source(here("functions/dna_db_functions.R"))
source(here("functions/join_dna_submission_sheets.R"))

pansolid_submissions <- join_pansolid_submission_sheets() |> 
  filter(!duplicated(labno)) |> 
  select(labno, stock_qubit)

ddpcr_cnv_samples <- validation_ddpcr_collated |> 
  filter(worksheet != "WS138419" & target_type == "Ch1Unknown") |> 
  filter(grepl(pattern = "\\d{8}",
               x = sample)) |> 
  filter(!duplicated(sample)) |> 
  select(sample, worksheet, experiment, cnv) |> 
  mutate(ddpcr_result = case_when(
    cnv >= 5 ~"amplified",
    cnv < 5 ~"not amplified"
  )) |> 
  left_join(pansolid_submissions, join_by("sample" == "labno")) |> 
  rename(ddpcr_worksheet = worksheet,
         ddpcr_copy_number = cnv)

sample_labnos <- ddpcr_cnv_samples$sample

ddpcr_cnv_extractions <- get_extraction_method(sample_vector = sample_labnos) |> 
  select(labno, method_name) |> 
  rename(dna_extraction_method = method_name)

ddpcr_cnv_with_extractions <- ddpcr_cnv_samples |> 
  left_join(ddpcr_cnv_extractions, join_by("sample" == "labno"))

write.csv(ddpcr_cnv_with_extractions,
          file = paste0(outputs_folder, "ddpcr_cnvs_with_extraction_type.csv"),
          row.names = FALSE)
