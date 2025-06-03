# EGFR version 3 

library(tidyverse)

source(here::here("scripts/del_val_load_processed_data.R"))

egfr_cohort <- c(
  "23034780",  # EGFR amp with v3
  "24027851",  # EGFR amp with deletion exon 25
  "24028153", # EGFR amp v3
  "24034555", # EGFR amp v3
  "24036195", # EGFR amp v3
  "23024556", # EGFR amp without v3
  "23036271", # EGFR amp without v3
  "24017321", # No EGFR amplification
  "24036644" # No EGFR amplification
)

egfr_wgs_ids <- del_val_wgs_html_ids |> 
  filter(labno %in% egfr_cohort) |> 
  select(filepath, wgs_r_no, wgs_p_no,
         wgs_pathno, nhsno,
         patient_name, labno) |> 
  rename(wgs_report_filepath = filepath,
         wgs_labno = labno)

egfr_pansolid_ids <- del_val_collated_stdev |> 
  left_join(del_val_sample_patient_info |> 
              select(labno, nhsno, extraction_method, pathno),
            by = "labno") |>  
  select(labno, suffix, worksheet, nhsno, extraction_method, 
         stdev_noise, pathno) |> 
  rename(pansolid_pathno = pathno,
         pansolid_labno = labno,
         pansolid_worksheet = worksheet)

egfr_joined_ids <- egfr_wgs_ids |> 
  left_join(egfr_pansolid_ids, 
            by = "nhsno") |> 
  filter(extraction_method %in% c("Fresh tissue", "QIAsymphony")) |> 
  filter(is.na(suffix)) |> 
  filter(!pansolid_labno %in% c("24030945", "23024556", "24017321")) |> 
  filter(!duplicated(pansolid_labno)) |> 
  select(pansolid_labno, pansolid_worksheet,
         wgs_labno, pansolid_pathno, wgs_pathno,
         patient_name, nhsno,
         wgs_r_no, wgs_p_no,
         extraction_method,
         wgs_report_filepath)

write_csv(egfr_joined_ids,
          "egfr_joined_ids.csv")


