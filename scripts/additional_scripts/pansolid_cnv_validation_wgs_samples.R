# PanSolid CNV Validation WGS Samples

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

# Functions -------------------------------------------------------------------------

source(here("functions/dna_database_connection.R"))
source(here("functions/dna_database_functions.R"))
source(here("functions/cnv_functions.R"))

# Load WGS Results ------------------------------------------------------------------

# Samples from the BRAIN MATRIX study.

brain_matrix_htmls_df <- data.frame(
  
  "filepath" = list.files(path = here::here("data/brain_matrix_htmls/"),
                          full.names = TRUE,
                          pattern = "*.html")) |> 
  mutate(lab_sample_id = str_extract(string = filepath,
                                     pattern = "(\\d{10})_p",
                                     group = 1))
  
# Samples from standard clinical WGS pathway which also have a DNA extraction 
# on the QIAsymphony method.

clinical_wgs_htmls_df <- data.frame(
  
  "filepath" = list.files(path = here::here("data/wgs_result_htmls/"),
             full.names = TRUE,
             pattern = "*.html")) |> 
    mutate(lab_sample_id = str_extract(string = filepath,
                                       pattern = "(\\d{10})_p",
                                       group = 1)) |> 
  # One sample is on both lists
  filter(!lab_sample_id %in% brain_matrix_htmls_df$lab_sample_id)

# Collate CNV Results ---------------------------------------------------------------

brain_matrix_fp <- brain_matrix_htmls_df$filepath

brain_matrix_fp |> 
  map(\(brain_matrix_fp) parse_wgs_html_table(
    brain_matrix_fp, 
    div_id = "d_svcnv_tier1")) |> 
  list_rbind()



brain_matrix_cnvs

clinical_wgs_cnvs

joined_cnvs



  
  
  



