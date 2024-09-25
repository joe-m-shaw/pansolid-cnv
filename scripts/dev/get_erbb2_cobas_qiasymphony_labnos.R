# Get lab numbers for ERBB2 Cobas samples

library(tidyverse)
library(here)

source(here("scripts/set_shared_drive_filepath.R"))
source(here("scripts/connect_to_dna_db.R"))
source(here("functions/dna_db_functions.R"))

erbb2_pansolid_labnos <- c("24007929",
                           "24007945",
                           "24007947",
                           "24007949",
                           "24007951",
                           "24007953",
                           "24008372",
                           "24008849",
                           "24008853",
                           "24009050",
                           "24010112")

erbb2_pansolid_nhsno_df <- sample_tbl |> 
  filter(labno %in% erbb2_pansolid_labnos) |> 
  select(labno, nhsno) |> 
  collect()

# 21015264 is 24007951
# 24008372 is 23019668

erbb2_pansolid_nhsnos <- erbb2_pansolid_nhsno_df$nhsno

erbb2_nhsno_df <- sample_tbl |> 
  filter(nhsno %in% erbb2_pansolid_nhsnos) |> 
  select(labno, nhsno) |> 
  collect()

erbb2_extractions <- get_extraction_method(sample_vector = erbb2_nhsno_df$labno) |> 
  filter(method_name %in% c("QIAsymphony_DNA_FFPE", "COBAS")) |> 
  select(labno, method_name) |>
  filter(!labno %in% c("22024049", "22029036")) |> 
  left_join(erbb2_nhsno_df, by = "labno") |> 
  pivot_wider(id_cols = nhsno,
              names_from = method_name,
              names_prefix = "labno_",
              values_from = labno)

write.csv(x = erbb2_extractions, file = str_c(data_folder, "erbb2_pansolid_cobas_labnos.csv"))
