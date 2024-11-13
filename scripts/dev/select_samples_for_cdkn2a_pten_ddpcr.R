# Selecting samples for CDKN2A and PTEN ddPCR

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)

source(here("scripts/set_shared_drive_filepath.R"))
source(here("scripts/connect_to_dna_db.R"))
source(here("functions/join_dna_submission_sheets.R"))
       
# Load samples -----------------------------------------------------------------------------------

folder_path <- "S:/central shared/Genetics/Mol_Shared/Development.Team/Pan Solid CLC Somatic Amplifications Validation/PanSolid live service CNVs/"

pansolid_cnv_list <- read_excel(path = str_c(folder_path, "pansolid_cnv_list.xlsx")) |> 
  mutate(labno = as.character(labno))

pansolid_controls <- read_csv(file = paste0(data_folder, "controls/qiasymphony_controls.csv"))

pansolid_control_vector <- pansolid_controls$labno

pansolid_control_details <- sample_tbl |> 
  filter(labno %in% pansolid_control_vector) |> 
  select(labno, firstname, surname) |> 
  collect() |> 
  mutate(category = "normal control") |> 
  filter(!duplicated(labno))

ckn2a_pten_df <- pansolid_cnv_list |> 
  filter(predicted_gene %in% c("CDKN2A", "PTEN"))

cdkn2a_pten_labnos <- ckn2a_pten_df$labno

cdkn2a_pten_details <- sample_tbl |> 
  filter(labno %in% cdkn2a_pten_labnos) |> 
  select(labno, firstname, surname) |> 
  collect() |> 
  mutate(category = "patient") |> 
  filter(!duplicated(labno))

# DNA concentrations -----------------------------------------------------------------------------

pansolid_submission_sheet <- join_pansolid_submission_sheets() |> 
  filter(!duplicated(labno))

all_sample_df <- pansolid_control_details |> 
  rbind(cdkn2a_pten_details) |> 
  left_join(pansolid_submission_sheet |> 
              select(labno, stock_qubit),
            relationship = "one-to-one",
            by = "labno") |> 
  arrange(labno) |>

write.csv(all_sample_df, file = "all_sample_df.csv",
          row.names = FALSE)
