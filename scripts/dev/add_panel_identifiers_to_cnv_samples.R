# Copy number variants seen on PanSolidv2

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

source(here("scripts/set_shared_drive_filepath.R"))

# Load Excel ------------------------------------------------------------------------

folder_path <- "S:/central shared/Genetics/Mol_Shared/Development.Team/Pan Solid CLC Somatic Amplifications Validation/PanSolid live service CNVs/"

pansolid_cnv_list <- read_excel(path = str_c(folder_path, "pansolid_cnv_list.xlsx")) |> 
  mutate(labno = as.character(labno))

# Summary ---------------------------------------------------------------------------

variant_summary <- pansolid_cnv_list |> 
  filter(variant_type %in% c("focal deletion", "focal amplification")) |> 
  group_by(variant_type, chromosome, chromosome_arm, predicted_gene) |> 
  summarise(total = n()) |> 
  arrange(desc(total))

variant_summary

# Add worksheet and panel -----------------------------------------------------------

worksheet_panel_info <- read_csv(paste0(data_folder, 
                                        "live_service_collated_data/pansolidv2_sample_worksheet_panel_information.csv"),
                                 col_types = c("ccccc")) |> 
  filter(!duplicated(labno)) 

pansolid_cnv_list_with_panels <- pansolid_cnv_list |> 
  left_join(worksheet_panel_info |> 
              select(labno, worksheet, panel), by = "labno",
            relationship = "many-to-one") |> 
  relocate(worksheet, panel) |> 
  arrange(worksheet)
  
setdiff(pansolid_cnv_list$labno, pansolid_cnv_list_with_panels$labno)

write.csv(pansolid_cnv_list_with_panels,
          str_c(folder_path, "pansolid_cnv_list_with_panels.csv"),
          row.names = FALSE)
