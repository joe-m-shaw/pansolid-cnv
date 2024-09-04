# Reformat PanSolid and WGS data ----------------------------------------------------

library(here)
source(here("functions/cnv_functions.R"))
source(here("scripts/set_shared_drive_filepath.R"))

# Load collated data ----------------------------------------------------------------

wgs_data_collated <- read_csv(file = paste0(data_folder, "validation/collated/",
                                                 "wgs_html_cnvs.csv"))

wgs_html_ids <- read_csv(file = paste0(data_folder, "validation/collated/",
                                       "wgs_html_ids.csv"),
                         col_types = list(
                           "filepath" = col_character(),
                           "wgs_r_no" = col_character(),
                           "wgs_p_no" = col_character(),
                           "wgs_version" = col_character(),
                           "wgs_analysis_date" = col_character(),
                           "patient_name" = col_character(),
                           "patient_dob" = col_character(),
                           "nhs_no" = col_character(),
                           "nhs_no_clean" = col_character(),
                           "labno" = col_character(),
                           "wgs_pathno" = col_character()
                         ))

# Reformat data ---------------------------------------------------------------------

wgs_htmls <- list.files(path = paste0(data_folder, "validation/raw/wgs/"),
                        full.names = TRUE,
                        pattern = "*.html")

wgs_del_data_reformatted <- wgs_htmls |> 
  map(\(wgs_htmls) reformat_wgs_cnv_result(filepath = wgs_htmls,
                                           cnv_type = "Deletions")) |> 
  list_rbind()

wgs_amp_data_reformatted <- wgs_htmls |> 
  map(\(wgs_htmls) reformat_wgs_cnv_result(filepath = wgs_htmls,
                                           cnv_type = "Amplifications")) |> 
  list_rbind()

wgs_data_reformatted_collated <- rbind(wgs_del_data_reformatted, 
                                       wgs_amp_data_reformatted)

# Add identifiers -------------------------------------------------------------------

wgs_reformatted_with_ids <- wgs_html_ids |> 
  left_join(wgs_data_reformatted_collated, by = "filepath",
            relationship = "one-to-many") |> 
  mutate(nhs_no_clean = as.character(nhs_no_clean),
         labno = as.character(labno))
  
