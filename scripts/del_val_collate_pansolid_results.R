# Collate deletions validation PanSolid results

# Packages ----------------------------------------------------------------

library(here)
library(tidyverse)
library(readxl)
library(janitor)

# Functions ---------------------------------------------------------------

source(here("functions/clc_raw_excel_functions.R"))
source(here("functions/pansolid_excel_functions.R"))

# Filepaths ---------------------------------------------------------------

data_folder <- config::get("data_filepath")

bio_cnv_folder <- "S:/central shared/Genetics/NGS/Bioinformatics/1_Pan-solid-Cancer/CNV/Deletions/"

# Load data ---------------------------------------------------------------

del_files <- list.files(path = bio_cnv_folder,
                        recursive = FALSE,
                        full.names = TRUE,
                        pattern  = "TSG\\s\\(Deleted\\).*.xlsx")

all_coarse_targets <- del_files |> 
  map(\(del_files) read_all_clc_targets(del_files, "Coarse Targets")) |> 
  list_rbind()

if(length(del_files) * 6174 != nrow(all_coarse_targets)){
  stop("Collated data has different number of rows than predicted")
}

all_fine_targets <- del_files |> 
  map(\(del_files) read_all_clc_targets(del_files, "Fine Targets")) |> 
  list_rbind()

if(length(del_files) * 6174 != nrow(all_fine_targets)){
  stop("Collated data has different number of rows than predicted")
}

all_targets <- rbind(all_coarse_targets, all_fine_targets)

stopifnot(nrow(all_targets) == ((6174*2)*length(del_files)))

stopifnot(anyNA.data.frame(x = all_targets |> 
                   select(-c(regional_consequence, 
                             regional_effect_size, comments))) == FALSE)

write_csv(all_targets, file = paste0(data_folder,
                                     "validation/DOC6567_deletions/",
                                     "processed/",
                                     "del_val_pansolid_targets_collated.csv"))

coarse_df <- del_files |> 
  map(\(del_files) read_del_raw_excel(filepath = del_files,
                                    sheet_no = 1)) |> 
  list_rbind() |> 
  mutate(graining = "coarse") |> 
  relocate(graining)

fine_df <- del_files |> 
  map(\(del_files) read_del_raw_excel(filepath = del_files,
                                      sheet_no = 2)) |> 
  list_rbind() |> 
  mutate(graining = "fine") |> 
  relocate(graining)

pansolid_ngs_cnvs <- rbind(coarse_df, fine_df)

write_csv(pansolid_ngs_cnvs, file = paste0(data_folder,
                                           "validation/DOC6567_deletions/",
                                           "processed/",
                                           "del_val_pansolid_ngs_collated.csv"))
