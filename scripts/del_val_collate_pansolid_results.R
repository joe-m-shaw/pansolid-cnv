# Collate deletions validation PanSolid results

# Packages ----------------------------------------------------------------

library(here)
library(tidyverse)
library(readxl)
library(janitor)

# Functions ---------------------------------------------------------------

source(here("functions/clc_raw_excel_functions.R"))

# Filepaths ---------------------------------------------------------------

data_folder <- config::get("data_filepath")

del_filepath <- paste0(data_folder, "validation/DOC6567_deletions/",
                       "raw/pansolid_ngs/raw_clc/")

# Load data ---------------------------------------------------------------

del_files <- list.files(path = del_filepath,
                        full.names = TRUE,
                        pattern  = "Results_TSG.*.xlsx")

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
