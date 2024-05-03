# Update PanSolidv2 Copy Number Variant Data

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

# Functions -------------------------------------------------------------------------

source(here::here("functions/cnv_functions.R"))

# S drive filepaths -----------------------------------------------------------------

pansolidv2_worksheets <- read_excel(here::here("data/pansolid_live_service_worksheets.xlsx"))

worksheet_list <- list(pansolidv2_worksheets$worksheet)

s_drive_filepaths <- worksheet_list |> 
  map(\(worksheet_list) get_annotated_filepaths(worksheet_list)) |> 
  flatten()

s_drive_file_df <- tibble(
  filepath = unlist(s_drive_filepaths)) |> 
  mutate(filename = str_extract(string = filepath, 
                                pattern = "Annotated_WS\\d{6}_.+.xlsx"))

# Local drive filepaths -------------------------------------------------------------

local_drive_file_df <- tibble(
  filepath = unlist(list.files(here::here("data/live_service_annotated_files/"),
                               full.names = TRUE))) |> 
  mutate(filename = str_extract(string = filepath, 
                                pattern = "Annotated_WS\\d{6}_.+.xlsx"))

local_filepaths <- list(local_drive_file_df$filepath) |> 
  flatten()

# Identify and copy new files -------------------------------------------------------

new_files <- s_drive_file_df |> 
  filter(!filename %in% local_drive_file_df$filename)

if (nrow(new_files) > 0) {
  
  file.copy(from = new_files$filepath, 
            to = here::here("data/live_service_annotated_files/"))
  
}

# Get new file local filepaths ------------------------------------------------------

local_drive_file_df <- tibble(
  filepath = unlist(list.files(here::here("data/live_service_annotated_files/"),
                               full.names = TRUE))) |> 
  mutate(filename = str_extract(string = filepath, 
                                pattern = "Annotated_WS\\d{6}_.+.xlsx"))

new_file_local_paths <- local_drive_file_df |> 
  filter(filename %in% new_files$filename)
