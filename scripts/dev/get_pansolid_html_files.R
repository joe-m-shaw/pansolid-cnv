
# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

# Functions and filepaths -----------------------------------------------------------

data_folder <- config::get("data_folderpath")

# S drive filepaths -----------------------------------------------------------------

html_ws_df <- read_csv(paste0(data_folder, 
                              "live_service/",
                              "pansolid_html_worksheets.csv"))

# Modify this list with the worksheets you want to check
worksheet_list <- html_ws_df$worksheet

get_html_filepaths <- function(
    repository_path = "S:/central shared/Genetics/Repository/WorksheetAnalysedData/",
    worksheet, full_names = TRUE) {
  
  html_regex <- "WS\\d{6}_\\d{8}_.*all_log2ratios.html"
  
  html_filepaths <- list.files(path = str_c(repository_path, {{ worksheet }},
                                                 "/"),
                                    recursive = TRUE, 
                                    pattern = html_regex,
                                    full.names = full_names)
  
  return(html_filepaths)
  
}

s_drive_filepaths <- worksheet_list |> 
  map(\(worksheet_list) get_html_filepaths(worksheet = worksheet_list)) |> 
  flatten()

s_drive_filepath_df <- tibble(
  "path" = unlist(s_drive_filepaths)) 

file.copy(from = s_drive_filepath_df$path,
          to = paste0(data_folder, "live_service/htmls/"))

