
# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

# Functions and filepaths -----------------------------------------------------------

data_folder <- config::get("data_folderpath")

# S drive filepaths -----------------------------------------------------------------

pansolid_ws_df <- read_excel(path = paste0(data_folder, 
                              "live_service/",
                              "pansolid_html_worksheets.xlsx"),
                         sheet = "pansolid_html_worksheets")

jbrca_ws_df <- read_excel(path = paste0(data_folder, 
                                           "live_service/",
                                           "pansolid_html_worksheets.xlsx"),
                             sheet = "jbrca_html_worksheets")

# New files to check ------------------------------------------------------

pansolid_ws_to_check <- pansolid_ws_df |> 
  filter(htmls_checked == "no")

jbrca_ws_to_check <- jbrca_ws_df |> 
  filter(htmls_checked == "no")

ws_to_check_bind <- rbind(pansolid_ws_to_check,
                          jbrca_ws_to_check)

worksheet_list <- list(ws_to_check_bind$worksheet)

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

# Files already checked ---------------------------------------------------

pansolid_ws_checked <- pansolid_ws_df |> 
  filter(htmls_checked == "yes")

jbrca_ws_checked <- jbrca_ws_df |> 
  filter(htmls_checked == "yes")

ws_checked_bind <- rbind(pansolid_ws_checked,
                         jbrca_ws_checked)

checked_worksheet_list <- list(ws_checked_bind$worksheet)

length(checked_worksheet_list[[1]])

checked_filepaths <- checked_worksheet_list |> 
  map(\(checked_worksheet_list) get_html_filepaths(worksheet = checked_worksheet_list)) |> 
  flatten()

length(checked_filepaths)

# There are 16 samples so far with evidence of contamination

(16 / length(checked_filepaths)) * 100

