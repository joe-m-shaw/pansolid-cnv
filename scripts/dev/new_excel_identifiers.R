# Get new Excel file identifiers to check

library(tidyverse)

source(here::here("functions/filename_functions.R"))

ws_files <- list.files(path = paste0(
  "S:/central shared/Genetics/Repository/WorksheetAnalysedData/",
  "WS152872"),
  pattern = "Annotated.*.xlsx",
  recursive = TRUE,
  full.names = TRUE)

parse_filename_panel <- function(filename) {
  
  output <- data.frame(
    worksheet = c(parse_filename(filename, 1)),
    labno = c(as.character(parse_filename(filename, 2))),
    patient_name = c(parse_filename(filename, 4))) |> 
    dplyr::mutate(
      panel = str_extract(string = filename,
                          pattern = "Annotated_(v2.*_PS)_WS.*.xlsx",
                          group = 1)) |> 
    relocate(panel, .after = worksheet)
  
  return(output)

}

ws_filename_df <- ws_files |> 
  map(\(ws_files) parse_filename_panel(ws_files)) |> 
  list_rbind()

write_csv(ws_filename_df,
          file = paste0(config::get("data_folderpath"),
                        "live_service/new_files_to_check.csv"))
        