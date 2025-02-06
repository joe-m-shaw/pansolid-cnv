# Get PanSolid QC metrics

library(here)
library(tidyverse)
library(readxl)

data_folder <- config::get("data_filepath")

pansolid_excels <- list.files(path = paste0(data_folder, 
                                            "live_service/raw/"),
                              pattern = ".*.xlsx",
                              full.names = TRUE)

read_qc_metrics <- function(file) {
  
  df <- read_excel(path = file, sheet = 1, skip = 4,
                   n_max = 5) |> 
    janitor::clean_names() |> 
    mutate(filepath = file,
           labno = str_extract(string = filepath,
                               pattern = "\\d{8}"),
           worksheet = str_extract(string = filepath,
                                   pattern = "WS\\d{6}")) |> 
    relocate(labno, worksheet)
  
  return(df)

}

# This took 12 minutes to run
pansolid_qc_tabs <- pansolid_excels |> 
  map(\(pansolid_excels) read_qc_metrics(pansolid_excels)) |>
  list_rbind()

pansolid_qc_tabs_mod <- pansolid_qc_tabs |> 
  select(labno, worksheet, filepath, summary_item, value) |> 
  pivot_wider(id_cols = c(labno, worksheet, filepath),
              names_from = summary_item,
              values_from = value) |> 
  janitor::clean_names()

write_csv(x = pansolid_qc_tabs_mod, 
          file = paste0(data_folder,
                      "live_service/collated/",
                      "live_service_qc_tabs_collated.csv"))

