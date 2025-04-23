# Collate ddPCR validation data (deletions validation)

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(here)

# Functions -------------------------------------------------------------------------

source(here("functions/ddpcr_functions.R"))

data_folder <- config::get("data_folderpath")

# ddPCR files -----------------------------------------------------------------------

ddpcr_validation_files <- list.files(
  path = paste0(data_folder, "/validation/DOC6567_deletions/raw/ddpcr/"),
  pattern = ".csv",
  full.names = TRUE)

if(length(ddpcr_validation_files) == 0) {
  stop("No ddPCR files in location")
}

# Collate results -------------------------------------------------------------------

ddpcr_validation_data <- ddpcr_validation_files |> 
  map(\(ddpcr_validation_files) read_biorad_ddpcr_csv(ddpcr_validation_files)) |> 
  list_rbind() |> 
  mutate(worksheet = str_extract(string = filepath,
                                 pattern = ".+(WS\\d{6}).+",
                                 group = 1),
         gene = str_extract(string = experiment,
                            pattern = "^(\\w{3,6})_Ex.+",
                            group = 1))

stopifnot(ncol(ddpcr_validation_data) == 66)

# Save collated results -------------------------------------------------------------

write.csv(x = ddpcr_validation_data,
          file = paste0(data_folder, 
                        "validation/DOC6567_deletions/processed/",
                        "del_val_ddpcr_collated.csv"),
          row.names = FALSE)
