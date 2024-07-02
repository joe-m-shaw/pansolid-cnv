# Get Pathology Lab Details for Louise

library(tidyverse)
library(here)

source("functions/dna_database_connection.R")

source("functions/dna_database_functions.R")

input_samples <- read_csv(file = here::here("data/2024-05-07_louise_path_lab_samples.csv"),
                          col_types = "c") |> 
  filter(!duplicated(labno))

samples_to_get <- input_samples$labno

length(samples_to_get)

sample_details <- sample_tbl |> 
  filter(labno %in% samples_to_get) |> 
  select(labno, consultant, consultant_address) |> 
  collect()

setdiff(samples_to_get, sample_details$labno)

write.csv(sample_details, here::here("data/sample_details.csv"),
          row.names = FALSE)
