# Organising PanSolid CNV Validation Testing

library(tidyverse)

source("functions/dna_database_connection.R")
source("functions/dna_database_functions.R")
source("functions/cnv_functions.R")

dev_samples <- sample_tbl |> 
  filter(disease == "293") |> 
  select(labno, firstname, surname, nhsno, date_in, pathno, concentration) |> 
  collect() |> 
  arrange(desc(date_in))

to_export <- dev_samples |> 
  mutate(sample_name = paste0(firstname, " ", surname),
         panel = "v2PANSOLID",
         enrichment = "PANSOLIDV2") |> 
  select(labno, sample_name, panel, enrichment, concentration)

write.csv(to_export, here("outputs/to_export.csv"),
          row.names = FALSE)

sample_tbl |> 
  filter(nhsno == "4426940699") |> 
  select(labno, firstname, surname, nhsno, date_in, pathno, concentration) |> 
  collect() 
