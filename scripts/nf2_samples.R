# NF2 Loss of Heterozygosity testing

# Eleanor's idea is to find samples which were tested using PanSolid and microsatellite
# testing for NF2 loss of heterozygosity.

library(tidyverse)

source("functions/dna_database_connection.R")
source("functions/dna_database_functions.R")

# PanSolidv2 was introduced in April 2024

recent_results <- results_tbl |> 
  filter(genodate > "2024-03-01 00:00:00") |> 
  select(labno, surname, test, exon, genotype, genotype2,
         genocomm) |> 
  collect() |> 
  filter(!labno %in% c("water", "Water", "WATER"))

nf2_loh_results <- recent_results |> 
  filter(grepl(pattern = "NF2 LOH", x = test, ignore.case = TRUE))

pansolid_results <- recent_results |> 
  filter(grepl(pattern = "pansol", x = test, 
               ignore.case = TRUE))

double_tested_samples <- intersect(nf2_loh_results$labno, pansolid_results$labno)

double_tested_extraction <- get_extraction_method(sample_vector = double_tested_samples) |> 
  select(labno, method_name)

nf2_loh_to_check <- nf2_loh_results |> 
  left_join(pansolid_results |> 
              select(labno, test), by = "labno") |> 
  left_join(double_tested_extraction, by = "labno") |> 
  filter(method_name == "QIAsymphony_DNA_FFPE") |> 
  select(-c(`test.x`, `test.y`)) |> 
  arrange(surname, exon)
