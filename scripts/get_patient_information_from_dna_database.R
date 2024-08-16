# Get Patient Information for PanSolid Gene Amplifications Cohort

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(here)

# Source scripts --------------------------------------------------------------------

source(here("scripts/set_shared_drive_filepath.R"))

source(here("functions/dna_database_connection.R"))

source(here("functions/dna_database_functions.R"))

# Get patient information -----------------------------------------------------------

validation_stdev_results_collated <- read_csv(paste0(data_folder, 
                                                     "validation/collated/",
                                                     "validation_stdev_results_collated.csv"))

sample_labnos <- validation_stdev_results_collated$labno

validation_sample_patient_info <- sample_tbl |> 
  filter(labno %in% sample_labnos) |> 
  select(labno, firstname, surname, nhsno, pathno,
         comments) |> 
  collect() |> 
  mutate(ncc = parse_ncc(comments))


# Extract cancer type from comments field -------------------------------------------

cancer_type_vector <- c("(?:O|o)varian",
                         "(?:C|c)olorectal",
                         "(?:L|l)ung",
                         "(?:B|b)ladder",
                         "CNS",
                         "(?:E|e)ndometrial",
                         "(?:M|m)elanoma",
                         "Oligodendroglioma",
                         "(?:P|p)rostate",
                         "Glioma",
                         "Pancreatic",
                         "Glioblastoma")

cancer_type_regex <- paste0("(",
                            paste(cancer_type_vector, collapse  = "|"),
                            ")")

validation_sample_patient_info_mod <- validation_sample_patient_info |> 
  mutate(cancer_comment = tolower(str_extract(pattern = cancer_type_regex,
                        string = comments,
                        group = 1)),
         
         cancer_group = case_when(
           
           cancer_comment %in% c("cns", "glioblastoma", "glioma",
                                 "oligodendroglioma") ~"central nervous system",
           
           TRUE ~cancer_comment))

cancer_type_summary <- validation_sample_patient_info_mod |> 
  filter(!is.na(cancer_group)) |> 
  count(cancer_group) |> 
  arrange(desc(n))


# Export results --------------------------------------------------------------------

write.csv(x = validation_sample_patient_info_mod,
          file = paste0(data_folder, "validation/collated/", 
                        "validation_sample_patient_info.csv"),
          row.names = FALSE)
