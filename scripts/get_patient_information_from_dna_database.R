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

sample_labnos <- unique(validation_stdev_results_collated$labno)

patient_info <- sample_tbl |> 
  filter(labno %in% sample_labnos) |> 
  select(labno, firstname, surname, date_in, tissue, nhsno, pathno,
         comments) |> 
  collect() 

# Add NCC and tissue type -----------------------------------------------------------

patient_info_ncc <- patient_info |> 
  mutate(ncc = parse_ncc(comments),
         tissue = as.numeric(tissue)) |> 
  left_join(tissue_types |> 
              select(tissue_type_id, tissue_type),
            join_by("tissue" == "tissue_type_id"))

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

patient_info_cancer_type <- patient_info_ncc |> 
  mutate(cancer_comment = tolower(str_extract(pattern = cancer_type_regex,
                        string = comments,
                        group = 1)),
         
         cancer_group = case_when(
           
           cancer_comment %in% c("cns", "glioblastoma", "glioma",
                                 "oligodendroglioma") ~"central nervous system",
           
           TRUE ~cancer_comment))

# Add DNA extraction method ---------------------------------------------------------

extraction_methods <- get_extraction_method(sample_vector = sample_labnos) |> 
  select(labno, method_name) |> 
  # Some samples were on a QIAsymphony run with an instrument failure
  # so they have multiple extraction entries for the same lab number
  filter(!duplicated(labno))

patient_info_extraction_method <- patient_info_cancer_type |> 
  left_join(extraction_methods, by = "labno")

# Export results --------------------------------------------------------------------

validation_sample_patient_info <- patient_info_extraction_method |> 
  select(-c(tissue))

write.csv(x = validation_sample_patient_info,
          file = paste0(data_folder, "validation/collated/", 
                        "validation_sample_patient_info.csv"),
          row.names = FALSE)
