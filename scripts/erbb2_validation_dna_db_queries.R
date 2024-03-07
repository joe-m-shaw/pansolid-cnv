## ERBB2 Validation: DNA Database Queries

# Packages --------------------------------------------------------------------------

library(here)

# Source functions ------------------------------------------------------------------

source(here::here("functions/dna_database_functions.R"))

source(here::here("functions/cnv_functions.R"))

# ERBB2 lab numbers -----------------------------------------------------------------

erbb2_labno_df <- read_csv(file = here::here("data/erbb2_validation_labnos.csv"),
                           col_types = "c")

erbb2_labnos <- erbb2_labno_df$labno

erbb2_ws_df <- read_csv(file = here::here("data/erbb2_validation_pansolid_worksheets.csv"),
                                          col_types = "c") |> 
  mutate(pcr_id = parse_number(x  = worksheet))

erbb2_ws <- erbb2_ws_df$pcr_id

# Worksheet details -----------------------------------------------------------------

erbb2_pansolid_ws_details <- dna_db_worksheets |> 
  filter(pcrid %in% erbb2_ws) |> 
  select(pcrid, date) |> 
  collect() |> 
  janitor::clean_names()

dna_db_export(erbb2_pansolid_ws_details)

# Extraction methods ----------------------------------------------------------------

erbb2_sample_extraction <- get_extraction_method(erbb2_labnos) |> 
  rename(extraction_method = method_name) |> 
  select(labno, extraction_method)

dna_db_export(erbb2_sample_extraction)

# Sample types ----------------------------------------------------------------------

erbb2_sample_types <- get_sample_tissue(erbb2_labnos) |> 
  select(labno, tissue_type)

dna_db_export(erbb2_sample_types)

# Neoplastic cell content -----------------------------------------------------------

erbb2_ncc <- sample_tbl |> 
  select(labno, comments) |> 
  filter(labno %in% erbb2_labnos) |> 
  collect() |> 
  janitor::clean_names() |> 
  mutate(ncc_db = parse_ncc(comments))

dna_db_export(erbb2_ncc)

# Tumour source ---------------------------------------------------------------------

erbb2_discodes <- sample_tbl |> 
  select(labno, disease, disease_2, disease_3, 
         disease_4) |> 
  filter(labno %in% erbb2_labnos) |> 
  collect() |> 
  janitor::clean_names()

erbb2_discodes_long <- erbb2_discodes |> 
  pivot_longer(cols = -c(labno),
               values_to = "discode") |> 
  left_join(discode |> 
              select(discode, disease), by = "discode") |> 
  filter(!is.na(discode))

erbb2_tumour_sources <- erbb2_discodes_long |> 
  
  # 168: Fusion NGS RNA Panel
  # 215: DNA Store
  filter(!discode %in% c(168, 215)) |>
  
  mutate(tumour_source = case_when(
    
    discode == 123 ~"Lung",
    
    discode %in% c(221, 222) ~"Central nervous system",
    
    discode %in% c(120, 209, 120,
                   174, 161) ~"Colorectal"
  ),
  
  tumour_source = case_when(
    
    labno == 21019092 ~"Lung",
    
    labno == 21011525 ~"Uterus",
    
    TRUE ~tumour_source)
  
  ) |> 
  filter(!duplicated(labno)) |> 
  select(-name)
  
dna_db_export(erbb2_tumour_sources)

# Qiaseq core panel results ---------------------------------------------------------

core_result_info <- results_tbl |> 
  select(labno, test, testtype, genotype, genotype2, genocomm) |> 
  filter(labno %in% erbb2_labnos) |> 
  collect() |> 
  janitor::clean_names() |> 
  filter(test %in% grep(pattern = "Q.{2,4}seq\\s(Core|NGS\\sCore|Merged)", 
                        x = test, 
                        ignore.case = TRUE,
                        value = TRUE)) |> 
  mutate(genotype = case_when(
    
    # Sample has "EGFR" instead of "ERBB2" written on DNA Database - confirmed on report
    labno == 23022389 ~"ERBB2 amplification detected (Mean DQ 25x)",
    
    # Sample has "ERRB2" instead of "ERBB2" written
    labno == 21015264 ~"No mutation identified; ERBB2 amplification detected (mean DQ 60.41x)",
    
    TRUE ~genotype)) |> 
  
  filter(!duplicated(labno))

dna_db_export(core_result_info)
