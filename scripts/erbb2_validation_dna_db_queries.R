## ERBB2 Validation: DNA Database Queries

# Packages --------------------------------------------------------------------------

library(here)

# Source functions ------------------------------------------------------------------

source(here::here("functions/dna_database_functions.R"))

source(here::here("functions/cnv_functions.R"))

# ERBB2 lab numbers -----------------------------------------------------------------

erbb2_labno_df <- read_csv(file = here::here("data/erbb2_validation_labnos_and_tissue_sources.csv"),
                           col_types = "c")

erbb2_labnos <- erbb2_labno_df$labno

# Worksheet details -----------------------------------------------------------------

worksheet_sample_list <- dna_db_pcr_records |> 
  select(pcrid, sample, name) |> 
  filter(sample %in% erbb2_labnos) |> 
  collect()

sample_worksheets <- unique(worksheet_sample_list$pcrid)

sample_worksheet_info <- dna_db_worksheets |> 
  filter(pcrid %in% sample_worksheets) |> 
  select(pcrid, date, description, disease, test_type) |> 
  collect()

pansolid_str_variants <- grep(pattern = "solid", x = sample_worksheet_info$description,
     ignore.case = TRUE, value = TRUE)

erbb2_pansolid_ws_details  <- sample_worksheet_info |> 
  filter(description %in% pansolid_str_variants & test_type == 19)

dna_db_export(erbb2_pansolid_ws_details)

# Extraction methods ----------------------------------------------------------------

erbb2_sample_extraction <- get_extraction_method(erbb2_labnos) |> 
  rename(extraction_method = method_name) |> 
  select(labno, extraction_method)

dna_db_export(erbb2_sample_extraction)

# Sample types ----------------------------------------------------------------------

erbb2_sample_types <- get_sample_tissue(erbb2_labnos) |> 
  select(labno, tissue_type)

if(length(setdiff(erbb2_sample_types$labno, erbb2_labnos)) != 0) {
  
  stop("Not all samples have sample types")

}

dna_db_export(erbb2_sample_types)

# Neoplastic cell content -----------------------------------------------------------

erbb2_ncc <- sample_tbl |> 
  select(labno, comments) |> 
  filter(labno %in% erbb2_labnos) |> 
  collect() |> 
  janitor::clean_names() |> 
  mutate(ncc_db = parse_ncc(comments))

dna_db_export(erbb2_ncc)

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
