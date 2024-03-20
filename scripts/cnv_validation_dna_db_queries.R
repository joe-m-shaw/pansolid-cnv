## PanSolid CNV Validation: DNA Database Queries

# Packages --------------------------------------------------------------------------

library(here)

# Source functions ------------------------------------------------------------------

source(here::here("functions/dna_database_functions.R"))

source(here::here("functions/cnv_functions.R"))

# ERBB2 lab numbers -----------------------------------------------------------------

labno_df <- read_csv(file = here::here("data/cnv_validation_labnos_and_tissue_sources.csv"),
                           col_types = "c")

sample_labnos <- labno_df$labno

# Worksheet details -----------------------------------------------------------------

worksheet_sample_list <- dna_db_pcr_records |> 
  select(pcrid, sample, name) |> 
  filter(sample %in% sample_labnos) |> 
  collect()

sample_worksheets <- unique(worksheet_sample_list$pcrid)

sample_worksheet_info <- dna_db_worksheets |> 
  filter(pcrid %in% sample_worksheets) |> 
  select(pcrid, date, description, disease, test_type) |> 
  collect()

pansolid_str_variants <- grep(pattern = "pan", x = sample_worksheet_info$description,
     ignore.case = TRUE, value = TRUE)

pansolid_ws_details  <- sample_worksheet_info |> 
  filter(description %in% pansolid_str_variants & test_type == 19)

dna_db_export(pansolid_ws_details)

# Extraction methods ----------------------------------------------------------------

# These samples were extracted by Cobas and QiaSymphony on the same lab number
# Only the Cobas extractions were used for PanSolid testing

dual_extracted_samples <- c("23037135", "23037279")

sample_extraction_details <- get_extraction_method(sample_labnos) |> 
  rename(extraction_method = method_name) |> 
  select(labno, extraction_method) |> 
  filter(!(labno %in% dual_extracted_samples & extraction_method == "QIAsymphony_DNA_FFPE"))

dna_db_export(sample_extraction_details)

# Sample types ----------------------------------------------------------------------

sample_types <- get_sample_tissue(sample_labnos) |> 
  select(labno, tissue_type)

if(length(setdiff(sample_types$labno, sample_labnos)) != 0) {
  
  stop("Not all samples have sample types")

}

dna_db_export(sample_types)

# Sample gender ---------------------------------------------------------------------

sample_gender <- get_sample_gender(sample_labnos) |> 
  select(labno, gender_string)

if(length(setdiff(sample_gender$labno, sample_labnos)) != 0) {
  
  stop("Not all samples have genders")
  
}

dna_db_export(sample_gender)


# NHS number ------------------------------------------------------------------------

sample_nhs_no <- get_sample_nhs_no(sample_labnos)

dna_db_export(sample_nhs_no)

# Neoplastic cell content -----------------------------------------------------------

sample_ncc <- get_sample_ncc(sample_labnos)

dna_db_export(sample_ncc)

# Qiaseq core panel results ---------------------------------------------------------

core_result_info <- results_tbl |> 
  select(labno, test, testtype, genotype, genotype2, genocomm) |> 
  filter(labno %in% sample_labnos) |> 
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
