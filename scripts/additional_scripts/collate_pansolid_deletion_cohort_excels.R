# Collate Deletion Cohort Results

# Packages --------------------------------------------------------------------------

library(here)
library(tidyverse)
library(readxl)

# Database connection ---------------------------------------------------------------

source(here("functions/dna_database_connection.R"))
source(here("functions/cnv_functions.R"))

# Functions -------------------------------------------------------------------------

filename_regex <- regex(
  r"[
  Results_SNVs_Indels-
  (WS\d{6})             # Worksheet
  _
  (\d{8})               # Lab number
  _
  (\D{5,30})            # Patient name
  (\d{10}_|_)           # Some names have numbers afterwards
  .*                    # Variable ending
  .xlsx
  ]",
  comments = TRUE
)

get_filename_identifiers <- function(filepath) {
  
  output <- data.frame(
    
    "worksheet" = str_extract(string = filepath, 
                                   pattern = filename_regex,
                                   group = 1),
    
    "labno" = str_extract(string = filepath, 
                               pattern = filename_regex,
                               group = 2),
    
    "patient_name" = str_extract(string = filepath, 
                                      pattern = filename_regex,
                                      group = 3),
    
    "filepath" = filepath) 
  
  return(output)
  
}

get_sheet_name <- function(excel_file, search_string) {
  
  sheet <- grep(x = excel_sheets(path = excel_file),
                pattern = search_string,
                value = TRUE)
  
  if(is.na(sheet)) {
    stop()
  }
  
  return(sheet)
  
}

make_table_for_empty_sheet <- function(table, sheet_string) {
  
  if(nrow(table) > 0 | ncol(table) > 0) {
    stop("Input table is not empty")
  }
  
  new_table <- data.frame(
    "chromosome" = 0,
    "region" = "0",
    "name" = str_c("No calls in ", sheet_string, " sheet"),
    "region_length" = 0,
    "type" = "0",
    "source" = "0",
    "id" = "0",
    "gene_id" = "0",
    "hgnc" = "0",
    "mim" = "0",
    "description" = "0",
    "gbkey" = "0",
    "gene" = "0",
    "gene_biotype" = "0",
    "gene_synonym" = "0",
    "cnv_region" = "0",
    "cnv_region_length" = 0,
    "consequence" = "0",         
    "fold_change_adjusted" = 0, 
    "p_value" = 0,            
    "number_of_targets" = 0,
    "targets" = "0")
  
  return(new_table)
  
}

read_clc_sheet <- function(excel_file, search_string) {
  
  sheet_name <- get_sheet_name(excel_file = {{ excel_file }},
                               search_string = {{ search_string }})
  
  sheet_table <- read_excel(path = excel_file, sheet = sheet_name) |> 
    janitor::clean_names()
  
  if(nrow(sheet_table) == 0) {
    
    sheet_table <- make_table_for_empty_sheet(sheet_table, search_string)
    
  }
  
  if("comments" %in% colnames(sheet_table)) {
    
    sheet_table <- select(sheet_table, -comments)
    
  }
  
  sheet_table <- sheet_table |> 
    mutate(filepath = excel_file)
  
  return(sheet_table)
  
}

# Collate data ----------------------------------------------------------------------

folder_path <- "S:/central shared/Genetics/NGS/Bioinformatics/1_Pan-solid-Cancer/CNV/Deletions/v1_Coarse_FC_1.33/"

pansolid_files <- list.files(path = folder_path, pattern = ".xlsx",
                             full.names = TRUE)

tsg_sheets_collated <- pansolid_files |> 
  map(\(pansolid_files) read_clc_sheet(excel_file = pansolid_files,
                                       search_string = "TSG")) |> 
  list_rbind()

amp_sheets_collated <- pansolid_files |> 
  map(\(pansolid_files) read_clc_sheet(excel_file = pansolid_files,
                                       search_string = "Oncogenes")) |> 
  list_rbind()

pansolid_results_collated <- extract_cnv_coordinates(
  rbind(tsg_sheets_collated, amp_sheets_collated),
  cnv_region)

# Get identifiers -------------------------------------------------------------------

pansolid_ids <- pansolid_files |> 
  map(\(pansolid_files) get_filename_identifiers(filepath = pansolid_files)) |> 
  list_rbind()

deletion_cohort_labnos <- unique(pansolid_ids$labno)

pansolid_nhs_path_nos <- sample_tbl |> 
  filter(labno %in% deletion_cohort_labnos) |> 
  select(labno, nhsno, pathno) |> 
  collect()

pansolid_all_ids <- pansolid_ids |> 
  left_join(pansolid_nhs_path_nos, by = "labno")

# Export results --------------------------------------------------------------------

write_csv(x = pansolid_all_ids,
          file = paste0(data_folder, "collated_validation_data/pansolid_ids.csv"))

write_csv(x = pansolid_results_collated,
          file = paste0(data_folder, "collated_validation_data/pansolid_results_collated.csv"))
