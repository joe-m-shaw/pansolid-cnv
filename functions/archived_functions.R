# Archived Functions

# The functions in this script are either no longer used or only apply to previous 
# output formats of the PanSolid analysis.

library(tidyverse)
library(readxl)
library(here)

coarse_tab <- "Oncogenes (Amplified) Coars..."
fine_tab <- "Oncogenes (Amplified) Fine-..."

read_clc_amp_calls <- function(file, input_sheet) {
  
  stopifnot(input_sheet %in% c(coarse_tab, fine_tab))
  
  x <- read_excel(path = file, sheet = input_sheet) 
  
  # Samples with no calls currently have a table with 0 rows
  if(nrow(x) == 0) {
    
    x <- tibble("Chromosome" = 0, 
                "Region" = "",
                "Name" = "", 
                "Region length" = 0, 
                "type" = "",
                "source" = "", 
                "ID"= "", 
                "GeneID"= "", 
                "HGNC"= "",
                "MIM"= "",
                "description"= "", 
                "gbkey"= "", 
                "gene"= "No calls",
                "gene_biotype"= "", 
                "gene_synonym"= "", 
                "CNV region"= "",
                "CNV region length"= 0,
                "Consequence"= "",
                "Fold-change (adjusted)"= 0,
                "p-value"= 0,
                "Number of targets"= 0,
                "Comments"= "", 
                "Targets"= "")
    
  }
  
  identifiers <- filename_to_df(file)
  
  results <- x |> 
    janitor::clean_names() |> 
    mutate(
      labno = as.character(parse_filename(file, 2)),
      setting = input_sheet,
      filename = file) |> 
    left_join(identifiers, by = "labno") |> 
    relocate(worksheet, labno, suffix, labno_suffix, patient_name,
             labno_suffix_worksheet, setting)
  
  return(results)
  
}


