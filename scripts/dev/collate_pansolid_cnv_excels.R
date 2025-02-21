# Collating new PanSolid CNV Excel layout

# Objective: collate information from 5 tables in the "Amplifications" tab
# of the new PanSolid CNV Excel layout, and annotate with sample identifiers.

# Packages ----------------------------------------------------------------

library(here)
library(tidyverse)
library(purrr)

# Functions ---------------------------------------------------------------

source(here::here("functions/pansolid_cnv_excel_functions.R"))

# Read data ---------------------------------------------------------------

data_path <- paste0("S:/central shared/Genetics/NGS/Bioinformatics/",
                    "1_Pan-solid-Cancer/CNV/Deletions/",
                    "Jan2025_withLOH_testing_del_visualisation/")

cnv_excel_files <- list.files(path = data_path, 
                         pattern = "Annotated_.*.xlsx",
                         full.names = TRUE)

file_cnv_tbl_list <- cnv_excel_files |> 
  map(\(cnv_excel_files) extract_cnv_tbls(cnv_excel_files))

# Isolate collated tables -------------------------------------------------

collated_stdev <- map(file_cnv_tbl_list, ~ .x[["stdev"]]) |> 
  list_rbind()

collated_138x <- map(file_cnv_tbl_list, ~ .x[["percent_138x"]]) |> 
  list_rbind()

collated_sig_cnvs <- map(file_cnv_tbl_list, ~ .x[["sig_cnvs"]]) |> 
  list_rbind()

collated_amp_genes <- map(file_cnv_tbl_list, ~ .x[["amp_genes"]]) |> 
  list_rbind()

collated_del_genes <- map(file_cnv_tbl_list, ~ .x[["del_genes"]]) |> 
  list_rbind()
