# Load Deletion Validation Processed Data

# Packages ---------------------------------------------------------------------------------------

library(tidyverse)
library(here)

# Scripts and functions --------------------------------------------------------------------------

source(here("functions/gene_table_functions.R"))

del_val_processed_folder <- config::get("del_val_processed_folder")

# Load data --------------------------------------------------------------------------------------

del_val_pansolid_ngs_collated <- read_csv(paste0(del_val_processed_folder,
                                               "del_val_pansolid_ngs_collated.csv"),
                                                col_types = list(
                                                  "graining" = col_character(),
                                                  "labno" = col_character(),
                                                  "worksheet" = col_character(),
                                                  "chromosome" = col_character(),
                                                  "region" = col_character(),
                                                  "name" = col_character(),
                                                  "region_length" = col_double(),
                                                  "type" = col_character(),
                                                  "source" = col_character(),
                                                  "id" = col_character(),
                                                  "gene_id" = col_character(),
                                                  "hgnc" = col_character(),
                                                  "mim" = col_character(),
                                                  "description" = col_character(),
                                                  "gbkey" = col_character(),
                                                  "gene" = col_character(),
                                                  "gene_biotype" = col_character(),
                                                  "gene_synonym" = col_character(),
                                                  "cnv_region" = col_character(),
                                                  "cnv_region_length" = col_double(),
                                                  "consequence" = col_character(),
                                                  "fold_change_adjusted" = col_double(),
                                                  "p_value" = col_double(),
                                                  "number_of_targets" = col_double(),
                                                  "comments" = col_character(),
                                                  "targets" = col_character(),
                                                  "filepath" = col_character(),
                                                  "suffix" = col_character()))

del_val_sample_patient_info <- read_csv(paste0(del_val_processed_folder, 
                                                  "del_val_sample_patient_info.csv"),
                                           col_types = list(
                                             "labno" = col_character(),
                                             "firstname" = col_character(),
                                             "surname" = col_character(),
                                             "date_in" = col_datetime(format = "%Y-%m-%d %H:%M:%OS"),
                                             "nhsno" = col_character(),
                                             "pathno" = col_character(),
                                             "comments" = col_character(),
                                             "ncc" = col_character(),
                                             "tissue_type" = col_character(),
                                             "cancer_comment" = col_character(),
                                             "cancer_group" = col_character(),
                                             "method_name" = col_character()
                                           ))

del_val_ddpcr_collated <- read_csv(paste0(del_val_processed_folder, 
                                             "del_val_ddpcr_collated.csv"),
                                      col_types = list(
                                        "sample" = col_character()
                                      ))


del_val_wgs_html_cnvs <- read_csv(paste0(del_val_processed_folder,
                                 "del_val_wgs_html_cnvs.csv"),
                          show_col_types = FALSE)

del_val_wgs_html_ids <- read_csv(paste0(del_val_processed_folder,
                                "del_val_wgs_html_ids.csv"),
                         col_types = list(
                           "nhsno" = col_character(),
                           "labno" = col_character()
                         ))
