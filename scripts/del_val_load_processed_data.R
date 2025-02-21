# Load Deletion Validation Processed Data

# Packages ---------------------------------------------------------------------

library(tidyverse)
library(here)

# Scripts and functions --------------------------------------------------------

del_val_processed_folder <- paste0(config::get("data_folderpath"),
                                   "validation/DOC6567_deletions/processed/")

# Load data --------------------------------------------------------------------

del_val_sample_patient_info <- read_csv(paste0(
  del_val_processed_folder, 
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

sample_id_col_list <- list("worksheet" = col_character(),
                           "labno" = col_character(),
                           "suffix" = col_character(),
                           "patient_name" = col_character(),
                           "labno_suffix" = col_character(),
                           "labno_suffix_worksheet" = col_character(),
                           "filepath" = col_character())

del_val_collated_138x <- read_csv(paste0(del_val_processed_folder,
                                         "del_val_collated_138x.csv"),
                                  col_types = 
                                    c(sample_id_col_list, 
                                      list("percent_138x" = col_double())))

del_val_collated_amp_genes <- read_csv(paste0(
  del_val_processed_folder,
  "del_val_collated_amp_genes.csv"),
  col_types = c(sample_id_col_list,
                list(
                 "gene" = col_character(),
                 "max_region_fold_change" = col_double(),
                 "min_region_fold_change" = col_double()
               )))

del_val_collated_del_genes <- read_csv(paste0(
  del_val_processed_folder,
  "del_val_collated_del_genes.csv"),
  col_types = c(sample_id_col_list,
                list(
                  "gene" = col_character(),
                  "max_region_fold_change" = col_double(),
                  "min_region_fold_change" = col_double()
                )))

del_val_collated_loh <- read_csv(paste0(
  del_val_processed_folder,
  "del_val_collated_loh.csv"),
  col_types = c(sample_id_col_list,
                list(
                  "chrom" = col_character(),
                  "gene" = col_character(),
                  "ploidy_state" = col_character(),
                  "loh_status"  = col_character(),
                  "no_targets_in_ploidy_region" = col_character()
              )))

del_val_collated_sig_cnvs <- read_csv(paste0(
  del_val_processed_folder,
  "del_val_collated_sig_cnvs.csv"),
  col_types = c(sample_id_col_list,
                list(
                  "gene" = col_character(),
                  "chromosome" = col_character(),
                  "cnv_co_ordinates" = col_character(),
                  "cnv_length" = col_double(),
                  "consequence" = col_character(),
                  "fold_change" = col_double(),
                  "p_value" = col_double(),
                  "no_targets" = col_double(),
                  "check_1" = col_character(),
                  "check_2" = col_character(),
                  "copy_number" = col_double(),
                  "start" = col_double(),
                  "end" = col_double()
                )))

del_val_collated_stdev <- read_csv(paste0(
  del_val_processed_folder,
  "del_val_collated_stdev.csv"),
  col_types = c(sample_id_col_list,
                list(
                  "stdev_noise" = col_double()
                )))
