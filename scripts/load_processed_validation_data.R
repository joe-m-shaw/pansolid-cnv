# Load PanSolid Processed Data

library(tidyverse)

source(here("scripts/set_shared_drive_filepath.R"))

validation_all_amp_gene_results_collated <- read_csv(paste0(data_folder, 
                                                            "validation/processed/",
                                                            "validation_all_amp_gene_results_collated.csv"),
                                                     col_types = list(
                                                       "worksheet" = col_character(),
                                                       "labno" = col_character(),
                                                       "suffix" = col_character(),
                                                       "patient_name" = col_character(),
                                                       "labno_suffix" = col_character(),
                                                       "labno_suffix_worksheet" = col_character(),
                                                       "gene" = col_character(),
                                                       "max_region_fold_change" = col_double(),
                                                       "min_region_fold_change" = col_double(),
                                                       "pansolid_call" = col_character()
                                                     ))

validation_pos_cnv_results_collated <- read_csv(paste0(data_folder, 
                                                       "validation/processed/",
                                                       "validation_pos_cnv_results_collated.csv"),
                                                col_types = list(
                                                  "worksheet" = col_character(),
                                                  "labno" = col_character(),
                                                  "suffix" = col_character(),
                                                  "patient_name" = col_character(),
                                                  "labno_suffix" = col_character(),
                                                  "labno_suffix_worksheet" = col_character(),
                                                  "filepath" = col_character(),
                                                  "gene" = col_character(),
                                                  "cnv_co_ordinates" = col_character(),
                                                  "cnv_length" = col_double(),
                                                  "consequence" = col_character(),
                                                  "fold_change" = col_double(),
                                                  "p_value" = col_double(),
                                                  "no_targets" = col_double(),
                                                  "start" = col_double(),
                                                  "end" = col_double()
                                                ))

validation_stdev_results_collated <- read_csv(paste0(data_folder, 
                                                     "validation/processed/",
                                                     "validation_stdev_results_collated.csv"),
                                              col_types = list(
                                                "worksheet" = col_character(),
                                                "labno" = col_character(),
                                                "suffix" = col_character(),
                                                "patient_name" = col_character(),
                                                "labno_suffix" = col_character(),
                                                "labno_suffix_worksheet" = col_character(),
                                                "filepath" = col_character(),
                                                "st_dev_signal_adjusted_log2_ratios" = col_double()))

validation_percent_138_collated <- read_csv(paste0(data_folder, 
                                                   "validation/processed/",
                                                   "validation_percent_138_collated.csv"),
                                            col_types = list(
                                              "worksheet" = col_character(),
                                              "labno" = col_character(),
                                              "suffix" = col_character(),
                                              "patient_name" = col_character(),
                                              "labno_suffix" = col_character(),
                                              "labno_suffix_worksheet" = col_character(),
                                              "filepath" = col_character(),
                                              "percent_whole_panel_covered_at_138x" = col_double()))

validation_sample_patient_info <- read_csv(paste0(data_folder, 
                                                  "validation/processed/",
                                                  "validation_sample_patient_info.csv"),
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

validation_ddpcr_collated <- read_csv(paste0(data_folder, 
                                             "validation/processed/",
                                             "validation_ddpcr_collated.csv"),
                                      col_types = list(
                                        "sample" = col_character()
                                      ))


wgs_html_cnvs <- read_csv(paste0(data_folder,
                                 "validation/processed/",
                                 "wgs_html_cnvs.csv"),
                          show_col_types = FALSE)

wgs_html_ids <- read_csv(paste0(data_folder,
                                 "validation/processed/",
                                 "wgs_html_ids.csv"),
                         col_types = list(
                           "nhs_no_clean" = col_character(),
                           "labno" = col_character()
                         ))
