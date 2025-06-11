# Load PanSolid Live Service Data

library(tidyverse)
library(here)

data_folder <- config::get("data_filepath")

live_service_amp_gene_results_collated <- read_csv(paste0(data_folder, 
                                                            "live_service/collated/",
                                                            "live_service_amp_gene_results_collated.csv"),
                                                     col_types = list(
                                                       "worksheet" = col_character(),
                                                       "labno" = col_character(),
                                                       "suffix" = col_character(),
                                                       "patient_name" = col_character(),
                                                       "labno_suffix" = col_character(),
                                                       "labno_suffix_worksheet" = col_character(),
                                                       "gene" = col_character(),
                                                       "max_region_fold_change" = col_double(),
                                                       "min_region_fold_change" = col_double()
                                                     ))

live_service_pos_cnv_results_collated <- read_csv(paste0(data_folder, 
                                                         "live_service/collated/",
                                                       "live_service_pos_cnv_results_collated.csv"),
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

live_service_std_dev_results_collated <- read_csv(paste0(data_folder, 
                                                     "live_service/collated/",
                                                     "live_service_std_dev_results_collated.csv"),
                                              col_types = list(
                                                "worksheet" = col_character(),
                                                "labno" = col_character(),
                                                "suffix" = col_character(),
                                                "patient_name" = col_character(),
                                                "labno_suffix" = col_character(),
                                                "labno_suffix_worksheet" = col_character(),
                                                "filepath" = col_character(),
                                                "st_dev_signal_adjusted_log2_ratios" = col_double()))

live_service_percent_138_results_collated <- read_csv(paste0(data_folder, 
                                                   "live_service/collated/",
                                                   "live_service_percent_138_results_collated.csv"),
                                            col_types = list(
                                              "worksheet" = col_character(),
                                              "labno" = col_character(),
                                              "suffix" = col_character(),
                                              "patient_name" = col_character(),
                                              "labno_suffix" = col_character(),
                                              "labno_suffix_worksheet" = col_character(),
                                              "filepath" = col_character(),
                                              "percent_whole_panel_covered_at_138x" = col_double()))

panel_info <- read_csv(paste0(data_folder, 
                              "live_service/collated/",
                             "/pansolidv2_sample_worksheet_panel_information.csv"),
                       col_types = "ccccc") 

