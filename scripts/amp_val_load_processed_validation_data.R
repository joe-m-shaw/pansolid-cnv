# Load PanSolid Processed Data

# Packages ---------------------------------------------------------------------------------------

library(tidyverse)
library(here)

# Scripts and functions --------------------------------------------------------------------------

source(here("functions/gene_table_functions.R"))

# Load data --------------------------------------------------------------------------------------

amp_genes <- load_pansolid_gene_table("Amplifications")

amp_validation_all_amp_gene_results_collated <- read_csv(paste0(config::get("data_folderpath"), 
                                                            "validation/",
                                                            "DOC6283_amplifications/",
                                                            "processed/",
                                                            "amp_validation_all_amp_gene_results_collated.csv"),
                                                     col_types = list(
                                                       "worksheet" = col_character(),
                                                       "labno" = col_character(),
                                                       "suffix" = col_character(),
                                                       "patient_name" = col_character(),
                                                       "labno_suffix" = col_character(),
                                                       "labno_suffix_worksheet" = col_character(),
                                                       "gene" = col_factor(levels = amp_genes$gene),
                                                       "max_region_fold_change" = col_double(),
                                                       "min_region_fold_change" = col_double(),
                                                       "pansolid_call" = col_character()
                                                     ))

amp_validation_pos_cnv_results_collated <- read_csv(paste0(config::get("data_folderpath"), 
                                                           "validation/",
                                                           "DOC6283_amplifications/",
                                                           "processed/",
                                                       "amp_validation_pos_cnv_results_collated.csv"),
                                                col_types = list(
                                                  "worksheet" = col_character(),
                                                  "labno" = col_character(),
                                                  "suffix" = col_character(),
                                                  "patient_name" = col_character(),
                                                  "labno_suffix" = col_character(),
                                                  "labno_suffix_worksheet" = col_character(),
                                                  "filepath" = col_character(),
                                                  "gene" = col_factor(levels = c(amp_genes$gene,
                                                                                 "no positive calls")),
                                                  "cnv_co_ordinates" = col_character(),
                                                  "cnv_length" = col_double(),
                                                  "consequence" = col_character(),
                                                  "fold_change" = col_double(),
                                                  "p_value" = col_double(),
                                                  "no_targets" = col_double(),
                                                  "start" = col_double(),
                                                  "end" = col_double()
                                                ))

amp_validation_stdev_results_collated <- read_csv(paste0(config::get("data_folderpath"), 
                                                         "validation/",
                                                         "DOC6283_amplifications/",
                                                         "processed/",
                                                     "amp_validation_stdev_results_collated.csv"),
                                              col_types = list(
                                                "worksheet" = col_character(),
                                                "labno" = col_character(),
                                                "suffix" = col_character(),
                                                "patient_name" = col_character(),
                                                "labno_suffix" = col_character(),
                                                "labno_suffix_worksheet" = col_character(),
                                                "filepath" = col_character(),
                                                "st_dev_signal_adjusted_log2_ratios" = col_double()))

amp_validation_percent_138_collated <- read_csv(paste0(config::get("data_folderpath"), 
                                                       "validation/",
                                                       "DOC6283_amplifications/",
                                                       "processed/",
                                                   "amp_validation_percent_138_collated.csv"),
                                            col_types = list(
                                              "worksheet" = col_character(),
                                              "labno" = col_character(),
                                              "suffix" = col_character(),
                                              "patient_name" = col_character(),
                                              "labno_suffix" = col_character(),
                                              "labno_suffix_worksheet" = col_character(),
                                              "filepath" = col_character(),
                                              "percent_whole_panel_covered_at_138x" = col_double()))

validation_sample_patient_info <- read_csv(paste0(config::get("data_folderpath"), 
                                                  "validation/",
                                                  "DOC6283_amplifications/",
                                                  "processed/",
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

validation_ddpcr_collated <- read_csv(paste0(config::get("data_folderpath"), 
                                             "validation/",
                                             "DOC6283_amplifications/",
                                             "processed/",
                                             "validation_ddpcr_collated.csv"),
                                      col_types = list(
                                        "sample" = col_character()
                                      ))


wgs_html_cnvs <- read_csv(paste0(config::get("data_folderpath"), 
                                 "validation/",
                                 "DOC6283_amplifications/",
                                 "processed/",
                                 "wgs_html_cnvs.csv"),
                          show_col_types = FALSE)

wgs_html_ids <- read_csv(paste0(config::get("data_folderpath"), 
                                "validation/",
                                "DOC6283_amplifications/",
                                "processed/",
                                 "wgs_html_ids.csv"),
                         col_types = list(
                           "nhsno" = col_character(),
                           "labno" = col_character()
                         ))

# Checks -----------------------------------------------------------------------------------------

if(length(dplyr::setdiff(amp_genes$gene, 
                  amp_validation_all_amp_gene_results_collated$gene)) != 0){
  
  stop("Gene names present which are not expected")
  
} else{
    
  message("All amps table contains expected genes")
  
  }
  
if(nrow(amp_validation_all_amp_gene_results_collated) != 
   (nrow(amp_validation_stdev_results_collated) * 9)) {
  
  stop("Different number of rows than expected in stdev and all_amp genes tables.")
  
}
