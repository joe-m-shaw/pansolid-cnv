# Packages --------------------------------------------------------------------------

library(shiny)
library(tidyverse)
library(here)

# Get CNV sample list for Leigh for her ThermoFisher project

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(here)

# Filepath --------------------------------------------------------------------------

source(here("scripts/set_shared_drive_filepath.R"))

# Load Data -------------------------------------------------------------------------

std_dev_results <- read_csv(str_c(data_folder, "live_service/collated/",
                                  "/live_service_std_dev_results_collated.csv"),
                            show_col_types = FALSE) |> 
  rename(noise = st_dev_signal_adjusted_log2_ratios)

pos_cnv_results <- read_csv(str_c(data_folder, "live_service/collated/",
                                  "/live_service_pos_cnv_results_collated.csv"),
                            show_col_types = FALSE)

percent_138_results <- read_csv(str_c(data_folder, "live_service/collated/",
                                      "/live_service_percent_138_results_collated.csv"),
                                show_col_types = FALSE)

# Join and filter data ---------------------------------------------------------------------------

pos_cnv_results_with_qc <- pos_cnv_results |> 
  left_join(std_dev_results |> 
              select(filepath, noise), by = "filepath",
            relationship = "many-to-one") |> 
  left_join(percent_138_results |> 
              select(filepath, percent_whole_panel_covered_at_138x), by = "filepath",
            relationship = "many-to-one") |> 
  filter(gene %in% c("EGFR", "ERBB2", "BRAF", "MYC", "MET") & !is.na(fold_change))

# Export results ---------------------------------------------------------------------------------

write.csv(pos_cnv_results_with_qc, file = str_c(outputs_folder, 
                                                     "pos_cnv_results_with_qc.csv"),
          row.names = FALSE)

