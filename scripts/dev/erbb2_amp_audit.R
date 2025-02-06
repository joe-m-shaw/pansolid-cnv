# ERBB2 Lung Cases

# 19/12/2024: Helene would like to know how many ERBB2 amplifications
# were detected on the lung panel.

library(tidyverse)
library(here)

source(here("scripts/load_processed_live_service_data.R"))
source(here("functions/pansolid_excel_functions.R"))
source(here("functions/join_dna_submission_sheets.R"))

pansolid_submission_2024 <- join_pansolid_submission_sheets() |> 
  filter(submission_sheet == "2024")

erbb2_amps <- live_service_pos_cnv_results_collated |> 
  filter(gene == "ERBB2" & !is.na(fold_change)) |> 
  mutate(worksheet = parse_filename(filepath, 1),
         labno = parse_filename(filepath, 2),
         suffix = parse_filename(filepath, 3)) |> 
  left_join(live_service_percent_138_results_collated |> 
              select(filepath, percent_whole_panel_covered_at_138x),
            by = "filepath") |> 
  left_join(live_service_std_dev_results_collated |> 
            select(filepath, st_dev_signal_adjusted_log2_ratios),
            by = "filepath") |> 
  # Take panel from submission sheets
  # Not all filenames contain panel information
  left_join(pansolid_submission_2024, by = "labno",
            relationship = "many-to-many") |> 
  filter(!duplicated(labno)) |> 
  select(labno, suffix, worksheet, panel, gene, fold_change, 
         st_dev_signal_adjusted_log2_ratios, 
         percent_whole_panel_covered_at_138x)

erbb2_amps_summary <- erbb2_amps |> 
  count(panel) |> 
  arrange(desc(n)) |> 
  janitor::adorn_totals()

write_csv(x = erbb2_amps, file = "erbb2_amps.csv")

write_csv(x = erbb2_amps_summary, file = "erbb2_amps_summary.csv")
