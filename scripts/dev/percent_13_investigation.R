

library(tidyverse)
library(readxl)
library(here)
library(janitor)

source(here("scripts/load_processed_live_service_data.R"))

ngs_progress_sheet <- read_excel(path = paste0("S:/central shared/Genetics/Repository/Technical Teams/NGS/",
                                               "NGS Technical Team Progress Sheet - 2024 August 2024.xlsm"),
                                 sheet = "QIAseq DNA PanSolid") |> 
  clean_names()

number_regex <- "(\\d{1,3}|\\d{1,3}.\\d{1,2})"

ngs_progress_sheet_clean <- ngs_progress_sheet |> 
  filter(!is.na(worksheet_number)) |> 
  mutate(percent_q30_number = as.numeric(str_extract(string = percent_q30,
                                                     pattern = number_regex)),
         cluster_pf_percent_number = as.numeric(str_extract(string = cluster_pf_percent,
                                                            pattern = number_regex)),
         indexing_percent_reads_identified_number = as.numeric(str_extract(string = indexing_percent_reads_identified,
                                                                           pattern = number_regex))) |> 
  rename(worksheet = worksheet_number) |> 
  mutate(worksheet_number = as.numeric(str_extract(string = worksheet,
                                                   pattern = "WS(\\d{6})",
                                                   group = 1)))

live_service_percent_138_results_collated |> 
  mutate(worksheet_number = as.numeric(str_extract(string = worksheet,
                                                   pattern = "WS(\\d{6})",
                                                   group = 1))) |> 
  filter(worksheet_number >= 146232) |> 
  ggplot(aes(x = worksheet, y = percent_whole_panel_covered_at_138x)) +
  geom_boxplot(outliers = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Percentage panel at 138X", x = "")

worksheet_cov_median <- live_service_percent_138_results_collated |> 
  group_by(worksheet) |> 
  summarise(median_138 = median(percent_whole_panel_covered_at_138x)) |> 
  mutate(coverage_category = case_when(
    median_138 >= 90 ~"good",
    median_138 < 90 ~"poor"
  ))

ngs_progress_sheet_with_138 <- worksheet_cov_median |> 
  left_join(ngs_progress_sheet_clean, by = "worksheet") |> 
  #filter(worksheet_number >= 146413) |> 
  select(worksheet, qi_aseq_library_prep_cp_worksheet_number,
         median_138, coverage_category, worksheet_picked_up_by, 
         frag_er_a_tail_completed_date, adapter_ligation_completed_date,
         target_enrichment_completed_date) 

"S:\central shared/Genetics/Repository/Reagents/Reagent Acceptance/NGS/QIAseq PanSolid/QIAseq PanSolid Batch records.xlsx"

qc_tabs <- read_csv(str_c(data_folder,
                          "live_service/collated/",
                          "live_service_qc_tabs_collated.csv"),
                    col_types = "cccddddd") |> 
  mutate(worksheet_number = as.numeric(str_extract(string = worksheet,
                                                   pattern = "WS(\\d{6})",
                                                   group = 1)))


qc_tabs |> 
  filter(worksheet_number >= 146413) |> 
  ggplot(aes(x = worksheet, y = average_number_of_reads_per_umi)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "")