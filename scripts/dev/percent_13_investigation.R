# Coverage at 138X investigation

# Packages ----------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)
library(janitor)
library(patchwork)

source(here("scripts/load_processed_live_service_data.R"))

source(here("scripts/connect_to_dna_db.R"))

# Data --------------------------------------------------------------------

ngs_progress_sheet24 <- read_excel(path = paste0("S:/central shared/Genetics/Repository/Technical Teams/NGS/",
                                               "NGS Technical Team Progress Sheet - 2024 August 2024.xlsm"),
                                 sheet = "QIAseq DNA PanSolid") |> 
  clean_names()

ngs_progress_sheet25 <- read_excel(path = paste0("S:/central shared/Genetics/Repository/Technical Teams/NGS/",
                                                 "NGS Technical Team Progress Sheet - 2025.xlsm"),
                                 sheet = "QIAseq DNA PanSolid") |> 
  clean_names()

ngs_progress_sheet <- ngs_progress_sheet24 |> 
  select(worksheet_number, worksheet_picked_up_by) |> 
  rbind(ngs_progress_sheet25 |> 
          select(worksheet_number, worksheet_picked_up_by))


colnames(live_service_std_dev_results_collated)

recent_data <- live_service_percent_138_results_collated |> 
  mutate(worksheet_number = as.numeric(str_extract(string = worksheet,
                                                   pattern = "WS(\\d{6})",
                                                   group = 1))) |> 
  filter(worksheet_number >= 146232) |> 
  left_join(live_service_std_dev_results_collated |> 
              select(filepath, st_dev_signal_adjusted_log2_ratios),
            by = "filepath") |> 
  left_join(ngs_progress_sheet, join_by("worksheet" == "worksheet_number"))

recent_data |> 
  filter(!is.na(worksheet_picked_up_by)) |> 
  ggplot(aes(x = worksheet, 
             y = percent_whole_panel_covered_at_138x)) +
  geom_boxplot(outliers = FALSE) +
  theme_bw() +
  labs(y = "Percentage panel at 138X", x = "",
       title = "PanSolid Worksheet Data") +
  facet_wrap(~worksheet_picked_up_by)








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

pansolid_labnos <- unique(live_service_percent_138_results_collated$labno)

pansolid_dna_concs <- sample_tbl |> 
  filter(labno %in% pansolid_labnos) |> 
  select(labno, concentration) |> 
  collect()

# Plots -------------------------------------------------------------------

cov_138_plot <- live_service_percent_138_results_collated |> 
  mutate(worksheet_number = as.numeric(str_extract(string = worksheet,
                                                   pattern = "WS(\\d{6})",
                                                   group = 1))) |> 
  filter(worksheet_number >= 146232) |> 
  ggplot(aes(x = worksheet, y = percent_whole_panel_covered_at_138x)) +
  geom_boxplot(outliers = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Percentage panel at 138X", x = "",
       title = "PanSolid Worksheet Data")

# DNA concentrations

dna_conc_plot <- live_service_percent_138_results_collated |> 
  left_join(pansolid_dna_concs, by = "labno") |> 
  mutate(dna_concentration = as.double(concentration)) |> 
  mutate(worksheet_number = as.numeric(str_extract(string = worksheet,
                                                   pattern = "WS(\\d{6})",
                                                   group = 1))) |> 
  filter(worksheet_number >= 146232) |> 
  ggplot(aes(x = worksheet, y = dna_concentration)) +
  geom_boxplot(outliers = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "DNA concentration (ng/ul)", x = "")

cov_138_plot + dna_conc_plot + plot_layout(ncol = 1)

# Panels

panel_plot <- panel_info |> 
  mutate(worksheet_number = as.numeric(str_extract(string = worksheet,
                                                   pattern = "WS(\\d{6})",
                                                   group = 1))) |> 
  filter(worksheet_number >= 146232 &!is.na(panel)) |> 
  ggplot(aes(x = worksheet, y =)) +
  geom_bar(position = "stack", aes(fill = panel)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90)) +
  labs(y = "Samples", x = "")

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

qc_tabs <- read_csv(str_c(data_folder,
                          "live_service/collated/",
                          "live_service_qc_tabs_collated.csv"),
                    col_types = "cccddddd") |> 
  mutate(worksheet_number = as.numeric(str_extract(string = worksheet,
                                                   pattern = "WS(\\d{6})",
                                                   group = 1)))

umi_plot <- qc_tabs |> 
  filter(worksheet_number >= 146413) |> 
  ggplot(aes(x = worksheet, y = average_number_of_reads_per_umi)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "")
