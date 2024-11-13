# Plots for NHS-R presentation

library(here)
library(tidyverse)

source(here("scripts/set_shared_drive_filepath.R"))
source(here("scripts/load_processed_validation_data.R"))
source(here("scripts/load_processed_live_service_data.R"))

amp_results_with_cancer_types <- amp_validation_all_amp_gene_results_collated |> 
  left_join(validation_sample_patient_info,
            by = "labno", 
            relationship = "many-to-one")

amp_results_with_cancer_types_filter <- amp_results_with_cancer_types |> 
  mutate(cancer_tissue_source = case_when(
    cancer_tissue_source == "central nervous system" ~"brain",
    cancer_tissue_source == "hepatocellular carcinoma" ~"liver",
    TRUE ~cancer_tissue_source
  ),
  max_region_fold_change = case_when(
    max_region_fold_change < 1 ~1,
    TRUE ~max_region_fold_change
  ),
  pansolid_call = case_when(
    pansolid_call == "normal result" ~"normal",
    pansolid_call == "amplification" ~"high dose"
  )) |> 
  filter(cancer_tissue_source %in% c("brain", "lung", "colorectal"))

dosage_cancer_type_plot <- ggplot(amp_results_with_cancer_types_filter, 
                                 aes(x = gene,
                                     y = max_region_fold_change)) +
  geom_jitter(shape = 21, aes(fill = pansolid_call), width = 0.2,
              alpha = 0.6, size = 3) +
  scale_fill_manual(values = c("#CC6677", "#FFFFFF")) +
  facet_wrap(~cancer_tissue_source, nrow = 1, ncol = 3, axes = "all_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.size = unit(2,"point"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  scale_y_log10(breaks = c(1, 10, 100)) +
  labs(x = "Gene", y = "Dosage",
       fill = "")

ggsave(
  filename =  "dosage_cancer_type_plot.png",
  plot = dosage_cancer_type_plot,
  device = "png",
  path = "scripts/dev/",
  units = "cm",
  width = 25,
  height = 14,
  dpi = 300
)

noise_worksheet_data <- live_service_std_dev_results_collated |> 
  mutate(ws_number = as.numeric(str_extract(string = worksheet,
                                            pattern = "WS(\\d{6})",
                                            group = 1))) |> 
  filter(ws_number >= 141915 &
           ws_number < 145671)

worksheet_noise_plot <- ggplot(noise_worksheet_data, aes(x = worksheet, 
                                 y = st_dev_signal_adjusted_log2_ratios)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12)) +
  coord_cartesian(ylim = c(0, 1.5)) +
  labs(x = "Batch", y = "Quality metric")

ggsave(
  filename =  "worksheet_noise_plot.png",
  plot = worksheet_noise_plot,
  device = "png",
  path = "scripts/dev/",
  units = "cm",
  width = 22,
  height = 10,
  dpi = 300
)


# MYC plot ---------------------------------------------------------------------------------------

myc_pansolid_vs_ddpcr <- validation_ddpcr_collated |> 
  mutate(sample_experiment = str_c(sample, "_", experiment)) |> 
  # Remove assay optimisation worksheet
  filter(worksheet != "WS138419") |> 
  filter(target == "ch1_target" & 
           experiment  == "MYC Ex3_1") |> 
  filter(!duplicated(sample_experiment)) |> 
  # Remove controls
  filter(grepl(pattern = "\\d{8}", x = sample)) |> 
  inner_join(amp_validation_all_amp_gene_results_collated |> 
               mutate(labno_suffix_gene = str_c(labno_suffix, "_", gene)) |> 
               # Some samples were run twice on PanSolid
               filter(!duplicated(labno_suffix_gene)),
             join_by(sample == labno, 
                     gene)) |> 
  mutate(
    pansolid_copy_number = max_region_fold_change * 2)

myc_dot_plot <- ggplot(myc_pansolid_vs_ddpcr |> 
         filter(pansolid_copy_number >= 0), 
       aes(x = cnv, y = pansolid_copy_number)) +
  geom_abline(linetype = "dashed", size = 1) +
  geom_point(shape = 21, size = 3) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)) +
  labs(x = "Gene dosage (different test)",
       y = "Gene dosage (new test)") +
  ylim(0, 220) +
  xlim(0, 220) 

ggsave(
  filename =  "myc_dot_plot.png",
  plot = myc_dot_plot,
  device = "png",
  path = "scripts/dev/",
  units = "cm",
  width = 14,
  height = 14,
  dpi = 300
)


