# SeraCare Control Results

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(here)
library(readxl)

# Source scripts --------------------------------------------------------------------

source(here("scripts/set_shared_drive_filepath.R"))
source(here("scripts/collate_pansolid_ddpcr_validation_data.R"))

# Load PanSolid data ----------------------------------------------------------------

validation_all_amp_gene_results_collated <- read_csv(paste0(data_folder, 
                                                            "validation/collated/",
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
                                                     )) |> 
  # The limit of detection normal control was listed as 24039973 (CNVMix12CopiesSERASEQ
  # at 0%) but can be changed to 24039975 (DNAWTmixSERASEQ at 100%)
  mutate(labno = case_when(
    
    labno_suffix_worksheet == "24039973d_WS144291" ~"24039975",
    TRUE ~labno),
    
    patient_name = case_when(
      labno_suffix_worksheet == "24039973d_WS144291" ~"DNAWTmixSERASEQ",
      TRUE ~patient_name),
    
    total_copies_ngs_pansolid = max_region_fold_change * 2
    )

# Load SeraCare results -------------------------------------------------------------

seracare_gene_copies <- read_excel(path = paste0(data_folder, 
                                                 "seracare_reference_materials/",
                                                 "seracare_gene_copies.xlsx")) |> 
  mutate(total_copies_ddpcr_seracare = additional_copies_ddpcr + 2,
         total_copies_ngs_seracare = additional_copies_ngs + 2)

# SeraCare lab numbers --------------------------------------------------------------

seracare_ids <- validation_all_amp_gene_results_collated |> 
  filter(!duplicated(labno)) |>
  select(labno, patient_name) |> 
  inner_join(seracare_gene_copies |> 
               filter(!duplicated(reference_material)) |> 
               select(reference_material),
             join_by("patient_name" == "reference_material"))

# NGS comparison --------------------------------------------------------------------

seracare_colours <- c(
  "#D55E00",
  "#E69F00", 
  "#F0E442",
  "#009E73"
)

seracare_ngs_comparison <- validation_all_amp_gene_results_collated |> 
  # Remove limit of detection samples
  filter(!(worksheet == "WS144291" & labno ==  "24039973")) |> 
  inner_join(seracare_gene_copies, 
             join_by("patient_name" == "reference_material",
                     "gene" == "gene")) |> 
  mutate(patient_name = factor(x = patient_name,
                               levels = c("CNVMix12CopiesSERASEQ",
                                          "CNVMix6CopiesSERASEQ",
                                          "CNVMix3copiesSERASEQ",
                                          "DNAWTmixSERASEQ")))

ngs_comparison_plot <- ggplot(seracare_ngs_comparison, 
                              aes(x = total_copies_ngs_pansolid,
                                    y = total_copies_ngs_seracare )) +
  geom_abline(linetype = "dashed") +
  geom_point(shape = 21, size = 3, alpha = 0.6, aes(fill = patient_name)) +
  scale_fill_manual(values = seracare_colours) +
  ggpubr::stat_cor(method = "pearson", label.x = 0, label.y = 30) +
  theme_bw() +
  scale_x_continuous(limits = c(-5, 35), 
                     breaks = c(0, 10, 20, 30)) +
  scale_y_continuous(limits = c(-5, 35), 
                    breaks = c(0, 10, 20, 30)) +
  facet_wrap(~gene) +
  labs(x = "Total copies NGS (Manchester)",
       y = "Total copies NGS (SeraCare)",
       title = "Comparison of NGS results for SeraCare reference materials",
       fill = "")
  
# ddPCR comparison ------------------------------------------------------------------

seracare_ddpcr_comparison <- ddpcr_validation_data |> 
  filter(target_type == "Ch1Unknown") |> 
  inner_join(seracare_ids, join_by("sample" == "labno")) |> 
  inner_join(seracare_gene_copies, join_by("patient_name" == "reference_material",
                                            "gene" == "gene")) |> 
  mutate(patient_name = factor(x = patient_name,
                               levels = c("CNVMix12CopiesSERASEQ",
                                          "CNVMix6CopiesSERASEQ",
                                          "CNVMix3copiesSERASEQ",
                                          "DNAWTmixSERASEQ")))

ddpcr_comparison_plot <- ggplot(seracare_ddpcr_comparison, 
                              aes(x = cnv,
                                  y = total_copies_ddpcr_seracare )) +
  geom_abline(linetype = "dashed") +
  geom_point(shape = 21, size = 3, alpha = 0.6, aes(fill = patient_name)) +
  scale_fill_manual(values = seracare_colours) +
  facet_wrap(~gene) +
  ggpubr::stat_cor(method = "pearson", label.x = 5, label.y = 2) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 25), 
                     breaks = c(0, 10, 20)) +
  scale_y_continuous(limits = c(0, 25), 
                     breaks = c(0, 10, 20)) +
  labs(x = "Total copies ddPCR (Manchester)",
       y = "Total copies ddPCR (SeraCare)",
       title = "Comparison of ddPCR results for SeraCare reference materials",
       subtitle = "Pearson's coefficient not calculated for genes with only 2 data points",
       fill = "")
