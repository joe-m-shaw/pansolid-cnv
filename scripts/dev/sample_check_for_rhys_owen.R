# Sample check for Rhys Owens

library(tidyverse)
library(patchwork)
library(readxl)
library(NHSRplotthedots)

source(here::here("scripts/connect_to_dna_db.R"))
source(here::here("functions/pansolid_cnv_excel_functions.R"))
source(here::here("functions/chromosome_arm_functions.R"))
source(here::here("functions/contamination_functions_INC10994.R"))

# Rhys wants to check that two samples from the same patient have matching
# SNP results.

WS156308_data <- collate_snp_data("WS156308")

WS156626_data <- collate_snp_data("WS156626")

all_data_25042104 <- rbind(WS156308_data, WS156626_data)

results_25048164 <- search_for_contaminant(all_data_25042104,
                                           "25048164_WS156626")

plot_25042104_25048164 <- make_contamination_snp_plot(all_data_25042104,
                            "25048164_WS156626",
                            all_data_25042104,
                            "25042104_WS156308",
                            hom_vaf_upper_threshold = 95,
                            het_vaf_upper_threshold = 90,
                            het_vaf_lower_threshold = 10,
                            hom_vaf_lower_threshold = 4)

ggsave("plot_25042104_25048164.png",
       plot_25042104_25048164)

