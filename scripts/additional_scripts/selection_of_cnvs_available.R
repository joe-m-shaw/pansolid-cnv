# Selection of Copy Number Variants Available

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

# Functions -------------------------------------------------------------------------

source(here::here("functions/cnv_functions.R"))

# Read in WGS reports ---------------------------------------------------------------

list.files(here("data/brain_matrix_htmls/"),
           full.names = TRUE)

sample1 <- parse_wgs_html_table_by_div_id(here("data/brain_matrix_htmls/4018657512_p02520587765_LP5100626-DNA_A01_LP5100650-DNA_C01.v2_47_1.supplementary.html"),
                                          div_id = "svcnv_tier1") |> 
  mutate(pid = "p02520587765") |> 
  relocate(pid) |> 
  select(pid, gene, impacted_transcript_region, variant_type)

sample2 <- parse_wgs_html_table_by_div_id(here("data/brain_matrix_htmls/4018664361_p94137667399_LP5100568-DNA_B02_LP5100591-DNA_E02.v2_47_1.supplementary.html"),
                                          div_id = "svcnv_tier1") |> 
  mutate(pid = "p94137667399") |> 
  relocate(pid)|> 
  select(pid, gene, impacted_transcript_region, variant_type)

sample3 <- parse_wgs_html_table_by_div_id(here("data/brain_matrix_htmls/4018665006_p96656846595_LP5100870-DNA_G01_LP5100839-DNA_F01.v2_47_1.supplementary.html"),
                                          div_id = "svcnv_tier1") |> 
  mutate(pid = "p96656846595") |> 
  relocate(pid)|> 
  select(pid, gene, impacted_transcript_region, variant_type)

sample4 <- parse_wgs_html_table_by_div_id(here("data/brain_matrix_htmls/4018664520_p12216122864_LP5100945-DNA_G01_LP5100946-DNA_G01.v2_47_1.supplementary.html"),
                                          div_id = "svcnv_tier1") |> 
  mutate(pid = "p12216122864") |> 
  relocate(pid)|> 
  select(pid, gene, impacted_transcript_region, variant_type)


sample5 <- parse_wgs_html_table_by_div_id(here("data/brain_matrix_htmls/4018660708_p99674341077_LP5100876-DNA_B01_LP5100844-DNA_B01.v2_47_1.supplementary.html"),
                                          div_id = "svcnv_tier1") |> 
  mutate(pid = "p99674341077") |> 
  relocate(pid)|> 
  select(pid, gene, impacted_transcript_region, variant_type)

sample6 <- parse_wgs_html_table_by_div_id(here("data/wgs_result_htmls/4034941694_p19236663685_LP5101208-DNA_B01_LP5101207-DNA_D01.v3_7.supplementary.html"),
                                          div_id = "svcnv_tier1") |> 
  mutate(pid = "p19236663685") |> 
  relocate(pid)|> 
  select(pid, gene, impacted_transcript_region, variant_type)

# There are 7 samples available for testing

joined_samples <- rbind(sample1, sample2, sample3, sample4, sample5, sample6)

variants_of_interest <- grep(pattern = "CDKN2A|EGFR|MDM2|PTEN|MSH6|MSH2|MLH1|MET",
                             x = joined_samples$gene, value = TRUE, ignore.case = TRUE)

joined_samples |> 
  filter(gene %in% variants_of_interest) |> 
  count(gene, variant_type) |>  view()


