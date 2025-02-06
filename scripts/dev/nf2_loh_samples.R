# Find NF2 LOH results

# Packages ---------------------------------------------------------------------

library(here)
library(tidyverse)
library(readxl)

# Scripts ----------------------------------------------------------------------

source(here("scripts/connect_to_dna_db.R"))
source(here("functions/dna_db_functions.R"))
source(here("functions/join_dna_submission_sheets.R"))

data_folder <- config::get("data_filepath")

# Loss of heterozygosity samples -----------------------------------------------

pansolidv2_worksheets <- read_excel(here(data_folder,
                                         "live_service/",
                                         "pansolid_live_service_worksheets.xlsx")) |> 
  mutate(pcrid = str_replace(string = worksheet,
                             pattern = "WS",
                             replacement = ""))

pcrid_list <- pansolidv2_worksheets$pcrid

pansolid_worksheet_samples <- dna_db_pcr_records |> 
  filter(pcrid %in% pcrid_list) |> 
  select(pcrid, sample, name, test_type) |> 
  collect() |> 
  rename(labno = sample)

dlims_results <- results_tbl |> 
  filter(genodate > "2024-03-01 00:00:00") |> 
  select(pcrid, labno, surname, test, exon, genotype, genotype2,
         genocomm) |> 
  collect() |> 
  filter(!labno %in% c("water", "Water", "WATER"))

loh_results <- dlims_results |> 
  filter(grepl(pattern = "LOH", x = test, ignore.case = TRUE))

loh_labnos <- loh_results$labno

loh_tissue_types <- get_sample_tissue(sample_vector = loh_labnos) |> 
  select(labno, tissue_type)

pansolid_submission_sheet <- join_pansolid_submission_sheets()

pansolid_samples_with_loh_results <- pansolid_worksheet_samples |> 
  filter(labno %in% loh_results$labno) |> 
  filter(!duplicated(name)) |> 
  mutate(worksheet = paste0("WS", pcrid)) |> 
  select(labno, worksheet, name) |> 
  left_join(pansolid_submission_sheet |> 
              select(labno, panel, enrichment),
            by = "labno") |> 
  mutate(category = "NF2 LOH testing") |> 
  select(-name)

pansolid_deletions_loh_cohort <- read_excel(
  path = paste0("S:/central shared/Genetics/NGS/Bioinformatics/",
                "1_Pan-solid-Cancer/CNV/Deletions/",
                "pansolid_deletions_loh_cohort.xlsx")) |> 
  janitor::clean_names() |> 
  filter(category == "NF2 LOH testing")

loh_regex <- regex(r"[
                   (\d{1,2}\.\d{1,2}|\d{1,2})   # 
                   %
                   (\sLOH|\sloh|\sloss)      #
                   .*                          # Text before or after
                   ]",
                   comments = TRUE)

loh_results_with_pansolid_samples <- loh_results |> 
  filter(labno %in% pansolid_deletions_loh_cohort$labno) |> 
  select(-c(genotype2, test)) |> 
  arrange(labno) |> 
  mutate(loh_percent = round(as.numeric(str_extract(
    string = genocomm,
    pattern = loh_regex,
    group = 1
  )), 1)) |> 
  left_join(loh_tissue_types,
            by = "labno")

loh_results_with_pansolid_samples_wide <- loh_results_with_pansolid_samples |> 
  select(-genotype) |> 
  filter(!exon %in% c("D22S446", "D22S303", "D22S1174", "D22S310")) |> 
  # Sample with multiple entries but no LOH results
  filter(labno != "24031830") |> 
  pivot_wider(id_cols = c(labno, surname, tissue_type),
              names_from = exon,
              values_from = c(genocomm, loh_percent)) |> 
  select(-c(tissue_type)) |> 
  rowwise() %>%
  mutate(median_loh = round(median(c_across(c(loh_percent_D22S268, 
                                              loh_percent_D22S275, 
                                              loh_percent_NF2CA3, 
                                              loh_percent_NF2intron10)), 
                               na.rm=TRUE), 1),
         loh_outcome = case_when(
           median_loh < 30 ~"No significant LOH",
           median_loh >= 30 ~"Significant LOH",
           TRUE ~NA
         ))

write_csv(loh_results_with_pansolid_samples_wide,
          file = paste0(data_folder,
                        "validation/DOC6567_deletions/raw/nf2_loh/",
                        "nf2_loh_dna_db_results.csv"))
