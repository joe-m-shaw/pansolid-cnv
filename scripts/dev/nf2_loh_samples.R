# Find NF2 LOH samples

# Packages ---------------------------------------------------------------------------------------

library(here)
library(tidyverse)
library(readxl)

# Scripts ----------------------------------------------------------------------------------------

source(here("scripts/set_shared_drive_filepath.R"))
source(here("scripts/connect_to_dna_db.R"))
source(here("functions/dna_db_functions.R"))
source(here("functions/join_dna_submission_sheets.R"))

# Loss of heterozygosity samples -----------------------------------------------------------------

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

loh_results_with_pansolid_samples <- loh_results |> 
  filter(labno %in% pansolid_worksheet_samples$labno) |> 
  select(-c(genotype2, test)) |> 
  arrange(labno) |> 
  left_join(loh_tissue_types,
            by = "labno")

loh_results_with_pansolid_samples_wide <- loh_results_with_pansolid_samples |> 
  select(-genotype) |> 
  pivot_wider(id_cols = c(labno, surname, tissue_type),
              names_from = exon,
              values_from = genocomm)

write.csv(pansolid_samples_with_loh_results,
          "pansolid_samples_with_loh_results.csv",
          row.names = FALSE)

# Get panel information for ddPCR samples --------------------------------------------------------

cdkn2a_ddpcr <- read_csv("cdkn2a_pten_ddpcr.csv",
                         col_types = "c")

cdkn2a_ddpcr_with_panels <- cdkn2a_ddpcr |> 
  left_join(pansolid_submission_sheet |> 
              select(labno, panel),
            by = "labno") |> 
  mutate(category = "CDKN2A and PTEN ddPCR") |> 
  left_join(pansolid_worksheet_samples,
            by = "labno") |> 
  mutate(worksheet = paste0("WS", pcrid)) |> 
  select(labno, worksheet, panel, category) |> 
  filter(!is.na(panel))

write.csv(cdkn2a_ddpcr_with_panels, "cdkn2a_ddpcr_with_panels.csv",
          row.names = FALSE)
