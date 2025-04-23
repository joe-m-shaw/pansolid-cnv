# Find NF2 LOH results

# Packages ---------------------------------------------------------------------

library(here)
library(tidyverse)
library(readxl)

# Scripts ----------------------------------------------------------------------

source(here("scripts/connect_to_dna_db.R"))

# Find LOH results --------------------------------------------------------

data_folder <- config::get("data_folderpath")

del_val_collated_stdev <- read_csv(paste0(config::get("data_folderpath"), 
                                          "validation/DOC6567_deletions/processed/",
                                          "del_val_collated_stdev.csv"),
                                   col_types = list(
                                     "labno" = col_character()))

dlims_results <- results_tbl |> 
  filter(genodate > "2024-03-01 00:00:00") |> 
  select(pcrid, labno, surname, test, exon, genotype, genotype2,
         genocomm) |> 
  collect() |> 
  filter(!labno %in% c("water", "Water", "WATER"))

loh_results <- dlims_results |> 
  filter(grepl(pattern = "LOH", x = test, ignore.case = TRUE))

loh_regex <- regex(r"[
                   (\d{1,2}\.\d{1,2}|\d{1,2})   # 
                   %
                   (\sLOH|\sloh|\sloss)      #
                   .*                          # Text before or after
                   ]",
                   comments = TRUE)

loh_results_with_pansolid_samples <- loh_results |> 
  filter(labno %in% del_val_collated_stdev$labno) |> 
  select(-c(genotype2, test)) |> 
  arrange(labno) |> 
  mutate(loh_percent = round(as.numeric(str_extract(
    string = genocomm,
    pattern = loh_regex,
    group = 1
  )), 1)) |> 
  filter(!exon %in% c("D22S446", "D22S303", "D22S1174", "D22S310"))

informative_markers_df <- loh_results_with_pansolid_samples |> 
  filter(!is.na(loh_percent)) |> 
  count(labno) |> 
  rename(num_informative_markers = n)

if(any(informative_markers_df$num_informative_markers > 4)){
  stop()
}

loh_results_with_pansolid_samples_wide <- loh_results_with_pansolid_samples |> 
  distinct(labno, exon, .keep_all = TRUE) |> 
  select(-genotype) |> 
  pivot_wider(id_cols = c(labno, surname),
              names_from = exon,
              values_from = c(genocomm, loh_percent)) |> 
  rowwise() %>%
  mutate(median_loh = round(median(c_across(c(loh_percent_D22S268, 
                                              loh_percent_D22S275, 
                                              loh_percent_NF2CA3, 
                                              loh_percent_NF2intron10)), 
                               na.rm=TRUE), 1),
         loh_outcome = case_when(
           loh_percent_D22S268 >= 30 | 
             loh_percent_D22S275 >= 30 |
             loh_percent_NF2CA3 >= 30 |
             loh_percent_NF2intron10 >= 30 ~"Significant LOH",
           TRUE ~"No significant LOH"
         )) |> 
  left_join(informative_markers_df, by = "labno")

write_csv(loh_results_with_pansolid_samples_wide,
          file = paste0(data_folder,
                        "validation/DOC6567_deletions/raw/nf2_loh/",
                        "nf2_loh_dna_db_results.csv"))
