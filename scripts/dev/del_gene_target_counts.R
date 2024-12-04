# Deletion gene targets

library(here)
library(tidyverse)
library(janitor)

source(here("scripts/set_shared_drive_filepath.R"))
source(here("functions/gene_table_functions.R"))

del_table <- load_pansolid_gene_table("Deletions")

target_df <- read_csv(paste0(data_folder,
                             "bed_files/PanSolidv2_GRCh38_noalt_BED.csv")) |> 
  clean_names() |> 
  mutate(target_type = case_when(
    
    grepl(x = name, pattern =  "chr(\\d{1,2}|X):\\d{1,3}.+") == TRUE ~"genomic backbone",
    TRUE ~"gene target"
  ))

del_gene_targets <- target_df |> 
  filter(name %in% del_table$gene) |> 
  count(name)

del_gene_targets
