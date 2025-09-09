# Checking PanSolid BED file

# 13/08/2025

library(readxl)
library(tidyverse)

# I took the BED file for the regions of interest on PanSolid v2 (CDHS-48608Z-11752)
# and converted the coordinates from GRCh37 to GRCh38 using the Ensembl converter at
# https://grch37.ensembl.org/Homo_sapiens/Tools/AssemblyConverter

grch37_bed <- read_delim(file = paste0(config::get("data_folderpath"),
                                       "validation/DOC6791_chromosome_arms/bed_files/",
                                       "QIAseq_DNA_panel.CDHS-48608Z-11752.roi.bed"),
                         skip = 1,
                         col_names = c("chrom", "start", "stop",
                                       "name"))

grch38_bed <- read_delim(file = paste0(config::get("data_folderpath"),
                         "validation/DOC6791_chromosome_arms/bed_files/",
                         "QIAseq_DNA_panel.CDHS-48608Z-11752.roi_GRCh38.bed"),
                         col_names = c("chrom", "start", "stop",
                                       "name"))
# There are 6173 rows

# Here is an example output from the PanSolid pipeline

example_sample <- read_excel(path = paste0(config::get("data_folderpath"),
                                           "validation/DOC6567_deletions/",
                                           "raw/pansolid_ngs/final_format/",
                                           "v2PANSOLID/Unannotated_Files/",
                                           "Results_SNVs_Indels-WS146593_24025266_S19_R1_001.xlsx"),
                             sheet = "CNV Targets Merged",
                             range = "A1:C6175",
                             col_types = c("text", "text", "text")) |> 
  janitor::clean_names()

# There are 6174 rows

example_sample_mod <- example_sample |> 
  mutate(
    chrom = str_replace(chromosome, pattern = "\\.0", replacement = ""),
    target = str_c("chr", chrom, "_",
                        region)) |> 
  select(chromosome, chrom, region, target)

grch38_bed_mod <- grch38_bed |> 
  mutate(
    # For some reason there is a difference of 1bp between the BED file and
    # the PanSolid pipeline output
    start_plus_one = start + 1,
    target = str_c(chrom, "_", 
                   start_plus_one, "..", stop))

test_df <- inner_join(example_sample_mod,
                      grch38_bed_mod,
                      by ="target",
                      relationship = "one-to-one")

# There is one target on the PanSolid output which is not in the BED file
# chr17:36.2Mb
# target 4901
setdiff(example_sample_mod$target, test_df$target)
# This doesn't get converted to GRCh38 for some reason
