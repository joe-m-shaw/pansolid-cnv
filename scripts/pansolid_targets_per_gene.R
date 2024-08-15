# QIAseq Targets per Gene on PanSolid

# Packages --------------------------------------------------------------------------

library(here)
library(tidyverse)
library(janitor)
library(readxl)

# Functions and filepath ------------------------------------------------------------

source(here("scripts/set_shared_drive_filepath.R"))
source(here("functions/cnv_functions.R"))

# Load tables -----------------------------------------------------------------------

primer_df <- read_csv(paste0(data_folder, 
                             "primers/QIAseq.CDHS-48608Z-11752.primers3 (Converted).csv")) |> 
  clean_names()

gene_tbl <- read_excel(path = paste0(data_folder, "transcripts/gene_labels.xlsx"),
                       col_types = c("text", "text", "text", "text", "numeric", "numeric"))

amp_genes <- load_gene_table("Amplifications")

target_df <- read_csv(paste0(data_folder,
                             "bed_files/PanSolidv2_GRCh38_noalt_BED.csv")) |> 
  clean_names()

# Get primer coordinates ------------------------------------------------------------

primer_df_mod <- extract_cnv_coordinates(df = primer_df, cnv_coord_col = region) |> 
  mutate(region_length = abs(start-end)) |> 
  rename(primer_start = start,
         primer_end = end)

# Count primers within gene regions -------------------------------------------------

gene_tbl_mod <- gene_tbl |> 
  mutate(gene_size_kb = round(abs(gene_start - gene_end) / 1000, 1)) |> 
  rowwise() |> 
  mutate(primers_in_gene = count_primers_in_region(chrom = chromosome,
                                                   coord1 = gene_start,
                                                   coord2 = gene_end,
                                                   df = primer_df_mod)) |> 
  filter(label %in% amp_genes$gene |
           label == "EGFRvIII") |> 
  arrange(label)

# Count targets for each gene -------------------------------------------------------

gene_target_counts <- target_df |> 
  filter(name %in% amp_genes$gene) |> 
  count(name)

gene_tbl_for_doc <- gene_tbl_mod |> 
  mutate(grch38_coordinates = str_c("chr", chromosome, ":", gene_start, "-", gene_end)) |> 
  rename(gene = label) |> 
  select(gene, grch38_coordinates, gene_size_kb, transcript_refseq,
         primers_in_gene) |> 
  left_join(gene_target_counts, 
            join_by("gene" == "name")) |> 
  rename(targets_in_gene = n)

csv_timestamp(gene_tbl_for_doc, paste0(outputs_folder, "tables"))
