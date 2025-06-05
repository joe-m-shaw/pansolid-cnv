# Collate gene transcript information

library(here)
library(tidyverse)

source(here("functions/transcript_functions.R"))

# This csv is manually curated
gene_labels <- read_csv(file = paste0(config::get("data_filepath"), 
                                      "validation/DOC6283_amplifications/",
                                      "gene_lists/",
                                      "gene_coordinates.csv"),
                        col_types = "ccccdd") |> 
  mutate(y_value = "Genes",
         # Place gene label half-way along gene locus
         start = pmin(gene_start, gene_end) + ((pmax(gene_start, gene_end) - 
                                                  pmin(gene_start, gene_end)) / 2))

write.csv(gene_labels, file = paste0(config::get("data_filepath"), 
                                     "validation/DOC6283_amplifications/",
                                    "transcripts/processed/",
                                     "gene_labels.csv"),
          row.names = FALSE)

# Exon tables (transcript files) are downloaded from Ensembl
# MANE transcript used for each gene
# Example: https://www.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000095002;r=2:47403067-47663146;t=ENST00000233146

transcript_files <- list.files(paste0(config::get("data_filepath"), 
                                      "validation/DOC6283_amplifications/",
                                      "transcripts/raw/"), full.names = TRUE,
                               pattern = ".csv")

all_transcripts <- transcript_files |>
  map(\(transcript_files) read_ensembl_exon_table(
    file = transcript_files
  )) |>
  list_rbind() |> 
  left_join(gene_labels |> 
              select(gene, chromosome, transcript_ensembl), 
            join_by(transcript == transcript_ensembl))

write.csv(all_transcripts, file = paste0(config::get("data_filepath"), 
                                         "validation/DOC6283_amplifications/",
                                         "transcripts/processed/",
                                         "collated_transcripts.csv"),
          row.names = FALSE)
