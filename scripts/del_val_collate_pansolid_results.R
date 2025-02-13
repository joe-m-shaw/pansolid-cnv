# Collate deletions validation PanSolid results

# Packages ----------------------------------------------------------------

library(here)
library(tidyverse)

# Functions ---------------------------------------------------------------

source(here("functions/pansolid_cnv_excel_functions.R"))

# Filepaths ---------------------------------------------------------------

bio_cnv_folder <- paste0("S:/central shared/Genetics/NGS/Bioinformatics/",
                         "1_Pan-solid-Cancer/CNV/Deletions/")

pansolid_results_folder <- paste0(bio_cnv_folder,
                                  "Jan2025_withLOH/",
                                  "v4_withoutNormalisation_MergeSize5_TF200/",
                                  "v2SchwannAll_PS/")

pansolid_files <- list.files(path = pansolid_results_folder,
                        recursive = FALSE,
                        full.names = TRUE,
                        pattern  = "Annotated.*.xlsx")

# LOH results -------------------------------------------------------------

collated_loh <- pansolid_files |> 
  map(\(pansolid_files) read_loh_table(filepath = pansolid_files)) |> 
  list_rbind() 

stopifnot(length(pansolid_files) * 7 == nrow(collated_loh))

stopifnot(setequal(unique(collated_loh$gene), 
          c("MSH2", "MSH6", "MLH1", "PMS2", "LZTR1", "SMARCB1", "NF2")))

stopifnot(anyNA.data.frame(collated_loh) == FALSE)

# CNV results -------------------------------------------------------------

file_cnv_tbl_list <- pansolid_files |> 
  map(\(pansolid_files) extract_cnv_tbls(pansolid_files, 
                                         sheet_regex = "CNVs_"))

collated_stdev <- map(file_cnv_tbl_list, ~ .x[["stdev"]]) |> 
  list_rbind()

stopifnot(length(setdiff(collated_stdev$filepath, pansolid_files)) == 0)

collated_138x <- map(file_cnv_tbl_list, ~ .x[["percent_138x"]]) |> 
  list_rbind()

stopifnot(length(setdiff(collated_138x$filepath, pansolid_files)) == 0)

collated_del_genes <- map(file_cnv_tbl_list, ~ .x[["del_genes"]]) |> 
  list_rbind()

stopifnot(nrow(collated_del_genes) == length(pansolid_files) * 37)

collated_amp_genes <- map(file_cnv_tbl_list, ~ .x[["amp_genes"]]) |> 
  list_rbind()

stopifnot(nrow(collated_amp_genes) == length(pansolid_files) * 9)

collated_sig_cnvs <- map(file_cnv_tbl_list, ~ .x[["sig_cnvs"]]) |> 
  list_rbind()

# Export collated data ----------------------------------------------------

collated_data_folder <- paste0("S:/central shared/Genetics/Mol_Shared/",
                               "Development.Team/PanSolid CNV/data/validation/",
                               "DOC6567_deletions/processed/")

export_del_val_data <- function(df, df_name) {
  
  filepath = paste0(collated_data_folder, 
                    "del_val_", 
                    df_name,
                    ".csv")
  
  write_csv(df, filepath)
  
}

df_list <- list(
  "collated_loh" = collated_loh,
  "collated_stdev" = collated_stdev,
  "collated_138x" = collated_138x,
  "collated_amp_genes" = collated_amp_genes,
  "collated_del_genes" = collated_del_genes,
  "collated_sig_cnvs" = collated_sig_cnvs
)

imap(df_list, export_del_val_data)
