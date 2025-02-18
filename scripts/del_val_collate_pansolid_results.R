# Collate deletions validation PanSolid results

# Packages ----------------------------------------------------------------

library(here)
library(tidyverse)

# Functions ---------------------------------------------------------------

source(here("functions/pansolid_cnv_excel_functions.R"))

# Filepaths ---------------------------------------------------------------

pansolid_files <- list.files(path = paste0(config::get("data_folderpath"),
                                           "validation/DOC6567_deletions/raw/",
                                           "pansolid_ngs/annotated_files/"),
                        recursive = TRUE,
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

export_del_val_data <- function(df, df_name) {
  
  filepath = paste0(config::get("data_folderpath"),
                    "validation/",
                    "DOC6567_deletions/processed/", 
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
