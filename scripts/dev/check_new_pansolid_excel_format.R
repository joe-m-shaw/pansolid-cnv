# Checking new PanSolid Excel format

library(tidyverse)
source(here::here("functions/pansolid_cnv_excel_functions.R"))

new_excel_files <- list.files(
  path = paste0(config::get("data_folderpath"),
                            "validation/DOC6567_deletions/raw/"),
  full.names = TRUE,
  recursive = FALSE,
  pattern = "Annotated.*.xlsx"
)

old_excel_files <- list.files(path = paste0(config::get("data_folderpath"),
                                            "validation/DOC6567_deletions/raw/",
                                            "pansolid_ngs/v2PANSOLID/"),
                              full.names = TRUE,
                              recursive = FALSE,
                              pattern = "Annotated.*.xlsx")

find_pred_ncc <- function(input_sheet, 
                          pred_ncc_string = "Predicted NCC (%)") {
  
  pred_ncc_start <- find_match(input_sheet, "a",
                               pred_ncc_string) + 1
  
  pred_ncc_df <- tibble::as_tibble(input_sheet[pred_ncc_start, 1]) |> 
    dplyr::rename(pred_ncc = a) |> 
    dplyr::mutate(pred_ncc = as.double(pred_ncc))
  
  if(pred_ncc_df[1,1] < 0 | pred_ncc_df[1,1] > 100) {
    stop("Predicted NCC must be between 0 and 100")
  }
  
  return(pred_ncc_df)
  
}

new_sheet <- read_cnv_sheet(filepath = new_excel_files[1], 
                            sheet_regex = "CNVs")

old_sheet <- read_cnv_sheet(filepath = old_excel_files[1],
                            sheet_regex = "CNVs")

del_genes_new <- find_del_genes(new_sheet, num_genes = 34)

del_genes_old <- find_del_genes(old_sheet, num_genes = 37)

amp_genes_new <- find_amp_genes(new_sheet, num_genes = 13)

amp_genes_old <- find_amp_genes(old_sheet, num_genes = 9)

new_ncc <- find_pred_ncc(new_sheet)

sig_cnvs_new <- find_sig_cnvs(new_sheet)

sig_cnvs_old <- find_sig_cnvs(old_sheet)
