## New way to read PanSolid Excels

library(here)
library(tidyverse)

source(here("functions/pansolid_excel_functions.R"))

ps_excels <- list.files(path = paste0(
  config::get("data_filepath"),
  "validation/DOC6283_amplifications/raw/pansolid_ngs_amplifications/"),
  full.names = TRUE)

multi_read <- function(file) {
  
  pos_cnv <- read_pos_cnv_results(file = file, 
                                  sheet = get_amp_sheetname(file))
  
  all_amp <- read_all_amp_genes_results(file = file,
                                        sheet = get_amp_sheetname(file))
  
  tbls <- list(
    "pos_cnv" = pos_cnv, 
    "all_amp" = all_amp)
  
  return(tbls)
  
}

ps_filelist <- ps_excels[1:20]

list_of_lists <- ps_filelist |> 
  map(\(ps_filelist) multi_read(ps_filelist))

joined_all_amp <- map(list_of_lists, ~ .x[["all_amp"]]) |> 
  list_rbind()

joined_pos_cnv <-  map(list_of_lists, ~ .x[["pos_cnv"]]) |> 
  list_rbind()
