# Collate WGS HTMLs for the deletions validation

library(tidyverse)
library(here)

source(here("functions/wgs_html_functions.R"))
source(here("scripts/set_shared_drive_filepath.R"))

wgs_del_htmls <- list.files(path = paste0(data_folder, "validation/raw/wgs_deletions/"),
                            full.names = TRUE)

wgs_amp_htmls <- list.files(path = paste0(data_folder, "validation/raw/wgs_amplifications/"),
                            full.names = TRUE)

wgs_htmls <- c(wgs_del_htmls, wgs_amp_htmls)

wgs_html_cnvs <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_table_by_div_id(wgs_htmls, 
                                                  "d_svcnv_tier1")) |> 
  list_rbind() |> 
  mutate(cnv_class = parse_wgs_cnv_class(x = variant_type),
         
         cnv_copy_number = parse_wgs_cnv_copy_number(x = variant_type),
         
         chromosome = parse_wgs_html_grch38_coordinates(x = variant_gr_ch38_coordinates,
                                                        group = 1),
         
         cnv_start = as.numeric(parse_wgs_html_grch38_coordinates(x = variant_gr_ch38_coordinates,
                                                                  group = 2)),
         
         cnv_end = as.numeric(parse_wgs_html_grch38_coordinates(x = variant_gr_ch38_coordinates,
                                                                group = 4)),
         gene = str_replace_all(string = gene, pattern = "\\*",
                                replacement = ""))

write_csv(x = wgs_html_cnvs,
          file = paste0(data_folder, "validation/processed/wgs_html_cnvs_del_validation.csv"))
