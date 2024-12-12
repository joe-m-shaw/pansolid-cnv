# Collate WGS HTMLs for the deletions validation

library(tidyverse)
library(here)
library(readxl)

source(here("functions/wgs_html_functions.R"))

data_folder <- config::get("data_filepath")

wgs_amp_htmls <- list.files(path = paste0(data_folder,
                                          "validation/DOC6283_amplifications/",
                                          "raw/wgs/"),
                            full.names = TRUE)

wgs_del_htmls <- list.files(path = paste0(data_folder, "validation/",
                                          "DOC6567_deletions/raw/wgs/"),
                            full.names = TRUE)

wgs_htmls <- c(wgs_del_htmls, wgs_amp_htmls)

# Read HTML CNV data ----------------------------------------------------------------

wgs_domain1_cnvs <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_table_by_div_id(
    html_filepath = wgs_htmls, 
    div_id = "d_svcnv_tier1")) |> 
  list_rbind() |> 
  mutate(domain = 1)

wgs_domain2_cnvs <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_table_by_div_id(
    html_filepath = wgs_htmls, 
    div_id = "d_svcnv_tier2")) |> 
  list_rbind() |> 
  mutate(domain = 2)

wgs_domain3_cnvs <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_table_by_div_id(
    html_filepath = wgs_htmls, 
    div_id = "d_svcnv_tier3")) |> 
  list_rbind() |> 
  mutate(domain = 3)

# Collate WGS CNV data --------------------------------------------------------------

wgs_html_cnvs <- rbind(wgs_domain1_cnvs  |> 
                         select(-c(population_germline_allele_frequency,
                                   gene_mode_of_action)), 
                       wgs_domain2_cnvs |> 
                         select(-c(population_germline_allele_frequency,
                                   gene_mode_of_action)), wgs_domain3_cnvs) |> 
  mutate(cnv_class = parse_wgs_cnv_class(x = variant_type),
         
         cnv_copy_number = parse_wgs_cnv_copy_number(x = variant_type),
         
         chromosome = parse_wgs_html_grch38_coordinates(x = variant_gr_ch38_coordinates,
                                                        group = 1),
         
         cnv_start = as.numeric(parse_wgs_html_grch38_coordinates(x = variant_gr_ch38_coordinates,
                                                                  group = 2)),
         
         cnv_end = as.numeric(parse_wgs_html_grch38_coordinates(x = variant_gr_ch38_coordinates,
                                                                group = 4)),
         
         # In the HTMLS,  clinical indication specific genes are annotated with "*".
         # Remove this asterisk for easier filtering later on.
         gene = str_replace_all(string = gene, pattern = "\\*",
                                replacement = ""),
         
         category = "WGS result")

# Read HTML identifiers -------------------------------------------------------------

wgs_pids <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_pid_text(wgs_htmls)) |> 
  list_rbind()

wgs_headers <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_header(wgs_htmls)) |> 
  list_rbind()

wgs_tumour_details <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_table_by_div_id(
    html_filepath = wgs_htmls, 
    div_id = "t_tumour_details")) |> 
  list_rbind() |> 
  rename(wgs_pathno = histopathology_or_sihmds_lab_id) |> 
  select(filepath, wgs_pathno)

wgs_tumour_samples <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_table_by_div_id(
    html_filepath = wgs_htmls, 
    div_id = "t_tumour_sample")) |> 
  list_rbind() |> 
  select(filepath, calculated_overall_ploidy, calculated_chromosome_count, calculated_tumour_content)

if(anyNA.data.frame(wgs_pids)) {
  stop("NA values in WGS patient IDs")
}

if(anyNA.data.frame(wgs_headers)) {
  stop("NA values in WGS headers")
}

# Get lab numbers -------------------------------------------------------------------

wgs_pathway_tracker <- read_excel(path = paste0(data_folder, 
                                                "validation/DOC6567_deletions/",
                                                "excel_spreadsheets/",
                                                "WGS pathway tracker_copy_2024-12-05.xlsx"),
                                  sheet = "Cancer") |> 
  janitor::clean_names()

wgs_pathway_tracker_dna_no_df <- wgs_pathway_tracker |> 
  filter(!is.na(mol_db_number) & 
           !is.na(ngis_referral_id) &
           sample_type == "Solid tumour" &
           !duplicated(ngis_referral_id)) |> 
  mutate(labno = str_extract(string = mol_db_number,
                             pattern = "\\d{8}")) |> 
  select(labno, ngis_referral_id)

# Collate identifiers ---------------------------------------------------------------

wgs_html_ids <- inner_join(x = wgs_headers, y = wgs_pids, by = "filepath") |> 
  left_join(wgs_pathway_tracker_dna_no_df, join_by("wgs_r_no" == "ngis_referral_id"),
            relationship = "one-to-one") |> 
  left_join(wgs_tumour_details, by = "filepath") |> 
  left_join(wgs_tumour_samples, by = "filepath")

# Export collated information -------------------------------------------------------

write_csv(x = wgs_html_ids,
          file = paste0(data_folder, "validation/DOC6567_deletions/",
                        "processed/del_val_wgs_html_ids.csv"))

write_csv(x = wgs_html_cnvs,
          file = paste0(data_folder, "validation/DOC6567_deletions/",
                        "processed/del_val_wgs_html_cnvs.csv"))
