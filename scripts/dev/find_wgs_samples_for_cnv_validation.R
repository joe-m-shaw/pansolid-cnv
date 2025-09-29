# Finding WGS Samples for the CNV Validation

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)

# Functions -------------------------------------------------------------------------

data_folder <- config::get("data_filepath")

source(here("scripts/connect_to_dna_db.R"))
source(here("functions/dna_db_functions.R"))
source(here("functions/wgs_html_functions.R"))

# Identify samples tested on PanSolid -----------------------------------------------

all_worksheets <- dna_db_worksheets |> 
  select(pcrid, date, description) |> 
  collect() |> 
  mutate(worksheet = paste0("WS", pcrid))

stopifnot(nrow(all_worksheets) > 0)

ps_ws_info <- all_worksheets |> 
  filter(pcrid >= 147437) |>
  filter(grepl(pattern = "pansolid|pan-solid|pan_solid|pan solid", 
               x = description,
               ignore.case = TRUE)) |> 
  mutate(ps_category = case_when(
    grepl(pattern = "jBRCA|j_BRCA|j-BRCA|jew",
          x = description,
          ignore.case = TRUE) ~"PanSolid Jewish BRCA",
    TRUE ~"PanSolid FFPE"
  ))

mlpa_ws_info <- all_worksheets |> 
  filter(pcrid >= 140721) |>
  filter(grepl(pattern = "P\\d{3}|PO\\d{2}", 
               x = description,
               ignore.case = TRUE)) 

mlpa_ws <- mlpa_ws_info$pcrid

stopifnot(anyNA.data.frame(ps_ws_info) == FALSE)

ps_ws_ffpe_only <- ps_ws_info |> 
  filter(ps_category == "PanSolid FFPE")

pcrid_list <- ps_ws_ffpe_only$pcrid

pansolid_worksheet_samples <- dna_db_pcr_records |> 
  filter(pcrid %in% pcrid_list) |> 
  select(pcrid, sample, name) |> 
  collect() |> 
  rename(labno = sample) |> 
  arrange(desc(labno))

mlpa_worksheet_samples <- dna_db_pcr_records |> 
  filter(pcrid %in% mlpa_ws) |> 
  select(pcrid, sample, name) |> 
  collect() |> 
  rename(labno = sample) |> 
  left_join(mlpa_ws_info |> 
              select(pcrid, description),
            by = "pcrid")

samples_on_ps_worksheets <- pansolid_worksheet_samples$labno

pansolid_sample_info <- sample_tbl |> 
  filter(labno %in% samples_on_ps_worksheets) |> 
  select(labno, firstname, surname, i_gene_r_no,
         pathno, nhsno) |> 
  collect()

# Finding WGS samples with QIAsymphony extractions ----------------------------------

# We want to identify samples which have had/are having whole genome sequencing
# and also have a DNA sample extracted by the QIASymphony method from the same or 
# similar pathology block to the whole genome sequencing sample.

## Find WGS samples -----------------------------------------------------------------

wgs_pathway_tracker <- read_excel(path = paste0(config::get("data_folderpath"),
                                                "validation/DOC6791_chromosome_arms/",
                                                "2025-07-16_WGS pathway tracker.xlsx"),
                                  sheet = "Cancer") |> 
  janitor::clean_names()

wgs_pathway_tracker_dna_no_df <- wgs_pathway_tracker |> 
  filter(!is.na(mol_db_number)) |> 
  mutate(labno = str_extract(string = mol_db_number,
                             pattern = "\\d{8}"),
         nhsno = str_replace_all(string = nhs_number,
                                 pattern = " ",
                                 replacement = ""))

wgs_labnos <- wgs_pathway_tracker_dna_no_df$mol_db_number

## Find WGS NHS numbers -------------------------------------------------------------

wgs_nhs_numbers_df <- sample_tbl |> 
  filter(labno %in% wgs_labnos) |> 
  select(labno, nhsno) |> 
  collect() |> 
  filter(!is.na(nhsno))

wgs_nhsnos <- wgs_nhs_numbers_df$nhsno

## Find WGS samples tested on PanSolid ----------------------------------------------

all_samples_from_wgs_patients_df <-  sample_tbl |> 
  filter(nhsno %in% wgs_nhsnos) |> 
  select(labno, nhsno, firstname, surname, pathno) |> 
  collect()

wgs_samples_tested_on_pansolid <- pansolid_worksheet_samples |> 
  filter(labno %in% all_samples_from_wgs_patients_df$labno) |> 
  left_join(pansolid_sample_info |> 
              filter(!is.na(nhsno)),
            by = "labno") |> 
  left_join(wgs_pathway_tracker_dna_no_df |> 
              rename(wgs_labno = labno) |> 
              select(wgs_labno, nhsno, clinical_indication_code_m_code) |> 
              filter(!is.na(nhsno)),
            by = "nhsno") |> 
  filter(!duplicated(labno)) |> 
  arrange(pcrid)

write_csv(wgs_samples_tested_on_pansolid,
          "wgs_samples_tested_on_pansolid.csv")

wgs_samples_tested_on_pansolid |> 
  group_by(clinical_indication_code_m_code) |> 
  count() |> 
  arrange(desc(n)) |> view()

mlpa_samples_tested_on_pansolid <- pansolid_worksheet_samples |> 
  filter(labno %in% mlpa_worksheet_samples$labno) |> 
  left_join(mlpa_worksheet_samples,
            by = "labno") |> 
  left_join
  filter(!duplicated(labno)) |> 
  mutate(mlpa_kit = str_extract(string = description,
                                pattern = "P\\d{3}"))

mlpa_samples_tested_on_pansolid |> 
  group_by(mlpa_kit) |> 
  count() |> 
  arrange(desc(n)) |>  view()


## Find WGS samples with QIAsymphony extractions ------------------------------------

all_samples_from_wgs_patients <- all_samples_from_wgs_patients_df$labno

wgs_samples_extraction_info <- get_extraction_method(sample_vector =
                                                       all_samples_from_wgs_patients)

wgs_other_files <- list.files(path = paste0(data_folder,
                                            "validation/DOC6567_deletions/raw/wgs_other/"),
                              full.names = TRUE)

wgs_amp_files <- list.files(path = paste0(data_folder,
                                            "validation/DOC6283_amplifications/raw/wgs/"),
                              full.names = TRUE)

wgs_del_files <- list.files(path = paste0(data_folder,
                                          "validation/DOC6567_deletions/raw/wgs/"),
                            full.names = TRUE)

previous_wgs_files <- c(wgs_other_files, wgs_amp_files, wgs_del_files)

previous_wgs_pids <- previous_wgs_files |> 
  map(\(previous_wgs_files) parse_wgs_html_pid_text(previous_wgs_files)) |> 
  list_rbind()

all_samples_from_wgs_patients_df_with_extraction <- all_samples_from_wgs_patients_df |> 
  left_join(wgs_samples_extraction_info |> 
              select(labno, method_name), by = "labno") |> 
  # Remove blood sample extractions
  filter(!method_name %in% c("Chemagen 360", "Saliva Chemagic 360-D", "FFPE RNA")) |> 
  arrange(nhsno) |> 
  filter(method_name != "COBAS" & pathno != "")

# Check HTMLs for deletions in relevant genes ----------------------------------------------------

source(here("functions/gene_table_functions.R"))
source(here("functions/wgs_html_functions.R"))

del_genes <- load_pansolid_gene_table("Deletions")

wgs_amp_htmls <- list.files(path = paste0(data_folder, 
                                      "validation/raw/wgs_amplifications/"),
                        full.names = TRUE,
                        pattern = "*.html")

wgs_del_htmls <- list.files(path = paste0(data_folder, 
                                          "validation/raw/wgs_deletions/"),
                            full.names = TRUE,
                            pattern = "*.html")

wgs_del_htmls_empty_tbls_removed <- wgs_del_htmls[-c(6, 9, 28)]

wgs_htmls <- c(wgs_amp_htmls, wgs_del_htmls_empty_tbls_removed)

wgs_domain1_cnvs <- wgs_del_htmls_empty_tbls_removed |> 
  map(\(wgs_del_htmls_empty_tbls_removed) parse_wgs_html_table_by_div_id(
    html_filepath = wgs_del_htmls_empty_tbls_removed, 
    div_id = "d_svcnv_tier1")) |> 
  list_rbind() |> 
  mutate(domain = 1) |> 
  mutate(cnv_class = parse_wgs_cnv_class(x = variant_type),
         cnv_copy_number = parse_wgs_cnv_copy_number(x = variant_type)) 

wgs_pids <- wgs_del_htmls_empty_tbls_removed |> 
  map(\(wgs_del_htmls_empty_tbls_removed) parse_wgs_html_pid_text(wgs_del_htmls_empty_tbls_removed)) |> 
  list_rbind()

wgs_tumour_samples <- wgs_del_htmls_empty_tbls_removed |> 
  map(\(wgs_del_htmls_empty_tbls_removed) parse_wgs_html_table_by_div_id(
    html_filepath = wgs_del_htmls_empty_tbls_removed, 
    div_id = "t_tumour_sample")) |> 
  list_rbind() |> 
  select(filepath, calculated_overall_ploidy, calculated_chromosome_count, calculated_tumour_content)

wgs_deletions_of_interest <- wgs_domain1_cnvs |> 
  mutate(gene = str_replace_all(string = gene, pattern = "\\*",
                         replacement = "")) |> 
  filter(cnv_class %in% c("DEL", "LOSS") &
           gene %in% del_genes$gene) |> 
  left_join(wgs_pids, by = "filepath") |> 
  left_join(wgs_tumour_samples, by = "filepath") |> 
  relocate(patient_name, nhsno, filepath, gene, cnv_class, cnv_copy_number, 
           calculated_overall_ploidy, calculated_tumour_content, confidence_support) |> 
  filter(grepl(pattern = "HC", x = confidence_support))

# HC - high confidence copy number variant (confidence score ≥ 10)
# LC - low confidence copy number variant (confidence score < 10)

wgs_deletions_of_interest |> 
  count(gene, cnv_class) |> 
  arrange(n)

wgs_deletions_of_interest |> 
  #filter(gene == "CDKN2A") |> 
  select(patient_name, gene, cnv_class, cnv_copy_number, calculated_tumour_content) |> 
  arrange(desc(calculated_tumour_content)) |>  view()


