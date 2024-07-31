# Reformat PanSolid and WGS data ----------------------------------------------------

library(here)
source(here("functions/cnv_functions.R"))
source(here("scripts/set_shared_drive_filepath.R"))

# Load collated data ----------------------------------------------------------------

pansolid_data_collated <- read_csv(file = paste0(data_folder, "collated_validation_data/",
                                                 "pansolid_results_collated.csv"))

wgs_data_collated <- read_csv(file = paste0(data_folder, "collated_validation_data/",
                                                 "wgs_html_cnvs.csv"))

pansolid_ids <- read_csv(file = paste0(data_folder, "collated_validation_data/",
                                       "pansolid_ids.csv"))

wgs_html_ids <- read_csv(file = paste0(data_folder, "collated_validation_data/",
                                       "wgs_html_ids.csv"))

# Reformatting functions ------------------------------------------------------------

reformat_wgs_cnv_result <- function(filepath) {
  
  del_gene_table <- read_excel(path = paste0(data_folder, 
                                             "pansolid_deletion_gene_list.xlsx"))
  
  sample_cnvs <- wgs_data_collated |> 
    filter(filepath == {{ filepath }}) |> 
    mutate(gene = str_replace_all(string = gene, pattern = "\\*",
                                  replacement = ""))
  
  sample_cnvs_filtered <- sample_cnvs |> 
    filter(cnv_class %in% c("LOSS", "DEL")) |> 
    filter(gene %in% del_gene_table$gene)
  
  output_table <- del_gene_table |> 
    mutate(wgs_result = case_when(
      gene %in% sample_cnvs_filtered$gene ~"Loss detected",
      !gene %in% sample_cnvs_filtered$gene ~"No loss detected"
    ))  |> 
    mutate(filepath = filepath) |> 
    relocate(filepath)
  
  return(output_table)
  
}

reformat_pansolid_cnv_result <- function(filepath) {
  
  del_gene_table <- read_excel(path = paste0(data_folder, 
                                             "pansolid_deletion_gene_list.xlsx"))
  
  sample_cnvs <- pansolid_data_collated |> 
    filter(filepath == {{ filepath }}) |> 
    filter(gene %in% del_gene_table$gene)
  
  output_table <- del_gene_table |> 
    mutate(pansolid_result = case_when(
      gene %in% sample_cnvs$gene ~"Loss detected",
      !gene %in% sample_cnvs$gene ~"No loss detected"
    )) |> 
    mutate(filepath = filepath) |> 
    relocate(filepath)
  
  return(output_table)
  
}

# Reformat data ---------------------------------------------------------------------

wgs_htmls <- list.files(path = paste0(data_folder, "wgs_result_htmls/"),
                        full.names = TRUE,
                        pattern = "*.html")

folder_path <- "S:/central shared/Genetics/NGS/Bioinformatics/1_Pan-solid-Cancer/CNV/Deletions/v1_Coarse_FC_1.33/"

pansolid_files <- list.files(path = folder_path, pattern = ".xlsx",
                             full.names = TRUE)

wgs_data_reformatted <- wgs_htmls |> 
  map(\(wgs_htmls) reformat_wgs_cnv_result(wgs_htmls)) |> 
  list_rbind()

pansolid_data_reformatted <- pansolid_files |> 
  map(\(pansolid_files) reformat_pansolid_cnv_result(pansolid_files)) |> 
  list_rbind()

# Add identifiers -------------------------------------------------------------------

ps_reformatted_with_ids <- pansolid_ids |> 
  left_join(pansolid_data_reformatted, by = "filepath",
            relationship = "one-to-many")

wgs_reformatted_with_ids <- wgs_html_ids |> 
  left_join(wgs_data_reformatted, by = "filepath",
            relationship = "one-to-many")

# Compare results -------------------------------------------------------------------

comparison <- ps_reformatted_with_ids |> 
  filter(nhsno %in% wgs_reformatted_with_ids$nhs_no_clean) |> 
  left_join(wgs_reformatted_with_ids |> 
              select(nhs_no_clean, gene, wgs_result), join_by(nhsno == nhs_no_clean,
                                              gene == gene)) |> 
  mutate(outcome = case_when(
    wgs_result == "Loss detected" & pansolid_result == "Loss detected" ~"True positive",
    
    wgs_result == "No loss detected" & pansolid_result == "No loss detected" ~"True negative",
    
    wgs_result == "Loss detected" & pansolid_result == "No loss detected" ~"False negative",
    
    wgs_result == "No loss detected" & pansolid_result == "Loss detected" ~"False positive"
    
  ))

comparison |> 
  group_by(gene, outcome) |> 
  count()
