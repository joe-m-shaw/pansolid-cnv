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

load_gene_table <- function(cnv_type) {
  
  if(!cnv_type %in% c("Deletions", "Amplifications")) {
    stop("cnv_type must be Deletions or Amplifications")
  }
  
  if(cnv_type == "Deletions") {
    
    gene_table <- read_excel(path = paste0(data_folder, 
                                           "pansolid_deletion_gene_list.xlsx"))
    
  }
  
  if(cnv_type == "Amplifications") {
    
    gene_table <- read_excel(path = paste0(data_folder, 
                                           "pansolid_amplification_gene_list.xlsx"))
    
  }
  
  return(gene_table)
  
}



reformat_wgs_cnv_result <- function(filepath, cnv_type) {
  
  gene_table <- load_gene_table({{ cnv_type }})
  
  sample_cnvs <- wgs_data_collated |> 
    filter(filepath == {{ filepath }}) |> 
    mutate(gene = str_replace_all(string = gene, pattern = "\\*",
                                  replacement = ""))
  
  if(cnv_type == "Deletions") {
    
    sample_cnvs_filtered <- sample_cnvs |> 
      filter(cnv_class %in% c("LOSS", "DEL")) |> 
      filter(gene %in% gene_table$gene)
    
  }
  
  if(cnv_type == "Amplifications") {
    
    sample_cnvs_filtered <- sample_cnvs |> 
      filter(cnv_class %in% c("DUP", "GAIN") 
             #& cnv_copy_number > 6
             ) |> 
      filter(gene %in% gene_table$gene)
    
  }
  
  if(cnv_type == "Deletions") {
    
    gene_result_table <- gene_table |> 
      mutate(wgs_result = case_when(
        gene %in% sample_cnvs_filtered$gene ~"Loss detected",
        !gene %in% sample_cnvs_filtered$gene ~"No loss detected"
      )) 
    
  }
  
  if(cnv_type == "Amplifications") {
    
    gene_result_table <- gene_table |> 
      mutate(wgs_result = case_when(
        gene %in% sample_cnvs_filtered$gene ~"Gain detected",
        !gene %in% sample_cnvs_filtered$gene ~"No gain detected"
      )) 
    
  }
  
  output_table <- gene_result_table |> 
    mutate(filepath = filepath) |> 
    relocate(filepath)
  
  return(output_table)
  
}

reformat_pansolid_cnv_result <- function(filepath, cnv_type) {
  
  gene_table <- load_gene_table({{ cnv_type }})
  
  sample_cnvs <- pansolid_data_collated |> 
    filter(filepath == {{ filepath }}) |> 
    filter(gene %in% gene_table$gene)
  
  if(cnv_type == "Deletions") {
    
    gene_result_table <- gene_table |> 
    mutate(pansolid_result = case_when(
      gene %in% sample_cnvs$gene ~"Loss detected",
      !gene %in% sample_cnvs$gene ~"No loss detected"
    )) 
    
  }
  
  if(cnv_type == "Amplifications") {
    
    gene_result_table <- gene_table |> 
      mutate(pansolid_result = case_when(
        gene %in% sample_cnvs$gene ~"Gain detected",
        !gene %in% sample_cnvs$gene ~"No gain detected"
      )) 
    
  }
  
  output_table <- gene_result_table |> 
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

wgs_del_data_reformatted <- wgs_htmls |> 
  map(\(wgs_htmls) reformat_wgs_cnv_result(filepath = wgs_htmls,
                                           cnv_type = "Deletions")) |> 
  list_rbind()

wgs_amp_data_reformatted <- wgs_htmls |> 
  map(\(wgs_htmls) reformat_wgs_cnv_result(filepath = wgs_htmls,
                                           cnv_type = "Amplifications")) |> 
  list_rbind()

pansolid_del_data_reformatted <- pansolid_files |> 
  map(\(pansolid_files) reformat_pansolid_cnv_result(filepath = pansolid_files,
                                                     cnv_type = "Deletions")) |> 
  list_rbind()

pansolid_amp_data_reformatted <- pansolid_files |> 
  map(\(pansolid_files) reformat_pansolid_cnv_result(filepath = pansolid_files,
                                                     cnv_type = "Amplifications")) |> 
  list_rbind()

wgs_data_reformatted_collated <- rbind(wgs_del_data_reformatted, 
                                       wgs_amp_data_reformatted)

pansolid_data_reformated_collated <- rbind(pansolid_del_data_reformatted,
                                           pansolid_amp_data_reformatted)

# Add identifiers -------------------------------------------------------------------

ps_reformatted_with_ids <- pansolid_ids |> 
  left_join(pansolid_data_reformated_collated, by = "filepath",
            relationship = "one-to-many")

wgs_reformatted_with_ids <- wgs_html_ids |> 
  left_join(wgs_data_reformatted_collated, by = "filepath",
            relationship = "one-to-many")
