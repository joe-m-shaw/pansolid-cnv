
source(here("functions/gene_table_functions.R"))

reformat_wgs_cnv_result <- function(filepath, cnv_type, 
                                    wgs_tbl = wgs_html_cnvs) {
  
  #' Reformat whole genome sequencing copy number variant results for comparison
  #' with PanSolid results
  #'
  #' @param filepath The filepath of the whole genome sequencing HTML file
  #' @param cnv_type The type of CNV, which must be "Deletions" or "Amplifications"
  #' @param wgs_tbl A dataframe of collated WGS CNV results, defaults to the 
  #' collated data from the collate_wgs_html.R script
  #'
  #' @return A dataframe of only the genes specified in the PanSolid gene list 
  #' along with the result on WGS
  #' @export
  #'
  #' @examples reformat_wgs_cnv_result(filepath = wgs_htmls[1], 
  #' cnv_type = "Amplifications",
  #' wgs_tbl = wgs_data_collated)
  
  gene_table <- load_pansolid_gene_table({{ cnv_type }})
  
  sample_cnvs <- wgs_tbl |> 
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
             & cnv_copy_number > 10
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
