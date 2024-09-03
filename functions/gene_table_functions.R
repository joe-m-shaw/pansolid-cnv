load_pansolid_gene_table <- function(cnv_type) {
  
  #' Load a PanSolid copy number variant gene table
  #'
  #' @param cnv_type The type of copy number variant, which must be "Deletions" or "Amplifications"
  #'
  #' @return Returns a dataframe of the current gene list.
  #' @export
  #'
  #' @examples amp_list <- load_gene_table("Amplifications")
  
  if(!cnv_type %in% c("Deletions", "Amplifications")) {
    stop("cnv_type must be Deletions or Amplifications")
  }
  
  if(cnv_type == "Deletions") {
    
    gene_table <- read_excel(path = paste0(data_folder, 
                                           "gene_lists/pansolid_deletion_gene_list.xlsx"))
    
  }
  
  if(cnv_type == "Amplifications") {
    
    gene_table <- read_excel(path = paste0(data_folder, 
                                           "gene_lists/pansolid_amplification_gene_list.xlsx"))
    
  }
  
  return(gene_table)
  
}
