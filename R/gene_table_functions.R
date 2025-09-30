load_pansolid_gene_table <- function(cnv_type) {
  
  #' Load a PanSolid copy number variant gene table
  #'
  #' @param cnv_type The type of copy number variant, which must be "Deletions" or "Amplifications"
  #'
  #' @return Returns a dataframe of the current gene list.
  #' 
  #' @export

  if(!cnv_type %in% c("Deletions", "Amplifications")) {
    stop("cnv_type must be Deletions or Amplifications")
  }
  
  if(cnv_type == "Deletions") {
    
    gene_table <- readr::read_csv(file = paste0(config::get("data_folderpath"), 
                                                "validation/",
                                                "DOC6283_amplifications/",
                                                "gene_lists/",
                                                "pansolid_deletion_gene_list.csv"),
                                  show_col_types = FALSE)
    
  }
  
  if(cnv_type == "Amplifications") {
    
    gene_table <- readr::read_csv(file = paste0(config::get("data_folderpath"), 
                                                "validation/",
                                                "DOC6283_amplifications/",
                                                "gene_lists/",
                                                "pansolid_amplification_gene_list.csv"),
                                  show_col_types = FALSE)
    
  }
  
  return(gene_table)
  
}
