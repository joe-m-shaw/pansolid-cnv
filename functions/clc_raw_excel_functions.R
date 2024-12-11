# CLC raw output functions

# This collection of functions are useful for working with Microsoft Excel
# files generated from the CLC Genomics Workbench (Qiagen). These files are
# referred to as "raw outputs" because no further formatting is added.
# This is in contrast to "annotated" Excel outputs, which are the final
# versions used by clinical scientists for interpretting copy number variants
# (CNVs).

make_empty_cnv_df <- function() {
  
  #' Create an empty dataframe for CLC outputs with no CNVs
  #' 
  #' When no copy number variants (CNVs) are detected with a 
  #' sample, the tab of the raw CLC output is empty,
  #' without column headings. This can make it difficult 
  #' when collating results from multiple files. This function 
  #' creates a placeholder dataframe with the correct
  #'  column headings to show that no CNVs have been identified.
  #'
  #' @return A dataframe with default columns for CLC exports
  #' @export
  #'
  #' @examples
  
  df <- data.frame(
    "chromosome" = "",
    "region" = "",
    "name" = "",
    "region_length" = "",
    "type" = "",
    "source" = "",
    "id" = "",
    "gene_id" = "",
    "hgnc" = "",
    "mim" = "",
    "description" = "",
    "gbkey" = "",
    "gene" = "No CNV detected",
    "gene_biotype" = "",
    "gene_synonym" = "",
    "cnv_region" = "",
    "cnv_region_length" = "",
    "consequence" = "",
    "fold_change_adjusted" = "",
    "p_value" = "",
    "number_of_targets" = "",
    "comments" = "",
    "targets" = "")
  
  return(df)
  
}

get_cnv_df_cols <- function() {
  
  column_names <- c("chromosome", "region", "name", "region_length", 
                    "type", "source", "id", "gene_id", "hgnc",
                    "mim", "description", "gbkey", "gene", "gene_biotype", 
                    "gene_synonym", "cnv_region", "cnv_region_length",
                    "consequence", "fold_change_adjusted", "p_value", 
                    "number_of_targets", "comments", "targets")
  
  return(column_names)
  
}

set_cnv_df_types <- function(df) {
  
  cols <- get_cnv_df_cols()
  
  if(length(setdiff(cols, colnames(df))) > 0) {
    stop("df does not have standard columns")
  }
  
  df <- df |> 
    dplyr::mutate(
      chromosome = as.character(chromosome),
      region = as.character(region),
      name = as.character(name),
      region_length = as.character(region_length),
      type = as.character(type),
      source = as.character(source),
      id = as.character(id),
      gene_id = as.character(gene_id),
      hgnc = as.character(hgnc),
      mim = as.character(mim),
      description = as.character(description),
      gbkey = as.character(gbkey),
      gene = as.character(gene),
      gene_biotype = as.character(gene_biotype),
      gene_synonym = as.character(gene_synonym),
      cnv_region = as.character(cnv_region),
      cnv_region_length = as.numeric(cnv_region_length),
      consequence = as.character(consequence),
      fold_change_adjusted = as.double(fold_change_adjusted),
      p_value = as.double(p_value),
      number_of_targets = as.numeric(number_of_targets),
      comments = as.character(comments),
      targets = as.character(targets)
    )
  
  return(df)
}

read_raw_cnv_excel <- function(filepath, sheet_no) {
  
  df <- readxl::read_excel(path = filepath,
                   sheet = sheet_no) |> 
    janitor::clean_names()
  
  if(nrow(df) == 0){
    
    df <- make_empty_cnv_df()
    
  }
  
  if(!"comments" %in% colnames(df)) {
    
    df <- df |> 
      dplyr::mutate(comments = "") |> 
      dplyr::relocate(comments, .after = number_of_targets)
    
  }
  
  df <- set_cnv_df_types(df) |> 
    dplyr::mutate(filepath = filepath)
  
  return(df)
  
}

read_del_raw_excel <- function(filepath, sheet_no) {
  
  df <- read_raw_cnv_excel(filepath = filepath, sheet_no = sheet_no) |> 
    dplyr::mutate(labno = str_extract(string = filepath, pattern = "\\d{8}"),
           suffix = str_extract(string = filepath, pattern = ".*\\d{8}(|a|b|c|d)_.*",
                                group = 1),
           worksheet = str_extract(string = filepath, pattern = "WS\\d{6}")) |> 
    dplyr::relocate(labno, worksheet)
  
  return(df)
  
}

