
get_columns <- function(table_input) {
  
  output <- odbc::odbcConnectionColumns(
    conn = dbi_con, 
    catalog_name = "MolecularDB",
    schema_name = "dbo",
    name = table_input)
  
  return(output)
  
}

get_extraction_method <- function(sample_vector) {
  
  extraction_tbl_samples <- extraction_tbl |> 
    dplyr::filter(lab_no %in% sample_vector) |> 
    dplyr::collect()
  
  batches <- unique(extraction_tbl_samples$extraction_batch_fk)
  
  extraction_batch_info <- extraction_batch_tbl |> 
    dplyr::filter(extraction_batch_id %in% batches) |> 
    dplyr::collect() |> 
    # Remove DNA dilutions
    dplyr::filter(extraction_method_fk != 11) |>
    dplyr::left_join(extraction_method_key, 
                     dplyr::join_by(extraction_method_fk == extraction_method_id))
  
  output <- extraction_tbl_samples |> 
    dplyr::left_join(extraction_batch_info, 
                     dplyr::join_by(extraction_batch_fk == extraction_batch_id)) |> 
    dplyr::filter(!is.na(method_name)) |> 
    dplyr::rename(labno = lab_no)
  
  return(output)
  
}

get_sample_tissue <- function(sample_vector) {
  
  output <- sample_tbl |> 
    dplyr::select(-c(status_comment, comments, consultant_address, address1)) |> 
    dplyr::filter(labno %in% sample_vector) |> 
    dplyr::collect() |> 
    dplyr::mutate(tissue = as.numeric(tissue)) |> 
    dplyr::left_join(tissue_types, dplyr::join_by(tissue == tissue_type_id))
  
  return(output)
  
}

get_sample_gender <- function(sample_vector) {
  
  output <- sample_tbl |> 
    dplyr::select(labno, gender) |> 
    dplyr::filter(labno %in% sample_vector) |> 
    dplyr::collect() |> 
    dplyr::mutate(gender_string = case_when(
      
      gender == "1" ~"Male",
      
      gender == "2" ~"Female",
      
      gender == "9" ~"Unknown"))
  
  return(output)
 
}

get_sample_nhs_no <- function(sample_vector) {
  
  output <- sample_tbl |> 
    dplyr::select(labno, nhsno) |> 
    dplyr::filter(labno %in% sample_vector) |> 
    dplyr::collect()
    
  return(output)
  
}

ncc_regex <- stringr::regex(
  r"[
  (>\d{2}% | \d{2}-\d{2}%)
  ]",
  comments = TRUE
)

parse_ncc <- function(input_string) {
  
  # Function for parsing neoplastic cell content values from the comments
  # column of DNA Database
  
  ncc <- stringr::str_extract(string = input_string, 
                     pattern = ncc_regex, 
                     group = 1)
  
  return(ncc)
  
}

get_sample_ncc <- function(sample_vector) {
  
  output <- sample_tbl |> 
    dplyr::select(labno, comments) |> 
    dplyr::filter(labno %in% sample_vector) |> 
    dplyr::collect() |> 
    dplyr::mutate(ncc_db = parse_ncc(comments))
  
  return(output)

}

get_pathno <- function(sample_vector) {
  
  output <- sample_tbl |> 
    dplyr::select(labno, pathno) |> 
    dplyr::filter(labno %in% sample_vector) |> 
    dplyr::collect()
  
  return(output)
  
}
