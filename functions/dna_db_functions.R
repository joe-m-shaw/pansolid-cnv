
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
    filter(lab_no %in% sample_vector) |> 
    collect()
  
  batches <- unique(extraction_tbl_samples$extraction_batch_fk)
  
  extraction_batch_info <- extraction_batch_tbl |> 
    filter(extraction_batch_id %in% batches) |> 
    collect() |> 
    # Remove DNA dilutions
    filter(extraction_method_fk != 11) |>
    left_join(extraction_method_key, join_by(extraction_method_fk == extraction_method_id))
  
  output <- extraction_tbl_samples |> 
    left_join(extraction_batch_info, join_by(extraction_batch_fk == extraction_batch_id)) |> 
    filter(!is.na(method_name)) |> 
    rename(labno = lab_no)
  
  return(output)
  
}

get_sample_tissue <- function(sample_vector) {
  
  output <- sample_tbl |> 
    select(-c(status_comment, comments, consultant_address, address1)) |> 
    filter(labno %in% sample_vector) |> 
    collect() |> 
    mutate(tissue = as.numeric(tissue)) |> 
    left_join(tissue_types, join_by(tissue == tissue_type_id))
  
  return(output)
  
}

get_sample_gender <- function(sample_vector) {
  
  output <- sample_tbl |> 
    select(labno, gender) |> 
    filter(labno %in% sample_vector) |> 
    collect() |> 
    mutate(gender_string = case_when(
      
      gender == "1" ~"Male",
      
      gender == "2" ~"Female",
      
      gender == "9" ~"Unknown"))
  
  return(output)
 
}

get_sample_nhs_no <- function(sample_vector) {
  
  output <- sample_tbl |> 
    select(labno, nhsno) |> 
    filter(labno %in% sample_vector) |> 
    collect()
    
  return(output)
  
}

ncc_regex <- regex(
  r"[
  (>\d{2}% | \d{2}-\d{2}%)
  ]",
  comments = TRUE
)

parse_ncc <- function(input_string) {
  
  # Function for parsing neoplastic cell content values from the comments
  # column of DNA Database
  
  ncc <- str_extract(string = input_string, 
                     pattern = ncc_regex, 
                     group = 1)
  
  return(ncc)
  
}

get_sample_ncc <- function(sample_vector) {
  
  output <- sample_tbl |> 
    select(labno, comments) |> 
    filter(labno %in% sample_vector) |> 
    collect() |> 
    mutate(ncc_db = parse_ncc(comments))
  
  return(output)

}

get_pathno <- function(sample_vector) {
  
  output <- sample_tbl |> 
    select(labno, pathno) |> 
    filter(labno %in% sample_vector) |> 
    collect()
  
  return(output)
  
}
