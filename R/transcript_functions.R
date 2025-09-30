read_ensembl_exon_table <- function(filepath) {
  
  transcript_regex <- stringr::regex(
    r"(
    .+
    (ENST\d{11})
    .csv
    )",
    comments = TRUE
  )
  
  transcript_id <- stringr::str_extract(string = filepath, 
                               pattern = transcript_regex,
                               group = 1)
  
  table <- readr::read_csv(file = filepath,
                    show_col_types = FALSE) |> 
    janitor::clean_names() |> 
    dplyr::filter(!is.na(no)) |> 
    dplyr::rename(exon = no) |> 
    dplyr::mutate(transcript = transcript_id) |> 
    dplyr::relocate(transcript) |> 
    dplyr::select(-sequence)
  
  return(table)
  
}
