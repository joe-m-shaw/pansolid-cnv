read_ensembl_exon_table <- function(filepath) {
  
  transcript_regex <- regex(
    r"(
    .+
    (ENST\d{11})
    .csv
    )",
    comments = TRUE
  )
  
  transcript_id <- str_extract(string = filepath, 
                               pattern = transcript_regex,
                               group = 1)
  
  table <- read_csv(file = filepath,
                    show_col_types = FALSE) |> 
    janitor::clean_names() |> 
    filter(!is.na(no)) |> 
    rename(exon = no) |> 
    mutate(transcript = transcript_id) |> 
    relocate(transcript) |> 
    select(-sequence)
  
  return(table)
  
}
