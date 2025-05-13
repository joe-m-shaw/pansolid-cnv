parse_filename <- function(input_file, input_group) {
  
  #' Parse an element from a PanSolid results Excel filename
  #'
  #' @param input_file The filename to parse
  #' @param input_group The number of the group to extract, defined by filename_regex
  #'
  #' @return The string of the selected regex group
  #' 
  #' @examples filename <- "Annotated_WS123456_12345678a_JoeShaw.xlsx"
  #' 
  #' worksheet <- parse_filename(filename, 1)
  
  filename_regex <- stringr::regex(
    r"[
    (WS\d{6})             # Worksheet number
    _
    (\d{5,8})             # Lab number
    (a|b|c|d|)            # Suffix
    _
    ([:alnum:]{2,30})     # Patient name (or alphanumeric identifier)
    (.*.xlsx|.xlsx)       # Variable ending                      
    ]",
    comments = TRUE)
  
  output <- stringr::str_extract(input_file, filename_regex,
                                 group = input_group)
  
  return(output)
  
}

filename_to_df <- function(file) {
  
  #' Convert PanSolid results Excel filename identifiers to a dataframe
  #'
  #' @param file The filename or full filepath
  #'
  #' @return A dataframe of patient filename identifiers
  #'
  #' @examples filename <- "Annotated_WS123456_12345678a_JoeShaw.xlsx"
  #' 
  #' id_df <- filename_to_df(filename)
  
  output <- data.frame(
    worksheet = c(parse_filename(file, 1)),
    labno = c(as.character(parse_filename(file, 2))),
    suffix = c(parse_filename(file, 3)),
    patient_name = c(parse_filename(file, 4))) |> 
    dplyr::mutate(
      labno_suffix = stringr::str_c(labno, suffix),
      labno_suffix_worksheet = stringr::str_c(labno_suffix, "_", worksheet))
  
  return(output)
  
}

add_identifiers <- function(file, tbl) {
  
  #' Add identifiers in a filename onto a dataframe
  #'
  #' @param file The filename with identifiers
  #' @param tbl The table to add identifiers to
  #'
  #' @return A dataframe of the initial table joined to the patient identifiers
  #'
  #' @examples filename <- "Annotated_WS123456_12345678a_JoeShaw.xlsx"
  #' 
  #' data <- data.frame(gene = c("ERBB2"), dosage = c(2))
  #' 
  #' data_with_ids <- add_identifiers(file = filename, tbl = data)
  
  identifiers <- filename_to_df(file)
  
  labno <- parse_filename(file, 2)
  
  output <- tbl |> 
    dplyr::mutate(labno = labno, 
                  filepath = file) |>  
    dplyr::left_join(identifiers, by = "labno") |> 
    dplyr::relocate(worksheet, labno, suffix, patient_name, 
                    labno_suffix, labno_suffix_worksheet, filepath)
  
  return(output)
  
}

source(here::here("tests/test_filename_functions.R"))
