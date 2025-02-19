parse_wgs_html_header <- function(html_filepath) {
  
  #' Parse the header element of a whole genome sequencing HTML
  #'
  #' @param html_filepath The full file-path of the HTML to parse
  #'
  #' @return The header identifiers as a data-frame
  #'
  #' @examples patient_ids <- parse_wgs_html_header(html_filepath)
  
  html <- rvest::read_html(x = html_filepath)
  
  header <- html |> 
    rvest::html_element("header") |> 
    rvest::html_text2()
  
  header_regex <- stringr::regex(
    r"[
    .
    Report\sversion:\sv(\d{1,2}.\d{1,2}.\d{1,2}|\d{1,2}.\d{1,2})
    .+
    Issue\sdate:\s(\d{2}-\d{2}-\d{4})
    .+
    (r\d{11}) # referral ID
    .+
    (p\d{11}) # patient ID
    .+
    ]",
    comments = TRUE
  )
  
  output <- data.frame( 
    "filepath" = html_filepath,
    "wgs_r_no" = stringr::str_extract(header, header_regex, group = 3),
    "wgs_p_no" = stringr::str_extract(header, header_regex, group = 4),
    "wgs_version" = stringr::str_extract(header, header_regex, group = 1),
    "wgs_analysis_date" = stringr::str_extract(header, header_regex, group = 2))
  
  return(output)
  
}

parse_wgs_html_pid_text <- function(html_filepath) {
  
  #' Parse the patient identifier text from a whole genome sequencing HTML
  #'
  #' @param html_filepath The full file-path of the HTML to parse
  #'
  #' @return The patient identifiers as a data-frame
  #'
  #' @note This text appears as 3 rows above the referral ID table on a whole
  #' genome sequencing HTML. The referral ID table itself needs to be parsed
  #' using the parse_wgs_html_table_by_number function.
  #'
  #' @examples pid_df <- parse_wgs_html_pid_text(html_filepath)
  
  html <- rvest::read_html(x = html_filepath)
  
  pid_text <- html |> 
    rvest::html_element("#pid") |> 
    rvest::html_text2()
  
  pid_regex <- stringr::regex(
    r"[
    ^Name:\s
    (.{4,30})
    \n
    Date\sof\sBirth:\s
    (\d{1,2}-\w{3}-\d{4})
    \n
    NHS\sNo.:\s
    (\d{3}\s\d{3}\s\d{4})
    ]",
    comments = TRUE
  )
  
  output <- data.frame( 
    "filepath" = html_filepath,
    "patient_name" = stringr::str_extract(string = pid_text, 
                                          pattern = pid_regex,
                                          group = 1),
    "patient_dob" = stringr::str_extract(string = pid_text, 
                                         pattern = pid_regex,
                                         group = 2),
    "nhsno_raw" = stringr::str_extract(string = pid_text, 
                                       pattern = pid_regex,
                                       group = 3)) |> 
    mutate(nhsno = stringr::str_replace_all(string = nhsno_raw,
                                            pattern = " ",
                                            replacement = ""))
  
  return(output)
  
}

parse_wgs_html_table_by_div_id <- function(html_filepath,
                                           div_id) {
  
  #' Parse whole genome sequencing HTML files using a CSS "div_id" identifier.
  #'
  #' @param html_filepath The filepath for the WGS HTML file
  #' @param div_id The CSS identifier for the table (examples below)
  #' 
  #' @return Returns the table from the HTML as a tibble with column names 
  #' in snake-case.
  #' 
  #' @section Useful Identifiers: 
  #' "t_tumour_details", "t_tumour_sample", 
  #' "t_germline_sample", "t_quality", "tier1" (Tier 1 sequence variants),
  #' "svcnv_tier1" (Tier 1 structural and copy number variants, for v2.28 
  #' and earlier), 
  #' "d_svcnv_tier1" (Tier 1 structural and copy number variants,later versions
  #'  than v2.28)
  #'
  #' @examples parse_wgs_html_table_by_div_id(filepath, "t_tumour_details")
  
  html <- rvest::read_html(x = html_filepath)
  
  output_table <- html |> 
    rvest::html_element( str_c("#", {{ div_id }} )) |> 
    rvest::html_table() |> 
    janitor::clean_names() |> 
    dplyr::mutate(filepath = html_filepath)
  
  if(nrow(output_table) == 0) {
    warning("No rows in table")
    return()
  }
  
  return(output_table)
  
}

parse_wgs_html_table_by_number <- function(html_filepath,
                                           table_number) {
  
  #' Parse whole genome sequencing HTML files using the table number
  #'
  #' @param html_filepath The filepath for the WGS HTML file
  #' @param table_number The number of the table to select
  #'
  #' @return Returns the table from the HTML as a tibble with column names 
  #' in snake-case.
  #' @note This function is a less specific version of 
  #' parse_wgs_html_table_by_div_id.
  #' It can be used for cases where the table does not have a div_id available.
  #' This is specifically relevant to the patient details table - see example.
  #' 
  #' @examples patient_ids <- parse_wgs_html_table_by_number(filepath, 1)
  
  html <- rvest::read_html(x = html_filepath)
  
  html_tables <- html |> 
    rvest::html_elements(".table") 
  
  output_table <- html_tables[[ {{ table_number }}]] |> 
    rvest::html_table() |> 
    janitor::clean_names()
  
  if(nrow(output_table) == 0) {
    warning("No rows in output table")
  }
  
  return(output_table)
  
}


wgs_html_variant_type_regex <- function() {
  
  #' Regular expression for whole genome sequencing "variant type" HTML column
  #'
  #' @return The regular expression for the "variant type" column 
  #' (example format: LOH(2); GAIN(3)).
  #' This output is then used for later functions for extracting different 
  #' parts of the regex.
  #' @export
  
  wgs_html_variant_type_regex <- stringr::regex(
    
    r"[
    (GAIN|LOSS|LOH|INV|DUP|DEL|  # Standard variant types
    .+)                          # Catch-all category
    (\(\d{1,3}\)                 # Dosage number
    |)                           # Value if absent
    ]",
    comments = TRUE
    
  )
  
  return(wgs_html_variant_type_regex)
  
}

parse_wgs_cnv_class <- function(x) {
  
  #' Parse the CNV class from the WGS variant type string
  #'
  #' @param x The WGS variant type string from a WGS HTML which has the 
  #' variant type and then may also have the copy number in brackets. 
  #' Examples: "LOH(2)", "GAIN(3)", "INV"
  #'
  #' @return The CNV class (examples: "LOH", "GAIN", "INV") as a string
  #' @export
  #'
  #' @examples cnv_class <- parse_wgs_cnv_class("LOH (2)")
  
  wgs_cnv_class <- stringr::str_extract(string = x,
                               pattern = wgs_html_variant_type_regex(),
                               group = 1)
  
  return(wgs_cnv_class)
  
}


parse_wgs_cnv_copy_number <- function(x) {
  
  #' Parse the CNV copy number from the WGS variant type string
  #'
  #' @param x The WGS variant type string from a WGS HTML which has the 
  #' variant type and then may also have the copy number in brackets. 
  #' Examples: "LOH(2)", "GAIN(3)", "INV"
  #'
  #' @return The CNV copy number as a number
  #' @export
  #'
  #' @examples cnv_number <- parse_wgs_cnv_copy_number("LOH (2)")
  
  wgs_cnv_copy_number <- readr::parse_number(stringr::str_extract(
                                                  string = x,
                                                  pattern = wgs_html_variant_type_regex(),
                                                  group = 2))
  
  return(wgs_cnv_copy_number)
  
}

wgs_html_grch38_coordinates_regex <- function() {
  
  #' Regular expression for the whole genome sequencing variant coordinates string
  #'
  #' @return The regular expression for the "Variant GRCh38 coordinates" 
  #' column in the WGS HTML file.
  #' The format gives the chromosome, start coordinate and end coordinate for a CNV.
  #' Example: "7:60912080-159335569"
  #' @export
  
  wgs_html_grch38_coordinates_regex <- stringr::regex(
    
      r"[
      (\d{1,2}|X|Y)        # Chromosome
      :
      (\d{1,10})           # First coordinate
      (-|:)
      (\d{1,10})           # Second coordinate
      ]",
      comments = TRUE
      
    )
  
  return(wgs_html_grch38_coordinates_regex)
  
}

parse_wgs_html_grch38_coordinates <- function(x, group) {
  
  #' Parse a group from the whole genome sequencing GRCh38 coordinate string
  #'
  #' @param x The coordinate string. Example: "7:60912080-159335569"
  #' @param group The group to be selected, which must be 1 for "chromosome", 
  #' 2 for "first coordinate" or 4 for "second coordinate".
  #'
  #' @return The selected regex group as a string.
  #' @export
  #'
  #' @examples cnv_chromosome <- parse_wgs_html_grch38_coordinates(x = "7:60912080-159335569",
  #' group = 1)
  
  output <- stringr::str_extract(string = x,
                        pattern = wgs_html_grch38_coordinates_regex(),
                        group = group)
  
  return(output)
  
}
