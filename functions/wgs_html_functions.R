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
  
  wgs_version <- stringr::str_extract(header, header_regex, group = 1)
  
  wgs_analysis_date <- stringr::str_extract(header, header_regex, group = 2)
  
  wgs_r_no <- stringr::str_extract(header, header_regex, group = 3)
  
  wgs_p_no <- stringr::str_extract(header, header_regex, group = 4)
  
  output <- data.frame( 
    "filepath" = html_filepath,
    "wgs_r_no" = wgs_r_no,
    "wgs_p_no" = wgs_p_no,
    "wgs_version" = wgs_version,
    "wgs_analysis_date" = wgs_analysis_date)
  
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
  
  name <- stringr::str_extract(string = pid_text, pattern = pid_regex,
                      group = 1)
  
  dob <- stringr::str_extract(string = pid_text, pattern = pid_regex,
                     group = 2)
  
  nhs_no <- stringr::str_extract(string = pid_text, pattern = pid_regex,
                        group = 3)
  
  nhs_no_clean <- stringr::str_replace_all(string = nhs_no,
                                  pattern = " ",
                                  replacement = "")
  
  output <- data.frame( 
    "filepath" = html_filepath,
    "patient_name" = name,
    "patient_dob" = dob,
    "nhs_no" = nhs_no,
    "nhs_no_clean" = nhs_no_clean)
  
  return(output)
  
}

parse_wgs_html_table_by_div_id <- function(html_filepath,
                                           div_id) {
  
  #' Parse whole genome sequencing HTML files using a CSS
  #' "div_id" identifier.
  #'
  #' @param html_filepath The filepath for the WGS HTML file
  #' @param div_id The CSS identifier for the table (examples below)
  #' 
  #' @return Returns the table from the HTML as a tibble with column names in snake-case.
  #' 
  #' @section Useful Identifiers: 
  #' "t_tumour_details", "t_tumour_sample", 
  #' "t_germline_sample", "t_quality", "tier1" (Tier 1 sequence variants),
  #' "svcnv_tier1" (Tier 1 structural and copy number variants, for v2.28 and earlier), 
  #' "d_svcnv_tier1" (Tier 1 structural and copy number variants,later versions than v2.28)
  #'
  #' @examples parse_wgs_html_table_by_div_id(filepath, "t_tumour_details")
  
  html <- rvest::read_html(x = html_filepath)
  
  output_table <- html |> 
    rvest::html_element( str_c("#", {{ div_id }} )) |> 
    rvest::html_table() |> 
    janitor::clean_names() |> 
    dplyr::mutate(filepath = html_filepath)
  
  return(output_table)
  
}

parse_wgs_html_table_by_number <- function(html_filepath,
                                           table_number) {
  
  #' Parse whole genome sequencing HTML files using the table number
  #'
  #' @param html_filepath The filepath for the WGS HTML file
  #' @param table_number The number of the table to select
  #'
  #' @return Returns the table from the HTML as a tibble with column names in snake-case.
  #' @note This function is a less specific version of parse_wgs_html_table_by_div_id.
  #' It can be used for cases where the table does not have a div_id available. This is 
  #' specifically relevant to the patient details table - see example.
  #' 
  #' @examples patient_ids <- parse_wgs_html_table_by_number(filepath, 1)
  
  html <- rvest::read_html(x = html_filepath)
  
  html_tables <- html |> 
    rvest::html_elements(".table") 
  
  output_table <- html_tables[[ {{ table_number }}]] |> 
    rvest::html_table() |> 
    janitor::clean_names()
  
  return(output_table)
  
}

wgs_html_variant_type_regex <- stringr::regex(
  r"[
  (GAIN|LOSS|LOH|INV|DUP|DEL|  # Standard variant types
  .+)                          # Catch-all category
  (\(\d{1,3}\)                 # Dosage number
  |)                           # Value if absent
  ]",
  comments = TRUE
)

parse_wgs_cnv_class <- function(col) {
  
  wgs_cnv_class <- stringr::str_extract(string = col,
                               pattern = wgs_html_variant_type_regex,
                               group = 1)
  
  if(anyNA(wgs_cnv_class, recursive = TRUE)) {
    stop("NA values present")
  }
  
  return(wgs_cnv_class)
  
}

parse_wgs_cnv_copy_number <- function(col) {
  
  wgs_cnv_copy_number <- readr::parse_number(stringr::str_extract(string = {{ col }} ,
                                                  pattern = wgs_html_variant_type_regex,
                                                  group = 2))
  
  return(wgs_cnv_copy_number)
  
}

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

parse_wgs_html_grch38_coordinates <- function(col, group) {
  
  if(!group %in% c("chromosome", "first coordinate", "second coordinate")) {
    stop("group must be either chromosome, first coordinate or second coordinate")
  }
  
  group_choice <- dplyr::case_when(
    
    group == "chromosome" ~1,
    group == "first coordinate" ~2,
    group == "second coordinate" ~4)
  
  output <- stringr::str_extract(string  = col,
                        pattern = wgs_html_grch38_coordinates_regex,
                        group = group_choice)
  
  return(output)
  
}