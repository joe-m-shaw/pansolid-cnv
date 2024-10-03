source(here("functions/extract_pansolid_cnv_coordinates.R"))

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
    (WS\d{6})                                   # Worksheet number
    _
    (\d{5,8})                                   # Lab number
    (a|b|c|d|)                                  # Suffix
    _
    ([:alnum:]{3,30})                           # Patient name - alphanumeric characters only
    (.xlsx|_S.+.xlsx|_S.+|_CNV_processed.xlsx)  # Ending varies between patients and controls
    ]",
    comments = TRUE)
  
  output <- stringr::str_extract(input_file, filename_regex,
                        group = input_group)
  
  return(output)
  
}

filename_to_df <- function(file) {
  
  #' Convert PanSolid results Excel filename identifiers to a dataframe
  #'
  #' @param file The filename
  #'
  #' @return A dataframe
  #'
  #' @examples filename <- "Annotated_WS123456_12345678a_JoeShaw.xlsx"
  #' 
  #' id_df <- filename_to_df(filename)
  
  output <- data.frame(
    worksheet = c(parse_filename(file, 1)),
    labno = c(as.character(parse_filename(file, 2))),
    suffix = c(parse_filename(file, 3)),
    patient_name = c(parse_filename(file, 4))) |> 
    mutate(
      labno_suffix = stringr::str_c(labno, suffix),
      labno_suffix_worksheet = stringr::str_c(labno_suffix, "_", worksheet))
  
  return(output)
  
}

get_full_amp_sheet <- function(file, sheet = "Amplifications") {
  
  #' Read the full sheet of the "Amplifications" tab of a PanSolid results Excel
  #'
  #' @param file The full filepath of the PanSolid results file
  #' @param sheet The name of the sheet to read
  #'
  #' @return All the information in the "Amplifications" tab as a data-frame.
  #'
  #' @examples x <- here::here("data/example_data/Annotated_WS123456_12345678a_JoeShaw.xlsx")
  #' 
  #' full_amp_sheet <- get_full_amp_sheet(file = x, sheet = get_amp_sheetname(x))
  
  full_tbl <- readxl::read_excel(path = file,
                         sheet = {{ sheet }},
                         col_names = FALSE) |> 
    janitor::clean_names()
  
  return(full_tbl)
  
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
  #' data <- data.frame(gene = c("ERBB2"), dq = c(2))
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

read_pos_cnv_results <- function(file, sheet = "Amplifications") {
  
  #' Read the "Positive CNV results" table from a PanSolid results Excel
  #'
  #' @param file The full filename of the Excel to read
  #' @param sheet The sheet to read the table from
  #'
  #' @return The positive CNV table as a dataframe with filename identifiers
  #'
  #' @examples x <- here::here("data/example_data/Annotated_WS123456_12345678a_JoeShaw.xlsx")
  #' 
  #' pos_cnv <- read_pos_cnv_results(file = x, sheet = get_amp_sheetname(x))
  
  full_tbl <- get_full_amp_sheet(file, sheet)
  
  pos_cnv_tbl_row <- match("Positive CNV results", full_tbl$x1)
  
  na_vector <- which(is.na(full_tbl$x1))
  
  first_na_after_pos_cnv_tbl <- min(na_vector[na_vector > pos_cnv_tbl_row])
  
  size_pos_cnv_tbl <- (first_na_after_pos_cnv_tbl - pos_cnv_tbl_row) - 1
  
  excel_range <- stringr::str_c("A", pos_cnv_tbl_row + 1, ":H", pos_cnv_tbl_row + size_pos_cnv_tbl)
  
  pos_cnv_tbl <- readxl::read_excel(path = file,
                            sheet = {{ sheet }},
                            range = excel_range,
                            col_types = c("text", "text", "text", 
                                          "numeric", "text", "numeric",
                                          "numeric","numeric")) |> 
    janitor::clean_names() 
  
  pos_cnv_coord <- extract_pansolid_cnv_coordinates(df = pos_cnv_tbl, 
                                           cnv_coord_col = cnv_co_ordinates)
  
  if (nrow(pos_cnv_coord) == 0) {
    
    pos_cnv_coord <- data.frame(
      "gene" = "no positive calls",
      "chromosome" = "",
      "cnv_co_ordinates" = "",
      "cnv_length" = 0,
      "consequence" = "no call",
      "fold_change" = 0,
      "p_value" = 0,
      "no_targets" = 0,
      "start" = 0,
      "end" = 0)
    
  }
  
  output <- add_identifiers(file, pos_cnv_coord)
  
  return(output)
  
}

read_all_amp_genes_results <- function(file, sheet = "Amplifications") {
  
  #' Read the "All amplification genes" table from a PanSolid results Excel
  #'
  #' @param file The full filename of the Excel to read
  #' @param sheet The sheet to read the table from
  #'
  #' @return The all amp gene CNV table as a dataframe with filename identifiers
  #'
  #' @examples x <- here::here("data/example_data/Annotated_WS123456_12345678a_JoeShaw.xlsx")
  #' 
  #' all_amp <- read_all_amp_genes_results(file = x, sheet = get_amp_sheetname(x))
  
  full_tbl <- get_full_amp_sheet(file, sheet)
  
  gene_table_row_start <- match("All amplification genes", full_tbl$x1)
  
  gene_tbl <- readxl::read_excel(path = file,
                         sheet = {{ sheet }},
                         skip = gene_table_row_start,
                         n_max = 9) |> 
    janitor::clean_names() 
  
  output <- add_identifiers(file, gene_tbl)
  
  return(output)
  
}

read_stdev_results <- function(file, sheet = "Amplifications") {
  
  #' Read the standard deviation of signal-adjusted noise ratios from a PanSolid Excel
  #'
  #' @param file The full filename of the Excel to read 
  #' @param sheet The sheet to read from
  #'
  #' @return A dataframe of the filename identifiers with the standard deviation value
  #'
  #' @examples x <- here::here("data/example_data/Annotated_WS123456_12345678a_JoeShaw.xlsx")
  #' 
  #' stdev <- read_stdev_results(file = x, sheet = get_amp_sheetname(x))
  
  full_tbl <- get_full_amp_sheet(file, sheet)
  
  stdev_start <- match("StDev Signal-adjusted Log2 Ratios", full_tbl$x1) - 1
  
  stdev <- readxl::read_excel(path = file,
                      sheet = {{ sheet }},
                      skip = stdev_start,
                      n_max = 1) |> 
    janitor::clean_names()
  
  output <- add_identifiers(file, stdev)
  
  return(output)
  
}

read_percent_138_results <- function(file, sheet = "Amplifications") {
  
  #' Read the percentage of the whole panel covered at 138X from a PanSolid Excel
  #'
  #' @param file The full filename of the Excel to read 
  #' @param sheet The sheet to read from
  #'
  #' @return A dataframe of the filename identifiers with the percentage of
  #' the whole panel covered at 138X value
  #'
  #' @examples  x <- here::here("data/example_data/Annotated_WS123456_12345678a_JoeShaw.xlsx")
  #' 
  #' percent_138 <- read_percent_138_results(file = x, sheet = get_amp_sheetname(x))
  
  full_tbl <- get_full_amp_sheet(file, sheet)
  
  percent_panel_start <- match("% Whole Panel Covered at 138X", full_tbl$x1) - 1
  
  percent_panel_tbl <- readxl::read_excel(path = file,
                                  sheet = {{ sheet }},
                                  skip = percent_panel_start,
                                  n_max = 1) |> 
    janitor::clean_names()
  
  output <- add_identifiers(file, percent_panel_tbl)
  
  return(output)
  
}

get_annotated_filepaths <- function(
    repository_path = "S:/central shared/Genetics/Repository/WorksheetAnalysedData/",
    worksheet, full_names = TRUE) {
  
  #' Get the filepaths of PanSolid results Excels from the S drive
  #'
  #' @param repository_path The filepath for the worksheet repository - defaults to 
  #' the standard S drive location.
  #' @param worksheet The PanSolid worksheet
  #' @param full_names TRUE or FALSE
  #'
  #' @return A list of filepaths
  #'
  #' @note PanSolid results Excels are automatically saved onto the S drive in a defined
  #' folder structure
  #'
  #' @examples get_annotated_filepaths(worksheet = "WS140721")
  
  pansolidv2_excel_regex <- "^Annotated(_|_v2.+_)WS\\d{6}_.+.xlsx"
  
  annotated_filepaths <- list.files(path = str_c(repository_path, {{ worksheet }},
                                                 "/"),
                                    recursive = TRUE, 
                                    pattern = pansolidv2_excel_regex,
                                    full.names = full_names)
  
  return(annotated_filepaths)
  
}

get_amp_sheetname <- function(filepath) {
  
  #' Get the sheet name of the amplifications tab for a PanSolid results Excel
  #'
  #' @param filepath The full filename of the Excel 
  #'
  #' @return A character string of the amplifications tab name
  #' 
  #' @note Prior to go-live, the results Excels defaulted to "Amplifications" as the 
  #'  name of the tab containing amplification information. After go live this changed
  #'  to include the sample DNA number (see example).
  #'
  #' @examples x <- here::here("data/example_data/Annotated_WS123456_12345678a_JoeShaw.xlsx")
  #' 
  #' get_amp_sheetname(x)
  
  sheets <- readxl::excel_sheets(filepath)
  
  amp_sheet_name <- grep(pattern = "Amplifications_", x = sheets, value = TRUE)
  
  return(amp_sheet_name)
  
}
