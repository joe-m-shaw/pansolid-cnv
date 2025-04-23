source(here::here("functions/extract_pansolid_cnv_coordinates.R"))
source(here::here("functions/filename_functions.R"))

get_sheetname <- function(filepath, sheet_regex = "CNVs_") {
  
  #' Get sheet names from an Excel file using regular expressions
  #'
  #' Excel files exported from the CLC pipeline have each sheet named with
  #' the 8 digit DNA laboratory number as a suffix.
  #' Example: "Amplifications_12345678"
  #' `get_sheetname` can be used to find the specific sheet name for a file.
  #'  
  #' @param filepath Full file path to an Excel file
  #' @param sheet_regex Regular expression for sheet name matching
  #'
  #' @returns The name of the sheet as a string
  #' @export
  #'
  #' @examples 
  #' 
  #' subfolder <- "validation/DOC6567_deletions/test_data/"
  #' 
  #' path <- paste0(config::get("data_folderpath"), subfolder)
  #' 
  #' files <- list.files(path, pattern = "*.xlsx", full.names = TRUE)
  #' 
  #' get_sheetname(files[2])

  sheets <- readxl::excel_sheets(filepath)
  
  sheetname <- grep(pattern = {{ sheet_regex }}, 
                    x = sheets, 
                    value = TRUE)
  
  if(length(sheetname) == 0) {
    stop("Sheet name not found")
  }
  
  if(length(sheetname) > 1) {
    stop("Sheet name must be unique")
  }
  
  return(sheetname)
  
}

read_cnv_sheet <- function(filepath, sheet_regex = "CNVs_") {
  
  #' Read the entire CNV sheet from a PanSolid CLC Excel file
  #' 
  #' The new PanSolid CLC Excel output includes multiple tables on the
  #' same tab. `read_cnv_sheet` reads the entire sheet as one table, which 
  #' allows individual tables to be subsequently extracted 
  #' using `extract_cnv_tbls` without the file being read multiple times. 
  #'
  #' @param filepath Full file path to an Excel file
  #' @param sheet_regex Regular expression for sheet name matching
  #'
  #' @returns The full CNV sheet as a tibble
  #' @export
  #'
  #' @examples
  
  sheet <- readxl::read_xlsx(path = filepath, 
                   sheet = get_sheetname(filepath = filepath, 
                                         sheet_regex = sheet_regex),
                   range = "A1:K100",
                   col_names = c("a", "b", "c", "d", "e", "f", "g", "h", 
                                 "i", "j", "k"),
                   col_types = c("text", "text", "text", "text", "text", 
                                 "text", "text", "text", "text", "text",
                                 "text"))
  
  return(sheet)
  
}

find_match <- function(input_sheet, input_col, match_string) {
  
  #' Find matches for strings within an Excel sheet
  #'
  #' @param input_sheet The Excel sheet to search
  #' @param input_col The column to search
  #' @param match_string The string to search for
  #'
  #' @returns The numeric position of the string match within the column
  #' @export
  #'
  #' @examples find_match(input_sheet, "a", "Amplification genes")
  
  output <- match(match_string, input_sheet[[ input_col ]])
  
  if(is.na(output)){
    stop(paste0(match_string, " not found"))
  }
  
  return(output)
  
}

find_stdev_ratios <- function(input_sheet, 
                              stdev_string = "StDev Signal-adjusted Log2 Ratios") {
  
  #' Find standard deviation of signal-adjusted log2 ratios within the CNV sheet
  #'
  #' @param input_sheet A tibble containing the stdev information. This is 
  #' intended to be the output from `read_cnv_sheet`.
  #' 
  #' @param stdev_string The string which identifies the position of the 
  #' standard deviation information. This defaults to the current header, but
  #' can be changed if necessary.
  #'
  #' @returns A tibble with the standard deviation of signal-adjusted 
  #' log2 ratios  
  #' @export
  #'
  #' @examples
  
  stdev_start <- find_match(input_sheet, "a", stdev_string) + 1
  
  stdev_df <- tibble::as_tibble(input_sheet[stdev_start, 1]) |> 
    dplyr::rename(stdev_noise = a) |> 
    dplyr::mutate(stdev_noise = as.double(stdev_noise))
  
  if(stdev_df[1,1] < 0) {
    stop("Standard deviation cannot be below 0")
  }
  
  return(stdev_df)
  
}


find_percent_138x <- function(input_sheet, 
                              percent_138x_string = "% Whole Panel Covered at 138X") {
  
  #' Find percentage coverage at 138X within the CNV sheet
  #'
  #' @param input_sheet A tibble containing the percent 138X information. 
  #' This is intended to be the output from `read_cnv_sheet`.
  #'
  #' @param percent_138x_string The string which identifies the position of the 
  #' percent coverage at 138x table. This defaults to the current header, but
  #' can be changed if necessary.
  #'
  #' @returns A tibble with the percent 138X information.
  #' @export
  #'
  #' @examples
  
  percent_138x_start <- find_match(input_sheet, "a", percent_138x_string) + 1
  
  percent_138x_df <- tibble::as_tibble(input_sheet[percent_138x_start, 1]) |> 
    dplyr::rename(percent_138x = a) |> 
    dplyr::mutate(percent_138x = as.double(percent_138x))
  
  if(percent_138x_df[1,1] < 0 | percent_138x_df[1,1] > 100) {
    stop("Percent 138X must be between 0 and 100")
  }
  
  return(percent_138x_df)
  
}

find_pred_ncc <- function(input_sheet, 
                          pred_ncc_string = "Predicted NCC (%)") {
  
  #' Find predicted neoplastic cell content (NCC) within the CNV sheet
  #'
  #' @param input_sheet A tibble containing the percent NCC information. 
  #' This is intended to be the output from `read_cnv_sheet`.
  #' @param pred_ncc_string The string which identifies the position of the 
  #' predicted NCC value. This defaults to the current header, but
  #' can be changed if necessary.
  #'
  #' @returns A tibble with the predicted NCC information.
  #' @export
  #'
  #' @examples
  
  pred_ncc_start <- find_match(input_sheet, "a", pred_ncc_string) + 1
  
  pred_ncc_df <- tibble::as_tibble(input_sheet[pred_ncc_start, 1]) |> 
    dplyr::rename(pred_ncc = a) |> 
    dplyr::mutate(pred_ncc = as.double(pred_ncc))
  
  if(pred_ncc_df[1,1] < 0 | pred_ncc_df[1,1] > 100) {
    stop("Predicted NCC must be between 0 and 100")
  }
  
  return(pred_ncc_df)
  
}

find_amp_genes <- function(input_sheet, num_genes = 9, 
                           amp_string = "Amplification genes") {
  
  #' Find fold change information for amplification genes with the CNV sheet
  #'
  #' @param input_sheet A tibble containing the amplified gene information. 
  #' This is intended to be the output from `read_cnv_sheet`.
  #' 
  #' @param num_genes The number of genes with amplification results returned
  #' by the CLC pipeline. This is currently 9 but may change in future.
  #'
  #' @param amp_string The string which identifies the position of the 
  #' amplification genes table. This defaults to the current header, but
  #' can be changed if necessary.
  #'
  #' @returns A tibble of the amplification genes results table
  #' @export
  #'
  #' @examples
  
  amp_tbl_header <- find_match(input_sheet, "a", amp_string)
  
  amp_tbl_colname_row <- amp_tbl_header + 1
  
  amp_tbl_start <- amp_tbl_header + 2
  
  amp_tbl_end <- amp_tbl_start + (num_genes-1)
  
  df <- tibble::as_tibble(input_sheet[amp_tbl_start:amp_tbl_end, 1:3])
  
  colnames(df) <- as.character(input_sheet[amp_tbl_colname_row, 1:3])          
  
  df_clean <- df |> 
    janitor::clean_names() |> 
    dplyr::mutate(max_region_fold_change = as.double(max_region_fold_change),
           min_region_fold_change = as.double(min_region_fold_change))
  
  return(df_clean)
  
}


find_del_genes <- function(input_sheet, num_genes = 37,
                           del_string = "Deletion genes",
                           col_length = 15) {
  
  #' For fold change information for deletion genes within the CNV sheet
  #'
  #' There are currently 37 genes identified for deletion analysis within the
  #' CLC pipeline. These are presented in 3 tables. `find_del_genes` identifies
  #' each table within the CNV sheet and collates them together.
  #'
  #' @param input_sheet A tibble containing the deleted gene information. 
  #' This is intended to be the output from `read_cnv_sheet`.
  #'
  #' @param num_genes The number of genes within the deletion genes table.
  #'
  #' @param del_string The string which identifies the position of the 
  #' deletion genes table. This defaults to the current header, but
  #' can be changed if necessary.
  #' 
  #' @param col_length The length in rows of the columns that genes are 
  #' presented in. Defaults to 15 and then blank space is removed 
  #' automatically. 
  #'
  #' @returns A tibble of the deletion genes results table
  #' @export
  #'
  #' @examples
  
  del_tbl_header <- find_match(input_sheet, "a", del_string)
  
  del_tbl_colname_row <- del_tbl_header + 1
  
  del_tbl_start <- del_tbl_header + 2
  
  del_tbl_end <- del_tbl_start + col_length
  
  del_tbl1 <- tibble::as_tibble(input_sheet[del_tbl_start:del_tbl_end, 1:3])
  
  colnames(del_tbl1) <- as.character(input_sheet[del_tbl_colname_row, 1:3])
  
  del_tbl2 <- tibble::as_tibble(input_sheet[del_tbl_start:del_tbl_end, 5:7])
  
  colnames(del_tbl2) <- as.character(input_sheet[del_tbl_colname_row, 5:7])
  
  del_tbl3 <- tibble::as_tibble(input_sheet[del_tbl_start:del_tbl_end, 9:11])
  
  colnames(del_tbl3) <- as.character(input_sheet[del_tbl_colname_row, 9:11])
  
  del_tbl <- rbind(del_tbl1, del_tbl2, del_tbl3) |> 
    janitor::clean_names() |> 
    dplyr::mutate(max_region_fold_change = as.double(max_region_fold_change),
           min_region_fold_change = as.double(min_region_fold_change)) |> 
    dplyr::filter(!is.na(gene))
  
  if(nrow(del_tbl) != num_genes){
    warning("Different number of genes found to expected")
  }
  
  return(del_tbl)
  
}

find_sig_cnvs <- function(input_sheet, 
                          sig_cnv_string = "Significant CNV results",
                          amp_gene_string = "Amplification genes") {
  
  #' Find information for significant CNV results within the CNV sheet
  #'
  #' @param input_sheet A tibble containing the significant CNV information. 
  #' This is intended to be the output from `read_cnv_sheet`.
  #'
  #' @param sig_cnv_string The string which identifies the position of the 
  #' significant CNVs table. This defaults to the current header, but
  #' can be changed if necessary.
  #' 
  #' @param amp_gene_string The string which identifiers the position
  #' of the amplification genes table.
  #'
  #' @returns A tibble of the significant CNV results table
  #' @export
  #'
  #' @examples
  
  sig_cnv_header <- find_match(input_sheet, "a", sig_cnv_string)
  
  sig_cnv_colname_row <- sig_cnv_header + 1
  
  amp_gene_header <- find_match(input_sheet, "a", amp_gene_string)
    
  na_row <- amp_gene_header - 1
  
  stopifnot(is.na(input_sheet[na_row, 1]))
  
  sig_cnv_tbl_start <- sig_cnv_colname_row + 1
  
  sig_cnv_tbl_end <- na_row - 1
  
  # Scenario where no significant CNVs are exported by the pipeline
  if(sig_cnv_tbl_start > sig_cnv_tbl_end) {
    
    sig_cnv_df <- tibble::tibble(
      "gene" = "no positive calls",
      "chromosome" = "",
      "cnv_co_ordinates" = "",
      "cnv_length" = 0,
      "consequence" = "no call",
      "fold_change" = 0,
      "p_value" = 0,
      "no_targets" = 0,
      "check_1" = "",
      "check_2" = "",
      "copy_number" = 0)
    
  }
  
  # Scenario where significant CNVs are exported by the pipeline
  if(sig_cnv_tbl_start <= sig_cnv_tbl_end) {
    
    sig_cnv_df <- tibble::as_tibble(input_sheet[sig_cnv_tbl_start:sig_cnv_tbl_end, 1:11])
    
    colnames(sig_cnv_df) <- as.character(input_sheet[sig_cnv_colname_row, 1:11])
    
    sig_cnv_df <- sig_cnv_df |> 
      janitor::clean_names() |> 
      dplyr::mutate(cnv_length = as.double(cnv_length),
                    fold_change = as.double(fold_change),
                    p_value = as.double(p_value),
                    no_targets = as.double(no_targets),
                    copy_number = as.double(copy_number))
    
  }
  
  sig_cnv_df_clean <- extract_pansolid_cnv_coordinates(sig_cnv_df,
                                                       cnv_co_ordinates)
  
  return(sig_cnv_df_clean)
  
}

extract_cnv_tbls <- function(filepath, 
                             sheet_regex = "CNVs_",
                             stdev_string = "StDev Signal-adjusted Log2 Ratios",
                             percent_138x_string = "% Whole Panel Covered at 138X",
                             pred_ncc_string = "Predicted NCC (%)",
                             num_amp_genes = 13, 
                             amp_string = "Amplification genes",
                             num_del_genes = 34,
                             del_string = "Deletion genes",
                             del_col_length = 15,
                             sig_cnv_string = "Significant CNV results") {
  
  #' Extract all information from the CNV sheet of a PanSolid CLC Excel output
  #'
  #' This function acts as a wrapper of the `find_*` functions to extract
  #' information from the 5 tables on the CLC pipeline Excel output CNV
  #' sheet, adds sample identifiers and stores them within a list. 
  #' 
  #' The defaults are set for the version of the CNV tab used in the live
  #' PanSolid service following addition of deletions and LOH.
  #'
  #' @param filepath Full file path to an Excel file
  #' @param sheet_regex Regular expression for sheet name matching
  #'
  #' @returns A named list of tables 
  #' @export
  #'
  #' @examples
  
  sheet <- read_cnv_sheet(filepath = filepath, 
                      sheet_regex = sheet_regex)
  
  tables <- list(
    
    "stdev" = add_identifiers(
      file = filepath,
      tbl = find_stdev_ratios(sheet,
                              stdev_string = stdev_string)),
    
    "percent_138x" = add_identifiers(
      file = filepath,
      tbl = find_percent_138x(sheet,
                              percent_138x_string = percent_138x_string)),
    
    "pred_ncc" = add_identifiers(
      file = filepath,
      tbl = find_pred_ncc(sheet,
                          pred_ncc_string = pred_ncc_string)),
    
    "sig_cnvs"  = add_identifiers(
      file = filepath,
      tbl = find_sig_cnvs(sheet,
                          sig_cnv_string = sig_cnv_string)),
    
    "amp_genes" = add_identifiers(
      file = filepath,
      tbl = find_amp_genes(sheet,
                           num_genes = num_amp_genes, 
                           amp_string = amp_string)),
    
    "del_genes" = add_identifiers(
      file = filepath,
      tbl = find_del_genes(sheet,
                           num_genes = num_del_genes,
                           del_string = del_string,
                           col_length = del_col_length))
  )
  
  return(tables)
  
}

read_loh_table <- function(filepath, 
                           sheet_regex = "LOH_",
                           sheet_range = "A1:G8",
                           loh_genes = c("MSH2", "MSH6", "MLH1",
                                         "PMS2", "LZTR1", "SMARCB1", "NF2")) {
  
  #' Read the loss of heterozygosity table from PanSolid CNV Excels
  #'
  #' @param filepath The full filepath to the Excel file
  #' @param sheet_regex Regular expression for sheet name matching, defaults to
  #' "LOH_".
  #' @param loh_genes Vector of genes expected to have LOH results.
  #'
  #' @returns The LOH results table as a tibble
  #' @export
  #'
  #' @examples
  
  loh_table <- readxl::read_xlsx(path = filepath, 
                                 sheet = get_sheetname(filepath,
                                                       sheet_regex = sheet_regex),
                                 range = sheet_range,
                                 col_types = c("text", "text","text", 
                                               "text", "text",
                                               "text", "text")) |> 
    janitor::clean_names()
  
  if(setequal(loh_table$gene, loh_genes) == FALSE){
    stop("LOH gene list not as expected")
  }
  
  loh_table_with_ids <- add_identifiers(filepath, loh_table)
  
  return(loh_table_with_ids)
}

source(here::here("tests/test_pansolid_cnv_excel_functions.R"))
