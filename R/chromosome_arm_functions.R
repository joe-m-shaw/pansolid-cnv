# Chromosome arm functions

define_autosomes <- function() {
  
  autosomes <- c("1", "2", "3", "4", "5", "6", "7", "8",
                 "9", "10", "11", "12", "13", "14", "15", 
                 "16", "17", "18", "19", "20", "21", "22")
  
  return(autosomes)
  
}

define_sex_chromosomes <- function() {
  
  sex_chromosomes <- c("X", "Y")
  
  return(sex_chromosomes)
  
}

define_chromosome_levels <- function() {
  
  #' Define factor levels for human chromosomes
  #'
  #' This function can be used to define the factor levels of a chromosome
  #' variable so that the observations will sort in chromosome order.  
  #'
  #' @returns A vector of human autosomes and sex chromosomes in numeric order 
  #' and then X and Y.
  #' 
  #' @export
    
  chr_levels <- c(define_autosomes(), define_sex_chromosomes())
  
  return(chr_levels)
  
}

format_chromosome_decimals <- function(df, col = chromosome) {
  
  #' Format chromosome character strings when formatted as decimal numbers
  #' 
  #' PanSolid Excel outputs contain columns with chromosome information.
  #' When these columns are read as text, the autosomes are converted into 
  #' strings with decimal places - example: "1" becomes "1.0".
  #' `format_chromosome_decimals` removes the decimal place to make these
  #' strings easier to work with
  #'
  #' @param df A dataframe  
  #'
  #' @param col The column of the dataframe containing 
  #' chromosome information formatted as a decimal number (i.e. "1.0").
  #'
  #' @returns The dataframe with an additional "chromosome_char" column with
  #' the decimal places removed.
  #' @export
  
  output <-  df |> 
    dplyr::mutate(chromosome_char = stringr::str_replace(string = {{ col }},
                                      pattern = "\\.0",
                                      replacement = ""))
  
  return(output)
  
}

factorise_chromosome <- function(df, col) {
  
  #' Factorise the chromosome column in a dataframe
  #'
  #' This function is designed to make sorting by chromosome easier.
  #' If chromosome information is formatted as a character then sorting does
  #' not work with autosomes (chromosome "10" is sorted after chromosome "1", not 
  #' after chromosome "9").
  #' If chromosome information is formatted as a numeric then sex chromosomes
  #' cannot be correctly formatted.
  #' Formatting as a factor allows easier by specifying the order of autosomes
  #' first, then sex chromosomes.
  #'
  #' @param df A dataframe 
  #' 
  #' @param col The column of the dataframe with chromosome information of 
  #' type character.
  #'
  #' @returns The input dataframe with a new column called "chromosome_fct" 
  #' containing the factorised version of the chromosome.
  #' 
  #' @export
  #'
  
  output <- df |> 
    dplyr::mutate(chromosome_fct = factor({{ col }}, 
                              levels = define_chromosome_levels()))
  
  return(output)
  
}

add_chromosome_arms <- function(df, chrom_col = chromosome_char, 
                                start_col = start, end_col = end) {
  
  #' Add chromosome arm information to GRCh38 region coordinates
  #' 
  #' Chromosome arm coordinates are taken from this github repo: 
  #' https://github.com/Neurosurgery-Brain-Tumor-Center-DiazLab/CONICS/blob/master/chromosome_arm_positions_grch38.txt
  #' 
  #' @param df A dataframe containing region information.
  #'
  #' @returns The input dataframe with an additional "chromosome_arm" column for
  #' autosome data only. Sex chromosome coordinates are removed.
  #' @export
  #'

  output <- df |> 
    dplyr::filter(!{{ chrom_col }} %in% c("X", "Y")) |> 
    dplyr::mutate(chromosome_arm = dplyr::case_when(
      # chr1
      {{ chrom_col }} == "1" & end < 122026459 ~"1p",
      {{ chrom_col }} == "1" & start >= 124932724 ~"1q",
      # chr2
      {{ chrom_col }} == "2" & end < 92188145 ~"2p",
      {{ chrom_col }} == "2" & start >= 94090557 ~"2q",
      # chr3
      {{ chrom_col }} == "3" & end < 90772458 ~"3p",
      {{ chrom_col }} == "3" & start >= 93655574 ~"3q",
      # chr4
      {{ chrom_col }} == "4" & end < 49712061 ~"4p",
      {{ chrom_col }} == "4" & start >= 51743951 ~"4q",
      # chr5
      {{ chrom_col }} == "5" & end < 46485900 ~"5p",
      {{ chrom_col }} == "5" & start >= 50059807 ~"5q",
      # chr6
      {{ chrom_col }} == "6" & end < 58553888 ~"6p",
      {{ chrom_col }} == "6" & start >= 59829934 ~"6q",
      # chr7
      {{ chrom_col }} == "7" & end < 58169653 ~"7p",
      # Coordinate change to reflect PanSolid coordinate
      {{ chrom_col }} == "7" & start >= 61091465 ~"7q",
      # chr8
      {{ chrom_col }} == "8" & end < 44033744 ~"8p",
      {{ chrom_col }} == "8" & start >= 45877265 ~"8q",
      # chr9
      {{ chrom_col }} == "9" & end < 43389635 ~"9p",
      {{ chrom_col }} == "9" & start >= 45518558 ~"9q",
      # chr10
      {{ chrom_col }} == "10" & end < 39686682 ~"10p",
      {{ chrom_col }} == "10" & start >= 41593521 ~"10q",
      # chr11
      {{ chrom_col }} == "11" & end < 51078348 ~"11p",
      {{ chrom_col }} == "11" & start >= 54425074 ~"11q",
      # chr12
      {{ chrom_col }} == "12" & end < 34769407 ~"12p",
      {{ chrom_col }} == "12" & start >= 37185252 ~"12q",
      # chr 13
      {{ chrom_col }} == "13" & start >= 18051248 ~"13q",
      # chr14
      {{ chrom_col }} == "14" & start >= 18173523 ~"14q",
      # chr15
      {{ chrom_col }} == "15" & start >= 19725254 ~"15q",
      # chr16
      {{ chrom_col }} == "16" & end < 36311158 ~"16p",
      {{ chrom_col }} == "16" & start >= 36334460 ~"16q",
      # chr17
      {{ chrom_col }} == "17" & end < 22813679 ~"17p",
      {{ chrom_col }} == "17" & start >= 26566633 ~"17q",
      # chr18
      {{ chrom_col }} == "18" & end < 15460899 ~"18p",
      {{ chrom_col }} == "18" & start >= 20861206 ~"18q",
      # chr19
      {{ chrom_col }} == "19" & end < 24498980 ~"19p",
      {{ chrom_col }} == "19" & start >= 27190874 ~"19q",
      # chr20
      {{ chrom_col }} == "20" & end < 26436232 ~"20p",
      {{ chrom_col }} == "20" & start >= 30038348 ~"20q",
      # chr21
      {{ chrom_col }} == "21" & end < 10864560 ~"21p",
      {{ chrom_col }} == "21" & start >= 12915808 ~"21q",
      # chr22
      {{ chrom_col }} == "22" & end < 12954788 ~"22p",
      {{ chrom_col }} == "22" & start >= 15054318 ~"22q"
    ))
  
  if(anyNA.data.frame(output) == TRUE) {
    message("There are NAs in the output")
  }
  
  return(output)

}

calculate_region_length <- function(df,  
                                  start_col = start,
                                  end_col = end) {
  
  output <- df |> 
    dplyr::mutate(region_length = {{ end_col }} - {{ start_col }}) 
  
  return(output)
  
}

add_chr_arm_region_percent <- function(df, 
                                       chr_arm_col = chromosome_arm,
                                       start_col = start,
                                       end_col = end) {
  
  # This dataframe is created by the pansolid_targets_per_chromosome_arm.R script
  chromosome_arm_positions_grch38 <- readr::read_delim(file = paste0(config::get("data_folderpath"),
                                                                     "validation/DOC6791_chromosome_arms/",
                                                                     "bed_files/",
                                                                     "chromosome_arm_positions_grch38.txt")) |> 
    janitor::clean_names()
  
  if(!"region_length" %in% colnames(df)) {
    
    df <- df |> 
      calculate_region_length(start_col = start_col,
                              end_col = end_col)
  }
  
  output <- df |>  
    dplyr::left_join(chromosome_arm_positions_grch38 |> 
                       dplyr::select(chromosome_arm, arm_length),
                     dplyr::join_by( {{ chr_arm_col }} == "chromosome_arm"),
                     relationship = "many-to-one") |> 
    dplyr::mutate(percent_chr_arm = round((region_length / arm_length) * 100, 1)) 
  
  return(output)
  
}

add_chr_arm_ps_target_region_percent <- function(df,
                                                 start_col = start,
                                                 end_col = end,
                                                 chr_arm_col = chromosome_arm) {
  
  pansolid_target_chr_arm_bed <- readr::read_csv(paste0(config::get("data_folderpath"),
                                                        "validation/DOC6791_chromosome_arms/",
                                                        "bed_files/",
                                                        "pansolid_target_region_chr_arm_bed.csv"))
  
  if(!"region_length" %in% colnames(df)) {
    
    df <- df |> 
      calculate_region_length(start_col = start_col,
                              end_col = end_col)
  }
  
  output <- df |> 
    dplyr::left_join(pansolid_target_chr_arm_bed |> 
                       dplyr::select(chromosome_arm, chr_arm_targets,
                                     chr_arm_target_covered_region_length),
                     dplyr::join_by({{ chr_arm_col }} == "chromosome_arm"),
                     relationship = "many-to-one") |> 
    dplyr::mutate(percent_chr_arm_target_region = 
             round((region_length / chr_arm_target_covered_region_length) * 100, 1))
  
  return(output)

}

add_cumulative_chr_coordinates <- function(df, col) {
  
  pansolid_chr_cumulative_coordinates <- readr::read_csv(paste0(config::get("data_folderpath"),
                                                                "validation/DOC6791_chromosome_arms/",
                                                                "bed_files/",
                                                                "pansolid_chr_cumulative_coordinates.csv"),
                                                         col_types = "ccdd")
  
  output <- df |> 
    dplyr::left_join(pansolid_chr_cumulative_coordinates |> 
              dplyr::select(chromosome_fct, cumulative_chr_coordinate),
              dplyr::join_by({{ col }} == "chromosome_fct")) 
  
  return(output)
  
}

read_snp_sheet <- function(filepath, 
                           sheetname = "Artefacts_removed_GIAB_filt...") {
  
  #' Read the SNP frequency sheet from PanSolid Excels
  #' 
  #' This function reads the "Artefacts_removed_GIAB_filt..." sheet from the
  #' Excel file which is saved in the "Unannotated_Files" folder of each
  #' PanSolid worksheet. This sheet includes data for the single nucleotide
  #' polymorphisms (SNPs) detected in the sample, which are presented in 
  #' the variant allele frequency (VAF) track of the interactive HTML file.
  #'
  #' @param filepath The filepath of the Excel to read from.
  #' @param sheetname The name of the sheet to read from.
  #'
  #' @returns A dataframe of the first 12 columns of the Excel sheet, which is
  #' annotated with the sample details from the filename.
  #' @export
  #'
  
  df <- readxl::read_excel(path = filepath,
                           sheet = sheetname,
                           range = readxl::cell_cols("A:L"),
                           col_types = c(
                             "text", "text", "text", "text", "text", 
                             "text", "numeric", "text", "text",
                             "numeric", "numeric", "numeric"
                           )) |> 
    janitor::clean_names()
  
  output <- add_identifiers(filepath, df)
  
  return(output)
  
}

format_snp_sheet_data <- function(df) {
  
  output <- df |> 
    format_chromosome_decimals(col = chromosome) |> 
    factorise_chromosome(col = chromosome_char) |> 
    dplyr::filter(type == "SNV") |> 
    dplyr::mutate(region_numeric = as.numeric(region)) |> 
    add_cumulative_chr_coordinates(col = chromosome_fct) |> 
    mutate(cumulative_region_coordinate = region_numeric + cumulative_chr_coordinate)
  
  return(output)
  
}

read_targets_merged <- function(filepath, sheetname = "CNV Targets Merged"){
  
  output <- readxl::read_excel(path = filepath,
                       sheet = sheetname,
                       range = "A1:N6175",
                       col_types = c("text", "text", "text", "text", "text",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric"))
  
  x <- add_identifiers(filepath, output) |> 
    janitor::clean_names()
  
  return(x)
  
}
