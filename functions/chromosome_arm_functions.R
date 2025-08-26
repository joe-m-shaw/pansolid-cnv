# Chromosome arm functions

chr_coords_grch38 <- readr::read_delim(file = paste0(config::get("data_folderpath"),
                                         "validation/DOC6791_chromosome_arms/",
                                         "chromosome_arm_positions_grch38.txt")) |> 
  janitor::clean_names()

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
  #' # and then X and Y.
  #' 
  #' @export
    
  chr_levels <- c(define_autosomes(), define_sex_chromosomes())
  
  return(chr_levels)
  
}


format_chromosome_decimals <- function(df) {
  
  #' Format chromosome character strings when formatted as decimal numbers
  #' 
  #' PanSolid Excel outputs contain columns with chromosome information.
  #' When these columns are read as text, the autosomes are converted into 
  #' strings with decimal places - example: "1" becomes "1.0".
  #' `format_chromosome_decimals` removes the decimal place to make these
  #' strings easier to work with
  #'
  #' @param df A dataframe with a column named "chromosome" containing 
  #' chromosome information formatted as a decimal number (i.e. "1.0").
  #'
  #' @returns The dataframe with an additional "chromosome_char" column with
  #' the decimal places removed.
  #' @export
  #'
  #' @examples
  
  stopifnot("chromosome" %in% colnames(df))
  
  output <-  df |> 
    dplyr::mutate(chromosome_char = stringr::str_replace(string = chromosome,
                                      pattern = "\\.0",
                                      replacement = ""))
  
  return(output)
  
}

factorise_chromosome <- function(df) {
  
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
  #' @param df A dataframe containing a column called "chromosome_char" with
  #' chromosome information of type character.
  #'
  #' @returns The input dataframe with a new column called "chromosome_fct" 
  #' containing the factorised version of the chromosome.
  #' @export
  #'
  #' @examples
  
  stopifnot("chromosome_char" %in% colnames(df))
  
  output <- df |> 
    mutate(chromosome_fct = factor(chromosome_char, 
                              levels = define_chromosome_levels()))
  
  return(output)
  
}

add_chromosome_arms <- function(df) {
  
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
  #' @examples
  
  stopifnot(c("chromosome_char", 
              "start", "end") %in% colnames(df))
  
  stopifnot(typeof(df$chromosome_char) == "character")
  
  stopifnot(length(setdiff(unique(df$chromosome_char),
                           define_chromosome_levels())) == 0)
  
  output <- df |> 
    filter(!chromosome_char %in% c("X", "Y")) |> 
    mutate(chromosome_arm = case_when(
      # chr1
      chromosome_char == "1" & end < 122026459 ~"1p",
      chromosome_char == "1" & start > 124932724 ~"1q",
      # chr2
      chromosome_char == "2" & end < 92188145 ~"2p",
      chromosome_char == "2" & start > 94090557 ~"2q",
      # chr3
      chromosome_char == "3" & end < 90772458 ~"3p",
      chromosome_char == "3" & start > 93655574 ~"3q",
      # chr4
      chromosome_char == "4" & end < 49712061 ~"4p",
      chromosome_char == "4" & start > 51743951 ~"4q",
      # chr5
      chromosome_char == "5" & end < 46485900 ~"5p",
      chromosome_char == "5" & start > 50059807 ~"5q",
      # chr6
      chromosome_char == "6" & end < 58553888 ~"6p",
      chromosome_char == "6" & start > 59829934 ~"6q",
      # chr7
      chromosome_char == "7" & end < 58169653 ~"7p",
      chromosome_char == "7" & start > 61528020 ~"7q",
      # chr8
      chromosome_char == "8" & end < 44033744 ~"8p",
      chromosome_char == "8" & start > 45877265 ~"8q",
      # chr9
      chromosome_char == "9" & end < 43389635 ~"9p",
      chromosome_char == "9" & start > 45518558 ~"9q",
      # chr10
      chromosome_char == "10" & end < 39686682 ~"10p",
      chromosome_char == "10" & start > 41593521 ~"10q",
      # chr11
      chromosome_char == "11" & end < 51078348 ~"11p",
      chromosome_char == "11" & start > 54425074 ~"11q",
      # chr12
      chromosome_char == "12" & end < 34769407 ~"12p",
      chromosome_char == "12" & start > 37185252 ~"12q",
      # chr 13
      chromosome_char == "13" & start > 18051248 ~"13q",
      # chr14
      chromosome_char == "14" & start > 18173523 ~"14q",
      # chr15
      chromosome_char == "15" & start > 19725254 ~"15q",
      # chr16
      chromosome_char == "16" & end < 36311158 ~"16p",
      chromosome_char == "16" & start > 36334460 ~"16q",
      # chr17
      chromosome_char == "17" & end < 22813679 ~"17p",
      chromosome_char == "17" & start > 26566633 ~"17q",
      # chr18
      chromosome_char == "18" & end < 15460899 ~"18p",
      chromosome_char == "18" & start > 20861206 ~"18q",
      # chr19
      chromosome_char == "19" & end < 24498980 ~"19p",
      chromosome_char == "19" & start > 27190874 ~"19q",
      # chr20
      chromosome_char == "20" & end < 26436232 ~"20p",
      chromosome_char == "20" & start > 30038348 ~"20q",
      # chr21
      chromosome_char == "21" & end < 10864560 ~"21p",
      chromosome_char == "21" & start > 12915808 ~"21q",
      # chr22
      chromosome_char == "22" & end < 12954788 ~"22p",
      chromosome_char == "22" & start > 15054318 ~"22q"
    ))
  
  if(anyNA.data.frame(output) == TRUE) {
    message("There are NAs in the output")
  }
  
  return(output)

}

calculate_region_size <- function(df) {
  
  output <- df |> 
    mutate(region_size = end-start) |> 
    left_join(chr_coords_grch38 |> 
                select(chrom_arm, length),
              join_by("chromosome_arm" == "chrom_arm")) |> 
    mutate(percent_chr_arm = round((region_size / length) * 100, 1),
           result_string = str_c(name, " ", percent_chr_arm, "%")) 
  
  return(output)
  
}

summarise_by_chromosome_arms <- function(df) {
  
  #' Summarise ploidy regions by chromosome arm
  #'
  #' @param df 
  #'
  #' @returns
  #' @export
  #'
  #' @examples
  
  stopifnot(c("labno_suffix_worksheet",
              "chromosome_char", 
              "chromosome_fct",
              "chromosome_arm", 
              "name", 
              "region_size") %in%
              colnames(df))
  
  output <- df |> 
    group_by(labno_suffix_worksheet, 
             chromosome_char, chromosome_fct, chromosome_arm, name) |> 
    summarise(ploidy_state_bp = sum(region_size),
              ploidy_state_percent = sum(percent_chr_arm)) 

  return(output)
  
}

add_chr_result_flag <- function(df, percent_threshold = 90) {
  
  #' Add a flag for chromosome ploidy states
  #'
  #' @param df 
  #' @param percent_threshold 
  #'
  #' @returns
  #' @export
  #'
  #' @examples
  
  output <- df |> 
    mutate(chr_flag = case_when(
      percent_chr_arm >= percent_threshold &
        name != "Normal diploid" ~str_c(chromosome_arm, " ", name),
      percent_chr_arm < percent_threshold |
        name == "Normal diploid" ~""
    ))
  
  return(output)
  
}

format_arm_ploidy_tbl <- function(df) {
  
  #' Title
  #'
  #' @param df 
  #'
  #' @returns
  #' @export
  #'
  #' @examples
  
  output <- df |> 
    arrange(labno_suffix_worksheet, 
            chromosome_fct, chromosome_arm, desc(percent_chr_arm)) |> 
    group_by(labno_suffix_worksheet, 
             chromosome_fct, chromosome_arm) |> 
    summarise(full_string = toString(result_string),
              full_flag = toString(chr_flag)) |> 
    mutate(result_flag = str_replace_all(full_flag,
                                           pattern = ",\\s",
                                           replacement = "")) |> 
    relocate(result_flag, .after = chromosome_arm) |> 
    select(-full_flag) |> 
    arrange(chromosome_fct) |> 
    ungroup()
  
  return(output)
  
}

make_chr_arm_ploidy_tbl <- function(df, input_threshold = 90) {
  
  #' Title
  #'
  #' @param df 
  #'
  #' @returns
  #' @export
  #'
  #' @examples
  
  output <- df |> 
    extract_pansolid_cnv_coordinates(cnv_coord_col = region) |> 
    format_chromosome_decimals() |> 
    factorise_chromosome() |> 
    add_chromosome_arms() |> 
    calculate_region_size() |> 
    summarise_by_chromosome_arms() |> 
    add_chr_result_flag(percent_threshold = input_threshold) |> 
    format_arm_ploidy_tbl()
  
  return(output)
  
}

source(here::here("tests/test_chromosome_arm_functions.R"))
