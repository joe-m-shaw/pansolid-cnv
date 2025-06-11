# Pan Solid CNV Functions

library(tidyverse)
library(readxl)
library(here)
library(rvest)
library(docstring)

source(here("functions/extract_pansolid_cnv_coordinates.R"))

# Export functions ------------------------------------------------------------------

csv_timestamp <- function(table, folder) {
  
  #' Save a table as a comma-separated file with a timestamp
  #'
  #' @param table The table to save
  #' @param folder The filepath of the folder to save the file into, without a 
  #' backlash at the end
  #'
  #' @return Saves the table in the desired folder with a timestamp
  #'
  #' @note This function is used for exporting tables from RStudio for inclusion
  #' in validation documentation.
  #'
  #' @examples csv_timestamp(gene_tbl_for_doc, paste0(outputs_folder, "tables"))
  
  write.csv(table,
            file = paste0(
              folder,
              "/",
              format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
              "_",
              deparse(substitute(table)), ".csv"
            ),
            row.names = FALSE
  )
}

dna_db_export <- function(input) {
  write_csv(input,
            file = paste0(config::get("data_folderpath"),
                          "validation/DOC6260_ERBB2/",
                          "dna_db_queries/",
                          deparse(substitute(input)), 
                          ".csv")
            )
}

plot_timestamp <- function(input_plot, 
                           input_width = 15, 
                           input_height = 12, 
                           dpi = 300,
                           folder) {
  
  #' Save a plot with a timestamp
  #'
  #' @param input_plot  Plot to save
  #' @param input_width Desired plot width
  #' @param input_height Desired plot height
  #' @param dpi Desired plot resolution in dots per inch
  #' @param folder The folder to save the plot into
  #'
  #' @return Saves the plot in the plots folder with a timestamp
  #'
  #' @note This function is used for exporting plots for inclusions in validation
  #' documentation. The default inputs are designed for presenting a plot as 
  #' half an A4 page.
  #'
  #' @examples plot_timestamp(erbb2_plot)
  
  ggsave(
    filename = paste0(
      format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
      "_",
      deparse(substitute(input_plot)), ".png"
    ),
    plot = input_plot,
    device = "png",
    path = paste0(folder, "/"),
    units = "cm",
    width = input_width,
    height = input_height,
    dpi = 300
  )
}

# Data functions --------------------------------------------------------------------

format_repeat_table <- function(df) {
  
  rpt_table <- df |> 
    arrange(labno) |> 
    select(labno_suffix, worksheet, gene, max_region_fold_change,
           st_dev_signal_adjusted_log2_ratios) |>  
    mutate("ERBB2 fold change" = round(max_region_fold_change, 1),
           "Signal adjusted noise" = round(st_dev_signal_adjusted_log2_ratios, 2)) |> 
    select(-c(gene, max_region_fold_change,
              st_dev_signal_adjusted_log2_ratios))

  return(rpt_table)
  
}

extract_cnv_calls <- function(df, input_gene) {
  
  #' Extract CNV calls from DNA Database genotype comments
  #'
  #' @param df A dataframe including a "genotype" column of QIAseq results from the 
  #' DNA Database.
  #' @param input_gene The gene of interest
  #'
  #' @return  A dataframe including the input dataframe with additional columns
  #' generated from the "genotype" column for the gene CNV result and dosage quotient.
  #' 
  #' @note Genotypes are entered onto the DNA Database as a text string with 
  #' a consistent structure. This function parses results for a 
  #'
  #' @examples data <- data.frame(labno = c(1),
  #' genotype = c("No SNVS. ERBB2 amplification detected (Mean DQ 11x)"))
  #' 
  #' cnvs <- extract_cnv_calls(df = data, input_gene = "ERBB2")
  
  stopifnot("genotype" %in% colnames(df))
  
  dq_regex <- regex(str_c(
    # Group - input gene
    "(",input_gene,")\\s",
    # Group - variable freetype,
    "(amplification\\sdetected|amplification)",
    # Use . for bracket
    "\\s.",
    # Group - variable freetype
    "(Mean\\sDQ|mean\\sDQ|DQ)",
    "\\s",
    # Group - dosage quotient regex
    "(\\d{1,3}|\\d{1,3}\\.\\d{2})",
    "x"))
  
  output <- df |> 
    select(labno, genotype) |> 
    mutate(gene_searched = input_gene,
           gene_match = str_extract(genotype, dq_regex, group = 1),
           gene_dq = as.numeric(str_extract(genotype, dq_regex, group = 4)),
           core_result = ifelse(!is.na(gene_dq), "Amplification", "No call"))
  
  return(output)
  
}

read_summary_tab <- function(file) {
  
  x <- read_excel(path = file,
                  sheet = "Whole Panel UMI Coverage Re...",
                  skip = 1,
                  n_max = 11) |> 
    dplyr::rename(value = "...2")
  
  x_wide <- x |> 
    pivot_wider(names_from = Summary,
                values_from = value) |> 
    # Renaming as >, < and ≥ are removed in clean_names
    # Names shortened for ease of use
    rename(number_target_regions_with_cov_lessthan_138 = `Number of target regions with coverage < 138`,
           total_length_target_regions_with_pos_cov_lessthan_138 = `Total length of target regions containing positions with coverage < 138`,
           total_length_target_region_pos_cov_lessthan_138 = `Total length of target region positions with coverage < 138`,
           total_length_target_region_pos_cov_greaterorequal_138 = `Total length of target region positions with coverage ≥ 138`,
           percent_target_region_pos_cov_greaterorequal_138 = `Percentage of target region positions with coverage ≥ 138 (%)`) |> 
    janitor::clean_names() 
  
  identifiers <- filename_to_df(file)
  
  output <- cbind(identifiers, x_wide)
  
  return(output)
  
}

format_chromosome <- function(df, input_col) {
  
  #' Format a chromosome string
  #'
  #' @param df Input dataframe
  #' @param input_col Column containing chromosome information
  #'
  #' @return The input dataframe with an additional column which presents chromosome
  #' information as a single character string. Autosomes are presented as characters
  #' without decimal places (i.e. "1" instead of "1.0")
  #'
  #' @examples formatted_data <- format_chromosome(df = data, input_col = "chromosome")
  
  output <- df |> 
    mutate(chrom_mod = case_when(
    
      {{ input_col }} %in% c("X", "Y") ~{{ input_col }},
      
      TRUE ~as.character(round(as.numeric({{ input_col }}), 0))),
    
    chromosome_formatted = fct(x = chrom_mod, levels = c("1", "2", "3", "4",
                                               "5",  "6",  "7",  "8",
                                               "9",  "10", "11", "12",
                                               "13", "14", "15", "16",
                                               "17", "18", "19", "20",
                                               "21", "22", "X", "Y")))
  
  return(output)
  
}

read_targeted_region_overview <- function(file) {
  
  # This function reads the table of reads mapped to each chromosome. 
  # The position of this table in the "Whole Panel UMI Coverage" tab varies in each file
  
  x <- read_excel(path = file,
                  sheet = "Whole Panel UMI Coverage Re...")
  
  num_skip <- match("Targeted region overview", x$`Target regions`) + 1
  
  targeted_region_overview <- read_excel(path = file,
                                         sheet = "Whole Panel UMI Coverage Re...",
                                         skip = num_skip,
                                         # 22 autosomes plus 2 sex chromosomes
                                         n_max = 24) |> 
    format_chromosome(input_col = Reference)

  identifiers <- filename_to_df(file)
  
  output <- cbind(identifiers, targeted_region_overview) |> 
    janitor::clean_names()
  
  return(output)
  
}

get_control_coverage <- function(file) {
  
  identifiers <- filename_to_df(file)
  
  df <- read_csv(file) |> 
    janitor::clean_names()
  
  cov <- data.frame(
    "median_coverage" = median(df$coverage),
    "mean_coverage" = mean(df$coverage))
  
  output <- cbind(identifiers, cov)
  
  return(output)
  
}

calculate_target_copies <- function(fold_change, ncc_percent) {
  
  ref_sample_copies <- 200
  
  ncc_fraction <- ncc_percent / 100
  
  sample_total_target_copies <- fold_change * ref_sample_copies
  
  sample_target_copies_in_tumour_cells <- sample_total_target_copies - ((100 * (1-ncc_fraction)) * 2)
  
  sample_target_copies_per_tumour_cell <- sample_target_copies_in_tumour_cells / (100 * ncc_fraction)
  
  return(sample_target_copies_per_tumour_cell)
  
}

true_pos <- "True positive"

true_neg <- "True negative"

false_pos <- "False positive"

false_neg <- "False negative"

classifiers = c(true_pos, true_neg, false_pos, false_neg)

make_confusion_matrix <- function(df, input_column = outcome,
                                  classifiers = c(true_pos, true_neg, false_pos, false_neg),
                                  initial_test,
                                  comparison_test,
                                  positive_state,
                                  negative_state) {
  
  # This function requires an input table with true and false positives and negatives already defined.
  
  true_positives <- nrow(df |> 
                           filter({{ input_column }} == classifiers[1]))
  
  true_negatives <- nrow(df |> 
                           filter({{ input_column }} == classifiers[2]))
  
  false_positives <- nrow(df |> 
                            filter({{ input_column }} == classifiers[3]))
  
  false_negatives <- nrow(df |> 
                            filter({{ input_column }} == classifiers[4]))
  
  tp_char <- as.character(true_positives)
  
  tn_char <- as.character(true_negatives)
  
  fp_char <- as.character(false_positives)
  
  fn_char <- as.character(false_negatives)
  
  conf_matrix <- tribble(
    ~"",               ~"",              ~"",                 ~"",         
    "",                "",               initial_test,        "", 
    "",                "",               positive_state,      negative_state, 
    comparison_test,   positive_state,   tp_char,             fn_char,
    "",                negative_state,   fp_char,             tn_char)
  
  # Overall percent agreement
  
  opa <- round((true_positives + true_negatives) / (true_positives + false_negatives +
                                                      false_positives + true_negatives) * 100, 1)
  
  # Positive percentage agreement
  
  ppa <- round((true_positives) / (true_positives + false_negatives) * 100, 1)
  
  # Negative percentage agreement
  
  npa <- round((true_negatives) / (true_negatives + false_positives) * 100, 1)
  
  return(list(conf_matrix, opa, ppa, npa))
  
}

read_clc_target_calls <- function(file) {
  
  identifiers <- filename_to_df(file)
  
  results <- read_excel(path = file, sheet = 2,
           col_types = c("text", "text", "text",
                         "numeric", "numeric", "numeric",
                         "numeric", "numeric", "numeric",
                         "numeric", "numeric", "numeric",
                         "text", "text",
                         "numeric", "numeric",
                         "text", "text", "text")) |> 
  janitor::clean_names() |> 
  mutate(
    labno = as.character(parse_filename(file, 2)),
    filename = file) |> 
  left_join(identifiers, by = "labno") |> 
  relocate(worksheet, labno, suffix, labno_suffix, patient_name,
           labno_suffix_worksheet)
  
  output <- extract_pansolid_cnv_coordinates(df = results,
                                    cnv_coord_col = region)
  
  return(output)
  
}

calculate_pooled_sd <- function(df, group = labno, target_col, round_places = 2) {
  
  output_table <- df |> 
    group_by( {{ group }}) |> 
    summarise(sd = sd( {{ target_col }} ),
              max = max( {{ target_col }} ),
              min = min( {{ target_col }} ),
              range = max - min,
              n = n(),
              z = (n-1)*sd^2)
  
  pooled_sd <- round(sqrt(sum(output_table$z) / 
                            (sum(output_table$n))), round_places)
  
  range <- str_c(round(min(output_table$range), round_places), 
                 "-", 
                 round(max(output_table$range), round_places))
  
  return(list(output_table, pooled_sd, range))
  
}

add_dna_db_info <- function(df, 
                            ps_version_df = pansolid_ws_details,
                            extraction_df = sample_extraction_details,
                            gender_df = sample_gender,
                            type_df = sample_types,
                            ncc_df = ncc_collated,
                            tissue_df = sample_tissue_sources_coded,
                            dna_conc_df = sample_dna_concentrations,
                            pathno_df = sample_pathnos,
                            nhsno_df = sample_nhs_no) {
  
  # This is a wrapper function that joins useful information from DNA database onto
  # the results.
  
  if (!"labno" %in% colnames(df) |
      !"worksheet" %in% colnames(df)) { stop("Join columns not present")}
  
  output <- df |> 
    left_join(ps_version_df, by = "worksheet") |> 
    left_join(extraction_df, by = "labno") |> 
    left_join(gender_df, by = "labno") |> 
    left_join(type_df, by = "labno") |> 
    left_join(ncc_df, by = "labno") |> 
    left_join(tissue_df, by = "labno") |> 
    left_join(dna_conc_df, by ="labno") |> 
    left_join(pathno_df, by = "labno") |> 
    left_join(nhsno_df, by = "labno")
  
  return(output)
  
}

add_case_group <- function(df) {
  
  stopifnot("patient_name" %in% colnames(df))
  
  output <- df |> 
    mutate(sample_group = "case",
           sample_subgroup = case_when(
             
             patient_name %in% grep(pattern = "seraseq", x = patient_name,
                                    ignore.case = TRUE, value = TRUE) ~"SeraCare reference material",
             
             TRUE ~"Patient FFPE sample"))
  
  return(output)
  
}

draw_lod_gene_plot <- function(df, chromosome, gene) {
  
  plot_limit_of_detection <- df |> 
    filter(chromosome == {{ chromosome }}) |> 
    ggplot(aes(x = start, y = fold_change_adjusted)) +
    geom_point(pch = 21) +
    geom_point(data = df |> 
                 filter(name == {{ gene }}), fill = safe_red, 
               pch = 21, size = 2) +
    facet_wrap(~ncc) +
    theme_bw() +
    scale_y_continuous(limits = c(-3, 6),
                       breaks = c(-3, -2, -1, 0, 1, 2, 2.8, 4, 5, 6)) +
    geom_hline(yintercept = 2.8, linetype = "dashed") +
    labs(x = str_c("Chromosome ", {{ chromosome }}),
         y = "Target fold change",
         title = str_c("Limit of detection results: ", {{ gene }}),
         caption = "Seracare +12 copies control spiked into Seracare wild type control",
         subtitle = str_c({{ gene }}, " in red"))
  
  return(plot_limit_of_detection)
  
}

calculate_pooled_sd <- function(df, group = labno, target_col, round_places = 2) {
  
  output_table <- df |> 
    group_by( {{ group }}) |> 
    summarise(sd = sd( {{ target_col }} ),
              max = max( {{ target_col }} ),
              min = min( {{ target_col }} ),
              range = max - min,
              n = n(),
              z = (n-1)*sd^2)
  
  pooled_sd <- round(sqrt(sum(output_table$z) / 
                            (sum(output_table$n))), round_places)
  
  range <- str_c(round(min(output_table$range), round_places), 
                 "-", 
                 round(max(output_table$range), round_places))
  
  return(list(output_table, pooled_sd, range))
  
}

make_noise_plot <- function(df = qc_data_with_ids, x_axis) {
  
  plot <- ggplot(qc_data_with_ids, aes(x = {{ x_axis }}, 
                                       y = st_dev_signal_adjusted_log2_ratios)) +
    geom_point(pch = 21) +
    theme_bw() +
    scale_y_continuous(limits = c(0, 1.25),
                       breaks = seq(0, 1.2, by = 0.2)) +
    labs(y = "Signal-adjusted noise") +
    geom_hline(yintercept = 1, linetype = "dashed")
  
  return(plot)
  
}
