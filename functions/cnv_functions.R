# Pan Solid CNV Functions

library(tidyverse)
library(readxl)
library(here)
library(rvest)
library(docstring)

source(here("scripts/set_shared_drive_filepath.R"))

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
  write.csv(input,
            file = here::here(paste0(
              "data/dna_db_queries/",
              deparse(substitute(input)), ".csv"
            )),
            row.names = FALSE
  )
}

plot_timestamp <- function(input_plot, input_width = 15, input_height = 12, dpi = 300,
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
  
  output <- extract_cnv_coordinates(df = results,
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

# Primers ---------------------------------------------------------------------------

grch38_primers <- read_csv(file = paste0(data_folder,
                                         "primers/CDHS-40079Z-11284.primer3_Converted.csv"),
                           show_col_types = FALSE) |> 
  janitor::clean_names()

grch38_primer_coordinates <- extract_cnv_coordinates(df = grch38_primers,
                                                     cnv_coord_col = region)

count_primers_in_region <- function(chrom, coord1, coord2, df) {
  #' Count the QIAseq primers within a region
  #'
  #' @param chrom The chromosome of the region of interest
  #' @param coord1 The first coordinate of the region of interest
  #' @param coord2 The second coordinate of the region of interest
  #' @param df The dataframe containing the primer coordinate information.
  #' This should be in the format of one primer per row with one column each for the 
  #' primer start and end coordinates. 
  #' 
  #' @return Returns the number of primers within the specified region
  #' 
  #' @examples count_primers_in_region(chrom = "17", coord1 = 39700064,
  #' coord2 = 39728658)
  
  if(typeof(chrom) != "character") {
    stop("Chromosome argument must be a character")
  }
  
  if (typeof(coord1) != "double" |
      typeof(coord2) != "double") {
    stop("Coordinates must be type double")
  }
  
  if(!("primer_start" %in% colnames(df) & "primer_end" %in% colnames(df))) {
    stop("Primer dataframe does not contain primer_start and primer_end columns")
  }
  
  region_start <- min(coord1, coord2)
  
  region_end <- max(coord1, coord2)
  
  df2 <- df |> 
    filter(chromosome == chrom) |> 
    mutate(in_region = case_when(
      
      (primer_start >= region_start &
         primer_start <= region_end) |
        (primer_end >= region_start &
           primer_end <= region_end) ~"Yes",
      TRUE ~"No"
      
    ))
  
  primers_in_region <- sum(df2$in_region == "Yes")
  
  return(primers_in_region)
  
}

# Genes and exons -------------------------------------------------------------------

transcript_regex <- regex(
  r"(
  .+
  (ENST\d{11})
  .csv
  )",
  comments = TRUE
)

read_ensembl_exon_table <- function(filename) {
  
  transcript_id <- str_extract(string = filename, 
                               pattern = transcript_regex,
                               group = 1)
  
  table <- read_csv(file = here::here(filename),
                    show_col_types = FALSE) |> 
    janitor::clean_names() |> 
    filter(!is.na(no)) |> 
    rename(exon = no) |> 
    mutate(transcript = transcript_id) |> 
    relocate(transcript) |> 
    select(-sequence)
  
  return(table)
  
}

gene_labels <- read_excel(path = paste0(data_folder, "transcripts/gene_labels.xlsx"),
                        col_types = c("text", "text", "text", "text",
                                      "numeric", "numeric")) |> 
  mutate(y_value = "Genes",
         # Place gene label half-way along gene locus
         start = pmin(gene_start, gene_end) + ((pmax(gene_start, gene_end) - pmin(gene_start, gene_end)) / 2))

transcript_files <- list.files(paste0(data_folder, "transcripts/"), full.names = TRUE,
                               pattern = ".csv")

all_transcripts <- transcript_files |>
  map(\(transcript_files) read_ensembl_exon_table(
    file = transcript_files
  )) |>
  list_rbind() |> 
  left_join(gene_labels |> 
              select(label, chromosome, transcript_ensembl), join_by(transcript == transcript_ensembl))
  
# Plot colours ----------------------------------------------------------------------

safe_blue <- "#88CCEE"
safe_red <- "#CC6677"
safe_grey <- "#888888"

# Plot functions --------------------------------------------------------------------

# CNV plots can be presented as triptychs: 
# Panel 1) The plot of CNV calls:
    # a) Either with fold change on the y axis (make_fold_change_plot)
    # b) Or with lab number on the y axis (make_labno_plot)
# Panel 2) A plot showing the locations of Qiaseq primers
# Panel 3) A plot showing annotated locations of gene exons

# The aim of these inter-related functions is to allow maximum flexibility and to keep 
# the x axes consistent between the different plots.

get_breaks <- function(interval, plot_xmin, plot_xmax) {
  
  breaks <- seq(plot_xmin, plot_xmax, by = interval)
  
  return(breaks)

}

get_data_for_plot <- function(df, 
                              gene) {
  
  data_for_plot <- df |> 
    filter(gene == {{ gene }})
  
  return(data_for_plot)
  
}

get_plot_xmin <- function(df, buffer) {
  
  plot_xmin <- min(df$start) - buffer
  
  return(plot_xmin)
  
}

get_plot_xmax <- function(df, buffer) {
  
  plot_xmax <- max(df$end) + buffer
  
  return(plot_xmax)
  
}

get_chromosome <- function(gene) {
  
  chromosome <- as.character(gene_labels[gene_labels$label == gene, 2])
  
  stopifnot(chromosome != "character(0)")
  
  return(chromosome)
  
}

make_fold_change_plot <- function(df = pos_cnv_results, 
                                  gene = "ERBB2",
                                  interval = 10000, 
                                  buffer = 5000, 
                                  ymin = 0,
                                  ymax = 40) {
  
  chromosome <- get_chromosome(gene = {{ gene }})
  
  data_for_plot <- get_data_for_plot(df = {{ df }}, 
                    gene = {{ gene }})
  
  plot_xmin <- get_plot_xmin(df = data_for_plot,
                                 buffer = buffer)
  
  plot_xmax <- get_plot_xmax(df = data_for_plot,
                                 buffer = buffer)
  
  fold_change_plot <- ggplot(data_for_plot, aes(x = start, y = fold_change)) +
    
    # Add theme
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    
    # Add CNV calls
    geom_segment(aes(x = start, xend = end, 
                     y = fold_change, yend = fold_change),
                 linewidth = 2,
                 colour = safe_red) +

    # Add x axes
    scale_x_continuous(breaks = get_breaks(interval = {{ interval}},
                                           plot_xmin = {{ plot_xmin }},
                                           plot_xmax = {{ plot_xmax }}),
                       minor_breaks = NULL,
                       limits = c({{ plot_xmin }}, {{ plot_xmax }} )) +
    
    scale_y_continuous(limits = c(ymin, ymax)) +
    
    # Add labels
    labs(
      y = "Fold change",
      x = "")
  
  return(list(plot_xmin, plot_xmax, interval, fold_change_plot, chromosome))
  
}

make_labno_plot <- function(df = pos_cnv_results, 
                            gene = "ERBB2",
                            interval = 10000, 
                            buffer = 5000, 
                            yaxis = labno) {
  
  chromosome <- get_chromosome(gene = {{ gene }})
  
  data_for_plot <- get_data_for_plot(df = {{ df }}, 
                                     gene = {{ gene }})
  
  max_fold_change <- max(data_for_plot$fold_change)
  
  min_fold_change <- min(data_for_plot$fold_change)
  
  plot_xmin <- get_plot_xmin(df = data_for_plot,
                             buffer = buffer)
  
  plot_xmax <- get_plot_xmax(df = data_for_plot,
                             buffer = buffer)
  
  labno_plot <- ggplot(data_for_plot, aes(x = start, y = {{ yaxis }},
                                                colour = fold_change)) +
    
    # Add theme
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    
    # Add CNV calls
    geom_segment(aes(x = start, xend = end, 
                     y = {{ yaxis }}, yend = {{ yaxis }}),
                 linewidth = 2) +
    
    scale_colour_gradient(low = "#FF9999", 
                          high = "#660000", 
                          limits = c(min_fold_change, max_fold_change),
                          n.breaks = 4) +
    
    # Add x axes
    scale_x_continuous(breaks = get_breaks(interval = {{ interval}},
                                           plot_xmin = {{ plot_xmin }},
                                           plot_xmax = {{ plot_xmax }}),
                       minor_breaks = NULL,
                       limits = c({{ plot_xmin }}, {{ plot_xmax }} )) +
    
    # Add labels
    labs(
      y = "Sample number",
      x = "",
      colour = "Fold change")
  
  return(list(plot_xmin, plot_xmax, interval, labno_plot, chromosome))
  
}

make_primer_plot <- function(plot_xmin, plot_xmax, interval, chromosome) {
  
  primers_filtered <- grch38_primer_coordinates |> 
    mutate(y_value = "Primers") |> 
    filter(chromosome == {{ chromosome }} ) |> 
    filter(start >= {{ plot_xmin }} & end <= {{ plot_xmax }} )
  
  output <-  ggplot(primers_filtered, aes(x = start, y = y_value)) +
    geom_point(pch = 21) +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    scale_x_continuous(breaks = get_breaks(interval = {{ interval}},
                                           plot_xmin = {{ plot_xmin }},
                                           plot_xmax = {{ plot_xmax }}),
                       minor_breaks = NULL,
                       limits = c({{ plot_xmin }}, {{ plot_xmax }} )) +
    labs (x = "", y = "")
  
  return(output)
  
}

make_exon_plot <- function(plot_xmin, plot_xmax, interval, chromosome) {
  
  exon_data_for_plot <- all_transcripts |> 
    mutate(y_value = "Exons") |> 
    filter(chromosome == {{ chromosome }}) |> 
    filter(start >= {{ plot_xmin }} & end <= {{ plot_xmax }})
  
  labels_for_plot <- gene_labels |> 
    filter(chromosome == {{ chromosome }} ) |> 
    filter(start >= {{ plot_xmin }} & start <= {{ plot_xmax }})
  
  output <- ggplot(exon_data_for_plot, 
                   aes(x = start, y = y_value)) +
    
    geom_segment(aes(x = start, xend = end, 
                     y = y_value, yend = y_value),
                 linewidth = 5) +
    
    theme_bw() +
    
    scale_x_continuous(breaks = get_breaks(interval = {{ interval}},
                                           plot_xmin = {{ plot_xmin }},
                                           plot_xmax = {{ plot_xmax }}),
                       minor_breaks = NULL,
                       limits = c({{ plot_xmin }}, {{ plot_xmax }} )) +
    
    scale_y_discrete(limits = rev) +
    
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    
    labs(y = "", x = str_c("Genome coordinate (GRCh38) Chr", 
                           chromosome)) +
    
    geom_label(data = labels_for_plot, label = labels_for_plot$label)
  
  return(output)
  
}

make_cnv_triptych <- function(input_plot) {
 
  # This function is a wrapper which takes the outputs of either the 
  # make_fold_change_plot or make_labno_plot functions
  
  plot_xmin <- input_plot[[1]]
  
  plot_xmax <- input_plot[[2]]
  
  interval <- input_plot[[3]]
  
  main_plot <- input_plot[[4]]
  
  chromosome <- input_plot[[5]]

  primer_plot <- make_primer_plot(plot_xmin = {{ plot_xmin }}, 
                                  plot_xmax = {{ plot_xmax }},
                                  interval = {{ interval }},
                                  chromosome = {{ chromosome }})
  
  exon_plot <- make_exon_plot(plot_xmin = {{ plot_xmin }}, 
                              plot_xmax = {{ plot_xmax }},
                              interval = {{ interval }},
                              chromosome = {{ chromosome }})
  
  triptych <- (main_plot / primer_plot / exon_plot) +
    plot_layout(
      heights = c(6, 1, 2)
    )
  
  return(triptych)
  
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

# PanSolid functions ----------------------------------------------------------------

join_pansolid_submission_sheets <- function() {
  
  #' Load and join PanSolid DNA submission sheets
  #'
  #' @return A tidy dataframe for all PanSolid submissions from 2022 to 2024, based on 
  #' files saved in the local drive.
  #' 
  #' @note Excel spreadsheets for coordinating testing on the PanSolid QIAseq 
  #' enrichment are manually curated by the tech team and stored either on the S
  #' drive or on the laboratory Sharepoint. This function uses versions of 
  #' these Excel spreadsheets copied onto a local drive.
  #' 
  #' The sample_id field has some surprises. Example: 23024772 has a degree sign (°)
  #' entered after it which is invisible in Excel and R.
  #'
  #' @examples pansolid_sheets <- join_pansolid_submission_sheets()
  #' 
  #' sample_info <- pansolid_sheets |> filter(labno == "12345678")
  
  pansolid_submission_2023 <- read_excel(path = paste0(data_folder, 
                                                       "dna_submission_sheets/DNA PanSolid QIAseq Submission Sheet 2023.xlsx")) |> 
    janitor::clean_names() |> 
    rename(stock_qubit = stock_qubit_ng_m_l) |> 
    mutate(submission_sheet = "2023",
           labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
    select(date_submitted, labno, sample_name,
           panel, enrichment, stock_qubit, submission_sheet)
  
  pansolid_submission_2024 <- read_excel(path = paste0(data_folder, "dna_submission_sheets/PanSolid Submission sheet 2024.xlsx"),
                                         sheet = "PanSolid samples") |> 
    janitor::clean_names()  |> 
    rename(stock_qubit = stock_qubit_ng_m_l) |> 
    mutate(submission_sheet = "2024",
           labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
    select(date_submitted, labno, sample_name,
           panel, enrichment, stock_qubit, submission_sheet)
  
  # Pansolid began in 2022 so the initial runs were recorded on the Qiaseq spreadsheet
  pansolid_submission_2022 <- read_excel(path = paste0(data_folder, 
                                                       "dna_submission_sheets/QIAseq DNA PanSolid Sample Submission 2022.xlsx")) |> 
    janitor::clean_names() |> 
    rename(date_submitted = date_sample_submitted,
           stock_qubit = stock_qubit_ng_m_l) |> 
    mutate(submission_sheet = "2022",
           labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
    select(date_submitted, labno, sample_name,
           panel, enrichment, stock_qubit, submission_sheet)
  
  output <- rbind(pansolid_submission_2024,
                  pansolid_submission_2023,
                  pansolid_submission_2022)
  
  return(output)
  
}



get_gene_list <-  function() {
  
  amp_gene_list <- read_excel(path = paste0(data_folder, 
                                            "gene_lists/pansolid_amplification_gene_list.xlsx"))
  
  del_gene_list <- read_excel(path = paste0(data_folder, 
                                            "gene_lists/pansolid_deletion_gene_list.xlsx"))
  
  gene_list <- rbind(amp_gene_list, del_gene_list)
  
  return(gene_list)
  
}

# Reformatting functions ------------------------------------------------------------

load_gene_table <- function(cnv_type) {
  
  if(!cnv_type %in% c("Deletions", "Amplifications")) {
    stop("cnv_type must be Deletions or Amplifications")
  }
  
  if(cnv_type == "Deletions") {
    
    gene_table <- read_excel(path = paste0(data_folder, 
                                           "gene_lists/pansolid_deletion_gene_list.xlsx"))
    
  }
  
  if(cnv_type == "Amplifications") {
    
    gene_table <- read_excel(path = paste0(data_folder, 
                                           "gene_lists/pansolid_amplification_gene_list.xlsx"))
    
  }
  
  return(gene_table)
  
}



reformat_wgs_cnv_result <- function(filepath, cnv_type) {
  
  gene_table <- load_gene_table({{ cnv_type }})
  
  sample_cnvs <- wgs_data_collated |> 
    filter(filepath == {{ filepath }}) |> 
    mutate(gene = str_replace_all(string = gene, pattern = "\\*",
                                  replacement = ""))
  
  if(cnv_type == "Deletions") {
    
    sample_cnvs_filtered <- sample_cnvs |> 
      filter(cnv_class %in% c("LOSS", "DEL")) |> 
      filter(gene %in% gene_table$gene)
    
  }
  
  if(cnv_type == "Amplifications") {
    
    sample_cnvs_filtered <- sample_cnvs |> 
      filter(cnv_class %in% c("DUP", "GAIN") 
             & cnv_copy_number > 10
      ) |> 
      filter(gene %in% gene_table$gene)
    
  }
  
  if(cnv_type == "Deletions") {
    
    gene_result_table <- gene_table |> 
      mutate(wgs_result = case_when(
        gene %in% sample_cnvs_filtered$gene ~"Loss detected",
        !gene %in% sample_cnvs_filtered$gene ~"No loss detected"
      )) 
    
  }
  
  if(cnv_type == "Amplifications") {
    
    gene_result_table <- gene_table |> 
      mutate(wgs_result = case_when(
        gene %in% sample_cnvs_filtered$gene ~"Gain detected",
        !gene %in% sample_cnvs_filtered$gene ~"No gain detected"
      )) 
    
  }
  
  output_table <- gene_result_table |> 
    mutate(filepath = filepath) |> 
    relocate(filepath)
  
  return(output_table)
  
}

reformat_pansolid_cnv_result <- function(filepath, cnv_type) {
  
  gene_table <- load_gene_table({{ cnv_type }})
  
  sample_cnvs <- pansolid_data_collated |> 
    filter(filepath == {{ filepath }}) |> 
    filter(gene %in% gene_table$gene)
  
  if(cnv_type == "Deletions") {
    
    gene_result_table <- gene_table |> 
      mutate(pansolid_result = case_when(
        gene %in% sample_cnvs$gene ~"Loss detected",
        !gene %in% sample_cnvs$gene ~"No loss detected"
      )) 
    
  }
  
  if(cnv_type == "Amplifications") {
    
    gene_result_table <- gene_table |> 
      mutate(pansolid_result = case_when(
        gene %in% sample_cnvs$gene ~"Gain detected",
        !gene %in% sample_cnvs$gene ~"No gain detected"
      )) 
    
  }
  
  output_table <- gene_result_table |> 
    mutate(filepath = filepath) |> 
    relocate(filepath)
  
  return(output_table)
  
}
