# Pan Solid CNV Functions

library(tidyverse)
library(here)

# DNA database connection -----------------------------------------------------------

dbi_con <- DBI::dbConnect(
  drv = odbc::odbc(),
  dsn = "moldb")

# Export functions ------------------------------------------------------------------

csv_timestamp <- function(input) {
  write.csv(input,
            file = here::here(paste0(
              "outputs/tables/",
              format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
              "_",
              deparse(substitute(input)), ".csv"
            )),
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

plot_timestamp <- function(input_plot, input_width = 15, input_height = 12, dpi = 300) {
  
  # Default inputs allow for presenting a plot as half an A4 page
  
  ggsave(
    filename = paste0(
      format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
      "_",
      deparse(substitute(input_plot)), ".png"
    ),
    plot = input_plot,
    device = "png",
    path = here::here("outputs/plots/"),
    units = "cm",
    width = input_width,
    height = input_height,
    dpi = 300
  )
}


# Data functions --------------------------------------------------------------------

filename_regex <- regex(
  r"[
  (WS\d{6})              # Worksheet number
  _
  (\d{8})                # Lab number
  (a|b|c|)               # Suffix
  _
  ([:alnum:]{5,30})      # Patient name - alphanumeric characters only
  (.xlsx|_S.+.xlsx|_S.+) # Ending varies between patients and controls
  ]",
  comments = TRUE
)

parse_filename <- function(input_file, input_group) {
  
  output <- str_extract(input_file, filename_regex,
                        group = input_group)
  
  return(output)
  
}

filename_to_df <- function(file) {
  
  output <- data.frame(
    worksheet = c(parse_filename(file, 1)),
    labno = c(as.character(parse_filename(file, 2))),
    suffix = c(parse_filename(file, 3)),
    patient_name = c(parse_filename(file, 4))) |> 
    mutate(
      labno_suffix = str_c(labno, suffix),
      labno_suffix_worksheet = str_c(labno_suffix, "_", worksheet))
  
  return(output)
  
}

coarse_tab <- "Oncogenes (Amplified) Coars..."
fine_tab <- "Oncogenes (Amplified) Fine-..."

read_clc_amp_calls <- function(file, input_sheet) {
  
  stopifnot(input_sheet %in% c(coarse_tab, fine_tab))
  
  x <- read_excel(path = file, sheet = input_sheet) 
  
  # Samples with no calls currently have a table with 0 rows
  if(nrow(x) == 0) {
    
    x <- tibble("Chromosome" = 0, 
                "Region" = "",
                "Name" = "", 
                "Region length" = 0, 
                "type" = "",
                "source" = "", 
                "ID"= "", 
                "GeneID"= "", 
                "HGNC"= "",
                "MIM"= "",
                "description"= "", 
                "gbkey"= "", 
                "gene"= "No calls",
                "gene_biotype"= "", 
                "gene_synonym"= "", 
                "CNV region"= "",
                "CNV region length"= 0,
                "Consequence"= "",
                "Fold-change (adjusted)"= 0,
                "p-value"= 0,
                "Number of targets"= 0,
                "Comments"= "", 
                "Targets"= "")
    
  }
  
  identifiers <- filename_to_df(file)
  
  results <- x |> 
    janitor::clean_names() |> 
    mutate(
      labno = as.character(parse_filename(file, 2)),
      setting = input_sheet,
      filename = file) |> 
    left_join(identifiers, by = "labno") |> 
    relocate(worksheet, labno, suffix, labno_suffix, patient_name,
             labno_suffix_worksheet, setting)
  
  return(results)
  
}

extract_cnv_coordinates <- function(df) {
  
  stopifnot("cnv_region" %in% colnames(df))
  
  cnv_coord_regex <- regex(
    r"[
        (|complement\()
        (\d{1,10})     # first coordinate number (1 to 10 digits)
        \.\.           # two full stops
        (\d{1,10})     # second coordinate number (1 to 10 digits)
        ]",
    comments = TRUE
  )
  
  output <- df |> 
    mutate(start = as.numeric(str_extract(string = cnv_region, 
                                          pattern = cnv_coord_regex, 
                                          group = 2)),
           end = as.numeric(str_extract(string = cnv_region, 
                                        pattern = cnv_coord_regex, 
                                        group = 3))) 
  
  return(output)
  
}

ncc_regex <- regex(
  r"[
  (>\d{2}% | \d{2}-\d{2}%)
  ]",
  comments = TRUE
)

parse_ncc <- function(input_string) {
  
  # Function for parsing neoplastic cell content values from the comments
  # column of DNA Database
  
  ncc <- str_extract(string = input_string, 
              pattern = ncc_regex, 
              group = 1)
  
  return(ncc)
  
}

# Plot functions --------------------------------------------------------------------

safe_blue <- "#88CCEE"
safe_red <- "#CC6677"
safe_grey <- "#888888"

plot_coarse_v_fine <- function(input_gene) {
  
  plot <- all_calls |> 
    filter(name == input_gene) |> 
    group_by(sample_id, setting) |> 
    summarise(calls = n()) |> 
    ggplot(aes(x = sample_id, y = calls)) +
    geom_col(aes(fill = setting), position = "dodge") +
    theme_bw() +
    labs(title = as.character(input_gene))
  
  return(plot)
  
}

# GRCh38 coordinates from NCBI
egfr_min <- 55019017
egfr_max <- 55211628

erbb2_min <- 39700064
erbb2_max <- 39728658

met_min <- 116672196
met_max <- 116798377

braf_min <- 140713328
braf_max <- 140924929

myc_min <- 127735434 
myc_max <- 127742951

by_n <- function(n) { 
  
  max_value <- max(egfr_min, egfr_max, erbb2_min, erbb2_max, 
                   met_min, met_max, braf_min, braf_max, myc_min, myc_max)
  
  seq(0, max_value + 10000000, by = n) 
  
  }


draw_cnv_plot <- function(df, input_gene, input_setting, 
                          interval = 200000,
                          ymax = 70,
                          gene_min,
                          gene_max,
                          buffer = 150000) {
  
  grch38_primers <- read_csv(file =
                               here::here("data/CDHS-40079Z-11284.primer3_Converted.csv"),
                             show_col_types = FALSE) |> 
    janitor::clean_names()
  
  grch38_primer_coordinates <- extract_cnv_coordinates(grch38_primers |> 
                                                         dplyr::rename(cnv_region = region))
  
  
  stopifnot(input_gene %in% c("EGFR", "ERBB2", "MET",
                              "BRAF", "MYC"))

  data_for_plot <- df |> 
    filter(setting == input_setting) |> 
    filter(gene == input_gene) 
  
  cnv_min <- min(data_for_plot$coordinate)
  
  cnv_max <- max(data_for_plot$coordinate)
  
  plot_min <- cnv_min - buffer
  
  plot_max <- cnv_max + buffer
  
  primers_for_plot <- grch38_primer_coordinates |> 
    dplyr::rename(coordinate = cnv_start) |> 
    mutate(fold_change_adjusted = 0) |> 
    filter(coordinate >= plot_min  & coordinate <= plot_max )
  
  cnv_plot <- ggplot(data_for_plot, aes(x = coordinate, y = fold_change_adjusted)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    geom_line(linewidth = 2, colour = safe_red) +
    geom_point(data = primers_for_plot, pch = 21) +
    geom_vline(xintercept = gene_min, linetype = "dashed") +
    geom_vline(xintercept = gene_max, linetype = "dashed") +
    facet_wrap(~sample_id) +
    ylim(0, ymax) +
    scale_x_continuous(breaks = by_n(interval),
                       minor_breaks = NULL,
                       limits = c(plot_min, plot_max)) +
    labs(title = str_c(input_gene, " CNV results"),
         subtitle = str_c("Setting: ", input_setting, ". Dashed lines show gene coordinates"),
         x = str_c("GRCh38 Genomic coordinates (", interval/1000, " kb intervals)"))
  
  return(cnv_plot)
  
}

draw_repeat_plot <- function(df, input_sample, input_setting, input_gene, input_ymin = 0,
                             input_ymax = 70) {
  
  plot <- df |> 
    filter(setting == input_setting) |> 
    filter(gene == input_gene) |> 
    filter(sample_id == input_sample) |> 
    ggplot(aes(x = coordinate, y = fold_change_adjusted)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    geom_line(linewidth = 3,
              colour = safe_blue,
              alpha = 1) +
    ylim(input_ymin, input_ymax) +
    geom_vline(xintercept = erbb2_min, linetype = "dashed") +
    geom_vline(xintercept = erbb2_max, linetype = "dashed") +
    facet_wrap(~sample_id_worksheet, ncol = 1) +
    labs(title = str_c("Repeat testing for ", input_sample),
         subtitle = str_c("Gene: ", input_gene))
  
  return(plot)
  
}

draw_repeat_coord_plot <- function(df, input_sample, input_setting, input_gene,
                                   bin_size) {
  
  plot <- df |> 
    filter(setting == input_setting) |> 
    filter(gene == input_gene) |> 
    filter(sample_id == input_sample) |> 
    ggplot(aes(x = coordinate, y = sample_id_worksheet)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    geom_line(linewidth = 3,
              aes(colour = fold_change_adjusted)) +
    scale_colour_binned(breaks = seq(0, 70, by = bin_size)) +
    geom_vline(xintercept = erbb2_min, linetype = "dashed") +
    geom_vline(xintercept = erbb2_max, linetype = "dashed") +
    labs(title = "Repeat testing",
         subtitle = str_c("Gene: ", input_gene, "  Sample: ", input_sample),
         caption = str_c("Colour bin size: ", bin_size))
  
  return(plot)
  
}

# Data wrangling functions ----------------------------------------------------------

extract_cnv_calls <- function(df, input_gene) {
  
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

get_excel_names <- function(filepath, folder) {
  
  output <- list.files(str_c(filepath, folder), pattern = "*.xlsx")
  
  return(output)
  
}



get_gene_result <- function(df, gene_name) {
  
  if (ncol(df) == 0) {
    
    output <- "No call"
    
  }
  
  if (ncol(df) > 1) {
    
    x <- df |> 
      filter(gene == gene_name)
    
    output <- case_when(
      nrow(x) < 1 ~"No call",
      nrow(x) >= 1 ~"Amplification")
    
  }
  
  return(output)
  
}

summarise_results <- function(file, input_sheet) {
  
  results <- read_excel(path = file,
                        sheet = input_sheet) |> 
    janitor::clean_names()
  
  worksheet <- parse_filename(file, 1)
  
  sample_id <- parse_filename(file, 2)
  
  qualifier <- parse_filename(file, 3)
  
  sample_id_suffix <- str_c(sample_id, qualifier)
  
  patient_name <- parse_filename(file, 4)
  
  egfr_result <- get_gene_result(df = results, gene_name = "EGFR")
  
  erbb2_result <- get_gene_result(df = results, gene_name = "ERBB2")
  
  met_result <- get_gene_result(df = results, gene_name = "MET")
  
  summary <- tribble(
    
    ~gene, ~result,
    "EGFR", egfr_result,
    "ERBB2", erbb2_result,
    "MET", met_result
  ) |> 
    mutate(suffix = qualifier,
           sample = sample_id,
           sample_suffix = sample_id_suffix,
           name = patient_name,
           worksheet = worksheet)
  
  return(summary)
  
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
    janitor::clean_names()
  
  identifiers <- filename_to_df(file)
  
  output <- cbind(identifiers, x_wide)
  
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

draw_confusion_matrix <- function(input_gene) {
  
  x <- joined |> 
    filter(gene == input_gene)
  
  true_positives <- nrow(x[x$outcome == "true positive", ])
  
  true_negatives <- nrow(x[x$outcome == "true negative", ])
  
  false_positives <- nrow(x[x$outcome == "false positive", ])
  
  false_negatives <- nrow(x[x$outcome == "false negative", ])
  
  confusion_matrix <- tribble(
    ~"",      ~"PanSolid CLC +",      ~"PanSolid CLC -", 
    "Core+",  true_positives,         false_negatives,  
    "Core-",  false_positives,        true_negatives
  )
  
  return(confusion_matrix)
  
}

calculate_target_copies <- function(fold_change, tcc_percent) {
  
  ref_sample_copies <- 200
  
  tcc_fraction <- tcc_percent / 100
  
  sample_total_target_copies <- fold_change * ref_sample_copies
  
  sample_target_copies_in_tumour_cells <- sample_total_target_copies - ((100 * (1-tcc_fraction)) * 2)
  
  sample_target_copies_per_tumour_cell <- sample_target_copies_in_tumour_cells / (100 * tcc_fraction)
  
  return(sample_target_copies_per_tumour_cell)
  
}

# Database tables -------------------------------------------------------------------

extraction_method_key <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                              schema = "dbo",
                                              table = "MOL_ExtractionMethods")) |> 
  # Have to remove large columns to avoid Invalid Descriptor Index error
  select(-c(Checks, Reagents)) |> 
  collect()

extraction_tbl <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                  schema = "dbo",
                                                  table = "MOL_Extractions"))

extraction_batch_tbl <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                        schema = "dbo",
                                                        table = "MOL_ExtractionBatches"))


sample_tbl <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                              schema = "dbo",
                                              table = "Samples"))

results_tbl <- tbl(dbi_con, 
                   dbplyr::in_catalog(
                     catalog = "MolecularDB",
                     schema = "dbo",
                     table = "ResultsAccess"))

tissue_types <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                schema = "dbo",
                                                table = "TissueTypes")) |> 
  collect() |> 
  janitor::clean_names()

discode <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                schema = "dbo",
                                                table = "Discode")) |> 
  select(-c("Description", "ReferralDetails")) |> 
  collect() |> 
  janitor::clean_names()

ngiscodes <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                             schema = "dbo",
                                             table = "NGISCodes")) |> 
  collect() |> 
  janitor::clean_names()


dna_db_worksheets <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                   schema = "dbo",
                                                   table = "PCR_New")) 

# Database functions ----------------------------------------------------------------

get_columns <- function(table_input) {
  
  output <- odbc::odbcConnectionColumns(
    conn = dbi_con, 
    catalog_name = "MolecularDB",
    schema_name = "dbo",
    name = table_input)
  
  return(output)
  
}

get_extraction_method <- function(sample_vector) {
  
  extraction_tbl_samples <- extraction_tbl |> 
    filter(LabNo %in% sample_vector) |> 
    collect()
  
  batches <- unique(extraction_tbl_samples$ExtractionBatchFK)
  
  extraction_batch_info <- extraction_batch_tbl |> 
    filter(ExtractionBatchId %in% batches) |> 
    collect() |> 
    # Remove DNA dilutions
    filter(ExtractionMethodFK != 11) |>
    left_join(extraction_method_key, join_by(ExtractionMethodFK == ExtractionMethodId))

  output <- extraction_tbl_samples |> 
    left_join(extraction_batch_info, join_by(ExtractionBatchFK == ExtractionBatchId)) |> 
    filter(!is.na(MethodName)) |> 
    janitor::clean_names() |> 
    rename(labno = lab_no)
  
  return(output)
  
}

get_sample_tissue <- function(sample_vector) {
  
  output <- sample_tbl |> 
    select(-c(StatusComment, COMMENTS, ConsultantAddress, ADDRESS1)) |> 
    filter(LABNO %in% sample_vector) |> 
    collect() |> 
    janitor::clean_names() |> 
    mutate(tissue = as.numeric(tissue)) |> 
    left_join(tissue_types, join_by(tissue == tissue_type_id))
    
  return(output)
  
}


# ddPCR functions -------------------------------------------------------------------

read_biorad_csv <- function(worksheet) {
  
  output <- read_csv(here::here(str_c("data/", worksheet)), 
              col_types = cols(
                "Well" = "c",
                "ExptType" = "c",
                "Experiment" = "c",
                "Sample" = "c",
                "TargetType" = "c",
                "Target" = "c",
                "Status" = "c",
                "Concentration" = "d",
                "Supermix" = "c",
                "CopiesPer20uLWell" = "d",
                "TotalConfMax" = "d",
                "TotalConfMin" = "d",
                "PoissonConfMax" = "d",
                "PoissonConfMin" = "d",
                "Positives" = "i",
                "Negatives" = "i",
                "Ch1+Ch2+" = "i",
                "Ch1+Ch2-" = "i",
                "Ch1-Ch2+" = "i",
                "Ch1-Ch2-" = "i",
                "Linkage"  = "d",
                "AcceptedDroplets" = "i",
                "CNV" = "d",
                "TotalCNVMax" = "d",
                "TotalCNVMin" = "d",
                "PoissonCNVMax" = "d",
                "PoissonCNVMin" = "d",
                "FractionalAbundance" = "d",
                "TotalFractionalAbundanceMax" = "d",
                "TotalFractionalAbundanceMin" = "d",
                "PoissonFractionalAbundanceMax" = "d",
                "PoissonFractionalAbundanceMin" = "d",
                "ReferenceAssayNumber" = "d",
                "TargetAssayNumber" = "d",
                "Threshold" = "d",
                "MeanAmplitudeofPositives" = "d",
                "MeanAmplitudeofNegatives" = "d",
                "MeanAmplitudeTotal" = "d",
                "ExperimentComments" = "c",
                "MergedWells" = "c",
                "TotalConfMax68" = "d",
                "TotalConfMin68" = "d",
                "PoissonConfMax68" = "d",
                "PoissonConfMin68" = "d",
                "TotalCNVMax68" = "d",
                "TotalCNVMin68" = "d",
                "PoissonCNVMax68" = "d",
                "PoissonCNVMin68" = "d",
                "PoissonCNVMin68" = "d",
                "PoissonRatioMax68" = "d",
                "TotalRatioMin68" = "d",
                "TotalFractionalAbundanceMax68" = "d",
                "TotalFractionalAbundanceMin68" = "d",
                "PoissonFractionalAbundanceMax68" = "d",                
                "PoissonFractionalAbundanceMin68" = "d")) |> 
  janitor::clean_names() |> 
  mutate(sample_well = str_c(sample, "_", well),
         worksheet = worksheet,
         labno = sample)
  
  return(output)
  
}     

draw_ddpcr_cnv_plot <- function(worksheet_df, y_max = 10) {
  
  output <- ggplot(worksheet_df, aes(x = sample_well, y = cnv)) +
    geom_point(aes(colour = gene), size = 3) +
    geom_errorbar(aes(ymin = poisson_cnv_min, max = poisson_cnv_max)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_y_continuous(breaks = c(1:y_max), limits = c(0, y_max),
                       minor_breaks = c(1:y_max)) +
    labs(y = "Copy number", x = "") +
    geom_hline(yintercept = 2, linetype = "dashed")

  return(output)
  
}
