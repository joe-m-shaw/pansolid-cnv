# INC10994 Contamination Functions

source(here::here("scripts/connect_to_dna_db.R"))
source(here::here("functions/dna_db_functions.R"))

library(patchwork)

compare_snp_results <- function(s1_df,
                                s2_df,
                                s1_id,
                                s2_id,
                                hom_vaf_upper_threshold = 90,
                                het_vaf_upper_threshold = 60,
                                het_vaf_lower_threshold = 40,
                                hom_vaf_lower_threshold = 10) {
  
  sample1 <- s1_df |> 
    filter(labno_suffix_worksheet == s1_id &
             type == "SNV") |> 
    select(chromosome, region, type, reference, reference_allele, 
           allele, frequency, cumulative_region_coordinate) |> 
    rename(sample1_frequency = frequency)
  
  sample2 <- s2_df |> 
    filter(labno_suffix_worksheet == s2_id &
             type == "SNV") |> 
    select(chromosome, region, type, 
           reference, reference_allele,
           allele, frequency) |> 
    rename(sample2_frequency = frequency)
  
  comparison_df <- sample1 |> 
    left_join(sample2,
              join_by(chromosome, region, type, 
                      reference, reference_allele,
                      allele)) |> 
    dplyr::mutate(sample2_frequency_category = case_when(
      # Above 90
      sample2_frequency >= hom_vaf_upper_threshold ~paste0("Above ", hom_vaf_upper_threshold, "%"),
      # 40-60
      sample2_frequency <= het_vaf_upper_threshold & 
        sample2_frequency >= het_vaf_lower_threshold ~paste0("Between ",
                                                             het_vaf_lower_threshold,
                                                             "-",
                                                             het_vaf_upper_threshold,
                                                             "%"),
      # Below 10
      sample2_frequency <= 10 ~paste0("Lower than ",
                                      hom_vaf_lower_threshold,
                                      "%"),
      # In between
      sample2_frequency < hom_vaf_upper_threshold & sample2_frequency > het_vaf_upper_threshold ~"In between",
      sample2_frequency < het_vaf_lower_threshold & sample2_frequency > hom_vaf_lower_threshold ~"In between",
      # No frequency
      is.na(sample2_frequency) ~"No frequency"),
      # Factorise
      sample2_frequency_category = factor(sample2_frequency_category,
                                          levels = c(
                                            paste0("Above ", hom_vaf_upper_threshold, "%"),
                                            paste0("Between ",
                                                     het_vaf_lower_threshold,
                                                     "-",
                                                     het_vaf_upper_threshold,
                                                     "%"),
                                          paste0("Lower than ",
                                                 hom_vaf_lower_threshold,
                                                 "%"),
                                          "In between",
                                          "No frequency")))
  
  return(comparison_df)
  
}

read_formatted_snp_data <- function(file, sheetname = "Artefacts_removed_GIAB_filt...") {
  
  pansolid_chr_cumulative_coordinates <- readr::read_csv(paste0(config::get("data_folderpath"),
                                                                "validation/DOC6791_chromosome_arms/",
                                                                "bed_files/",
                                                                "pansolid_chr_cumulative_coordinates.csv"),
                                                         col_types = "ccdd")
  
  snp_data <- readxl::read_excel(path = file,
                                 sheet = sheetname,
                                 range = readxl::cell_cols("A:L"),
                                 col_types = c(
                                   "text", "numeric", "text", "text", "text", 
                                   "text", "numeric", "text", "text",
                                   "numeric", "numeric", "numeric"
                                 )) |> 
    janitor::clean_names() |> 
    format_chromosome_decimals(col = chromosome) |> 
    factorise_chromosome(col = chromosome_char) |> 
    filter(!is.na(region)) |> 
    left_join(pansolid_chr_cumulative_coordinates |> 
                select(chromosome_fct, cumulative_chr_coordinate),
              by = "chromosome_fct") |> 
    mutate(cumulative_region_coordinate = region + cumulative_chr_coordinate)
  
  output_df <- add_identifiers(file = file,
                               tbl = snp_data)
  
  return(output_df)
  
}

collate_snp_data <- function(worksheet) {
  
  files <- get_worksheet_filepaths(worksheet = worksheet,
                                   file_regex  = "Results_SNVs_Indels.*.xlsx")
  
  collated_data <- files |> 
    map(\(files) read_snp_sheet(file = files)) |> 
    list_rbind()
  
  formatted_data <- format_snp_sheet_data(collated_data)
  
  return(formatted_data)
  
}

check_hom_snps <- function(s1_df,
                           s2_df,
                           s1_id,
                           s2_id) {
  
  compdf <- compare_snp_results(s1_df = s1_df,
                                s2_df = s2_df,
                                s1_id = s1_id,
                                s2_id = s2_id) 
  
  high_snps <- compdf |> 
    filter(sample1_frequency > 98)
  
  above90 <- high_snps[high_snps$sample2_frequency_category == "Above 90%",]
  
  prop_above90 <- (nrow(above90) / nrow(high_snps))*100
  
  output_df <- tibble(
    "sample1" = c(s1_id),
    "sample2" = c(s2_id),
    "proportion_matching_hom_snps" = c(prop_above90),
    "number_matching_hom_snps" = c(nrow(above90))
  )
  
  return(output_df)
  
}

search_for_contaminant <- function(df, s1_id) {
  
  query_df <- df |> 
    filter(labno_suffix_worksheet != s1_id)
  
  query_ids <- unique(query_df$labno_suffix_worksheet)
  
  final_df <- query_ids |> 
    map(\(query_ids) check_hom_snps(s1_df = df,
                                    s1_id = s1_id,
                                    s2_df = df,
                                    s2_id = query_ids)) |> 
    list_rbind() |> 
    arrange(desc(proportion_matching_hom_snps))
  
  return(final_df)
  
  arr}

make_contamination_snp_plot <- function(s1_df,
                                        s1_id,
                                        s2_df,
                                        s2_id,
                                        hom_vaf_upper_threshold = 90,
                                        het_vaf_upper_threshold = 60,
                                        het_vaf_lower_threshold = 40,
                                        hom_vaf_lower_threshold = 10) {
  
  comparison_df <- compare_snp_results(s1_df = s1_df,
                                       s1_id = s1_id,
                                       s2_df = s2_df,
                                       s2_id = s2_id,
                                       hom_vaf_upper_threshold = hom_vaf_upper_threshold,
                                       het_vaf_upper_threshold = het_vaf_upper_threshold,
                                       het_vaf_lower_threshold = het_vaf_lower_threshold,
                                       hom_vaf_lower_threshold = hom_vaf_lower_threshold)
  
  comparison_df_no_na <- comparison_df |> 
    dplyr::filter(!is.na(sample2_frequency))
  
  s1_plot <- comparison_df |> 
    filter(reference_allele == "No") |> 
    ggplot2::ggplot(aes(x = cumulative_region_coordinate, y = sample1_frequency)) +
    geom_point(shape = 21, aes(fill = sample2_frequency_category),
               alpha = 0.8) +
    scale_fill_manual(values = c("#56B4E9", 
                                 "#D55E00", 
                                 "#33FFFF",
                                 "#FFFFFF",
                                 "#999999")) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(size=8),
          plot.subtitle = element_text(size=6)) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    labs(x = "SNPs", y = "SNP Allele Frequency",
         fill = "",
         title = paste0("SNP data from ", s1_id),
         subtitle = paste0("SNPs coloured by frequency in ", s2_id),
         caption = paste0(nrow(comparison_df_no_na), " SNPs match"))
  
  s2_plot <- comparison_df |> 
    filter(reference_allele == "No") |> 
    ggplot2::ggplot(aes(x = cumulative_region_coordinate, y = sample2_frequency)) +
    geom_point(shape = 21, aes(fill = sample2_frequency_category),
               alpha = 0.8) +
    scale_fill_manual(values = c("#56B4E9", 
                                 "#D55E00", 
                                 "#33FFFF",
                                 "#FFFFFF",
                                 "#999999")) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          plot.title = element_text(size=8)) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    labs(x = "SNPs", y = "SNP Allele Frequency",
         title = paste0("SNP data from ", s2_id))
  
  output <- s2_plot + s1_plot +
    plot_layout(ncol = 2)
  
  return(output)
  
}

# DNA Database Functions


define_pansolid_worksheets <- function() {
  
  all_worksheets <- dna_db_worksheets |> 
    select(pcrid, date, description) |> 
    collect() |> 
    mutate(worksheet = paste0("WS", pcrid))
  
  stopifnot(nrow(all_worksheets) > 0)
  
  ps_string_vars <- paste(c("pansolid",
                            "pan-solid", 
                            "pan_solid", 
                            "pan\\ssolid", 
                            "PnaSolid", 
                            "Pandolid", 
                            "PamSolid"), collapse = "|")
  
  ps_ws_info <- all_worksheets |> 
    filter(grepl(pattern = ps_string_vars, 
                 x = description,
                 ignore.case = TRUE)) |> 
    filter(!grepl(pattern = c("Limit of detection|cobas|ddpcr|confs|RNA"),
                  x = description,
                  ignore.case = TRUE)) |> 
    mutate(ps_category = case_when(
      grepl(pattern = "jBRCA|j_BRCA|j-BRCA|jew",
            x = description,
            ignore.case = TRUE) ~"PanSolid Jewish BRCA",
      TRUE ~"PanSolid FFPE"
    ))
  
  return(ps_ws_info)
  
}

define_pansolid_samples <- function() {
  
  ps_ws_info <- define_pansolid_worksheets()
  
  pansolid_worksheets <- unique(ps_ws_info$pcrid)
  
  all_ps_samples <- dna_db_pcr_records |> 
    filter(pcrid %in% pansolid_worksheets) |> 
    select(pcrid, sample) |> 
    collect()
  
  return(all_ps_samples)
  
}

define_plate_coordinates_96 <- function() {
  
  plate_coordinates_96 <- data.frame(
    "position" = seq(1, 96, by = 1),
    "y_coordinate" = rep(LETTERS[1:8], 12),
    "x_coordinate" = rep(1:12, each = 8)) |> 
    mutate(y_factor = factor(y_coordinate,
                             levels = rev(LETTERS[1:8])),
           x_string = case_when(
             str_length(x_coordinate) == 1 ~paste0("0",
                                                   as.character(x_coordinate)),
             str_length(x_coordinate) == 2 ~as.character(x_coordinate)
           ),
           coordinate_string = paste0(y_coordinate, x_string))
  
  return(plate_coordinates_96)
  
}

plate_position_to_well <- function(position_input) {
  
  stopifnot(position_input %in% seq(1, 96, by = 1))
  
  plate_coordinates_96 <- define_plate_coordinates_96()
  
  output <- plate_coordinates_96 |> 
    filter(position == position_input)
  
  return(output)
  
}

get_plate_location <- function(labno, worksheet) {
  
  worksheet_number <- parse_number(worksheet)
  
  location_df <- dna_db_pcr_records |> 
    filter(sample == labno &
             pcrid == worksheet_number) |> 
    select(sample, pcrid, position) |> 
    collect()
  
  plate_coordinates_96 <- define_plate_coordinates_96()
  
  output <- location_df |> 
    left_join(plate_coordinates_96, by = "position")
  
  return(output)
  
}


ps_samples_per_qs_batch <- function(qs_batch) {
  
  all_ps_samples <- define_pansolid_samples()
  
  qs_batch_df <- extraction_tbl |> 
    filter(extraction_batch_fk %in% qs_batch) |> 
    collect()
  
  tested_on_pansolid <- qs_batch_df |> 
    filter(lab_no %in% all_ps_samples$sample)
  
  qs_ps_samples <- all_ps_samples |> 
    filter(sample %in% tested_on_pansolid$lab_no)
  
  output <- data.frame(
    "qs_batch" = c(qs_batch),
    "qs_batch_total_samples" = c(nrow(qs_batch_df)),
    "qs_batch_samples_on_ps" = c(nrow(tested_on_pansolid)),
    "ps_worksheets" = c(paste(unique(qs_ps_samples$pcrid), collapse = ","))
  ) |> 
    mutate(percent_on_ps = round((qs_batch_samples_on_ps / qs_batch_total_samples)*100, 1))
  
  return(output)
  
}

make_contamination_summary_table <- function(labno, worksheet, input_text) {
  
  pansolid_df <- get_plate_location(labno = labno, worksheet = worksheet) |> 
    select(sample, pcrid, coordinate_string) |> 
    rename(labno = sample,
           ps_ws = pcrid,
           ps_pos = coordinate_string) |> 
    mutate(category = input_text) |> 
    relocate(category, .before = labno)
  
  extraction_df <- get_extraction_method(c(labno)) |>  
    select(extraction_batch_fk, sort_order.x) |> 
    rename(qs_batch = extraction_batch_fk,
           qs_pos = sort_order.x)
  
  output <- cbind(pansolid_df, extraction_df)
  
  return(output)
  
}
