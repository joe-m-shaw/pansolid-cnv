# Pan Solid CNV Functions

# CSV timestamp ---------------------------------------------------------------------

export_timestamp <- function(input) {
  write.csv(input,
            file = here::here(paste0(
              "outputs/",
              format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
              "_",
              deparse(substitute(input)), ".csv"
            )),
            row.names = FALSE
  )
}

# Plot functions --------------------------------------------------------------------

safe_blue <- "#88CCEE"
safe_red <- "#CC6677"
safe_grey <- "#888888"

save_plot <- function(input_plot, input_width = 15, input_height = 12, dpi = 300) {
  
  # Default inputs allow for presenting a plot as half an A4 page
  
  ggsave(
    filename = paste0(
      format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
      "_",
      deparse(substitute(input_plot)), ".png"
    ),
    plot = input_plot,
    device = "png",
    path = "plots/",
    units = "cm",
    width = input_width,
    height = input_height,
    dpi = 300
  )
}

plot_coarse_v_fine <- function(input_gene) {
  
  plot <- all_calls |> 
    filter(name == input_gene) |> 
    group_by(sample, setting) |> 
    summarise(calls = n()) |> 
    ggplot(aes(x = sample, y = calls)) +
    geom_col(aes(fill = setting), position = "dodge") +
    theme_bw() +
    labs(title = as.character(input_gene))
  
  return(plot)
  
}

draw_cnv_plot <- function(df, input_gene, input_setting) {
  
  stopifnot(input_gene %in% c("EGFR", "ERBB2", "MET"))
  
  if(input_gene == "EGFR") {
    
    gene_min <- egfr_min
    gene_max <- egfr_max
    
  }
  
  if(input_gene == "ERBB2") {
    
    gene_min <- erbb2_min
    gene_max <- erbb2_max
    
  }
  
  if(input_gene == "MET") {
    
    gene_min <- met_min
    gene_max <- met_max
    
  }
  
  cnv_plot <- df |> 
    filter(setting == input_setting) |> 
    filter(gene == input_gene) |> 
    ggplot(aes(x = coordinate, y = fold_change_adjusted)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    geom_line(linewidth = 2, colour = safe_red) +
    geom_vline(xintercept = gene_min, linetype = "dashed") +
    geom_vline(xintercept = gene_max, linetype = "dashed") +
    facet_wrap(~sample) +
    ylim(0, 70) +
    labs(title = str_c(input_gene, " CNV results"),
         subtitle = str_c("Setting: ", input_setting, ". Dashed lines show gene coordinates"))
  
  return(cnv_plot)
  
}

# Data wrangling functions ----------------------------------------------------------

extract_cnv_calls <- function(df, input_gene) {
  
  stopifnot("Genotype" %in% colnames(df))
  
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
    select(LABNO, Genotype) |> 
    mutate(gene_searched = input_gene,
           gene_match = str_extract(Genotype, dq_regex, group = 1),
           gene_dq = as.numeric(str_extract(Genotype, dq_regex, group = 4)),
           core_result = ifelse(!is.na(gene_dq), "Amplification", "No call"))
  
  return(output)
  
}

get_excel_names <- function(filepath, folder) {
  
  output <- list.files(str_c(filepath, folder), pattern = "*.xlsx")
  
  return(output)
  
}

filename_regex <- regex(
  r"[
  (WS\d{6})             # Worksheet number
  _
  (\d{8})               # Sample number
  (a_|b_|c_|_)
  ([A-z]+)              # Patient name
  (.xlsx|_.+)
  ]",
  comments = TRUE
)

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
  
  worksheet <- str_extract(file, filename_regex,
                          group = 1)
  
  sample_id <- str_extract(file, filename_regex,
                           group = 2)
  
  qualifier <- str_extract(file, filename_regex,
                           group = 3)
  
  sample_id_suffix <- str_c(sample_id, qualifier)
  
  patient_name <- str_extract(file, filename_regex,
                              group = 4)
  
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
           name = patient_name)
  
  return(summary)
  
}

read_clc_excel <- function(file, input_sheet) {
  
  sample_id <- str_extract(file, filename_regex,
                           group = 2)
  
  patient_name <- str_extract(file, filename_regex,
                              group = 4)
  
  results <- read_excel(path = file,
                        sheet = input_sheet) |> 
    janitor::clean_names() |> 
    mutate(sample = sample_id,
           patient = patient_name)
  
  return(results)
  
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

extract_cnv_coordinates <- function(df) {
  
  stopifnot("cnv_region" %in% colnames(df))
  
  cnv_coord_regex <- regex(
      r"[
        \D{0,11}       # Either 0 or 11 non-digit characters. 11 is "complement("
        (\d{1,10})     # first coordinate number (1 to 10 digits)
        \.\.           # two full stops
        (\d{1,10})     # second coordinate number (1 to 10 digits)
        ]",
        comments = TRUE
        )
  
  output <- df |> 
    mutate(cnv_start = as.numeric(str_extract(string = cnv_region, 
                                              pattern = cnv_coord_regex, 
                                              group = 1)),
           cnv_end = as.numeric(str_extract(string = cnv_region, 
                                            pattern = cnv_coord_regex, 
                                            group = 2))) 
  
  return(output)
  
}



