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

