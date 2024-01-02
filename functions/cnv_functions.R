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

extract_cnv_calls <- function(df) {
  
  stopifnot("Genotype" %in% colnames(df))
  
  dq_regex <- regex(
        r"[
        (E.{2}B2|EGFR|MYC|MET|ARID1A|SUFU|)\s  # Gene names
        (amplification\sdetected|amplification)
        \s.                                 # Use . for bracket
        (Mean\sDQ|mean\sDQ|DQ)
        \s
        (\d{1,3}|\d{1,3}\.\d{2})            # Dosage quotient
        x
        ]",
        comments = TRUE
      )

  cnv_calls <- unique(grep(pattern = "DQ", x = df$Genotype, ignore.case = TRUE,
                           value = TRUE))
  
  no_cnv_calls <- unique(grep(pattern = "No evidence for a clinically relevant CNV", 
                              x = df$Genotype, ignore.case = TRUE,
                              value = TRUE))
  
  no_erbb2_calls <- unique(grep(pattern = "No ERBB2 amplification detected", 
                                x = df$Genotype, ignore.case = TRUE,
                                value = TRUE))
  
  output <- df |> 
    mutate(
      core_cnv_status = case_when(
        Genotype %in% cnv_calls ~"CNVs detected",
        Genotype %in% no_cnv_calls | Genotype %in% no_erbb2_calls ~"No CNVs detected",
        TRUE ~"Not specified"
      ),
      target_gene = str_extract(Genotype, dq_regex, group = 1), 
      target_gene = ifelse(target_gene == "ERRB2", "ERBB2", target_gene),
      target_gene_dq = as.numeric(str_extract(Genotype, dq_regex, group = 4)))
  
  return(output)
}

get_excel_names <- function(filepath, folder) {
  
  output <- list.files(str_c(filepath, folder), pattern = "*.xlsx")
  
  return(output)
  
}
