# How many samples have been tested so far?

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(odbc)
library(DBI)
library(dbplyr)
library(ggpubr)

# Database connection ---------------------------------------------------------------

dbi_con <- DBI::dbConnect(
  drv = odbc::odbc(),
  dsn = "moldb")

# Filepaths -------------------------------------------------------------------------

cnv_path <- "S:/central shared/Genetics/NGS/Bioinformatics/1_Pan-solid-Cancer/CNV/"

f1 <- "Added BioBank controls"

f2 <- "Adjusted_FC_and_additional_controls_100323"

f3 <- "cnv_analysis_R_outputs"

f4 <- "ERBB2_amplifications"

f5 <- "Misc_QC_reports"

f6 <- "No_single_exon_sensitivity" 

f7 <- "Oncogene_FC_3_120923"

f8 <- "Reduced_Onc_TSG_genes_220323"

# Functions -------------------------------------------------------------------------

get_excel_names <- function(filepath, folder) {
  
  output <- list.files(str_c(filepath, folder), pattern = "*.xlsx")
  
  return(output)
  
}

# Read data -------------------------------------------------------------------------

filename_regex <- regex(
  r"[
  (\d{8})     # Sample number
  ]"
)


all_files <- data.frame(
  
  "file" = c(get_excel_names(cnv_path, f1),
               get_excel_names(cnv_path, f2),
               get_excel_names(cnv_path, f3),
               get_excel_names(cnv_path, f4),
               get_excel_names(cnv_path, f5),
               get_excel_names(cnv_path, f6),
               get_excel_names(cnv_path, f7),
               get_excel_names(cnv_path, f8))) |> 
  mutate(sample = str_extract(file, pattern = "\\d{8}"))

summary <- all_files |> 
  filter(!is.na(sample)) |> 
  count(sample) |> 
  arrange(desc(n))

# How many samples have been run through CLC?
length(unique(summary$sample))

# Check sample types ----------------------------------------------------------------

sample_tbl <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                              schema = "dbo",
                                              table = "Samples"))

tissue_types <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                schema = "dbo",
                                table = "TissueTypes")) |> 
  collect()

disease_codes <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                schema = "dbo",
                                                table = "Discode")) |> 
  select(-c(Description, ReferralDetails)) |> 
  collect()

cnv_sample_ids <- summary$sample

cnv_sample_info <- sample_tbl |> 
  select(-c(StatusComment, COMMENTS, ConsultantAddress, ADDRESS1)) |>
  filter(LABNO %in% cnv_sample_ids) |> 
  collect() 

cnv_info_disease_codes <- cnv_sample_info |> 
  pivot_longer(cols = c("DISEASE", "DISEASE 2", "DISEASE 3", "DISEASE 4"),
               names_to = "original_column",
               values_to = "disease_code") |> 
  relocate(disease_code) |> 
  filter(!is.na(disease_code)) |> 
  left_join(disease_codes,
            join_by(disease_code == DISCODE)) |> 
  relocate(DISEASE)

disease_code_summary <- cnv_info_disease_codes |> 
  count(DISEASE)

disease_code_plot <- ggplot(disease_code_summary, aes(x = reorder(DISEASE, desc(n)), 
                                                      y = n)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "", title = "Referral Disease Code")

cnv_info_sample_types <- cnv_sample_info |>
  left_join(tissue_types |> 
              mutate(TissueTypeId = as.character(TissueTypeId)),
            join_by(TISSUE == "TissueTypeId")) |> 
  relocate(TISSUE, TissueType)

sample_types_summary <- cnv_info_sample_types |> 
  count(TissueType)

sample_type_plot <- ggplot(sample_types_summary, aes(x = reorder(TissueType, desc(n)), 
                                                     y = n)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "", title = "Sample Types")

ggarrange(disease_code_plot, sample_type_plot, ncol = 2, nrow = 1)


# Check results ---------------------------------------------------------------------

results_tbl <- tbl(dbi_con, 
                   dbplyr::in_catalog(
                     catalog = "MolecularDB",
                     schema = "dbo",
                     table = "ResultsAccess"))

all_results <- results_tbl |> 
  select(-c(GENOCOMM, QUALCOMM, WithdrawReason, EditReason)) |> 
  collect()

myc_results <- all_results |> 
  filter(Genotype %in% grep(pattern = "MYC", x = Genotype, value = TRUE))

result_info <- results_tbl |> 
  select(LABNO, TEST, TESTTYPE, Genotype, Genotype2, GENOCOMM) |> 
  filter(LABNO %in% cnv_sample_ids) |> 
  collect()

# New idea - get all Core results with positive results
testtype_19 <- results_tbl |> 
  select(LABNO, TEST, TESTTYPE, Genotype, Genotype2, GENOCOMM) |> 
  filter(TESTTYPE == 19) |> 
  collect()

dq_regex <- regex(
  r"[
  (ERBB2|EGFR|MYC|MET|ARID1A|SUFU)\s  # Gene names
  (amplification\sdetected|amplification)
  \s.                                 # Use . for bracket
  (Mean\sDQ|mean\sDQ|DQ)
  \s
  (\d{1,3}|\d{1,3}\.\d{2})            # Dosage quotient
  x
  ]",
  comments = TRUE
)

positive_results <- testtype_19 |> 
  filter(Genotype %in% grep(pattern = "DQ", x = Genotype, value = TRUE)) |> 
  filter(Genotype %in% grep(pattern = "ERBB2|EGFR|MYC|MET|ARID1A|SUFU",
                            x = Genotype,
                            value = TRUE)) |> 
  mutate(target_gene = case_when(
    
    Genotype %in% grep(pattern = "ERBB2",
                       x = Genotype,
                       value = TRUE) ~"ERBB2",
    
    Genotype %in% grep(pattern = "EGFR",
                       x = Genotype,
                       value = TRUE) ~"EGFR",
    
    Genotype %in% grep(pattern = "MYC",
                       x = Genotype,
                       value = TRUE) ~"MYC",
    
    Genotype %in% grep(pattern = "MET",
                       x = Genotype,
                       value = TRUE) ~"MET",
    
    Genotype %in% grep(pattern = "ARID1A",
                       x = Genotype,
                       value = TRUE) ~"ARID1A",
    
    Genotype %in% grep(pattern = "SUFU",
                       x = Genotype,
                       value = TRUE) ~"SUFU"),
    target_gene_dq = str_extract(Genotype, dq_regex, group = 4)) 

positive_results |> 
  count(target_gene)

positive_results |> 
  count(target_gene, DISEASE)

# Qiaseq Primers --------------------------------------------------------------------

pan_solid_only <- read_excel("data/Primer and Gene Comparison.xlsx",
                            sheet = "Primer Overlap",
                            range = "A2:D11024",
                            col_names = c("chromosome", "coordinates", "sequence",
                                          "gene"),
                            col_types = c("text", "text", "text", "text")) |> 
  mutate(category = "Unique to 44038Z-11379",
         text = "PanSolid only")

core_only <- read_excel("data/Primer and Gene Comparison.xlsx",
                             sheet = "Primer Overlap",
                             range = "G2:J595",
                             col_names = c("chromosome", "coordinates", "sequence",
                                           "gene"),
                             col_types = c("text", "text", "text", "text")) |> 
  mutate(category = "Unique to 17500Z-950",
         text = "Core only")

intersect(pan_solid_only$sequence, core_only$sequence)

both_panels <- read_excel("data/Primer and Gene Comparison.xlsx",
                          sheet = "Primer Overlap",
                          range = "L2:O357",
                          col_names = c("chromosome", "coordinates", "sequence",
                                        "gene"),
                          col_types = c("text", "text", "text", "text")) |> 
  mutate(category = "common to both 44038Z-11379 and 17500Z-950",
         text = "Both")

all_primers <- rbind(pan_solid_only, core_only, both_panels)

all_primers |> 
  filter(gene %in% c("ERBB2", "ERBB2, MIR4728",
                     "EGFR", "EGFR, EGFR-AS1",
                     "MYC", "MET", "ARID1A", "SUFU")) |> 
  mutate(gene_clean = case_when(
    
    gene == "ERBB2, MIR4728" ~"ERBB2",
    gene == "EGFR, EGFR-AS1" ~"EGFR",
    TRUE ~gene)) |> 
  ggplot(aes(x = gene_clean, y = , fill = text)) +
  geom_bar() +
  theme(legend.title = element_blank()) +
  facet_wrap(~text) +
  labs(x = "", y = "Number of primers", 
       title = "EGFR, ERBB2 and MET share primers between Pan Solid and Core")
