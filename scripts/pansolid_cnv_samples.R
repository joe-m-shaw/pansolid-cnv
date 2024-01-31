# How many samples have been tested so far?

rm(list=ls())

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(odbc)
library(DBI)
library(dbplyr)
library(ggpubr)
library(here)

source(here::here("functions/cnv_functions.R"))

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

f9 <- "00_Amplifications_Fine_vs_Coarse"

# Read data -------------------------------------------------------------------------

filename_regex <- regex(
  r"[
  (WS\d{6})             # Worksheet number
  (_|Reload_|reload_)
  (\d{8})               # Sample number
  (a|b|c|Reload|)       # Variable sample descriptors
  _
  ([A-z]+)              # Patient name
  _
  ]",
  comments = TRUE
)

all_files <- data.frame(
  
  "file" = c(get_excel_names(cnv_path, f1),
             get_excel_names(cnv_path, f2),
             get_excel_names(cnv_path, f3),
             get_excel_names(cnv_path, f4),
             get_excel_names(cnv_path, f5),
             get_excel_names(cnv_path, f6),
             get_excel_names(cnv_path, f7),
             get_excel_names(cnv_path, f8),
             get_excel_names(cnv_path, f9))) |> 
  mutate(
    worksheet = str_extract(file, pattern = filename_regex, group = 1),
    sample = str_extract(file, pattern = filename_regex, group = 3),
    suffix = str_extract(file, pattern = filename_regex, group = 4),
    name = str_extract(file, pattern = filename_regex, group = 5),
    worsheet_sample_id = str_c(worksheet, "_", sample, suffix))

summary <- all_files |> 
  filter(!is.na(sample)) |> 
  count(sample) |> 
  arrange(desc(n))

# How many samples have been run through CLC?
length(unique(summary$sample))

# Repeated samples ------------------------------------------------------------------

repeats <- all_files |> 
  filter(!duplicated(worsheet_sample_id)) |> 
  filter(duplicated(sample, fromLast = TRUE) | duplicated(sample, fromLast = FALSE)) |> 
  arrange(sample)

# 8 samples for intra-run variability
intra_run_samples <- repeats |> 
  filter(suffix %in% c("a", "b", "c"))

inter_run_samples <- repeats |> 
  filter(!suffix %in% c("b", "c"))

# 7 samples for inter-run variability (17 data points)
inter_run_summary <- inter_run_samples |> 
  group_by(name) |> 
  count() |> 
  filter(n > 1)

# Core results of samples already tested on PanSolid --------------------------------

results_tbl <- tbl(dbi_con, 
                   dbplyr::in_catalog(
                     catalog = "MolecularDB",
                     schema = "dbo",
                     table = "ResultsAccess"))

cnv_sample_ids <- c(summary$sample)

result_info <- results_tbl |> 
  select(LABNO, TEST, TESTTYPE, Genotype, Genotype2, GENOCOMM) |> 
  filter(LABNO %in% cnv_sample_ids) |> 
  collect()

# All samples with Core results -----------------------------------------------------

all_results <- results_tbl |> 
  select(LABNO, TEST, TESTTYPE, Genotype, Genotype2, GENOCOMM) |> 
  collect()

core_results <- all_results |> 
  filter(TEST %in% grep(pattern = "Q.{2,4}seq.+core", x = TEST, ignore.case = TRUE,
                        value = TRUE))

pansolid_results <- all_results |> 
  filter(TEST %in% grep(pattern = "pansolid", x = TEST, ignore.case = TRUE,
                        value = TRUE))


# All samples tested on PanSolid ----------------------------------------------------

standard_columns <- c("sample_id", "sample_name", "panel", 
                      "enrichment", "qi_aseq_worksheet")

pansolid_samples_2022 <- read_excel(path = "data/QIAseq DNA PanSolid Sample Submission 2022.xlsx",
                                    sheet = "PanSolid Samples") |> 
  janitor::clean_names() |> 
  fill(qi_aseq_worksheet)

pansolid_samples_2023 <- read_excel(path = "data/DNA PanSolid QIAseq Submission Sheet 2023.xlsx",
                               sheet = "QIAseq samples") |> 
  janitor::clean_names()

pansolid_samples <- pansolid_samples_2022 |> 
  select(all_of(standard_columns)) |> 
  rbind(pansolid_samples_2023 |> 
          select(all_of(standard_columns))) |> 
  filter(enrichment == "PANSOLID") |> 
  filter(!duplicated(sample_id))

# Samples tested on PanSolid with Core results --------------------------------------

# Method 1: use PanSolid Excel spreadsheet

method1_samples <- pansolid_samples |> 
  inner_join(core_results, join_by(sample_id == LABNO))

length(unique(method1_samples$sample_id))

# Method 2: use outputs from CLC CNV pipeline

method2_samples <- all_files |> 
  filter(!duplicated(sample)) |> 
  inner_join(core_results, join_by(sample == LABNO))

all_ids <- c(method2_samples$sample, method1_samples$sample_id)

all_unique_ids <- length(unique(all_ids))

both_lists <- length(intersect(method1_samples$sample_id, method2_samples$sample))

method1_only <- length(setdiff(method1_samples$sample_id, method2_samples$sample))

method2_only <- length(setdiff(method2_samples$sample, method1_samples$sample_id))

# Sense check
all_unique_ids == both_lists + method1_only + method2_only

all_unique_ids
both_lists
method1_only
method2_only

# Join tables together --------------------------------------------------------------

x <- method1_samples |> 
    select(-c(TESTTYPE, Genotype2, GENOCOMM)) |> 
    rename(pansolid_worksheet = qi_aseq_worksheet) |> 
    select(sample_id, sample_name, panel, enrichment, 
           pansolid_worksheet, TEST, Genotype)

y <- method2_samples |> 
    select(-file, TESTTYPE, Genotype2, GENOCOMM) |> 
    rename(sample_id = sample,
           pansolid_worksheet = worksheet, 
           sample_name = name) |> 
    mutate(panel = "",
           enrichment = "") |> 
    select(sample_id, sample_name, panel, enrichment, 
           pansolid_worksheet, TEST, Genotype)
         
samples_pansolid_core <- rbind(x, y) |> 
  #filter(!duplicated(sample_id)) |> 
  #extract_cnv_calls() |> 
  rename(core_panel_genotype = Genotype,
         core_panel_description = TEST) 


repeated_samples <- samples_pansolid_core |> 
  filter(duplicated(sample_id, fromLast = TRUE) | 
           duplicated(sample_id, fromLast = FALSE))


export_timestamp(samples_pansolid_core)

samples_pansolid_core_filtered <- samples_pansolid_core |>
  filter(panel != "") |> 
  filter((core_cnv_status == "CNVs detected" & !is.na(target_gene)) |
           core_cnv_status == "No CNVs detected" & is.na(target_gene))

export_timestamp(samples_pansolid_core_filtered)

samples_pansolid_core_filtered |> 
  count(core_cnv_status)

samples_pansolid_core_filtered |> 
  count(target_gene)

samples_pansolid_core |> 
  filter(!is.na(target_gene)) |> 
  ggplot(aes(x = reorder(sample_id, target_gene_dq), y = target_gene_dq,
             colour = target_gene)) +
  geom_point(size = 2) +
  ylim(0, 120) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept = 5, linetype = "dashed")


# Samples from Dosage Quotient validation -------------------------------------------

dq_samples <- read_excel(path = here::here("data/dq_validation_samples.xlsx")) |> 
  janitor::clean_names() |> 
  mutate(lab_no = as.character(lab_no))

sample_tbl <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                              schema = "dbo",
                                              table = "Samples")) 
  
dq_sample_list <- dq_samples$lab_no

dlms_info <- sample_tbl |> 
  filter(LABNO %in% dq_sample_list) |> 
  select(LABNO, FIRSTNAME, SURNAME) |> 
  collect()

dq_list_with_names <- dq_samples |> 
  left_join(dlms_info, join_by(lab_no == LABNO))

for_export <- dq_list_with_names |> 
  arrange(lab_no) |> 
  mutate("DNA volume" = "",
         name = str_c(FIRSTNAME, " ", SURNAME)) |> 
  select(lab_no, name, "DNA volume")

export_timestamp(for_export)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
dna_volumes <- read_csv(file = here::here("data/2024_01_05_09_06_52_dna_volumes.csv")) |> 
  filter(!duplicated(lab_no)) |> 
  mutate(lab_no = as.character(lab_no))

dna_concentrations <- sample_tbl |> 
  filter(LABNO %in% dq_sample_list) |> 
  select(LABNO, CONCENTRATION) |> 
  collect() |> 
  mutate(CONCENTRATION = as.numeric(CONCENTRATION))

all_sample_info <- dq_list_with_names |> 
  left_join(dna_volumes, by = "lab_no") |> 
  left_join(dna_concentrations, join_by(lab_no == LABNO)) |> 
  janitor::clean_names() |> 
  mutate(dna_available = dna_volume * concentration,
         enough_for_pansolid = ifelse(dna_available >= 100, "Yes", "No")) |> 
  relocate(enough_for_pansolid, .after = orthogonal_result_methodology)


# MET amplifications ----------------------------------------------------------------

all_results <- results_tbl |> 
  select(LABNO, TEST, TESTTYPE, Genotype, Genotype2, GENOCOMM) |> 
  collect()

get_amp_calls <- function(df, gene) {
  
  output <- extract_cnv_calls(df = df, input_gene = gene) |> 
    filter(!is.na(gene_match)) |> 
    filter(!duplicated(LABNO))
  
  return(output)
  
}

met_calls <- get_amp_calls(all_results, "MET")

erbb2_calls <- get_amp_calls(all_results, "ERBB2")

egfr_calls <- get_amp_calls(all_results, "EGFR")

plot_gene_results <- function(df, gene_title) {
  
  ggplot(df, aes(reorder(LABNO, gene_dq), y = gene_dq)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    labs(x = "Sample", y = "Dosage Quotient",
         title = str_c(gene_title, ": ", nrow(df), " samples"))
  
}

met_p <- plot_gene_results(met_calls, "MET")

erbb2_p <- plot_gene_results(erbb2_calls, "ERBB2")

egfr_p <- plot_gene_results(egfr_calls, "EGFR")

ggarrange(met_p, erbb2_p, egfr_p, nrow = 2, ncol = 2)






# MDM2 amplifications ---------------------------------------------------------------

mdm2_samples <- c(23045472, 23055487, 23041652)

mdm2_details <- get_extraction_method(mdm2_samples) |> 
  janitor::clean_names()

mdm2_sample_details <- sample_tbl |> 
  select(LABNO, FIRSTNAME, SURNAME, iGeneRNo, iGeneSNo) |> 
  filter(LABNO %in% mdm2_samples) |> 
  collect() |> 
  janitor::clean_names()

all_mdm2_details <- mdm2_details |> 
  left_join(mdm2_sample_details, join_by ("lab_no" == "labno")) |> 
  select(lab_no, extraction_batch_fk, run_date, firstname, 
         surname, i_gene_r_no, i_gene_s_no)

write.csv(x = all_mdm2_details, file = here::here("outputs/mdm2_samples.csv"),
          row.names = FALSE)





samples_pansolid_core |> 
  filter(sample_id %in% repeats$sample) |>  view()



