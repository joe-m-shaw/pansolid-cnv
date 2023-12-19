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

source("functions/cnv_functions.R")

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
               get_excel_names(cnv_path, f8))) |> 
  mutate(
    worksheet = str_extract(file, pattern = filename_regex, group = 1),
    sample = str_extract(file, pattern = filename_regex, group = 3),
    name = str_extract(file, pattern = filename_regex, group = 5))

summary <- all_files |> 
  filter(!is.na(sample)) |> 
  count(sample) |> 
  arrange(desc(n))

# How many samples have been run through CLC?
length(unique(summary$sample))

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
  filter(!duplicated(sample_id)) |> 
  extract_cnv_calls() |> 
  rename(core_panel_genotype = Genotype,
         core_panel_description = TEST) |> 
  arrange(target_gene)

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
  ggplot(aes(x = reorder(sample_id, target_gene_dq), y = target_gene_dq)) +
  geom_jitter() +
  facet_wrap(~target_gene) +
  ylim(0, 120) +
  theme(axis.text.x = element_blank())




  