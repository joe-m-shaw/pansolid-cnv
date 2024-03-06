# DNA Database Connection

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(odbc)
library(DBI)
library(dbplyr)

# DNA database connection -----------------------------------------------------------

dbi_con <- DBI::dbConnect(
  drv = odbc::odbc(),
  dsn = "moldb")

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

