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
  collect() |> 
  janitor::clean_names()

extraction_tbl <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                  schema = "dbo",
                                                  table = "MOL_Extractions")) |> 
  janitor::clean_names()

extraction_batch_tbl <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                        schema = "dbo",
                                                        table = "MOL_ExtractionBatches")) |> 
  janitor::clean_names()

sample_tbl <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                              schema = "dbo",
                                              table = "Samples")) |> 
  janitor::clean_names()

results_tbl <- tbl(dbi_con, 
                   dbplyr::in_catalog(
                     catalog = "MolecularDB",
                     schema = "dbo",
                     table = "ResultsAccess")) |> 
  janitor::clean_names()

tissue_types <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                schema = "dbo",
                                                table = "TissueTypes")) |> 
  janitor::clean_names() |> 
  collect() 

discode <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                           schema = "dbo",
                                           table = "Discode")) |> 
  select(-c("Description", "ReferralDetails")) |> 
  janitor::clean_names() |> 
  collect() 
  

ngiscodes <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                             schema = "dbo",
                                             table = "NGISCodes")) |> 
  janitor::clean_names() |> 
  collect()
  

dna_db_worksheets <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                     schema = "dbo",
                                                     table = "PCR_New"))|> 
  janitor::clean_names()

