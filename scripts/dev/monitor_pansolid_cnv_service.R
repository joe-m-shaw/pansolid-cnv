# Monitor PanSolid CNV service

# Connect to DLIMS
source(here::here("scripts/connect_to_dna_db.R"))

source(here::here("functions/pansolid_excel_functions.R"))

source(here::here("functions/pansolid_cnv_excel_functions.R"))

# Get a list of PanSolid worksheets ---------------------------------------

all_worksheets <- dna_db_worksheets |> 
  select(pcrid, date, description) |> 
  collect() |> 
  mutate(worksheet = paste0("WS", pcrid))

ps_ws_info <- all_worksheets |> 
  # New CNV Excel layout started with WS152758
  filter(pcrid >= 152758) |> 
  filter(grepl(pattern = "pansolid", 
               x = description,
               ignore.case = TRUE)) |> 
  mutate(ps_category = case_when(
    grepl(pattern = "jBRCA|j_BRCA|j-BRCA|jew",
      x = description,
      ignore.case = TRUE) ~"PanSolid Jewish BRCA",
    TRUE ~"PanSolid FFPE"
  ))

stopifnot(anyNA.data.frame(ps_ws_info) == FALSE)

ps_ws_ffpe_only <- ps_ws_info |> 
  filter(ps_category == "PanSolid FFPE")

ps_worksheets <- ps_ws_ffpe_only$worksheet

if(length(ps_worksheets) == 0) {
  stop("No PanSolid worksheets found on DNA Database")
}

# Find filepaths for PanSolid Excels --------------------------------------

ps_filepaths <- ps_worksheets |> 
  map(\(ps_worksheets) get_annotated_filepaths(worksheet = ps_worksheets)) |> 
  flatten()

stopifnot(length(ps_filepaths) != 0)

filename_regex <- "Annotated_.*_WS\\d{6}_\\d{8}.*.xlsx"

ps_filepath_df <- tibble(
  filepath = unlist(ps_filepaths)) |> 
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex))

stopifnot(anyNA.data.frame(ps_filepath_df) == FALSE)

# Load current collated data ----------------------------------------------

collated_data_folder <- paste0(config::get("data_folderpath"),
                               "live_service/collated/")

stdev_live <- read_csv(paste0(collated_data_folder,
                              "stdev_live.csv"))

percent_138x_live <- read_csv(paste0(collated_data_folder,
                                     "percent_138x_live.csv"))

pred_ncc_live <- read_csv(paste0(collated_data_folder,
                                 "pred_ncc_live.csv"))

sig_cnvs_live <- read_csv(paste0(collated_data_folder,
                                 "sig_cnvs_live.csv"))

amp_genes_live <- read_csv(paste0(collated_data_folder,
                                  "amp_genes_live.csv"))

del_genes_live <- read_csv(paste0(collated_data_folder,
                                  "del_genes_live.csv"))

loh_live <- read_csv(paste0(collated_data_folder,
                                  "loh_live.csv"))

# Identify files not already collated -------------------------------------

ps_new_filepath_df <- ps_filepath_df |> 
  filter(!filename %in% stdev_live$filename)

message(paste0(length(ps_new_filepath_df$filepath),
               " new files identified"))

# Prepare raw data folder -------------------------------------------------

raw_folder_path <- paste0(config::get("data_folderpath"),
                          "live_service/raw/")

if(length(list.files(raw_folder_path)) != 0){
  stop("Raw file folder is not empty")
}

# Copy new files to raw data folder ---------------------------------------

file.copy(from = ps_new_filepath_df$filepath,
          to = raw_folder_path)

new_filepaths <- list.files(path = raw_folder_path,
                            full.names = TRUE,
                            pattern = "Annotated.*.xlsx")

message(paste0(length(new_filepaths), " files copied into raw data folder"))

# Collate data from new files ---------------------------------------------

cnv_new <- new_filepaths |> 
  map(\(new_filepaths) extract_cnv_tbls(new_filepaths, 
                                        sheet_regex = "CNVs_"))

loh_new <- new_filepaths |> 
  map(\(new_filepaths) read_loh_table(filepath = new_filepaths)) |> 
  list_rbind() 

message(paste0(length(new_filepaths), " CNV and LOH results collated"))

stdev_new <- map(cnv_new, ~ .x[["stdev"]]) |> 
  list_rbind() |>
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex)) |> 
  relocate(filename, .after = filepath)

percent_138x_new <- map(cnv_new, ~ .x[["percent_138x"]]) |> 
  list_rbind()|>
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex)) |> 
  relocate(filename, .after = filepath)

pred_ncc_new <-  map(cnv_new, ~ .x[["pred_ncc"]]) |> 
  list_rbind()|>
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex)) |> 
  relocate(filename, .after = filepath)

sig_cnvs_new <- map(cnv_new, ~ .x[["sig_cnvs"]]) |> 
  list_rbind()|>
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex)) |> 
  relocate(filename, .after = filepath)

amp_genes_new <- map(cnv_new, ~ .x[["amp_genes"]]) |> 
  list_rbind()|>
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex)) |> 
  relocate(filename, .after = filepath)

del_genes_new <-  map(cnv_new, ~ .x[["del_genes"]]) |> 
  list_rbind()|>
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex)) |> 
  relocate(filename, .after = filepath)

# Perform checks ----------------------------------------------------------

# Signal-adjusted noise

stopifnot(anyNA.data.frame(stdev_new) == FALSE)

stopifnot(min(stdev_new$stdev_noise) >= 0)

# Percentage coverage at 138X

stopifnot(anyNA.data.frame(percent_138x_new) == FALSE)

stopifnot(min(percent_138x_new$percent_138x) >= 0)

stopifnot(max(percent_138x_new$percent_138x) <= 100)

# Predicted NCC

stopifnot(anyNA.data.frame(pred_ncc_new) == FALSE)

stopifnot(min(pred_ncc_new$pred_ncc) >= 20)

stopifnot(max(pred_ncc_new$pred_ncc) <= 100)

# Significant CNVs

stopifnot(anyNA.data.frame(sig_cnvs_new |> 
                             select(-c(check_1,	check_2,	
                                    copy_number, copy_number,
                                    start, end))) == FALSE)

# Amplification genes

stopifnot(anyNA.data.frame(amp_genes_new) == FALSE)

# Deletion genes

stopifnot(anyNA.data.frame(del_genes_new) == FALSE)

# LOH

stopifnot(anyNA.data.frame(loh_new |> 
                             select(-c(check_1, check_2))) == FALSE)

# Check all tables have been read from each file

stopifnot(length(setdiff(stdev_new$filepath, percent_138x_new$filepath)) == 0)

stopifnot(length(setdiff(stdev_new$filepath, pred_ncc_new$filepath)) == 0)

stopifnot(length(setdiff(stdev_new$filepath, unique(sig_cnvs_new$filepath))) == 0)

stopifnot(length(setdiff(stdev_new$filepath, unique(amp_genes_new$filepath))) == 0)

stopifnot(length(setdiff(stdev_new$filepath, unique(del_genes_new$filepath))) == 0)

stopifnot(length(setdiff(stdev_new$filepath, unique(loh_new$filepath))) == 0)

# Add new results to existing data ----------------------------------------

stdev_df <- rbind(stdev_live, stdev_new)

percent_138x_df <- rbind(percent_138x_live, percent_138x_new)

pred_ncc_df <- rbind(pred_ncc_live, pred_ncc_new)

sig_cnvs_df <- rbind(sig_cnvs_live, sig_cnvs_new)

amp_genes_df <- rbind(amp_genes_live, amp_genes_new)

del_genes_df <- rbind(del_genes_live, del_genes_new)

loh_df <- rbind(loh_live, loh_new)

# Export collated data ----------------------------------------------------

export_collated_data <- function(df, df_name) {
  
  filepath = paste0(config::get("data_folderpath"),
                    "live_service/collated/", 
                    df_name, 
                    "_live",
                    ".csv")
  
  write_csv(df, filepath)
  
}

df_list <- list(
  "stdev" = stdev_df,
  "percent_138x" = percent_138x_df,
  "pred_ncc" = pred_ncc_df,
  "sig_cnvs" = sig_cnvs_df,
  "amp_genes" = amp_genes_df,
  "del_genes" = del_genes_df,
  "loh" = loh_df
)

imap(df_list, export_collated_data)

# Delete copied files -----------------------------------------------------

message("Deleting new raw files")

file.remove(new_filepaths)

if(length(list.files(paste0(config::get("data_folderpath"),
                            "live_service/raw/"))) == 0){
  message("Raw file folder is empty")
}

# Clear environment -------------------------------------------------------

rm(list=ls())
