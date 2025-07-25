# Monitor PanSolid CNV service

# Connect to DLIMS
source(here::here("scripts/connect_to_dna_db.R"))

source(here::here("functions/pansolid_cnv_excel_functions.R"))

# Get a list of PanSolid worksheets ---------------------------------------

message("Finding list of all PanSolid worksheets")

all_worksheets <- dna_db_worksheets |> 
  select(pcrid, date, description) |> 
  collect() |> 
  mutate(worksheet = paste0("WS", pcrid))

stopifnot(nrow(all_worksheets) > 0)

ps_ws_info <- all_worksheets |> 
  # New CNV Excel layout started with WS152758
  filter(pcrid >= 152758) |> 
  filter(grepl(pattern = "pansolid|pan-solid|pan_solid|pan solid", 
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

message("Finding filepaths for all PanSolid Excels")

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

# Identify files without CNV tabs -----------------------------------------

# Some samples have such low coverage that a CNV tab is not included in 
# the output Excel. These files will cause an error message at the data 
# collation stage if they are included.

files_without_cnv_tabs <- c("WS152872_25022765",
                            "WS153248_25026440",
                            "WS153962_25031207",
                            "WS154357_25030473",
                            "WS154624_25035774")

# Identify files not already collated -------------------------------------

message("Identifying new files for collation")

ps_new_filepath_df <- ps_filepath_df |> 
  mutate(worksheet = parse_filename(filename, 1),
         labno = parse_filename(filename, 2),
         worksheet_labno = str_c(worksheet, "_", labno)) |> 
  # Filename (not filepath) is used to check if a file has been collated
  # before. This is to avoid duplication of data when files have been 
  # moved into different folders.
  filter(!filename %in% stdev_live$filename) |> 
  filter(!worksheet_labno %in% files_without_cnv_tabs)

if(length(ps_new_filepath_df) > 0) {
  message(paste0(length(ps_new_filepath_df$filepath),
               " new files identified"))
} else {
  stop("No new files identified")
}

# Prepare raw data folder -------------------------------------------------

raw_folder_path <- paste0(config::get("data_folderpath"),
                          "live_service/raw/")

if(length(list.files(raw_folder_path)) != 0){
  stop("Raw file folder is not empty")
} else {
  message("Raw data folder is empty")
}

# Copy new files to raw data folder ---------------------------------------

# Files have to be copied prior to data collation, as read_excel() cannot be
# used when someone has the file open
file.copy(from = ps_new_filepath_df$filepath,
          to = raw_folder_path)

new_filepaths <- list.files(path = raw_folder_path,
                            full.names = TRUE,
                            pattern = "Annotated.*.xlsx")

message(paste0(length(new_filepaths), " files copied into raw data folder"))

# Collate data from new files ---------------------------------------------

message("Collating new data files")

cnv_new <- new_filepaths |> 
  map(\(new_filepaths) extract_cnv_tbls(new_filepaths, 
                                        sheet_regex = "CNVs_"))

loh_new <- new_filepaths |> 
  map(\(new_filepaths) read_loh_table(filepath = new_filepaths)) |> 
  list_rbind() |> 
  select(worksheet,	labno,	suffix,	patient_name,	labno_suffix,
         labno_suffix_worksheet,	filepath,	chrom,	gene,
         ploidy_state,	loh_status,	no_targets_in_ploidy_region,
         check_1,	check_2)

message(paste0(length(new_filepaths), " CNV and LOH results collated"))

# Split into separate dataframes ------------------------------------------

# Filename column added to allow later rbind step with previously 
# collated data.

stdev_new <- map(cnv_new, ~ .x[["stdev"]]) |> 
  list_rbind() |>
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex)) |> 
  relocate(filename, .after = filepath) |> 
  # Specify columns to remove any random columns added by scientists
  select(worksheet, labno,	suffix, patient_name, labno_suffix, 
         labno_suffix_worksheet,	filepath,	filename,	stdev_noise)

percent_138x_new <- map(cnv_new, ~ .x[["percent_138x"]]) |> 
  list_rbind()|>
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex)) |> 
  relocate(filename, .after = filepath) |> 
  select(worksheet,	labno,	suffix,	patient_name,	labno_suffix,
         labno_suffix_worksheet,	filepath,	filename,	percent_138x)

pred_ncc_new <-  map(cnv_new, ~ .x[["pred_ncc"]]) |> 
  list_rbind()|>
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex)) |> 
  relocate(filename, .after = filepath) |> 
  select(worksheet,	labno,	suffix,	patient_name,	labno_suffix,
         labno_suffix_worksheet,	filepath,	filename,	pred_ncc)

sig_cnvs_new <- map(cnv_new, ~ .x[["sig_cnvs"]]) |> 
  list_rbind()|>
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex)) |> 
  relocate(filename, .after = filepath) |> 
  select(worksheet,	labno,	suffix,	patient_name,	labno_suffix,
         labno_suffix_worksheet,	filepath,	filename,	gene,
         chromosome,	cnv_co_ordinates,	cnv_length,	consequence,
         fold_change,	p_value,	no_targets,	check_1,	check_2,
         copy_number,	start,	end)

amp_genes_new <- map(cnv_new, ~ .x[["amp_genes"]]) |> 
  list_rbind()|>
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex)) |> 
  relocate(filename, .after = filepath) |> 
  select(worksheet,	labno,	suffix,	patient_name,	labno_suffix,
         labno_suffix_worksheet,	filepath,	filename,	gene,
         max_region_fold_change,	min_region_fold_change)

del_genes_new <-  map(cnv_new, ~ .x[["del_genes"]]) |> 
  list_rbind()|>
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex)) |> 
  relocate(filename, .after = filepath) |> 
  select(worksheet,	labno,	suffix,	patient_name,	labno_suffix,
         labno_suffix_worksheet,	filepath,	filename,	gene,
         max_region_fold_change,	min_region_fold_change)

# Perform checks ----------------------------------------------------------

message("Performing checks")

# Signal-adjusted noise

stopifnot(anyNA.data.frame(stdev_new) == FALSE)

stopifnot(min(stdev_new$stdev_noise) >= 0)

stopifnot(typeof(stdev_new$stdev_noise) == "double")

# Percentage coverage at 138X

stopifnot(anyNA.data.frame(percent_138x_new) == FALSE)

stopifnot(min(percent_138x_new$percent_138x) >= 0)

stopifnot(max(percent_138x_new$percent_138x) <= 100)

stopifnot(typeof(percent_138x_new$percent_138x) == "double")

# Predicted NCC

stopifnot(anyNA.data.frame(pred_ncc_new) == FALSE)

stopifnot(min(pred_ncc_new$pred_ncc) >= 20)

stopifnot(max(pred_ncc_new$pred_ncc) <= 100)

stopifnot(typeof(pred_ncc_new$pred_ncc) == "double")

# Significant CNVs

stopifnot(anyNA.data.frame(sig_cnvs_new |> 
                             select(-c(check_1,	check_2,	
                                    copy_number, copy_number,
                                    start, end))) == FALSE)

# Amplification genes

stopifnot(anyNA.data.frame(amp_genes_new) == FALSE)

stopifnot(typeof(amp_genes_new$max_region_fold_change) == "double" &
  typeof(amp_genes_new$min_region_fold_change) == "double")
 
# Deletion genes

stopifnot(anyNA.data.frame(del_genes_new) == FALSE)

stopifnot(typeof(del_genes_new$max_region_fold_change) == "double" &
            typeof(del_genes_new$min_region_fold_change) == "double")

# LOH

stopifnot(anyNA.data.frame(loh_new |> 
                             select(-c(check_1, check_2))) == FALSE)

# Check all tables have been read from each file

stopifnot(anyDuplicated(stdev_new$filepath) == FALSE)

stopifnot(length(setdiff(stdev_new$filepath, percent_138x_new$filepath)) == 0)

stopifnot(length(setdiff(stdev_new$filepath, pred_ncc_new$filepath)) == 0)

stopifnot(length(setdiff(stdev_new$filepath, unique(sig_cnvs_new$filepath))) == 0)

stopifnot(length(setdiff(stdev_new$filepath, unique(amp_genes_new$filepath))) == 0)

stopifnot(length(setdiff(stdev_new$filepath, unique(del_genes_new$filepath))) == 0)

stopifnot(length(setdiff(stdev_new$filepath, unique(loh_new$filepath))) == 0)

message("Checks passed")

# Add new results to existing data ----------------------------------------

message("Adding new results to existing data")

stdev_df <- rbind(stdev_live, stdev_new)

percent_138x_df <- rbind(percent_138x_live, percent_138x_new)

pred_ncc_df <- rbind(pred_ncc_live, pred_ncc_new)

sig_cnvs_df <- rbind(sig_cnvs_live, sig_cnvs_new)

amp_genes_df <- rbind(amp_genes_live, amp_genes_new)

del_genes_df <- rbind(del_genes_live, del_genes_new)

loh_df <- rbind(loh_live, loh_new)

# Export collated data ----------------------------------------------------

message("Exporting collated data")

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

if(length(list.files(raw_folder_path)) == 0){
  message("Raw file folder is empty")
}

# Clear environment -------------------------------------------------------

rm(list=ls())
