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
  filter(pcrid == 152758) |> 
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

ps_filepath_df <- tibble(
  filepath = unlist(ps_filepaths)) |> 
  mutate(filename = str_extract(string = filepath,
                                pattern = "Annotated_.*_WS\\d{6}_\\d{8}.*.xlsx"))

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

# Identify files not already collated -------------------------------------


ps_filepath_df |> 
  filter(!filename %in% stdev_live$filename)


# Stop if there are already files in the raw folder


# Check whether the filenames have already been collated

# Copy any new results into a local folder

# Read in the results of those new files and collate them

# CNV tab, LOH, SVs

# Perform checks to make sure the new data is the right format

# Add those results to the collated data





# Copy new files into raw folder ------------------------------------------

file.copy(from = ps_filepath_df$filepath,
          to = paste0(config::get("data_folderpath"),
                      "live_service/raw/"))

new_filepaths <- list.files(path = paste0(config::get("data_folderpath"),
                                          "live_service/raw/"),
                            full.names = TRUE,
                            pattern = "Annotated.*.xlsx")


time_start <- Sys.time()

new_data <- new_filepaths |> 
  map(\(new_filepaths) extract_cnv_tbls(new_filepaths, 
                                         sheet_regex = "CNVs_"))
time_end <- Sys.time()

message(paste0(length(new_filepaths), " CNV results collated"))

stdev_df <- map(new_data, ~ .x[["stdev"]]) |> 
  list_rbind()

percent_138x_df <- map(new_data, ~ .x[["percent_138x"]]) |> 
  list_rbind()

pred_ncc_df <-  map(new_data, ~ .x[["pred_ncc"]]) |> 
  list_rbind()

sig_cnvs_df <- map(new_data, ~ .x[["sig_cnvs"]]) |> 
  list_rbind()

amp_genes_df <- map(new_data, ~ .x[["amp_genes"]]) |> 
  list_rbind()

del_genes_df <-  map(new_data, ~ .x[["del_genes"]]) |> 
  list_rbind()


loh <- new_filepaths |> 
  map(\(new_filepaths) read_loh_table(filepath = new_filepaths)) |> 
  list_rbind() 



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
  "loh" = loh
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
