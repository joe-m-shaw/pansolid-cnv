# Monitor PanSolid CNV service

# Connect to DLIMS
source(here::here("scripts/connect_to_dna_db.R"))

# Find a list of worksheets with pansolid in the description, excluding jBRCA

all_worksheets <- dna_db_worksheets |> 
  select(pcrid, date, description) |> 
  collect() |> 
  mutate(worksheet = paste0("WS", pcrid))

ps_ws_info <- all_worksheets |> 
  filter(pcrid >= 152758) |> 
  filter(grepl(pattern = "pansolid", 
               x = description,
               ignore.case = TRUE)) |> 
  filter(str_detect(string = description,
                    pattern = "jBRCA",
                    negate = TRUE))

ps_worksheets <- ps_ws_info$worksheet

# Look for results Excels from those worksheets on the S drive

source(here::here("functions/pansolid_excel_functions.R"))

ps_filepaths <- ps_worksheets |> 
  map(\(ps_worksheets) get_annotated_filepaths(worksheet = ps_worksheets)) |> 
  flatten()

# Copy new results into a local folder, 
# Read in the results of those new files and collate them
# Add those results to the collated data
# Delete the copied results in the local folder

