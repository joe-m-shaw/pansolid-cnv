# Get PanSolid DNA concentrations

library(here)

source("scripts/connect_to_dna_db.R")

source("scripts/load_processed_live_service_data.R")

source("scripts/set_shared_drive_filepath.R")

pansolid_samples <- unique(live_service_std_dev_results_collated$labno)

pansolid_sample_dna_conc <- sample_tbl |> 
  filter(labno %in% pansolid_samples) |> 
  select(labno, concentration) |> 
  collect()

write_csv(x = pansolid_sample_dna_conc, 
          file = paste0(data_folder, "live_service/collated/pansolid_sample_dna_conc.csv"))


DBI::dbDisconnect(dbi_con)
