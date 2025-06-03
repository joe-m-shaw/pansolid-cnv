# jBRCA noise monitoring

# Connect to DLIMS
source(here::here("scripts/connect_to_dna_db.R"))

source(here::here("functions/pansolid_excel_functions.R"))

# Get a list of PanSolid worksheets ---------------------------------------

all_worksheets <- dna_db_worksheets |> 
  select(pcrid, date, description) |> 
  collect() |> 
  mutate(worksheet = paste0("WS", pcrid))

ps_ws_info <- all_worksheets |> 
  # New CNV Excel layout started with WS152758
  #filter(pcrid >= 141200) |> 
  filter(grepl(pattern = "pansolid", 
               x = description,
               ignore.case = TRUE)) |> 
  mutate(ps_category = case_when(
    grepl(pattern = "jBRCA|j_BRCA|j-BRCA|jew",
          x = description,
          ignore.case = TRUE) ~"PanSolid Jewish BRCA",
    TRUE ~"PanSolid FFPE"
  ))

jbrca_ws_info <- ps_ws_info |> 
  filter(ps_category == "PanSolid Jewish BRCA") |> 
  filter(pcrid < 152828)

jbrca_worksheets <- jbrca_ws_info$worksheet

jbrca_filepaths <- jbrca_worksheets |> 
  map(\(jbrca_worksheets) get_annotated_filepaths(worksheet = jbrca_worksheets)) |> 
  flatten()

brca_stdev <- jbrca_filepaths |> 
  map(\(jbrca_filepaths) 
      read_stdev_results(file = jbrca_filepaths,
                         sheet = get_amp_sheetname(jbrca_filepaths))) |> 
  list_rbind()

ggplot(brca_stdev, aes(x = worksheet, y = st_dev_signal_adjusted_log2_ratios)) +
  geom_boxplot(outliers = FALSE) +
  theme(axis.text.x = element_text(angle = 90))



write_csv(brca_stdev,
          file = paste0(config::get("data_folderpath"),
                        "live_service/collated/",
                        "jbrca_noise.csv"))

nrow(brca_stdev) / 42

length(unique(brca_stdev$worksheet))


# Compare with FFPE noise -------------------------------------------------

data_folder <- config::get("data_folderpath")

collated_data_folder <- paste0(data_folder, "live_service/collated/")

old_stdev <- read_csv(paste0(
  collated_data_folder,
  "pansolid_amplifications_live_service/",
  "live_service_std_dev_results_collated.csv")) |> 
  rename(stdev_noise =
           st_dev_signal_adjusted_log2_ratios)  |> 
  mutate(ps_category = "FFPE")

new_stdev <- read_csv(paste0(collated_data_folder,
                             "stdev_live.csv")) |> 
  select(-filename) |> 
  mutate(ps_category = "FFPE")

jbrca_stdev <- read_csv(paste0(collated_data_folder,
                               "jbrca_noise.csv")) |> 
  mutate(ps_category = "jBRCA") |> 
  rename(stdev_noise = st_dev_signal_adjusted_log2_ratios)

colnames(jbrca_stdev)

all_stdev <- rbind(old_stdev, new_stdev, jbrca_stdev) |> 
  mutate(worksheet_number = parse_number(worksheet)) |> 
  filter(!is.na(worksheet))

all_stdev |> 
  #filter(worksheet_number > 147000 &
          # worksheet_number < 150000) |> 
  ggplot(aes(x = worksheet, y = stdev_noise)) +
  geom_boxplot(outliers = FALSE,
               aes(colour = ps_category)) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(y = "Signal-adjusted noise") 

all_stdev |> 
  filter(ps_category == "jBRCA") |> 
  ggplot(aes(x = worksheet, y = stdev_noise)) +
  geom_boxplot(outliers = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Signal-adjusted noise")



?geom_boxplot()

ggplot(all_stdev, aes(x = worksheet, y = stdev_noise)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))


jbrca_mean <- jbrca_stdev |> 
  group_by(worksheet) |> 
  summarise(mean = round(mean(stdev_noise), 1)) |> 
  arrange(worksheet)

write_csv(jbrca_mean, "jbrca_mean.csv")


