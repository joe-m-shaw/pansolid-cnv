# Collate whole genome sequencing HTML data

# Libraries and functions -----------------------------------------------------------

library(here)

source(here("scripts/set_shared_drive_filepath.R"))
source(here("functions/cnv_functions.R"))
source(here("functions/dna_database_connection.R"))

# WGS HTML filepaths ----------------------------------------------------------------

wgs_htmls <- list.files(path = paste0(data_folder, "wgs_result_htmls/"),
                        full.names = TRUE,
                        pattern = "*.html")

# Read HTML data --------------------------------------------------------------------

wgs_tier1_cnvs <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_table_by_div_id(
    html_filepath = wgs_htmls, 
    div_id = "d_svcnv_tier1")) |> 
  list_rbind()

wgs_pids <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_pid_text(wgs_htmls)) |> 
  list_rbind()

wgs_headers <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_header(wgs_htmls)) |> 
  list_rbind()

wgs_tumour_details <- wgs_htmls |> 
  map(\(wgs_htmls) parse_wgs_html_table_by_div_id(
    html_filepath = wgs_htmls, 
    div_id = "t_tumour_details")) |> 
  list_rbind() |> 
  rename(wgs_pathno = histopathology_or_sihmds_lab_id) |> 
  select(filepath, wgs_pathno)

# Get lab numbers -------------------------------------------------------------------

wgs_pathway_tracker <- read_excel(path = paste0(data_folder, 
                                                "WGS pathway tracker_copy_2024-07-01.xlsx"),
                                  sheet = "Cancer") |> 
  janitor::clean_names()

wgs_pathway_tracker_dna_no_df <- wgs_pathway_tracker |> 
  filter(!is.na(mol_db_number) & 
           !is.na(ngis_referral_id) &
           sample_type == "Solid tumour" &
           !duplicated(ngis_referral_id)) |> 
  mutate(labno = str_extract(string = mol_db_number,
                             pattern = "\\d{8}")) |> 
  select(labno, ngis_referral_id)

# Collate data ----------------------------------------------------------------------

wgs_html_data <- inner_join(x = wgs_headers, y = wgs_pids, by = "filepath") |> 
  left_join(wgs_pathway_tracker_dna_no_df, join_by("wgs_r_no" == "ngis_referral_id"),
            relationship = "one-to-one") |> 
  left_join(wgs_tumour_details, by = "filepath") |> 
  left_join(wgs_tier1_cnvs, by = "filepath", relationship = "one-to-many") |> 
  mutate(cnv_class = parse_wgs_cnv_class(col = variant_type),
         cnv_copy_number = parse_wgs_cnv_copy_number(col = variant_type),
         chromosome = parse_wgs_html_grch38_coordinates(col = variant_gr_ch38_coordinates,
                                                        group = "chromosome"),
         cnv_start = as.numeric(parse_wgs_html_grch38_coordinates(col = variant_gr_ch38_coordinates,
                                                              group = "first coordinate")),
         cnv_end = as.numeric(parse_wgs_html_grch38_coordinates(col = variant_gr_ch38_coordinates,
                                                               group = "second coordinate")))

# Plots -----------------------------------------------------------------------------

plot <- wgs_html_data |> 
  filter(grepl(pattern = "pten", x = gene, ignore.case = TRUE) &
           cnv_class == "LOSS" &
           cnv_copy_number == 0) |> 
  ggplot(aes(x = cnv_start, y = wgs_p_no,
             colour = cnv_copy_number)) +
  geom_segment(aes(x = cnv_start, xend = cnv_end,
               y = wgs_p_no, yend = wgs_p_no),
               colour = "#660000") +
  theme_bw() +
  geom_vline(xintercept = 87863113, linetype = "dashed") +
  geom_vline(xintercept = 87971930, linetype = "dashed") +
  labs(title = "Homozygous PTEN deletions by coordinate",
       subtitle = "Dashed lines show gene coordinates")

