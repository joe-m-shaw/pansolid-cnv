
library(tidyverse)
source(here::here("functions/extract_pansolid_cnv_coordinates.R"))
source(here::here("functions/chromosome_arm_functions.R"))

pansolid_bed <- read_csv(paste0(config::get("data_folderpath"),
                                "validation/DOC6283_amplifications/bed_files/",
                                "PanSolidv2_GRCh38_noalt_BED.csv"),
                         col_types = "cccd") |> 
  janitor::clean_names() 

# Find max and min target coordinates for each arm ------------------------

pansolid_bed_mod <- pansolid_bed |> 
  extract_pansolid_cnv_coordinates(cnv_coord_col = region) |> 
  rename(chromosome_char = chromosome) |> 
  factorise_chromosome(col = chromosome_char) |> 
  add_chromosome_arms()

pansolid_target_region_chr_arm_bed <- pansolid_bed_mod |> 
  group_by(chromosome_fct, chromosome_arm) |> 
  summarise(chr_arm_min_target_coordinate = min(start),
            chr_arm_max_target_coordinate = max(end),
            chr_arm_targets = n()) |> 
  mutate(chr_arm_target_covered_region_length = chr_arm_max_target_coordinate - chr_arm_min_target_coordinate)

write_csv(pansolid_target_region_chr_arm_bed |> 
            select(chromosome_fct, 
                   chromosome_arm,
                   chr_arm_min_target_coordinate, 
                   chr_arm_max_target_coordinate, 
                   chr_arm_targets, 
                   chr_arm_target_covered_region_length),
          paste0(config::get("data_folderpath"),
                 "validation/DOC6791_chromosome_arms/bed_files/",
                 "pansolid_target_region_chr_arm_bed.csv"))

# Get minimum coordinates -------------------------------------------------

chr_coords_long <- pansolid_bed_mod |> 
  group_by(chromosome_fct) |> 
  summarise(min_coord = min(start),
            max_coord = max(end)) |> 
  pivot_longer(cols = c(min_coord, max_coord),
               names_to  = "category",
               values_to = "coordinate") |>  
  mutate(cumulative_coordinate = case_when(
    chromosome_fct == "1" ~coordinate,
    chromosome_fct != "1" ~cumsum(coordinate)
  ))

pansolid_chr_cumulative_coordinates <- chr_coords_long |> 
  filter(category == "min_coord")
  
write_csv(pansolid_chr_cumulative_coordinates,
          paste0(config::get("data_folderpath"),
                 "validation/DOC6791_chromosome_arms/",
                 "bed_files/",
                 "pansolid_chr_cumulative_coordinates.csv"))
