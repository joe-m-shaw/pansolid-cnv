# LOH limit of detection

library(tidyverse)
library(here)

source(here("scripts/del_val_load_processed_data.R"))

# PanSolid results --------------------------------------------------------

analyst_comment_regex <- "(.*);\\s(.*)"

del_val_collated_loh <- del_val_collated_loh |> 
  mutate(comment_present = str_detect(string = check_1,
                                      pattern = ";"),
         analyst_result = case_when(
           comment_present == TRUE ~str_extract(string = check_1,
                                                pattern = analyst_comment_regex,
                                                group = 1),
           comment_present == FALSE ~check_1
         ),
         analyst_comment = str_extract(string = check_1,
                                       pattern = analyst_comment_regex,
                                       group = 2))

# NF2 PCR results ---------------------------------------------------------

nf2_fragpcr_results <- read_csv(file = paste0(config::get("data_folderpath"),
                                              "validation/DOC6567_deletions/",
                                              "raw/nf2_loh/",
                                              "nf2_loh_dna_db_results.csv"),
                                col_types = "ccccccdddddcd") |> 
  rename(pcr_median_loh = median_loh,
         pcr_loh = loh_outcome)

# Compare results ---------------------------------------------------------

pansolid_vs_nf2_loh <- del_val_collated_loh |> 
  filter(gene == "NF2") |> 
  inner_join(nf2_fragpcr_results, by = "labno") |> 
  left_join(del_val_collated_stdev |> 
              select(filepath, stdev_noise), by = "filepath",
            relationship = "one-to-one") |> 
  left_join(del_val_collated_138x |> 
              select(filepath, percent_138x), by = "filepath",
            relationship = "one-to-one") |> 
  left_join(del_val_collated_del_genes |> 
              filter(gene == "NF2") |> 
              select(filepath, min_region_fold_change, 
                     max_region_fold_change),
            by = "filepath",
            relationship = "one-to-one") |> 
  mutate(
    outcome_pipeline = case_when(
      loh_status %in% c("Yes", "Yes, No") &
        pcr_loh == "Significant LOH" ~"true_positive",
      loh_status == "No" &
        pcr_loh == "No significant LOH" ~"true_negative",
      loh_status %in% c("Yes", "Yes, No") &
        pcr_loh == "No significant LOH" ~"false_positive",
      loh_status == "No" &
        pcr_loh == "Significant LOH" ~"false_negative"),
    
    outcome_pipeline = factor(outcome_pipeline, levels = c("true_positive",
                                                           "false_positive",
                                                           "false_negative",
                                                           "true_negative")),
    
    outcome_combined = case_when(
      analyst_result == "LOH" &
        pcr_loh == "Significant LOH" ~"true_positive",
      analyst_result == "no LOH" &
        pcr_loh == "No significant LOH" ~"true_negative",
      analyst_result == "LOH" &
        pcr_loh == "No significant LOH" ~"false_positive",
      analyst_result == "no LOH" &
        pcr_loh == "Significant LOH" ~"false_negative"),
    
    outcome_combined = factor(outcome_combined, 
                              levels = c("true_positive",
                                         "false_positive",
                                         "false_negative",
                                         "true_negative")))

# Pivot longer ------------------------------------------------------------

pansolid_vs_nf2_loh_long <- pansolid_vs_nf2_loh |> 
  select(labno, gene, ploidy_state,
         loh_status, no_targets_in_ploidy_region,
         check_1, analyst_result,
         loh_percent_D22S268,
         loh_percent_D22S275,
         loh_percent_NF2CA3,
         loh_percent_NF2intron10,
         pcr_median_loh,
         stdev_noise,
         percent_138x) |> 
  pivot_longer(cols = c(loh_percent_D22S268,
                        loh_percent_D22S275,
                        loh_percent_NF2CA3,
                        loh_percent_NF2intron10),
               names_to = "marker",
               values_to = "percent_loh")

# Draw plots --------------------------------------------------------------

loh_lod_plot <- pansolid_vs_nf2_loh_long |> 
  filter(percent_138x >= 75 & stdev_noise <= 0.7 &
           # Remove sample with contamination
           labno != "24054291") |> 
  filter(!is.na(percent_loh)) |> 
  ggplot(aes(x = reorder(labno, pcr_median_loh),
                                y = percent_loh)) +
  geom_hline(yintercept = 30, linetype = "dashed") +
  geom_jitter(shape = 21, aes(fill = analyst_result),
              alpha = 0.7,
              width = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  labs(x = "Sample", 
       y = "Relative loss of peak height by PCR (%)",
       fill = "PanSolid analyst result",
       title = "LOH testing: PCR testing vs PanSolid", 
       subtitle = "Up to 4 dincucleotide marker results for each sample")


# Theoretical limit of detection ------------------------------------------

loh_df <- data.frame(
  "ncc" = seq(0, 100, by = 1)) |> 
  mutate(stroma_cell_content = 100 - ncc,
         a_count_stroma = stroma_cell_content / 2,
         b_count_stroma = stroma_cell_content / 2,
         # Copy neutral LOH
         a_count_tumour_copy_neutral = ncc,
         b_count_tumour_copy_neutral = 0,
         a_count_total_copy_neutral = a_count_stroma + a_count_tumour_copy_neutral,
         b_count_total_copy_neutral = b_count_stroma + b_count_tumour_copy_neutral,
         b_percentage_copy_neutral = (b_count_total_copy_neutral / 
                                        (b_count_total_copy_neutral + 
                                           a_count_total_copy_neutral)) * 100,
         # Copy loss LOH
         a_count_tumour_copy_loss = ncc / 2,
         b_count_tumour_copy_loss = 0,
         a_count_total_copy_loss = a_count_stroma + a_count_tumour_copy_loss,
         b_count_total_copy_loss = b_count_stroma + b_count_tumour_copy_loss,
         b_percentage_copy_loss = (b_count_total_copy_loss / 
                                        (b_count_total_copy_loss + 
                                           a_count_total_copy_loss)) * 100)

loh_theoretical_lod_plot <- loh_df |>
  select(ncc, b_percentage_copy_neutral, b_percentage_copy_loss) |> 
  pivot_longer(cols = c(b_percentage_copy_neutral, b_percentage_copy_loss),
               names_to = "scenario",
               values_to = "b_percentage") |> 
  ggplot(aes(x = ncc, y = b_percentage)) +
    geom_point(shape = 21, aes(fill = scenario)) +
  scale_y_continuous(breaks = seq(0, 100, by = 5))+
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 35, linetype = "dashed") +
  labs(x = "Neoplastic cell content (%)",
       y = "B allele percentage",
       title = "Copy neutral LOH can be detected at a lower NCC than copy loss LOH")

