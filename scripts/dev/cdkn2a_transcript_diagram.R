

pms2_results <- del_val_collated_loh |> 
  filter(gene == "PMS2")

colnames(del_val_collated_loh)

unique(del_val_collated_loh$ploidy_state)



pms2_results |> 
  group_by(ploidy_state) |> 
  count() |> 
  arrange(desc(n))

pms2_results |> 
  filter(ploidy_state %in% c("Deletion"))


del_val_collated_sig_cnvs |> 
  filter(gene == "PMS2")

standard_loh_calls <- c("Normal diploid", "Uniparental disomy",
                        "Duplication", "Bi-allelic deletion", "Deletion",
                        "Whole genome duplication")

del_val_collated_loh_mod <- del_val_collated_loh |> 
  left_join(del_val_collated_138x |> 
              select(filepath, percent_138x),
            by = "filepath") |> 
  left_join(del_val_collated_stdev |> 
              select(filepath, stdev_noise),
            by = "filepath") |> 
  mutate(loh_number = str_extract(pattern = "^\\(\\d{1},\\d{1}\\)$",
               string = ploidy_state),
         loh_simplify = case_when(
          ploidy_state %in% standard_loh_calls ~ploidy_state,
          !is.na(loh_number) ~ploidy_state,
          TRUE ~"Multiple calls"
         ),
         gene = factor(gene, levels = c("MSH2",
                                        "MSH6",
                                        "MLH1",
                                        "PMS2",
                                        "LZTR1",
                                        "SMARCB1",
                                        "NF2")))
  
del_val_collated_loh_mod |> 
  filter(percent_138x >= 75 &
           stdev_noise <= 0.7) |> 
  group_by(gene, loh_simplify) |> 
  count() |> 
  ggplot(aes(x = n, y = reorder(loh_simplify, n))) +
  geom_col() +
  facet_wrap(~gene, nrow = 2) +
  theme_bw() +
  labs(x = "Samples", y = "Ploidy state")

grch38_primer_coordinates |> 
  filter(chromosome  == "9") |>  view()


transcript_folder <- paste0(config::get("data_folderpath"),
                                        "validation/DOC6567_deletions/transcripts/")

transcript_files <- list.files(transcript_folder, full.names = TRUE)


cdkn2ab_df <- transcript_files |> 
  map(\(transcript_files) read_ensembl_exon_table(transcript_files)) |> 
  list_rbind()

cdkn2ab_df_mod <- cdkn2ab_df |> 
  mutate(
    gene = case_when(
      transcript == "ENST00000304494" ~"CDKN2A",
      transcript == "ENST00000579755" ~"CDKN2A",
      transcript == "ENST00000276925" ~"CDKN2B"
    ),
    refseq = case_when(
      transcript == "ENST00000304494" ~"NM_000077.5",
      transcript == "ENST00000579755" ~"NM_058195.4",
      transcript == "ENST00000276925" ~"NM_004936.4"
    ),
    refseq = factor(refseq, levels = c("NM_004936.4",
                                       "NM_058195.4",
                                       "NM_000077.5")),
    exon_char = case_when(
      transcript == "ENST00000304494" & exon == 1 ~"1a",
      transcript == "ENST00000579755" & exon == 1 ~"1b",
      TRUE ~as.character(exon)
    ))

str(cdkn2ab_df_mod$refseq)

plot_buffer <- 5000

plot_start <- 21967752 - plot_buffer

plot_end <- 22009305 + plot_buffer

cdkn2a_targets <- target_df_with_coordinates |> 
  filter(chromosome == "9") |> 
  filter(start > plot_start &
           end < plot_end) |> 
  mutate(y_value = "Targets")

cdkn2a_target_plot <- ggplot(cdkn2a_targets, aes(x = start, y = y_value)) +
  geom_point(shape = 21, size = 3) +
  theme_bw() +
  labs (x = "GRCh38 coordinates (chr9)", y = "") +
  scale_x_continuous(limits = c(plot_start, plot_end))

cdkn2a_transcripts_plot <- ggplot(cdkn2ab_df_mod, aes(x = start, y = refseq,
                           label = exon_char)) +
  geom_segment(aes(x = start, xend = end,
                   y = refseq, yend = refseq),
               linewidth = 4) +
  geom_segment(aes(
    x = 21974857, xend = 21967752,
    y = "NM_000077.5", yend = "NM_000077.5")) +
  geom_segment(aes(
    x = 21994392, xend = 21967752,
    y = "NM_058195.4", yend = "NM_058195.4")) +
  geom_segment(aes(
    x = 22002903, xend = 22009305,
    y = "NM_004936.4", yend = "NM_004936.4")) +
  geom_label(nudge_y = 0.2) +
  theme_bw() +
  labs(title = "The CDKN2A/B locus", y = "", x = "") +
  theme(axis.text.x = element_blank()) +
  scale_x_continuous(limits = c(plot_start, plot_end)) +
  geom_text(x = 21980000, y = "Gene",
            label = "CDKN2A", nudge_y = 0.4) +
  geom_text(x = 22005000, y = "Gene",
            label = "CDKN2B", nudge_y = 0.4)

(cdkn2a_transcripts_plot / cdkn2a_target_plot) +
  plot_layout(
    heights = c(3, 1)
  )






  
