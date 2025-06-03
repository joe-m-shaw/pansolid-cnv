# QIAsymphony issue

source(here::here("scripts/connect_to_dna_db.R"))

batches_since_friday <- extraction_batch_tbl |> 
  filter(run_date >= "2025-05-30 00:00:00") |> 
  collect() |> 
  filter(extraction_method_fk %in% c(33))

batches_to_check <- unique(batches_since_friday$extraction_batch_id)

samples_since_friday <- extraction_tbl |> 
  filter(extraction_batch_fk %in% batches_to_check) |> 
  collect() |> 
  select(extraction_batch_fk, lab_no)

write_csv(samples_since_friday,
          "qiasymphony_dna_samples_since_2025_05_30.csv")

