# Qiaseq Primers --------------------------------------------------------------------

library(tidyverse)
library(readxl)

pan_solid_only <- read_excel("data/Primer and Gene Comparison.xlsx",
                             sheet = "Primer Overlap",
                             range = "A2:D11024",
                             col_names = c("chromosome", "coordinates", "sequence",
                                           "gene"),
                             col_types = c("text", "text", "text", "text")) |> 
  mutate(category = "Unique to 44038Z-11379",
         text = "PanSolid only")

core_only <- read_excel("data/Primer and Gene Comparison.xlsx",
                        sheet = "Primer Overlap",
                        range = "G2:J595",
                        col_names = c("chromosome", "coordinates", "sequence",
                                      "gene"),
                        col_types = c("text", "text", "text", "text")) |> 
  mutate(category = "Unique to 17500Z-950",
         text = "Core only")

intersect(pan_solid_only$sequence, core_only$sequence)

both_panels <- read_excel("data/Primer and Gene Comparison.xlsx",
                          sheet = "Primer Overlap",
                          range = "L2:O357",
                          col_names = c("chromosome", "coordinates", "sequence",
                                        "gene"),
                          col_types = c("text", "text", "text", "text")) |> 
  mutate(category = "common to both 44038Z-11379 and 17500Z-950",
         text = "Both")

all_primers <- rbind(pan_solid_only, core_only, both_panels) |> 
  mutate(gene_clean = case_when(
    
    gene == "ERBB2, MIR4728" ~"ERBB2",
    gene == "EGFR, EGFR-AS1" ~"EGFR",
    TRUE ~gene)) 

all_primers |> 
  filter(gene_clean %in% c("ERBB2", "EGFR", 
                           "MYC", "MET", "ARID1A", "SUFU")) |> 
  ggplot(aes(x = gene_clean, y = , fill = text)) +
  geom_bar() +
  theme(legend.title = element_blank()) +
  facet_wrap(~text) +
  labs(x = "", y = "Number of primers", 
       title = "EGFR, ERBB2 and MET share primers between Pan Solid and Core")

all_primers |> 
  filter(gene %in% c("ERBB2", "EGFR", "MET")) |>
  ggplot(aes(x = gene_clean, y = , fill = text)) +
  geom_bar()

primer_table <- all_primers |> 
  filter(gene %in% c("ERBB2", "EGFR", "MET")) |>
  filter(text != "Core only") |> 
  group_by(gene_clean, text) |> 
  summarise(total = n()) |> 
  mutate(prop = round(total/sum(total) * 100, 1)) |> 
  ungroup() |> 
  pivot_wider(names_from = c(text),
              values_from = c(total, prop)) |> 
  arrange(desc(total_Both))

primer_table_formatted <- primer_table |> 
  rename(Gene = gene_clean,
         "Primers shared with Core panel" = total_Both,
         "Primers shared with Core panel (%)" = prop_Both,
         "Primers unique to PanSolid" = "total_PanSolid only",
         "Primers unique to PanSolid (%)" = "prop_PanSolid only") |> 
  select(Gene, "Primers shared with Core panel",
         "Primers shared with Core panel (%)",
         "Primers unique to PanSolid",
         "Primers unique to PanSolid (%)")

export_timestamp(primer_table_formatted)
