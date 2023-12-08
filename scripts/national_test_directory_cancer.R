# National Test Directory - Cancer

# Packages --------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)

# Test Directory --------------------------------------------------------------------

read_excel("data/cancer-national-genomic-test-directory-v7.3.xlsx")

read_directory_tab <- function(input_sheet, input_skip = 0) {
  
  output <- read_excel("data/cancer-national-genomic-test-directory-v7.3.xlsx",
                       sheet = input_sheet,
                       skip = input_skip) |> 
    janitor::clean_names() |> 
    mutate(sheet = input_sheet)
  
  return(output)
  
}

solid <- read_directory_tab(input_sheet = "Solid tumours")

neuro <- read_directory_tab(input_sheet = "Neurological tumours", input_skip = 1) |> 
  dplyr::rename(specialist_test_group = test_group)

sarcoma <- read_directory_tab(input_sheet = "Sarcoma", input_skip = 1)

haem <- read_directory_tab(input_sheet = "Haematological", input_skip = 1)

paed <- read_directory_tab(input_sheet = "Paediatric", input_skip = 1)

all_cancer <- rbind(solid, neuro, sarcoma, haem, paed) |> 
  mutate(gene = str_split(string = target_gene_s_essential, pattern = ", ")) |> 
  unnest(gene) |> 
  mutate(gene = trimws(x = gene, which = "both"))


# CNV analysis ----------------------------------------------------------------------

cnv_variants <- unique(grep(pattern = "copy number variant", x = all_cancer$test_name,
                       value = TRUE))

cnv_targets <- all_cancer |> 
  filter(sheet != "Haematological") |> 
  filter(test_name %in% cnv_variants) |> 
  mutate(gene_clean = ifelse(
    gene %in% c("1p19q codel", "1p19q", "1p19qcodel"), "1p19q codel", gene)) |> 
  arrange(test_code)

cnv_target_summary <- cnv_targets |> 
  group_by(gene_clean) |> 
  summarise(total = n()) |> 
  arrange(desc(total)) 

length(unique(cnv_targets$test_code))
length(unique(cnv_targets$gene_clean))

ggplot(cnv_target_summary, aes(x = reorder(gene_clean, total,
                         decreasing = TRUE), y = total)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "", title = "Number of times genes are on test directory codes for CNV analysis")

ggplot(cnv_targets, aes(x = test_code, y = gene_clean)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "", y = "Test directory code")

# Checking gene frequency -----------------------------------------------------------

count_gene_occurences <- function(input_gene) {
  
  output <- all_cancer |> 
    filter(sheet != "Haematological") |> 
    filter(gene == input_gene) |> 
    select(test_code, test_name)
  
  return(output)
  
}

count_gene_occurences("19q") |>  view()

count_gene_occurences("CDKN2A") |>  view()

