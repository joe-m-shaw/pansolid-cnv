# Contamination functions

make_contamination_snp_plot <- function(df,
                                        sample1_ws, 
                                        sample1_labno,
                                        sample2_ws, 
                                        sample2_labno) {
  
  #' Make a SNP plot to investigate contamination
  #'
  #' @param df The dataframe containing all the SNP data. This can be 
  #' created using the `read_snp_sheet` function.
  #' @param sample1_ws The worksheet number of the sample in which
  #' contamination is suspected.
  #' @param sample1_labno The sample number (lab number) of the sample in which
  #' contamination is suspected.
  #' @param sample2_ws The worksheet number of the sample which is suspected to
  #' be contaminating sample 1.
  #' @param sample2_labno The sample number (lab number)  of the sample 
  #' which is suspected to be contaminating sample 1.
  #'
  #' @returns A list of a plot and a dataframe of the SNP comparisons.
  #' @export
  #'
  #' @examples contamination_plot <- make_contamination_snp_plot(
  #' df = all_snp_data, "WS150465", "24026628b",
  #' "WS150529", "24053299")[1]
  
  sample1_df <- df |> 
    filter(worksheet == sample1_ws &
             labno_suffix == sample1_labno)
  
  sample2_df <- df |> 
    filter(worksheet == sample2_ws &
             labno_suffix == sample2_labno)
  
  comparison_df <- sample1_df |> 
    filter(type == "SNV") |> 
    select(chromosome, region, type, reference, reference_allele, allele, frequency) |> 
    left_join(sample2_df |> 
                filter(type == "SNV") |> 
                select(chromosome, region, type, reference, reference_allele,
                       allele, frequency) |> 
                rename(sample2_frequency = frequency),
              join_by(chromosome, region, type, reference, reference_allele,
                      allele)) |> 
    mutate(sample2_frequency_category = case_when(
      sample2_frequency > 90 ~"Above 90%",
      sample2_frequency < 60 &
        sample2_frequency > 40 ~"Between 40-60%",
      sample2_frequency < 10 ~"Lower than 10%",
      TRUE ~"Not present")) |> 
    rename(sample1_frequency = frequency)
  
  comparison_df_no_na <- comparison_df |> 
    filter(!is.na(sample2_frequency))
  
  comparison_plot <- comparison_df |> 
    ggplot(aes(x = region, y = sample1_frequency)) +
    geom_point(shape = 21, aes(fill = sample2_frequency_category),
               alpha = 0.8) +
    scale_fill_manual(values = c("#56B4E9", "#D55E00", "#999999", "white")) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank()) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    labs(x = "SNPs", y = "SNP Allele Frequency",
         fill = paste0("Percentage in ", sample2_labno),
         title = paste0("SNP data from ", sample1_labno),
         subtitle = paste0("SNPs coloured by frequency in ", sample2_labno),
         caption = paste0(nrow(comparison_df_no_na), " SNPs match"))
  
  return(list(comparison_plot, comparison_df))
  
}
