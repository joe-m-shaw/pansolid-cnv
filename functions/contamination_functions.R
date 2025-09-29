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
  #' @returns A plot of the SNP comparisons.
  #' @export
  #'
  #' @examples contamination_plot <- make_contamination_snp_plot(
  #' df = all_snp_data, "WS150465", "24026628b",
  #' "WS150529", "24053299")
  
  sample1_df <- df |> 
    dplyr::filter(worksheet == sample1_ws &
             labno_suffix == sample1_labno)
  
  sample2_df <- df |> 
    dplyr::filter(worksheet == sample2_ws &
             labno_suffix == sample2_labno)
  
  comparison_df <- sample1_df |> 
    dplyr::filter(type == "SNV") |> 
    dplyr::select(chromosome, region, type, reference, reference_allele, 
                  allele, frequency) |> 
    dplyr::left_join(sample2_df |> 
                       dplyr::filter(type == "SNV") |> 
                       dplyr::select(chromosome, region, type, 
                                     reference, reference_allele,
                                     allele, frequency) |> 
                       dplyr::rename(sample2_frequency = frequency),
                     dplyr::join_by(chromosome, region, type, 
                                    reference, reference_allele,
                                    allele)) |> 
    dplyr::mutate(sample2_frequency_category = case_when(
      sample2_frequency >= 90 ~"Above 90%",
      sample2_frequency <= 60 &
        sample2_frequency >= 40 ~"Between 40-60%",
      sample2_frequency <= 10 ~"Lower than 10%",
      TRUE ~"Not present")) |> 
    dplyr::rename(sample1_frequency = frequency)
  
  comparison_df_no_na <- comparison_df |> 
    dplyr::filter(!is.na(sample2_frequency))
  
  comparison_plot <- comparison_df |> 
    ggplot2::ggplot(aes(x = region, y = sample1_frequency)) +
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
  
  return(comparison_plot)
  
}

make_contamination_snp_plotlist <- function(df,
                                            sample1_labno,
                                            sample1_ws,
                                            sample2_ws,
                                            plotlist_title) {
  
  #' Generate a plotlist of comparison SNP plots for a single sample
  #' 
  #' This function allows you to quickly investigate whether a sample shares a 
  #' similar SNP profile (and therefore is likely contaminated with) another
  #' sample on the same PanSolid worksheet.
  #'
  #' @param df The dataframe containing all the SNP data.
  #' @param sample1_labno The sample number (lab number) of the sample in which
  #' contamination is suspected.
  #' @param ws The worksheet that `sample1 labno` is on.
  #' @param plotlist_title The title for the final plot file, which should end 
  #' in an appropriate suffix such as ".pdf"
  #'
  #' @returns This functions exports the plotlist as a multi-page file to your
  #' local directory.
  #' @export
  #'
  #' @examples 
  #' 
  #' # Here is an example of using the function to investigate potential 
  #' # contamination of sample 25048537 on WS156496
  #' 
  #' WS156496_files <- list.files(path = paste0("S:/central shared/",
  #'                               "Genetics/Repository/",
  #'                                "WorksheetAnalysedData/",
  #'                                 "WS156496"),
  #'                                 pattern = "Results_SNVs_Indels.*.xlsx",
  #'                                 full.names = TRUE,
  #'                                 recursive = TRUE)
  #'                                 
  #' WS156496_snp_data <- WS156496_files |> 
  #'   map(\(WS156496_files) read_snp_sheet(WS156496_files)) |> 
  #'   list_rbind()
  #'   
  #' make_contamination_snp_plotlist(df = WS156496_snp_data,
  #'                                 sample1_labno = "25048537",
  #'                                 ws = "WS156496",
  #'                                 plotlist_title = "WS156496_plots_25048537.pdf")

  query_sample_df <- df |> 
    dplyr::filter(labno != sample1_labno)
  
  query_labnos <- unique(query_sample_df$labno)
  
  plotlist <- list()
  
  for(i in query_labnos) {
    
    new_plot <- make_contamination_snp_plot(df = df,
                                            sample1_ws = sample1_ws, 
                                            sample1_labno = sample1_labno,
                                            sample2_ws = sample2_ws, 
                                            sample2_labno = i)
    
    plotlist <- list(plotlist, new_plot)
    
    rm(new_plot)
    
  }
  
  ggpubr::ggexport(plotlist = plotlist,
                   filename = plotlist_title,
                   width=15, 
                   height=8, 
                   res=300)
  
}
