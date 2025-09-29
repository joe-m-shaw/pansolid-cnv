
read_biorad_ddpcr_csv <- function(filepath) {
  
  #' Read a CSV file from a BioRad ddPCR experiment as a data-frame
  #'
  #' @param filepath The full filepath of the CSV
  #'
  #' @return The data as a data-frame with column names in snakecase and 
  #' appropriate column types.
  #' 
  #' @note This function is designed to handle droplet digital polymerase chain reaction 
  #' (ddPCR) csv files exported from Quantasoft v1.7.4 (BioRad).
  
  ddpcr_df <- readr::read_csv(filepath, 
                       col_types = cols(
                         "Well" = "c",
                         "ExptType" = "c",
                         "Experiment" = "c",
                         "Sample" = "c",
                         "TargetType" = "c",
                         "Target" = "c",
                         "Status" = "c",
                         "Concentration" = "d",
                         "Supermix" = "c",
                         "CopiesPer20uLWell" = "d",
                         "TotalConfMax" = "d",
                         "TotalConfMin" = "d",
                         "PoissonConfMax" = "d",
                         "PoissonConfMin" = "d",
                         "Positives" = "i",
                         "Negatives" = "i",
                         "Ch1+Ch2+" = "i",
                         "Ch1+Ch2-" = "i",
                         "Ch1-Ch2+" = "i",
                         "Ch1-Ch2-" = "i",
                         "Linkage"  = "d",
                         "AcceptedDroplets" = "i",
                         "CNV" = "d",
                         "TotalCNVMax" = "d",
                         "TotalCNVMin" = "d",
                         "PoissonCNVMax" = "d",
                         "PoissonCNVMin" = "d",
                         "FractionalAbundance" = "d",
                         "TotalFractionalAbundanceMax" = "d",
                         "TotalFractionalAbundanceMin" = "d",
                         "PoissonFractionalAbundanceMax" = "d",
                         "PoissonFractionalAbundanceMin" = "d",
                         "ReferenceAssayNumber" = "d",
                         "TargetAssayNumber" = "d",
                         "Threshold" = "d",
                         "MeanAmplitudeofPositives" = "d",
                         "MeanAmplitudeofNegatives" = "d",
                         "MeanAmplitudeTotal" = "d",
                         "ExperimentComments" = "c",
                         "MergedWells" = "c",
                         "TotalConfMax68" = "d",
                         "TotalConfMin68" = "d",
                         "PoissonConfMax68" = "d",
                         "PoissonConfMin68" = "d",
                         "TotalCNVMax68" = "d",
                         "TotalCNVMin68" = "d",
                         "PoissonCNVMax68" = "d",
                         "PoissonCNVMin68" = "d",
                         "PoissonCNVMin68" = "d",
                         "PoissonRatioMax68" = "d",
                         "TotalRatioMin68" = "d",
                         "TotalFractionalAbundanceMax68" = "d",
                         "TotalFractionalAbundanceMin68" = "d",
                         "PoissonFractionalAbundanceMax68" = "d",                
                         "PoissonFractionalAbundanceMin68" = "d")) |> 
    janitor::clean_names() 
  
  if(ncol(ddpcr_df) != 63) {
    stop("Incorrect number of columns in ddPCR dataframe")
  }
  
  if(nrow(ddpcr_df) == 0) {
    stop("No rows in ddPCR dataframe")
  }
  
  output <- ddpcr_df |> 
    dplyr::mutate(filepath = filepath)
  
  return(output)
  
}   