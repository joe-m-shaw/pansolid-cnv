group_sd <- function(df, group_variable1, group_variable2, measurement_variable) {
  
  #' Calculate standard deviations with grouped variables
  #'
  #' This function performs calculations described in 
  #' [NISTIR 6969](https://nvlpubs.nist.gov/nistpubs/ir/2019/NIST.IR.6969-2019.pdf) 
  #' Selected Laboratory and Measurement Practices and Procedures to Support 
  #' Basic Mass Calibrations (2019 Ed)]
  #' (page 211, section 8.4).
  #' 
  #' @param df A dataframe containing measurement data
  #' @param group_variable1 First variable to group by
  #' @param group_variable2 Second variable to group by
  #' @param measurement_variable The variable to calculate the standard deviation for
  #'
  #' @return A grouped data frame showing standard deviation
  #' @export
  #'
  #' @examples group_sd(repeat_data_fold_change_range, labno, gene, fold_change_for_uom)
  
  if(is.data.frame(df) == FALSE){
    stop("df must be a data frame")
  }
  
  if(nrow(df) == 0) {
    stop("df must not be empty")
  }
  
  output_df <- df |> 
    group_by( {{ group_variable1 }}, {{ group_variable2 }} ) |> 
    dplyr::summarise(sd = sd( {{ measurement_variable }} ),
              n = n(),
              n_minus_1 = n-1,
              z = (n_minus_1)*sd^2)
  
  return(output_df)
  
}

pool_sd <- function(df) {
  
  #' Calculate the pooled standard deviation (pooled variance) for grouped standard deviation data
  #'
  #' @param df A dataframe returned from the [group_sd] function.
  #'
  #' @return The pooled standard deviation value
  #' @export
  #'
  #' @examples pool_sd(df = group_sd(repeat_data_fold_change_range, labno, gene, fold_change_for_uom))
  
  pooled_sd <- sqrt(sum(df$z) / (sum(df$n_minus_1)))
  
  return(pooled_sd)
  
}

define_uom <- function(pooled_sd, num_sd) {
  
  #' Define the uncertainty of measurement at a given confidence level
  #'
  #' When a normal (Gaussian) distribution is assumed, a defined number of results will fall within
  #' different distances from the mean, which can be described by [standard deviations](https://en.wikipedia.org/wiki/Normal_distribution#/media/File:Standard_deviation_diagram_micro.svg).
  #' 
  #' - ± 1 standard deviation: 68.2% of results
  #' - ± 2 standard deviations: 95.4% of results
  #' - ± 3 standard deviations: 99.6% of results
  #'
  #' @param pooled_sd The pooled standard deviation
  #' @param num_sd The number of standard deviations from the mean, which defines the confidence level
  #'
  #' @return The uncertainty of measurement value.
  #' @export
  #'
  #' @examples define_uom(0.0781, 2)
  
  uom <- pooled_sd * num_sd
  
  return(uom)
  
}
