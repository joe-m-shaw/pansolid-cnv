
calc_sensitivity <- function(tp, tn, fp, fn) {
  
  #' Calculate the sensitivity for a test in comparison to another test
  #'
  #' @param tp True positives
  #' @param tn True negatives
  #' @param fp False positives
  #' @param fn False negatives
  #'
  #' @return The sensitivity (true positive rate) of the test as a percentage
  #' @export
  #'
  #' @examples calc_sensitivity(tp = 40, tn = 44, fn = 3, fp = 0)
  
  sensitivity <- (tp / (tp + fn)) * 100
  
  return(sensitivity)
  
}

calc_specificity <- function(tp, tn, fp, fn) {
  
  #' Calculate the specifity for a test in comparison to another test
  #'
  #' @param tp True positives
  #' @param tn True negatives
  #' @param fp False positives
  #' @param fn False negatives
  #'
  #' @return The specificity (true negative rate) as a percentage
  #' @export
  #'
  #' @examples calc_specificity(tp = 40, tn = 44, fn = 3, fp = 0)
  
  specificity <- (tn / (tn + fp)) * 100
  
  return(specificity)
  
}


get_epir_metrics <- function(outcome_vector) {
  
  #' Get test metrics using epiR
  #'
  #' @param outcome_vector Numeric vector in the order: true positives, false positives,
  #' false negatives, true negatives
  #'
  #' @return A dataframe with test metrics formatted for publication.
  #' @export
  #'
  #' @examples
  
  epi_test <- epiR::epi.tests(dat = outcome_vector,
                   digits = 2, conf.level = 0.95, method = "exact")

  epi_test_df <- tibble::tibble(epi_test$detail) |> 
    dplyr::filter(statistic %in% c("se", "sp", "pv.pos", "pv.neg",
                                   "diag.ac")) |> 
    dplyr::mutate(
      `Percentage (%)` = round(est*100, 1),
      `Confidence interval (95%)` = str_c(round(lower*100, 1), "-", round(upper*100, 1)),
      Metric = case_when(
        statistic == "se" ~"Sensitivity",
        statistic == "sp" ~"Specificity",
        statistic == "pv.pos" ~"Positive predictive value",
        statistic == "pv.neg" ~"Negative predictive value",
        statistic == "diag.ac" ~"Correctly classified proportion"
      )) |> 
    dplyr::select(Metric, `Percentage (%)`, `Confidence interval (95%)`)
  
  return(epi_test_df)
  
}

make_se_sp_table <- function(outcome_vector) {
  
  epi_test <- epiR::epi.tests(dat = outcome_vector,
                              digits = 2, conf.level = 0.95, method = "exact")
  
  se_sp_df <- tibble::tibble(epi_test$detail) |> 
    dplyr::filter(statistic %in% c("se", "sp")) |> 
    dplyr::mutate(
      percent = round(est*100, 1),
      ci = str_c(round(lower*100, 1), "-", round(upper*100, 1)),
      percent_ci = str_c(percent, " (", ci, ")")) |> 
    dplyr::select(statistic, percent_ci) |> 
    tidyr::pivot_wider(names_from = statistic,
                values_from = percent_ci) |> 
    dplyr::mutate(`True positives` = outcome_vector[1],
           `False positives` = outcome_vector[2],
           `False negatives` = outcome_vector[3],
           `True negatives` = outcome_vector[4]) |> 
    dplyr::relocate(se, .after = `True negatives`) |> 
    dplyr::relocate(sp, .after = se) |> 
    dplyr::rename(`Sensitivity (%)` = se,
           `Specificity (%)` = sp)
  
  return(se_sp_df)
  
}
