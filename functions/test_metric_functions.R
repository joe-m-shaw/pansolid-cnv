
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

