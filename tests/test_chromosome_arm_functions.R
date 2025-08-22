library(testthat)

test_that("format_chromosome_decimals works with correct input", {
  
  input_df <- data.frame(
    "chromosome" = c("1.0", "2.0", "3.0", "10.0", "X")
  )
  
  expected_df <- data.frame(
    "chromosome" = c("1.0", "2.0", "3.0", "10.0", "X"),
    "chromosome_char" = c("1", "2", "3", "10", "X")
  )
  
  expect_equal(format_chromosome_decimals(input_df),
               expected_df)
  
})



