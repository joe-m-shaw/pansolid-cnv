library(testthat)

# format_chromosome_decimals

test_that("format_chromosome_decimals works with correct input", {
  
  input_df <- data.frame(
    "chromosome" = c("1.0", "2.0", "3.0", "10.0", "X")
  )
  
  expected_df <- data.frame(
    "chromosome" = c("1.0", "2.0", "3.0", "10.0", "X"),
    "chromosome_char" = c("1", "2", "3", "10", "X")
  )
  
  expect_equal(format_chromosome_decimals(df = input_df,
                                          col = chromosome),
               expected_df)
  
})

test_that("format_chromosome_decimals works without decimal input", {
  
  input_df <- data.frame(
    "chromosome" = c("1", "2", "3", "10", "10", "X")
  )
  
  expected_df <- data.frame(
    "chromosome" = c("1", "2", "3", "10", "10", "X"),
    "chromosome_char" = c("1", "2", "3", "10", "10", "X")
  )
  
  expect_equal(format_chromosome_decimals(df = input_df,
                                          col = chromosome),
               expected_df)
  
})

test_that("format_chromosome_decimals works with mixed input", {
  
  input_df <- data.frame(
    "chromosome" = c("1", "2.0", "3", "10.0", "10", "X")
  )
  
  expected_df <- data.frame(
    "chromosome" = c("1", "2.0", "3", "10.0", "10", "X"),
    "chromosome_char" = c("1", "2", "3", "10", "10", "X")
  )
  
  expect_equal(format_chromosome_decimals(df = input_df,
                                          col = chromosome),
               expected_df)
  
})

test_that("format_chromosome_decimals works with numeric input", {
  
  input_df <- data.frame(
    "chromosome" = c(1.0, 2.0, 3.0, 10.0, 11.0))
  
  expected_df <- data.frame(
    "chromosome" =      c(1.0, 2.0, 3.0, 10.0, 11.0),
    "chromosome_char" = c("1", "2", "3", "10", "11"))
  
  expect_equal(format_chromosome_decimals(input_df,
                                          col = chromosome),
               expected_df)
  
})

# factorise_chromosome

test_that("factorise_chromosome works with character input", {
  
  input_df <- data.frame(
    "chromosome" = c("1", "2", "10", "X")
  )
  
  output_df <- factorise_chromosome(input_df, 
                                    col = chromosome)
  
  expect_equal(class(output_df$chromosome_fct),
               "factor")
  
  expect_equal(levels(output_df$chromosome_fct),
               c("1", "2", "3", "4", "5", "6", "7", "8",
                 "9", "10", "11", "12", "13", "14", "15", 
                 "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
  
})

test_that("factorise_chromosome works with numeric input", {
  
  input_df <- data.frame(
    "chromosome" = c(1, 2, 3.0, 4, 10, 11))
  
  output_df <- factorise_chromosome(input_df, 
                                    col = chromosome)
  
  expect_equal(class(output_df$chromosome_fct),
               "factor")
  
  expect_equal(levels(output_df$chromosome_fct),
               c("1", "2", "3", "4", "5", "6", "7", "8",
                 "9", "10", "11", "12", "13", "14", "15", 
                 "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
  
})


test_that("factorise_chromosome and format_chromosome_decimals works together", {
  
  input_df <- data.frame(
      "chromosome" = c("1.0", "2.0", "3.0", "10.0", "10", "X")
    )
  
  expected_df <- data.frame(
    "chromosome" =      c("1.0", "2.0", "3.0", "10.0", "10", "X"),
    "chromosome_char" = c("1", "2", "3", "10", "10", "X"),
    "chromosome_fct" = factor(x = c("1", "2", "3", "10", "10", "X"),
                              levels = c("1", "2", "3", "4", "5", 
                                         "6", "7", "8", "9", "10", 
                                         "11", "12", "13", "14", "15", 
                                         "16", "17", "18", "19", "20", 
                                         "21", "22", "X", "Y")))
    
    expect_equal(input_df |> 
                   format_chromosome_decimals(col = chromosome) |> 
                   factorise_chromosome(col = chromosome_char),
                 expected_df)
  
})
