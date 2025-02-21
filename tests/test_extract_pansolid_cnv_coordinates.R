library(testthat)

test_that("extract_pansolid_cnv_coordinates handles standard format", {
  
  df <- data.frame(
    "cnv_coord" = c("55174776..55174793"))
  
  df_expected <- data.frame(
    "cnv_coord" = c("55174776..55174793"),
    "start" = c(55174776),
    "end" = c(55174793)) 
  
  expect_equal(extract_pansolid_cnv_coordinates(df, cnv_coord),
               df_expected)
  
})

test_that("extract_pansolid_cnv_coordinates handles complement format", {
  
  df <- data.frame(
    "cnv_coord" = c("complement(21968224..21995324)"))
  
  df_expected  <- data.frame(
    "cnv_coord" = c("complement(21968224..21995324)"),
    "start" = c(21968224),
    "end" = c(21995324))
  
  expect_equal(extract_pansolid_cnv_coordinates(df, cnv_coord),
               df_expected)
  
})

test_that("small and large coordinates are handled", {
  
  df <- data.frame(
    "cnv_coord" = c("1..3",
                    "123456789..1234567890"))
  
  df_expected <- data.frame(
    "cnv_coord" = c("1..3",
                    "123456789..1234567890"),
    "start" = c(1, 123456789),
    "end" = c(3, 1234567890))
  
  expect_equal(extract_pansolid_cnv_coordinates(df, cnv_coord),
               df_expected)
  
})
