library(testthat)

test_that("parse_filename handles original filename format",{
  
  x <- "Annotated_WS140721_24018287_NatashaROSTOVA.xlsx"
  
  expect_equal(parse_filename(x, 1),
               "WS140721")
  
  expect_equal(parse_filename(x, 2),
               "24018287")
  
  expect_equal(parse_filename(x, 3),
               "")
  
  expect_equal(parse_filename(x, 4),
               "NatashaROSTOVA")
  
})

test_that("parse_filename handles patient names with numbers",{
  
  x <- "Annotated_WS140954_24017042_AndreiBolkonsky0133841080.xlsx"
  
  expect_equal(parse_filename(x, 1),
               "WS140954")
  
  expect_equal(parse_filename(x, 2),
               "24017042")
  
  expect_equal(parse_filename(x, 3),
               "")
  
  expect_equal(parse_filename(x, 4),
               "AndreiBolkonsky0133841080")
  
})

test_that("parse_filename handles new format with panel and patient names", {
  
  x <- "Annotated_v2M7_MELA_PS_WS147208_24060851_PierreBEZUKHOV.xlsx"
  
  expect_equal(parse_filename(x,1),
               "WS147208")
  
  expect_equal(parse_filename(x,2),
               "24060851")
  
  expect_equal(parse_filename(x,3),
               "")
  
  expect_equal(parse_filename(x,4),
               "PierreBEZUKHOV")
  
})

test_that("parse_filename handles lab numbers with suffixes", {
  
  x <- "Annotated_v2SchwannAll_PS_WS142945_24032975c_PetyaROSTOV.xlsx"
  
  expect_equal(parse_filename(x, 1),
               "WS142945")
  
  expect_equal(parse_filename(x, 2),
               "24032975")
  
  expect_equal(parse_filename(x, 3),
               "c")
  
  expect_equal(parse_filename(x, 4),
               "PetyaROSTOV")
  
})

test_that("parse_filename handles new format with panel and without patient names", {
  
  x <- "Annotated_v2aM4_LUNG_PS_WS148086_24066935_S14.xlsx"
  
  expect_equal(parse_filename(x, 1),
               "WS148086")
  
  expect_equal(parse_filename(x, 2),
               "24066935")
  
  expect_equal(parse_filename(x, 3),
               "")
  
  expect_equal(parse_filename(x, 4),
               "S14")
  
})

test_that("filename_to_df handles filename input", {
  
  x <- "Annotated_v2aM7_MELA_PS_WS123456_12345678_AnnaKARENINA.xlsx"
  
  expected_df <- data.frame(
    "worksheet" = c("WS123456"),
    "labno" = c("12345678"),
    "suffix" = c(""),
    "patient_name" = c("AnnaKARENINA"),
    "labno_suffix" = c("12345678"),
    "labno_suffix_worksheet" = c("12345678_WS123456")
  )
  
  expect_equal(filename_to_df(x),
               expected_df)
  
})

test_that("filename_to_df handles filepath input", {
  
  s_drive_path <- "S:/central shared/Genetics/Repository/WorksheetAnalysedData/"
  
  ws_path <- "WS123456/v2aM7_MELA_PS/"
  
  file <- "Annotated_v2aM7_MELA_PS_WS123456_12345678_AnnaKARENINA.xlsx"
  
  filepath <- paste0(s_drive_path, ws_path, file)
  
  expected_df <- data.frame(
    "worksheet" = c("WS123456"),
    "labno" = c("12345678"),
    "suffix" = c(""),
    "patient_name" = c("AnnaKARENINA"),
    "labno_suffix" = c("12345678"),
    "labno_suffix_worksheet" = c("12345678_WS123456")
  )
  
  expect_equal(filename_to_df(filepath),
               expected_df)
  
})

test_that("filename_to_df handles edge cases with low DNA numbers", {
  
  s_drive_path <- "S:/central shared/Genetics/Repository/WorksheetAnalysedData/"
  
  ws_path <- "WS141272/v2R215_CDH1_PS/"
  
  file <- "Annotated_WS141272_088727_IvanILYICH.xlsx"
  
  filepath <- paste0(s_drive_path, ws_path, file)
  
  expected_df <- data.frame(
    "worksheet" = c("WS141272"),
    "labno" = c("088727"),
    "suffix" = c(""),
    "patient_name" = c("IvanILYICH"),
    "labno_suffix" = c("088727"),
    "labno_suffix_worksheet" = c("088727_WS141272")
  )
  
  expect_equal(filename_to_df(filepath),
               expected_df)
  
})
