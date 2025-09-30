# parse_filename

config_file <- file.path("..", "..", "config.yml")

data_folderpath <- config::get(file = config_file,
                               config = "data_folderpath")

test_datapath <- paste0(data_folderpath,
                          "validation/DOC6567_deletions/test_data/")

test_that("parse_filename handles original filename format",{
  
  x <- paste0(test_datapath,
              "Annotated_v2PANSOLID_WS143574_24030966_NatashaROSTOVA.xlsx")
    
  expect_equal(parse_filename(x, 1),
               "WS143574")
  
  expect_equal(parse_filename(x, 2),
               "24030966")
  
  expect_equal(parse_filename(x, 3),
               "")
  
  expect_equal(parse_filename(x, 4),
               "NatashaROSTOVA")
  
})

test_that("parse_filename handles patient names with numbers",{
  
  x <- "Annotated_WS140954_24017042_BorisDrubetskoy0133841080.xlsx"
  
  expect_equal(parse_filename(x, 1),
               "WS140954")
  
  expect_equal(parse_filename(x, 2),
               "24017042")
  
  expect_equal(parse_filename(x, 3),
               "")
  
  expect_equal(parse_filename(x, 4),
               "BorisDrubetskoy0133841080")
  
})

test_that("parse_filename handles new format with panel and patient names", {
  
  x <- paste0(test_datapath,
              "Annotated_v2PANSOLID_WS143415_24030946_PierreBEZUKHOV.xlsx")
  
  expect_equal(parse_filename(x,1),
               "WS143415")
  
  expect_equal(parse_filename(x,2),
               "24030946")
  
  expect_equal(parse_filename(x,3),
               "")
  
  expect_equal(parse_filename(x,4),
               "PierreBEZUKHOV")
  
})

test_that("parse_filename handles lab numbers with suffixes", {
  
  x <- paste0(test_datapath,
              "Annotated_v2PANSOLID_WS150465_24026628a_PetyaROSTOV.xlsx")
  
  expect_equal(parse_filename(x, 1),
               "WS150465")
  
  expect_equal(parse_filename(x, 2),
               "24026628")
  
  expect_equal(parse_filename(x, 3),
               "a")
  
  expect_equal(parse_filename(x, 4),
               "PetyaROSTOV")
  
})

test_that("parse_filename handles new format with panel and without patient names", {
  
  x <- paste0(test_datapath,
              "Annotated_v2M4_LUNG_PS_WS143513_24037405_S32.xlsx")
    
  expect_equal(parse_filename(x, 1),
               "WS143513")
  
  expect_equal(parse_filename(x, 2),
               "24037405")
  
  expect_equal(parse_filename(x, 3),
               "")
  
  expect_equal(parse_filename(x, 4),
               "S32")
  
})

test_that("parse_filename handles .csv filetype", {
  
  x <- paste0(test_datapath,
              "All ROI Coverage Table-UMI Read Mapping-WS123790_22000209_NikolaiLEVIN_S25_R1_001.csv")
  
  expect_equal(parse_filename(x, 1),
               "WS123790")
  
  expect_equal(parse_filename(x, 2),
               "22000209")
  
  expect_equal(parse_filename(x, 3),
               "")
  
  expect_equal(parse_filename(x, 4),
               "NikolaiLEVIN")
  
})

# filename_to_df

test_that("filename_to_df handles filename input", {
  
  x <- "Annotated_v2aSchwannAll_PS_WS140775_24018518_AnnaKARENINA.xlsx"
  
  expected_df <- data.frame(
    "worksheet" = c("WS140775"),
    "labno" = c("24018518"),
    "suffix" = c(""),
    "patient_name" = c("AnnaKARENINA"),
    "labno_suffix" = c("24018518"),
    "labno_suffix_worksheet" = c("24018518_WS140775")
  )
  
  expect_equal(filename_to_df(x),
               expected_df)
  
})

test_that("filename_to_df handles filepath input", {
  
  s_drive_path <- "S:/central shared/Genetics/Repository/WorksheetAnalysedData/"
  
  ws_path <- "WS140775/v2aSchwannAll_PS/"
  
  file <- "Annotated_v2aSchwannAll_PS_WS140775_24018518_AnnaKARENINA.xlsx"
  
  filepath <- paste0(s_drive_path, ws_path, file)
  
  expected_df <- data.frame(
    "worksheet" = c("WS140775"),
    "labno" = c("24018518"),
    "suffix" = c(""),
    "patient_name" = c("AnnaKARENINA"),
    "labno_suffix" = c("24018518"),
    "labno_suffix_worksheet" = c("24018518_WS140775")
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
