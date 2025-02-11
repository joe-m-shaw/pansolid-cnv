
library(testthat)

test_datapath <- paste0(config::get("data_filepath"),
                        "validation/DOC6567_deletions/test_data/")

test_filepath <- paste0(test_datapath,
                        "Annotated_Jan2025_withLOH_testing_del_visualisation_WS123456_12345678_PierreBEZUKHOV.xlsx")

test_sheet <- read_cnv_sheet(test_filepath)

# get_sheetnames

test_that("get_sheetnames works with standard input", {
  
  expect_equal(get_sheetname(test_filepath),
               "Amplifications_12345678")
  
})

test_that("get_sheetnames fails when sheet name is absent", {

  expect_error(get_sheetname(test_filepath, sheet_regex = "new_sheet"),
               "Sheet name not found")
  
})

test_that("get_sheetnames fail when sheet name is not unique", {
  
  expect_error(get_sheetname(test_filepath, sheet_regex = "Hotspots"),
               "Sheet name must be unique")
  
})

# find_match

test_that("find_match fails when the string is wrong", {
  
  test_sheet_missing_titles <- read_cnv_sheet(paste0(test_datapath,
                                                     "Annotated_Jan2025_withLOH_testing_del_visualisation_WS123456_12345678_PlatonKARATAEV.xlsx"))
  
  expect_error(find_match(test_sheet_missing_titles, "a",
                          "StDev Signal-adjusted Log2 Ratios"),
               "StDev Signal-adjusted Log2 Ratios not found")
  
  
})

# find_stdev_ratios

test_that("find_stdev_ratios works with standard format", {
  
  df_expected <- tibble(
    "stdev_noise" = c(0.300228451587477)
  )
  
  expect_equal(find_stdev_ratios(test_sheet), 
               df_expected)
  
})

# find_percent_138x

test_that("find_percent_138x works with standard format", {
  
  df_expected <- tibble(
    "percent_138x" = c(99.2912455702848)
  )
  
  expect_equal(find_percent_138x(test_sheet),
               df_expected)
  
})

# find_sig_cnvs

test_that("find_sig_cnvs handles multiple significant cnvs", {
  
  tribble_expected <- tibble::tribble(
    ~gene, ~chromosome, ~cnv_co_ordinates, ~cnv_length, ~consequence, ~fold_change, ~p_value, ~no_targets, ~check_1, ~check_2, ~copy_number, ~start, ~end,
    "MLH1",	"3",	"10094284..39962881",	29868598,	"Loss",	-1.69,	0,	112, as.character(NA), as.character(NA), as.numeric(NA), 10094284, 39962881,
    "FGFR3",	"4",	"72152..47007452",	46935301,	"Loss",	-1.47,	0,	64, as.character(NA), as.character(NA), as.numeric(NA), 72152, 47007452,
    "CDKN2A",	"9",	"21242158..28208471",	6966314,	"Loss",	-20.83,	0,	24, as.character(NA), as.character(NA), as.numeric(NA), 21242158, 28208471,
    "CDKN2B",	"9",	"21242158..23207860",	1965703,	"Loss",	-22.29,	0,	19, as.character(NA), as.character(NA), as.numeric(NA), 21242158, 23207860)
  
  expect_equal(find_sig_cnvs(test_sheet),
               tribble_expected)
    
})

test_that("find_sig_cnvs handles no significant CNVs",{
  
  no_sig_cnv_path <- paste0(test_datapath, 
                           "Annotated_Jan2025_withLOH_testing_del_visualisation_WS123456_12345678_AnnaKARENINA.xlsx")
  
  no_sig_cnv_sheet <- read_cnv_sheet(no_sig_cnv_path)
  
  df_expected <- tibble(
    "gene" = "no positive calls",
    "chromosome" = "",
    "cnv_co_ordinates" = "",
    "cnv_length" = 0,
    "consequence" = "no call",
    "fold_change" = 0,
    "p_value" = 0,
    "no_targets" = 0,
    "check_1" = "",
    "check_2" = "",
    "copy_number" = 0,
    "start" = as.numeric(NA),
    "end" = as.numeric(NA))
  
  expect_equal(find_sig_cnvs(no_sig_cnv_sheet),
               df_expected)
  
})

test_that("find_sig_cnvs handles one significant CNV", {
  
  one_sig_cnv_path <- paste0(test_datapath, 
                            "Annotated_Jan2025_withLOH_testing_del_visualisation_WS123456_12345678_KonstantinLEVIN.xlsx")
  
  one_sig_cnv_sheet <- read_cnv_sheet(one_sig_cnv_path)
  
  df_expected <- tibble(
    "gene" = "MLH1",
    "chromosome" = "3",
    "cnv_co_ordinates" = "10094284..39962881",
    "cnv_length" = 29868598,
    "consequence" = "Loss",
    "fold_change" = -1.69,
    "p_value" = 0,
    "no_targets" = 112,
    "check_1" = as.character(NA),
    "check_2" = as.character(NA),
    "copy_number" = as.numeric(NA),
    "start" = 10094284,
    "end" = 39962881)
  
  expect_equal(find_sig_cnvs(one_sig_cnv_sheet),
               df_expected)
  
})

# find_amp_genes

test_that("find_amp_genes loads table correctly", {
  
  tribble_expected <- tibble::tribble(
    ~gene,   ~max_region_fold_change, ~min_region_fold_change,
    "ALK",	 1.108169324,	 1.108169324,
    "BRAF",	 1.45630905,	 1.45630905,
    "EGFR",	 1.101747319,	 1.101747319,
    "ERBB2", 1.090836841,	 1.090836841,
    "MDM2",	 1.357472306,	 1.114555468,
    "MET",	 1.45630905,	 1.45630905,
    "MYC",	 1.741610741,	 1.741610741,
    "MYCN",	 -1.06841134,	 -1.06841134,
    "SMO",	 1.45630905,	 1.45630905)
    
  expect_equal(find_amp_genes(test_sheet), 
               tribble_expected)
  
})

# find_del_genes

test_that("find_del_genes loads table correctly", {
  
  tribble_expected <- tibble::tribble(
    ~gene, ~max_region_fold_change, ~min_region_fold_change,
    "APC", -1.116238365, -1.116238365,
    "ATM", -1.004510532, -1.004510532,
    "BMPR1A", -1.107152631, -1.107152631,
    "BRCA1", 1.090836841, 1.090836841,
    "BRCA2", -1.594467871, -1.594467871,
    "BRIP1", 1.616280083, 1.616280083,
    "CDH1", 1.0397996, 1.0397996,
    "CDK12", 1.090836841, 1.090836841,
    "CDKN2A", -20.82574205, -20.82574205,
    "CDKN2B", -22.28795334, -22.28795334,
    "CHEK2", 1.495115834, 1.495115834,
    "DICER1", 1.208210707, 1.208210707,
    "FGFR1", 1.410163351, 1.410163351,
    "FGFR2", 1.509195128, 1.509195128,
    "FGFR3", -1.465906486, -1.465906486,
    "FH", 1.323745509, 1.323745509,
    "GREM1", -1.210006223, -1.210006223,
    "MLH1", -1.692864938, -1.692864938,
    "MSH2", 1.446246255, 1.174131547,
    "MSH6", 1.174131547, 1.174131547,
    "MUTYH", 1.020898587, 1.020898587,
    "NF1", 1.090836841, 1.090836841,
    "NF2", -1.273819864, -1.273819864,
    "NTHL1", 1.191047918, 1.191047918,
    "PALB2", 1.0397996, 1.0397996,
    "PMS2", 1.409828034, 1.409828034,
    "PTCH1", -1.186303008, -1.186303008,
    "PTEN", -1.182256748, -1.182256748,
    "RAD51C", 1.097850552, 1.097850552,
    "RAD51D", 1.090836841, 1.090836841,
    "RB1", -1.594467871, -1.594467871,
    "RNF43", 1.097850552, 1.097850552,
    "SMAD4", -1.575295864, -1.575295864,
    "SMARCA4", -1.487787655, -1.487787655,
    "STK11", -1.481281658, -1.481281658,
    "SUFU", -1.107152631, -1.107152631,
    "TP53", 1.064460108, 1.064460108
  )

  expect_equal(find_del_genes(test_sheet),
               tribble_expected)
  
})

rm(test_datapath, test_filepath, test_sheet)