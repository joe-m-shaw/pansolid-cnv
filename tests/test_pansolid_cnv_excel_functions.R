library(testthat)

test_datapath <- paste0(config::get("data_folderpath"),
                        "validation/DOC6567_deletions/test_data/")

test_filepath <- paste0(test_datapath,
                        "Annotated_v2PANSOLID_WS143415_24030946_PierreBEZUKHOV.xlsx")

test_sheet <- read_cnv_sheet(test_filepath)

old_format_path <- paste0(test_datapath,
                          "Annotated_v2PANSOLID_WS143574_24030966_NatashaROSTOVA.xlsx")

old_format_sheet <- read_cnv_sheet(old_format_path)

# get_sheetnames

test_that("get_sheetnames works with standard input", {
  
  expect_equal(get_sheetname(filepath = test_filepath),
               "CNVs_24030946")
  
})

test_that("get_sheetnames fails when sheet name is absent", {

  expect_error(get_sheetname(filepath = test_filepath, 
                             sheet_regex = "new_sheet"),
               "Sheet name not found")
  
})

test_that("get_sheetnames fails when sheet name is not unique", {
  
  expect_error(get_sheetname(filepath = test_filepath, 
                             sheet_regex = "Hotspots"),
               "Sheet name must be unique")
  
})

# find_match

test_that("find_match fails when the string is wrong", {
  
  test_sheet_missing_titles <- read_cnv_sheet(paste0(test_datapath,
                                                     "Annotated_v2PANSOLID_WS143415_24030946_PlatonKARATAEV.xlsx"))
  
  expect_error(find_match(test_sheet_missing_titles, "a",
                          "StDev Signal-adjusted Log2 Ratios"),
               "StDev Signal-adjusted Log2 Ratios not found")
  
  
})

# find_stdev_ratios

test_that("find_stdev_ratios works with standard format", {
  
  df_expected <- tibble::tibble(
    "stdev_noise" = c(0.44043682299167)
  )
  
  expect_equal(find_stdev_ratios(test_sheet), 
               df_expected)
  
})

# find_percent_138x

test_that("find_percent_138x works with standard format", {
  
  df_expected <- tibble::tibble(
    "percent_138x" = c(98.1529430013483)
  )
  
  expect_equal(find_percent_138x(test_sheet),
               df_expected)
  
})

# find_pred_ncc

test_that("find_pred_ncc works with standard format", {
  
  df_expected <- tibble::tibble(
    "pred_ncc" = c(93)
  )
  expect_equal(find_pred_ncc(test_sheet),
               df_expected)
  
})

# find_sig_cnvs

test_that("find_sig_cnvs handles multiple significant cnvs", {
  
  tribble_expected <- tibble::tribble(
   
    ~gene,	~chromosome, ~cnv_co_ordinates, ~cnv_length, ~consequence, ~fold_change, ~p_value, ~no_targets, ~check_1,          ~check_2,           ~copy_number,  ~start,   ~end,
    "EGFR",	  "7",	   "54185293..55205622",	  1020330,	"Gain",	      111.75,	     0,           	42,			as.character(NA),  as.character(NA),  as.numeric(NA), 54185293, 55205622,
    "FH",	    "1",	   "241497823..247750014",	6252192,	"Loss",	      -1.45,	     1.03741E-09, 	17,			as.character(NA),  as.character(NA),  as.numeric(NA), 241497823, 247750014,
    "CDKN2A",	"9",	   "21968224..22008958",	  40735,	  "Loss",	      -14.05,	     0,	            16,			as.character(NA),  as.character(NA),  as.numeric(NA), 21968224, 22008958,
    "CDKN2B",	"9",	   "21968224..22008958",	  40735,	  "Loss",	      -14.05,	     0,	            16,			as.character(NA),  as.character(NA),  as.numeric(NA), 21968224, 22008958,
    "BMPR1A",	"10",	   "87329..86923724",	    86836396,	  "Loss",	      -1.75,	     0,	            123,		as.character(NA),  as.character(NA),  as.numeric(NA), 87329,    86923724,
    "PTEN",	  "10",	   "87864465..87965477",	  101013,	  "Loss",	      -15.36,	     0,	            9,			as.character(NA),  as.character(NA),  as.numeric(NA), 87864465, 87965477,
    "SUFU",	  "10",	   "90436067..132677850",	42241784,	  "Loss",	      -1.87,	     0,	            75,			as.character(NA),  as.character(NA),  as.numeric(NA), 90436067, 132677850,
    "RB1",	  "13",	   "48367489..48465373",	  97885,	  "Loss",	      -1.26,	     0,	            14,			as.character(NA),  as.character(NA),  as.numeric(NA), 48367489, 48465373,
    "GREM1",	"15",	   "28928393..54918513",	  25990121,	"Loss",	      -1.81,	     0,	            74,			as.character(NA),  as.character(NA),  as.numeric(NA), 28928393, 54918513
    
    )
  
  expect_equal(find_sig_cnvs(test_sheet),
               tribble_expected)
    
})

test_that("find_sig_cnvs handles no significant CNVs",{
  
  no_sig_cnv_path <- paste0(test_datapath, 
                           "Annotated_v2aSchwannAll_PS_WS140775_24018518_AnnaKARENINA.xlsx")
  
  no_sig_cnv_sheet <- read_cnv_sheet(no_sig_cnv_path)
  
  df_expected <- tibble::tibble(
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
                            "Annotated_v2M119_PAEDEXT_PS_WS140815_24018922_KonstantinLEVIN.xlsx")
  
  one_sig_cnv_sheet <- read_cnv_sheet(one_sig_cnv_path)
  
  df_expected <- tibble::tibble(
    "gene" = "SMAD4",
    "chromosome" = "18",
    "cnv_co_ordinates" = "28984159..51067192",
    "cnv_length" = 22083034,
    "consequence" = "Loss",
    "fold_change" = -1.49,
    "p_value" = 1.87627691161652E-13,
    "no_targets" = 29,
    "check_1" = as.character(NA),
    "check_2" = as.character(NA),
    "copy_number" = as.numeric(NA),
    "start" = 28984159,
    "end" = 51067192)
  
  expect_equal(find_sig_cnvs(one_sig_cnv_sheet),
               df_expected)
  
})

# find_amp_genes

test_that("find_amp_genes loads table correctly", {
  
  tribble_expected <- tibble::tribble(
    ~gene,   ~max_region_fold_change, ~min_region_fold_change,
    "ALK",	1.059691798,	1.059691798,
    "BRAF",	1.057454359,	1.057454359,
    "CDK12",	-1.013880552,	-1.013880552,
    "EGFR",	111.7454691,	111.7454691,
    "ERBB2",	1.053855281,	-1.225449632,
    "FGFR1",	1.121534495,	1.121534495,
    "FGFR2",	-1.867858861,	-1.867858861,
    "FGFR3",	1.087045245,	1.087045245,
    "MDM2",	1.118774909,	1.118774909,
    "MET",	1.631485437,	1.631485437,
    "MYC",	-1.072607912,	-1.072607912,
    "MYCN",	1.075771694,	1.075771694,
    "SMO",	1.057454359,	1.057454359
    
    )
    
  expect_equal(find_amp_genes(test_sheet, num_genes = 13), 
               tribble_expected)
  
})

test_that("find_amp_genes loads with old gene list", {
  
  tribble_expected <- tibble::tribble(
    ~gene, ~max_region_fold_change, ~min_region_fold_change,
    "ALK",	-1.162531526,	-1.162531526,
    "BRAF",	1.628365713,	1.628365713,
    "EGFR",	21.77182049,	3.718416458,
    "ERBB2",	1.09505526,	-1.258704461,
    "MDM2",	1.187791379,	1.187791379,
    "MET",	1.628365713,	1.628365713,
    "MYC",	-1.005257589,	-1.005257589,
    "MYCN",	1.149774419,	1.149774419,
    "SMO",	1.628365713,	1.628365713,
  )
  
  expect_equal(find_amp_genes(old_format_sheet),
               tribble_expected)
  
})

# find_del_genes

test_that("find_del_genes loads table correctly", {
  
  tribble_expected <- tibble::tribble(
    ~gene, ~max_region_fold_change, ~min_region_fold_change,
    "APC",	1.312340678,	1.312340678,
    "ATM",	1.161092492,	1.047549833,
    "BMPR1A",	-1.746078662,	-1.746078662,
    "BRCA1",	1.053855281,	1.053855281,
    "BRCA2",	1.077306237,	1.077306237,
    "BRIP1",	1.053855281,	1.053855281,
    "CDH1",	-1.149986369,	-1.149986369,
    "CDK12",	-1.013880552,	-1.013880552,
    "CDKN2A",	-14.04557779,	-14.04557779,
    "CDKN2B",	-14.04557779,	-14.04557779,
    "CHEK2",	-1.038818917,	-1.038818917,
    "DICER1",	1.016807806,	1.016807806,
    "FH",	-1.454001098,	-1.454001098,
    "GREM1",	-1.812118953,	-1.812118953,
    "MLH1",	1.0067468,	1.0067468,
    "MSH2",	1.059691798,	1.059691798,
    "MSH6",	1.059691798,	1.059691798,
    "MUTYH",	-1.107864821,	-1.107864821,
    "NF1",	-1.013880552,	-1.013880552,
    "NF2",	-1.110092924,	-1.110092924,
    "NTHL1",	-1.030602896,	-1.030602896,
    "PALB2",	-1.149986369,	-1.149986369,
    "PMS2",	1.434916931,	1.434916931,
    "PTCH1",	-1.045946119,	-1.045946119,
    "PTEN",	-15.3636289,	-15.3636289,
    "RAD51C",	1.057301334,	1.057301334,
    "RAD51D",	-1.013880552,	-1.013880552,
    "RB1",	1.163907907,	-1.256292757,
    "RNF43",	1.057301334,	1.057301334,
    "SMAD4",	1.116294864,	1.116294864,
    "SMARCA4",	1.39451318,	1.39451318,
    "STK11",	1.393082686,	1.393082686,
    "SUFU",	-1.867858861,	-1.867858861,
    "TP53",	1.371084737,	-1.013880552
  )
  
  expect_equal(find_del_genes(test_sheet,
                              num_genes = 34),
               tribble_expected)
  
})

test_that("find_del_genes throws warning when gene number is wrong", {
  
  expect_warning(find_del_genes(test_sheet,
                                num_genes = 37),
                 "Different number of genes found to expected")
  
})

test_that("find_del_genes works with old gene list", {
  
  tribble_expected <- tibble::tribble(
    ~gene, ~max_region_fold_change, ~min_region_fold_change, 
    "APC",	1.24689262,	1.24689262,
    "ATM",	1.173120823,	1.173120823,
    "BMPR1A",	1.163514339,	1.163514339,
    "BRCA1",	1.464356208,	1.464356208,
    "BRCA2",	1.221956278,	1.221956278,
    "BRIP1",	1.229714675,	1.229714675,
    "CDH1",	-1.031850967,	-1.031850967,
    "CDK12",	1.09505526,	1.09505526,
    "CDKN2A",	-1.567332902,	-1.567332902,
    "CDKN2B",	-13.60924907,	-13.60924907,
    "CHEK2",	1.163849616,	-1.023158645,
    "DICER1",	1.100415838,	1.100415838,
    "FGFR1",	1.206991917,	-1.221353538,
    "FGFR2",	-2.049500436,	-2.049500436,
    "FGFR3",	-1.231258367,	-1.231258367,
    "FH",	1.114605365,	1.114605365,
    "GREM1",	-1.189153779,	-1.189153779,
    "MLH1",	-1.019437087,	-1.019437087,
    "MSH2",	1.230490141,	-1.162531526,
    "MSH6",	1.230490141,	1.230490141,
    "MUTYH",	-1.242798691,	-1.242798691,
    "NF1",	1.205077509,	1.205077509,
    "NF2",	-1.023158645,	-1.023158645,
    "NTHL1",	-1.025359579,	-1.025359579,
    "PALB2",	-1.031850967,	-1.031850967,
    "PMS2",	1.646078388,	1.646078388,
    "PTCH1",	-1.004195022,	-1.004195022,
    "PTEN",	1.06189978,	1.06189978,
    "RAD51C",	1.220344713,	1.220344713,
    "RAD51D",	-1.229422986,	-1.229422986,
    "RB1",	1.221956278,	1.221956278,
    "RNF43",	-1.489376831,	-1.489376831,
    "SMAD4",	1.126827043,	1.126827043,
    "SMARCA4",	-1.215548051,	-1.215548051,
    "STK11",	1.386967793,	1.386967793,
    "SUFU",	-1.034262707,	-1.034262707,
    "TP53",	2.127024825,	2.127024825
  )
  
  expect_equal(find_del_genes(old_format_sheet),
               tribble_expected)
  
})

# read_loh_table

test_that("read_loh_table loads table correctly", {
  
  tribble_expected <- tibble::tribble(
    ~chrom,	~gene,	~ploidy_state,	~loh_status,	~no_targets_in_ploidy_region, ~check_1,         ~check_2,
    "2",	"MSH2",	"Normal diploid",	    "No",	     "162",                      as.character(NA),  as.character(NA),
    "2",	"MSH6",	"Normal diploid",	    "No",	     "162",                      as.character(NA),  as.character(NA),
    "3",	"MLH1",	"Normal diploid",	    "No",	     "247",                      as.character(NA),  as.character(NA),
    "7",	"PMS2",	"Duplication",	      "No",	     "59",                       as.character(NA),  as.character(NA),
    "22",	"LZTR1",	"Normal diploid", 	"No",	     "182",                      as.character(NA),  as.character(NA),
    "22",	"SMARCB1",	"Normal diploid",	"No",	     "182",                      as.character(NA),  as.character(NA),
    "22",	"NF2",	"Normal diploid",	    "No",	     "182",                      as.character(NA),  as.character(NA)
  )
  
  expect_equal(read_loh_table(test_filepath)[,8:14],
               tribble_expected)
  
})

rm(test_datapath, test_filepath, test_sheet, old_format_path, 
   old_format_sheet)