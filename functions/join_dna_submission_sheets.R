library(here)
source(here("scripts/set_shared_drive_filepath.R"))

join_pansolid_submission_sheets <- function() {
  
  #' Load and join PanSolid DNA submission sheets
  #'
  #' @return A tidy dataframe for all PanSolid submissions from 2022 to 2024, based on 
  #' files saved in the local drive.
  #' 
  #' @note Excel spreadsheets for coordinating testing on the PanSolid QIAseq 
  #' enrichment are manually curated by the tech team and stored either on the S
  #' drive or on the laboratory Sharepoint. This function uses versions of 
  #' these Excel spreadsheets copied onto a local drive.
  #' 
  #' The sample_id field has some surprises. Example: 23024772 has a degree sign (°)
  #' entered after it which is invisible in Excel and R.
  #'
  #' @examples pansolid_sheets <- join_pansolid_submission_sheets()
  #' 
  #' sample_info <- pansolid_sheets |> filter(labno == "12345678")
  
  pansolid_submission_2023 <- read_excel(path = paste0(data_folder, 
                                                       "excel_spreadsheets/",
                                                       "DNA PanSolid QIAseq Submission Sheet 2023.xlsx")) |> 
    janitor::clean_names() |> 
    rename(stock_qubit = stock_qubit_ng_m_l) |> 
    mutate(submission_sheet = "2023",
           labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
    select(date_submitted, labno, sample_name,
           panel, enrichment, stock_qubit, submission_sheet)
  
  pansolid_submission_2024 <- read_excel(path = paste0(data_folder, 
                                                       "excel_spreadsheets/",
                                                       "PanSolid Submission sheet 2024.xlsx"),
                                         sheet = "PanSolid samples") |> 
    janitor::clean_names()  |> 
    rename(stock_qubit = stock_qubit_ng_m_l) |> 
    mutate(submission_sheet = "2024",
           labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
    select(date_submitted, labno, sample_name,
           panel, enrichment, stock_qubit, submission_sheet)
  
  # Pansolid began in 2022 so the initial runs were recorded on the Qiaseq spreadsheet
  pansolid_submission_2022 <- read_excel(path = paste0(data_folder, 
                                                       "excel_spreadsheets/",
                                                       "QIAseq DNA PanSolid Sample Submission 2022.xlsx")) |> 
    janitor::clean_names() |> 
    rename(date_submitted = date_sample_submitted,
           stock_qubit = stock_qubit_ng_m_l) |> 
    mutate(submission_sheet = "2022",
           labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
    select(date_submitted, labno, sample_name,
           panel, enrichment, stock_qubit, submission_sheet)
  
  output <- rbind(pansolid_submission_2024,
                  pansolid_submission_2023,
                  pansolid_submission_2022)
  
  return(output)
  
}

join_qiaseq_core_submission_sheets <- function() {
  
  core_qiaseq_submission_2023 <- readxl::read_excel(path = paste0(data_folder,
                                                                "excel_spreadsheets/",
                                                                "DNA QIAseq Submission Sheet 2023.xlsx"),
                                                  sheet = "QIAseq samples") |> 
  janitor::clean_names() |> 
  mutate(submission_sheet = "2023",
         worksheet = paste0("WS", qi_aseq),
         labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
  filter(!is.na(labno)) |> 
  select(date_submitted, labno, sample_name,
         panel, enrichment, stock_qubit, worksheet, submission_sheet)

  core_qiaseq_submission_2024 <- readxl::read_excel(path = paste0(data_folder,
                                                                  "excel_spreadsheets/",
                                                                  "DNA QIAseq Submission Sheet 2024.xlsx"),
                                                    sheet = "QIAseq samples") |> 
    janitor::clean_names() |> 
    rename(stock_qubit = stock_qubit_ng_m_l) |> 
    mutate(submission_sheet = "2024",
           worksheet = paste0("WS", qi_aseq_worksheet),
           labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
    filter(!is.na(labno)) |> 
    select(date_submitted, labno, sample_name,
           panel, enrichment, stock_qubit, worksheet, submission_sheet)
  
  output <- rbind(core_qiaseq_submission_2023,
                  core_qiaseq_submission_2024)
  
  return(output)
  
}
