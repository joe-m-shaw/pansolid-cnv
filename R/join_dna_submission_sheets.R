data_folder <- config::get("data_filepath")

spreadsheet_folder <- paste0(data_folder,
                             "validation/DOC6283_amplifications/excel_spreadsheets/")

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

  pansolid_submission_2023 <- readxl::read_excel(path = paste0(spreadsheet_folder, 
                                                       "DNA PanSolid QIAseq Submission Sheet 2023.xlsx")) |> 
    janitor::clean_names() |> 
    dplyr::rename(stock_qubit = stock_qubit_ng_m_l) |> 
    dplyr::mutate(submission_sheet = "2023",
           labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
    dplyr::select(date_submitted, labno, sample_name,
           panel, enrichment, stock_qubit, submission_sheet)
  
  pansolid_submission_2024 <- readxl::read_excel(path = paste0(spreadsheet_folder, 
                                                       "PanSolid Submission sheet 2024.xlsx"),
                                         sheet = "PanSolid samples") |> 
    janitor::clean_names()  |> 
    dplyr::rename(stock_qubit = stock_qubit_ng_m_l) |> 
    dplyr::mutate(submission_sheet = "2024",
           labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
    dplyr::select(date_submitted, labno, sample_name,
           panel, enrichment, stock_qubit, submission_sheet)
  
  # Pansolid began in 2022 so the initial runs were recorded on the Qiaseq spreadsheet
  pansolid_submission_2022 <- readxl::read_excel(path = paste0(spreadsheet_folder, 
                                                       "QIAseq DNA PanSolid Sample Submission 2022.xlsx")) |> 
    janitor::clean_names() |> 
    dplyr::rename(date_submitted = date_sample_submitted,
           stock_qubit = stock_qubit_ng_m_l) |> 
    dplyr::mutate(submission_sheet = "2022",
           labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
    dplyr::select(date_submitted, labno, sample_name,
           panel, enrichment, stock_qubit, submission_sheet)
  
  output <- rbind(pansolid_submission_2024,
                  pansolid_submission_2023,
                  pansolid_submission_2022)
  
  return(output)
  
}

join_qiaseq_core_submission_sheets <- function() {
  
  core_qiaseq_submission_2023 <- readxl::read_excel(path = paste0(spreadsheet_folder,
                                                                "DNA QIAseq Submission Sheet 2023.xlsx"),
                                                  sheet = "QIAseq samples") |> 
  janitor::clean_names() |> 
    dplyr::mutate(submission_sheet = "2023",
         worksheet = paste0("WS", qi_aseq),
         labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
    dplyr::filter(!is.na(labno)) |> 
    dplyr::select(date_submitted, labno, sample_name,
         panel, enrichment, stock_qubit, worksheet, submission_sheet)

  core_qiaseq_submission_2024 <- readxl::read_excel(path = paste0(spreadsheet_folder,
                                                                  "DNA QIAseq Submission Sheet 2024.xlsx"),
                                                    sheet = "QIAseq samples") |> 
    janitor::clean_names() |> 
    dplyr::rename(stock_qubit = stock_qubit_ng_m_l) |> 
    dplyr::mutate(submission_sheet = "2024",
           worksheet = paste0("WS", qi_aseq_worksheet),
           labno = str_extract(string = sample_id, pattern = "\\d{8}")) |> 
    dplyr::filter(!is.na(labno)) |> 
    dplyr::select(date_submitted, labno, sample_name,
           panel, enrichment, stock_qubit, worksheet, submission_sheet)
  
  output <- rbind(core_qiaseq_submission_2023,
                  core_qiaseq_submission_2024)
  
  return(output)
  
}
