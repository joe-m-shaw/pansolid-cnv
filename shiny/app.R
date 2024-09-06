# PanSolidv2 Audit Shiny

# Packages --------------------------------------------------------------------------

library(shiny)
library(tidyverse)
library(here)

# Filepath --------------------------------------------------------------------------

source(here("scripts/set_shared_drive_filepath.R"))

pansolidv2_excel_regex <- "^Annotated(_|_v2.+_)WS\\d{6}_.+.xlsx"

# Load Data -------------------------------------------------------------------------

amp_gene_results <- read_csv(str_c(data_folder, "live_service/collated/",
                                   "/live_service_amp_gene_results_collated.csv"),
                             show_col_types = FALSE) |> 
  mutate(filename = str_extract(string = filepath, 
                         pattern = str_replace(string = pansolidv2_excel_regex, 
                                               pattern = "\\^", 
                                               replacement = "")))

std_dev_results <- read_csv(str_c(data_folder, "live_service/collated/",
                                  "/live_service_std_dev_results_collated.csv"),
                            show_col_types = FALSE) |> 
  rename(noise = st_dev_signal_adjusted_log2_ratios)

percent_138_results <- read_csv(str_c(data_folder, "live_service/collated/",
                                      "/live_service_percent_138_results_collated.csv"),
                                show_col_types = FALSE)

pos_cnv_results <- read_csv(str_c(data_folder, "live_service/collated/",
                                  "/live_service_pos_cnv_results_collated.csv"),
                            show_col_types = FALSE) |> 
  mutate(gene = factor(gene, levels = unique(amp_gene_results$gene)),
         filename = str_extract(string = filepath, 
                                pattern = str_replace(string = pansolidv2_excel_regex, 
                                                      pattern = "\\^", 
                                                      replacement = "")))

panel_info <- read_csv(str_c(data_folder, "live_service/collated/",
                             "/pansolidv2_sample_worksheet_panel_information.csv"),
                       show_col_types = FALSE)
  
pos_cnv_results_with_qc <- pos_cnv_results |> 
  left_join(std_dev_results |> 
              select(filepath, noise),
            by = "filepath") |> 
  left_join(percent_138_results |> 
              select(filepath, percent_whole_panel_covered_at_138x),
            by = "filepath") |> 
  left_join(panel_info |> 
              select(filename, panel),
            by = "filename",
            # Some samples are run on multiple panels so have multiple rows
            relationship = "many-to-many")

fold_change_threshold <- 2.8

crc_erbb2_results <- amp_gene_results |> 
  left_join(panel_info |>
              select(filename, panel),
            by = "filename") |> 
  filter(panel == "v2M1_CRC_PS" & gene == "ERBB2") |> 
  mutate(result = case_when(
    max_region_fold_change >= fold_change_threshold ~"ERBB2 Amplification",
    max_region_fold_change < fold_change_threshold ~"No ERBB2 amplification"
  ))

qc_joined <- std_dev_results |> 
  left_join(percent_138_results |> 
              select(filepath, percent_whole_panel_covered_at_138x),
            by = "filepath")

num_samples <- length(unique(amp_gene_results$filepath))

gene_options <- unique(amp_gene_results$gene)

ws_options <- unique(std_dev_results$worksheet)

weeks_live <- round(as.numeric(difftime(time1 = today(), 
                                        time2 = "2024-04-08", 
                                        units = "weeks")), 1)

total_samples <- length(unique(std_dev_results$filepath))

samples_per_week <- round(total_samples / weeks_live, 0)

worksheet_quality_outcomes <- std_dev_results |> 
  mutate(quality_category = case_when(
    noise >= 1 ~"poor",
    noise < 1 & noise >= 0.7 ~"sub-optimal",
    noise < 0.7 ~"good"
  ),
  quality_category = factor(quality_category, levels = c("poor", "sub-optimal", "good"))) |> 
  group_by(quality_category, worksheet) |> 
  summarise(total = n())
  
# User Interface --------------------------------------------------------------------

ui <- fluidPage(
  
  theme = bslib::bs_theme(bootswatch = "sandstone"),
  
  titlePanel("PanSolidv2 CNV Audit"),
  
  str_c("Date: ", format(today(), "%d-%m-%Y")),
  
  str_c("Total samples: ", total_samples),
  
  str_c("Weeks live: ", weeks_live),
  
  str_c("Samples per week: ", samples_per_week),

  fluidRow(
    h1("Summary"),
    
    column(3,
           h2("All amplification results"),
           tableOutput("summary_gene_table")),
    column(3,
           h2("QIAseq PanSolid Colorectal Results"),
           tableOutput("crc_summary_table"))
  ),
  
  fluidRow(
    h1("Quality Tracking"),
    
    column(10,
           plotOutput("noise_by_ws_plot")),
    
    column(10,
           plotOutput("percent_138_plot")),
    
    column(10,
           plotOutput("noise_category_by_ws_plot"))

  ),
  
  fluidRow(
    h1("Interactive plots and tables"),
    
    column(5,
           checkboxGroupInput(inputId = "gene_multiple", 
                              label = "Genes to plot", 
                              choices = gene_options,
                              selected = gene_options, inline = TRUE),
           sliderInput("y_axis_range", "Fold change range", min = -5, 
                       max = 200,
                       value = c(-5, 200)),
           plotOutput("gene_amp_plot")),
    
    column(5,
           selectInput(inputId = "ws_select", label = "Worksheet to highlight",
                       choices = ws_options, selected = "WS140721"),
           
           plotOutput("ws_noise_plot"))
    
    
  ),
  fluidRow(
    
    column(5,
           h2("Sample details"),
           
           selectInput(inputId = "gene_select", label = "Gene", 
                       choices = gene_options, selected = "ERBB2"),
           
           numericInput("noise_qc", "Signal-adjusted noise must be less than", 
                        value = 1, min = 0, max = 6),
           
           tableOutput("amp_tbl")
    )
  ),
  fluidRow(
    column(3,
           h2("Panels"),
           tableOutput("summary_panel_table")),
  )
)

# Server ----------------------------------------------------------------------------

server <- function(input, output) {
  
  output$summary_gene_table <- renderTable({
    
    pos_cnv_results_with_qc |> 
      filter(gene != "no positive calls") |> 
      count(gene, .drop = FALSE) |> 
      arrange(desc(n)) |> 
      rename(Cases = n,
             `Amplification Result` = gene)
    
  })

  output$summary_panel_table <- renderTable({
    panel_info |> 
      count(panel) |> 
      arrange(desc(n)) |> 
      rename("Panel" = panel,
             "Cases" = n)
  })
  
  output$crc_summary_table <- renderTable({
    
    crc_erbb2_results |> 
      count(result) |> 
      arrange(desc(n)) |> 
      mutate(Percentage = round(n/sum(n) * 100, 1)) |> 
      rename(Cases = n,
             `Amplification Result` = result)
    
  })
  
  data_filtered <- reactive({
    
    pos_cnv_results_with_qc |> 
      filter(!is.na(fold_change)) |> 
      filter(gene == input$gene_select & 
               noise < input$noise_qc) |> 
      select(worksheet, panel, labno,  
             fold_change, noise)

  })
  
  all_amp_data_filtered <- reactive({
    
    amp_gene_results |> 
      filter(gene %in% input$gene_multiple)
    
  })
  
  noise_data_filtered <- reactive({
    
    qc_joined |> 
      filter(worksheet == input$ws_select)
    
  })
  
  output$gene_amp_plot <- renderPlot({
    
    ggplot(all_amp_data_filtered(), aes(x = gene, y = max_region_fold_change)) +
      geom_jitter(pch = 21, size = 2) +
      geom_hline(yintercept = 2.8, linetype = "dashed") +
      theme_bw() +
      ylim(input$y_axis_range) +
      labs(y = "Fold change", x = "")
    
  })
  
  output$amp_tbl <- renderTable({
    
    data_filtered() 
    
  })
  
  output$qc_plot <- renderPlot({
    
    ggplot(std_dev_results, aes(x = worksheet, y = noise)) +
      geom_boxplot() +
      theme_bw() +
      labs(x = "", y = "Signal-adjusted noise")
    
  })
  
  output$ws_noise_plot <- renderPlot({
    
    ggplot(qc_joined, aes(x = noise,
                                    y = percent_whole_panel_covered_at_138x)) +
      geom_point(shape = 21, size = 3, alpha = 0.5) +
      theme_bw() +
      labs(x = "Signal-adjusted noise",
           y = "Percent panel covered > 138X") +
      ylim(0, 100) +
      xlim(0, max(qc_joined$noise) + 0.1) +
      geom_point(data = noise_data_filtered(), size = 3,
                 shape = 21, fill = "red")
    
  })
  
  output$noise_by_ws_plot <- renderPlot({
    
    ggplot(std_dev_results, aes(x = worksheet, y = noise)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "", y = "Signal-adjusted noise")
    
  })
  
  output$percent_138_plot <- renderPlot({
    
    ggplot(percent_138_results, aes(x = worksheet, 
                                    y = percent_whole_panel_covered_at_138x)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90)) +
      labs(x = "", y = "Percentage panel covered at 138x")
    
  })
  
  output$noise_category_by_ws_plot <- renderPlot({
    
    ggplot(worksheet_quality_outcomes, aes(x = worksheet, y = total)) +
      geom_col(aes(fill = quality_category)) +
      scale_fill_manual(values = c("#D55E00", "#E69F00","#009E73")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90),
            legend.position = "bottom") +
      labs(x = "Worksheet", y = "Number of samples",
           title = "Quality categories of noise by worksheet",
           subtitle = "Sub-optimal noise threshold: 0.7; poor noise threshold: 1",
           fill = "")
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
