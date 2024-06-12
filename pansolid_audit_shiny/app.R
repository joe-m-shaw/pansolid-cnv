# PanSolidv2 Audit Shiny

# Packages --------------------------------------------------------------------------

library(shiny)
library(tidyverse)

source(here::here("functions/cnv_functions.R"))

# Load Data -------------------------------------------------------------------------

amp_gene_results <- read_csv(here::here("data/live_service_collated_data/live_service_amp_gene_results_collated.csv")) |> 
  mutate(filename = str_extract(string = filepath, 
                         pattern = str_replace(string = pansolidv2_excel_regex, 
                                               pattern = "\\^", 
                                               replacement = "")))

std_dev_results <- read_csv(here::here("data/live_service_collated_data/live_service_std_dev_results_collated.csv")) |> 
  rename(noise = st_dev_signal_adjusted_log2_ratios)

percent_138_results <- read_csv(here::here("data/live_service_collated_data/live_service_percent_138_results_collated.csv"))

pos_cnv_results <- read_csv(here::here("data/live_service_collated_data/live_service_pos_cnv_results_collated.csv")) |> 
  mutate(gene = factor(gene, levels = unique(amp_gene_results$gene)),
         filename = str_extract(string = filepath, 
                                pattern = str_replace(string = pansolidv2_excel_regex, 
                                                      pattern = "\\^", 
                                                      replacement = "")))

panel_info <- read_csv(here::here("data/live_service_collated_data/pansolidv2_sample_worksheet_panel_information.csv"))

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
              filter(panel == "v2M1_CRC_PS") |> 
              select(filename, panel),
            by = "filename") |> 
  filter(!is.na(panel) & gene == "ERBB2") |> 
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
    
    column(5,
           h2("All amplification results"),
           tableOutput("summary_gene_table")),
    column(5,
           h2("Colorectal referrals"),
           tableOutput("crc_summary_table")),
  ),
  
  fluidRow(
    h1("Quality Tracking"),
    
    column(10,
           plotOutput("noise_by_ws_plot")),
    
    column(10,
           plotOutput("percent_138_plot"))
    
  ),
  
  fluidRow(
    h1("Interactive plots and tables"),
    
    column(5,
           checkboxGroupInput(inputId = "gene_multiple", 
                              label = "Genes to plot", 
                              choices = gene_options,
                              selected = gene_options, inline = TRUE),
           sliderInput("y_axis_range", "Fold change range", min = -5, 
                       max = max(pos_cnv_results$fold_change + 5),
                       value = c(-5, max(pos_cnv_results$fold_change + 5))),
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

  output$crc_summary_table <- renderTable({
    
    crc_erbb2_results |> 
      count(result) |> 
      arrange(desc(n)) |> 
      rename(Cases = n,
             `Amplification Result` = result)
    
  })
  
  data_filtered <- reactive({
    
    pos_cnv_results_with_qc |> 
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
  
}

# Run the application 
shinyApp(ui = ui, server = server)
