# PanSolidv2 Audit Shiny

# Packages --------------------------------------------------------------------------

library(shiny)
library(tidyverse)

# Load Data -------------------------------------------------------------------------

amp_gene_results <- read_csv(here::here("data/live_service_collated_data/live_service_amp_gene_results_collated.csv"))

std_dev_results <- read_csv(here::here("data/live_service_collated_data/live_service_std_dev_results_collated.csv"))

percent_138_results <- read_csv(here::here("data/live_service_collated_data/live_service_percent_138_results_collated.csv"))

pos_cnv_results <- read_csv(here::here("data/live_service_collated_data/live_service_pos_cnv_results_collated.csv"))

pos_cnv_results_with_qc <- pos_cnv_results |> 
  left_join(std_dev_results |> 
              select(filepath, st_dev_signal_adjusted_log2_ratios),
            by = "filepath") |> 
  left_join(percent_138_results |> 
              select(filepath, percent_whole_panel_covered_at_138x),
            by = "filepath")

num_samples <- length(unique(amp_gene_results$filepath))

gene_options <- unique(amp_gene_results$gene)

ws_options <- unique(std_dev_results$worksheet)

# User Interface --------------------------------------------------------------------

ui <- fluidPage(
  
  titlePanel("Welcome to the PanSolidv2 audit Shiny!"),
  
  str_c("There are ", num_samples, " files in the dataset for this audit."),

  fluidRow(
    column(10,
           plotOutput("qc_plot")
    )
  ),
  
  fluidRow(
    column(5,
           "Amplifications detected",
           
           selectInput(inputId = "gene_select", label = "Gene", 
              choices = gene_options, selected = "ERBB2"),
  
          numericInput("noise_qc", "Signal-adjusted noise must be less than", 
               value = 1, min = 0, max = 6),
  
  tableOutput("amp_tbl")
  ),
  column(5,
         selectInput(inputId = "ws_select", label = "Worksheet",
                     choices = ws_options, selected = "WS140721"),
         
         plotOutput("ws_noise_plot")
         )
  )
)

# Server ----------------------------------------------------------------------------

server <- function(input, output) {

  data_filtered <- reactive({
    
    pos_cnv_results_with_qc |> 
      filter(gene == input$gene_select & 
               st_dev_signal_adjusted_log2_ratios < input$noise_qc) |> 
      select(worksheet, labno, patient_name, 
             gene, fold_change)

  })
  
  noise_data_filtered <- reactive({
    
    std_dev_results |> 
      filter(worksheet == input$ws_select)
    
  })
  
  
  output$amp_tbl <- renderTable({
    
    data_filtered() 
    
  })
  
  output$qc_plot <- renderPlot({
    
    ggplot(std_dev_results, aes(x = worksheet, y = st_dev_signal_adjusted_log2_ratios)) +
      geom_boxplot() +
      theme_bw() +
      labs(x = "", y = "Signal-adjusted noise")
    
  })
  
  output$ws_noise_plot <- renderPlot({
    
    ggplot(noise_data_filtered(), aes(x = st_dev_signal_adjusted_log2_ratios,
                                    y = )) +
      geom_histogram(binwidth = 0.1) +
      theme_bw() +
      labs(x = "Signal-adjusted noise",
           y = "Number of samples") +
      ylim(0, 40) +
      xlim(0, 2)
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
