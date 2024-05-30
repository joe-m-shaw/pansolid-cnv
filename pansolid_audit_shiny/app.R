# PanSolidv2 Audit Shiny

# Packages --------------------------------------------------------------------------

library(shiny)
library(tidyverse)

# Load Data -------------------------------------------------------------------------

amp_gene_results <- read_csv(here::here("data/live_service_collated_data/live_service_amp_gene_results_collated.csv"))

std_dev_results <- read_csv(here::here("data/live_service_collated_data/live_service_std_dev_results_collated.csv")) |> 
  rename(noise = st_dev_signal_adjusted_log2_ratios)

percent_138_results <- read_csv(here::here("data/live_service_collated_data/live_service_percent_138_results_collated.csv"))

pos_cnv_results <- read_csv(here::here("data/live_service_collated_data/live_service_pos_cnv_results_collated.csv")) |> 
  mutate(gene = factor(gene, levels = unique(amp_gene_results$gene)))

pos_cnv_results_with_qc <- pos_cnv_results |> 
  left_join(std_dev_results |> 
              select(filepath, noise),
            by = "filepath") |> 
  left_join(percent_138_results |> 
              select(filepath, percent_whole_panel_covered_at_138x),
            by = "filepath")

qc_joined <- std_dev_results |> 
  left_join(percent_138_results |> 
              select(filepath, percent_whole_panel_covered_at_138x),
            by = "filepath")

num_samples <- length(unique(amp_gene_results$filepath))

gene_options <- unique(amp_gene_results$gene)

ws_options <- unique(std_dev_results$worksheet)

# User Interface --------------------------------------------------------------------

ui <- fluidPage(
  
  theme = bslib::bs_theme(bootswatch = "sandstone"),
  
  titlePanel("Welcome to the PanSolidv2 audit Shiny!"),

  fluidRow(
  column(5,
         h2("Fun gene plot"),
         checkboxGroupInput(inputId = "gene_multiple", 
                            label = "Genes to plot", 
                            choices = gene_options,
                            selected = gene_options, inline = TRUE),
         sliderInput("y_axis_range", "Fold change range", min = -5, 
                     max = max(pos_cnv_results$fold_change + 5),
                     value = c(-5, max(pos_cnv_results$fold_change + 5))),
         plotOutput("gene_amp_plot")),
  
  column(5,
         h2("Quality control plot"),
         selectInput(inputId = "ws_select", label = "Worksheet to highlight",
                     choices = ws_options, selected = "WS140721"),
         
         plotOutput("ws_noise_plot"))
  ),
  fluidRow(
    column(5,
           h2("Summary table"),
           tableOutput("summary_gene_table")
    ),
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

?checkboxGroupInput()

# Server ----------------------------------------------------------------------------

server <- function(input, output) {
  
  output$summary_gene_table <- renderTable({
    
    pos_cnv_results_with_qc |> 
      filter(gene != "no positive calls") |> 
      count(gene, .drop = FALSE) |> 
      arrange(desc(n)) |> 
      rename(Cases = n,
             Result = gene)
    
  })

  data_filtered <- reactive({
    
    pos_cnv_results_with_qc |> 
      filter(gene == input$gene_select & 
               noise < input$noise_qc) |> 
      select(worksheet, labno,  
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
  
}

# Run the application 
shinyApp(ui = ui, server = server)
