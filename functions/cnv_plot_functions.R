source(here("functions/extract_cnv_coordinates.R"))
source(here("scripts/set_shared_drive_filepath.R"))

get_cnv_plot_x_breaks <- function(interval, plot_xmin, plot_xmax) {
  #' Get consistent X axis breaks for creating CNV plots
  #'
  #' @param interval The desired interval between breaks
  #' @param plot_xmin The desired x axis minimum value
  #' @param plot_xmax The desired x axis maximum value
  #'
  #' @return A list of breaks
  #' @export
  #'
  #' @examples get_cnv_plot_x_breaks(interval = 100, plot_xmin = 0, plot_xmax = 1000)
  
  if(typeof(interval) != "double" |
     typeof(plot_xmin) != "double" |
     typeof(plot_xmax) != "double") {
    stop("All inputs must be numbers")
  }
  
  if(plot_xmin >= plot_xmax) {
    stop("plot_xmin must be smaller than plot_xmax")
  }
  
  breaks <- seq(plot_xmin, plot_xmax, by = interval)
  
  return(breaks)
  
}


get_data_for_cnv_plot <- function(df, gene) {
  
  #' Get data for making a copy number variant plot
  #'
  #' @param df The dataframe containing CNV information to be plotted. This should be based
  #' on the "Positive CN results" table from the PanSolid Excel files.
  #' @param gene The gene of interest
  #'
  #' @return Returns the dataframe filtered to include CNVs in the gene of interest
  #' @export
  #'
  #' @examples erbb2_cnvs <- get_data_for_cnv_plot(df = validation_pos_cnv_results_collated,
  #' gene = "ERBB2")
  
  data_for_plot <- df |> 
    filter(gene == {{ gene }})
  
  return(data_for_plot)
  
}

get_cnv_plot_xmin <- function(df, buffer) {
  
  #' Get the X axis minimum value for a CNV plot
  #'
  #' @param df A dataframe containing CNVs
  #' @param buffer The desired buffer to add onto the width of the X axis. This can be
  #' adjusted to preference based on how close or far from the X axis minimum the CNV
  #' calls should be.
  #'
  #' @return The value for the X axis minimum.
  #' @export
  #'
  #' @examples 
  
  if(typeof(buffer) != "double") {
    stop("buffer must be a number")
  }
  
  if(!"start" %in% colnames(df)) {
    stop("Input dataframe must have a column named start")
  }
  
  plot_xmin <- min(df$start) - buffer
  
  return(plot_xmin)
  
}

get_plot_xmax <- function(df, buffer) {
  
  plot_xmax <- max(df$end) + buffer
  
  return(plot_xmax)
  
}

get_chromosome <- function(gene) {
  
  gene_coordinates <- read_csv(file = paste0(data_folder,
                                             "gene_lists/",
                                             "gene_coordinates.csv"))
  
  chromosome <- as.character(gene_coordinates[gene_coordinates$gene == gene, 2])
  
  stopifnot(chromosome != "character(0)")
  
  return(chromosome)
  
}

make_fold_change_plot <- function(df, 
                                  gene,
                                  interval = 10000, 
                                  buffer = 5000, 
                                  ymin = 0,
                                  ymax = 40) {
  
  chromosome <- get_chromosome(gene = {{ gene }})
  
  data_for_plot <- get_data_for_cnv_plot(df = {{ df }}, 
                                     gene = {{ gene }})
  
  plot_xmin <- get_cnv_plot_xmin(df = data_for_plot,
                             buffer = buffer)
  
  plot_xmax <- get_plot_xmax(df = data_for_plot,
                             buffer = buffer)
  
  fold_change_plot <- ggplot(data_for_plot, aes(x = start, y = fold_change)) +
    
    # Add theme
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    
    # Add CNV calls
    geom_segment(aes(x = start, xend = end, 
                     y = fold_change, yend = fold_change),
                 linewidth = 2,
                 colour = safe_red) +
    
    # Add x axes
    scale_x_continuous(breaks = get_cnv_plot_x_breaks(interval = {{ interval}},
                                           plot_xmin = {{ plot_xmin }},
                                           plot_xmax = {{ plot_xmax }}),
                       minor_breaks = NULL,
                       limits = c({{ plot_xmin }}, {{ plot_xmax }} )) +
    
    scale_y_continuous(limits = c(ymin, ymax)) +
    
    # Add labels
    labs(
      y = "Fold change",
      x = "")
  
  return(list(plot_xmin, plot_xmax, interval, fold_change_plot, chromosome))
  
}

make_labno_plot <- function(df, 
                            gene,
                            interval = 10000, 
                            buffer = 5000, 
                            yaxis = labno) {
  
  chromosome <- get_chromosome(gene = {{ gene }})
  
  data_for_plot <- get_data_for_cnv_plot(df = {{ df }}, 
                                     gene = {{ gene }})
  
  max_fold_change <- max(data_for_plot$fold_change)
  
  min_fold_change <- min(data_for_plot$fold_change)
  
  plot_xmin <- get_cnv_plot_xmin(df = data_for_plot,
                             buffer = buffer)
  
  plot_xmax <- get_plot_xmax(df = data_for_plot,
                             buffer = buffer)
  
  labno_plot <- ggplot(data_for_plot, aes(x = start, y = {{ yaxis }},
                                          colour = fold_change)) +
    
    # Add theme
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    
    # Add CNV calls
    geom_segment(aes(x = start, xend = end, 
                     y = {{ yaxis }}, yend = {{ yaxis }}),
                 linewidth = 2) +
    
    scale_colour_gradient(low = "#FF9999", 
                          high = "#660000", 
                          limits = c(min_fold_change, max_fold_change),
                          n.breaks = 4) +
    
    # Add x axes
    scale_x_continuous(breaks = get_cnv_plot_x_breaks(interval = {{ interval}},
                                           plot_xmin = {{ plot_xmin }},
                                           plot_xmax = {{ plot_xmax }}),
                       minor_breaks = NULL,
                       limits = c({{ plot_xmin }}, {{ plot_xmax }} )) +
    
    # Add labels
    labs(
      y = "Sample number",
      x = "",
      colour = "Fold change")
  
  return(list(plot_xmin, plot_xmax, interval, labno_plot, chromosome))
  
}

make_primer_plot <- function(plot_xmin, plot_xmax, interval, chromosome) {
  
  grch38_primers <- read_csv(file = paste0(data_folder,
                                           "primers/CDHS-40079Z-11284.primer3_Converted.csv"),
                             show_col_types = FALSE) |> 
    janitor::clean_names()
  
  grch38_primer_coordinates <- extract_cnv_coordinates(df = grch38_primers,
                                                       cnv_coord_col = region)
  
  primers_filtered <- grch38_primer_coordinates |> 
    mutate(y_value = "Primers") |> 
    filter(chromosome == {{ chromosome }} ) |> 
    filter(start >= {{ plot_xmin }} & end <= {{ plot_xmax }} )
  
  output <-  ggplot(primers_filtered, aes(x = start, y = y_value)) +
    geom_point(pch = 21) +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    scale_x_continuous(breaks = get_cnv_plot_x_breaks(interval = {{ interval}},
                                           plot_xmin = {{ plot_xmin }},
                                           plot_xmax = {{ plot_xmax }}),
                       minor_breaks = NULL,
                       limits = c({{ plot_xmin }}, {{ plot_xmax }} )) +
    labs (x = "", y = "")
  
  return(output)
  
}

make_exon_plot <- function(plot_xmin, plot_xmax, interval, chromosome) {
  
  exon_data_for_plot <- all_transcripts |> 
    mutate(y_value = "Exons") |> 
    filter(chromosome == {{ chromosome }}) |> 
    filter(start >= {{ plot_xmin }} & end <= {{ plot_xmax }})
  
  gene_labels <- read_csv(file = paste0(data_folder, 
                                        "transcripts/processed/",
                                        "gene_labels.csv"))
  
  labels_for_plot <- gene_labels |> 
    filter(chromosome == {{ chromosome }} ) |> 
    filter(start >= {{ plot_xmin }} & start <= {{ plot_xmax }})
  
  output <- ggplot(exon_data_for_plot, 
                   aes(x = start, y = y_value)) +
    
    geom_segment(aes(x = start, xend = end, 
                     y = y_value, yend = y_value),
                 linewidth = 5) +
    
    theme_bw() +
    
    scale_x_continuous(breaks = get_cnv_plot_x_breaks(interval = {{ interval}},
                                           plot_xmin = {{ plot_xmin }},
                                           plot_xmax = {{ plot_xmax }}),
                       minor_breaks = NULL,
                       limits = c({{ plot_xmin }}, {{ plot_xmax }} )) +
    
    scale_y_discrete(limits = rev) +
    
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    
    labs(y = "", x = str_c("Genome coordinate (GRCh38) Chr", 
                           chromosome)) +
    
    geom_label(data = labels_for_plot, label = labels_for_plot$label)
  
  return(output)
  
}

make_cnv_triptych <- function(input_plot) {
  
  # This function is a wrapper which takes the outputs of either the 
  # make_fold_change_plot or make_labno_plot functions
  
  plot_xmin <- input_plot[[1]]
  
  plot_xmax <- input_plot[[2]]
  
  interval <- input_plot[[3]]
  
  main_plot <- input_plot[[4]]
  
  chromosome <- input_plot[[5]]
  
  primer_plot <- make_primer_plot(plot_xmin = {{ plot_xmin }}, 
                                  plot_xmax = {{ plot_xmax }},
                                  interval = {{ interval }},
                                  chromosome = {{ chromosome }})
  
  exon_plot <- make_exon_plot(plot_xmin = {{ plot_xmin }}, 
                              plot_xmax = {{ plot_xmax }},
                              interval = {{ interval }},
                              chromosome = {{ chromosome }})
  
  triptych <- (main_plot / primer_plot / exon_plot) +
    plot_layout(
      heights = c(6, 1, 2)
    )
  
  return(triptych)
  
}