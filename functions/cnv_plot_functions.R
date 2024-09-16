source(here("functions/extract_pansolid_cnv_coordinates.R"))
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
  
  breaks <- plyr::round_any(seq(plot_xmin, plot_xmax, by = interval),
                            accuracy = interval)
  
  return(breaks)
  
}

get_data_for_cnv_plot <- function(df, gene) {
  
  #' Get data for making a copy number variant plot
  #'
  #' @param df The dataframe containing CNV information to be plotted. This should be based
  #' on the "Positive CNV results" table from the PanSolid Excel files.
  #' @param gene The gene of interest
  #'
  #' @return Returns the dataframe filtered to include CNVs in the gene of interest
  #' @export
  #'
  #' @examples erbb2_cnvs <- get_data_for_cnv_plot(df = validation_pos_cnv_results_collated,
  #' gene = "ERBB2")
  
  data_for_plot <- df |> 
    dplyr::filter(gene == {{ gene }})
  
  if(nrow(data_for_plot) == 0){
    stop("Data for plot has 0 rows")
  }
  
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


get_cnv_plot_xmax <- function(df, buffer) {
  
  #' Get the X axis maximum value for a CNV plot
  #'
  #' @param df A dataframe containing CNVs
  #' @param buffer The desired buffer to add onto the width of the X axis.
  #'
  #' @return The value for the X axis maximum
  #' @export
  #'
  #' @examples 
  
  if(typeof(buffer) != "double") {
    stop("buffer must be a number")
  }
  
  if(!"end" %in% colnames(df)) {
    stop("Input dataframe must have a column named end")
  }
  
  plot_xmax <- max(df$end) + buffer
  
  return(plot_xmax)
  
}

get_gene_chromosome <- function(gene) {
  
  #' Get the chromosome for a gene
  #'
  #' @param gene The gene of interest
  #'
  #' @return The chromosome which the gene is on, formatted as a string.
  #' @export
  #'
  #' @examples chrom <- get_gene_chromosome("ERBB2")
  
  gene_coordinates <- readr::read_csv(file = paste0(data_folder,
                                             "gene_lists/",
                                             "gene_coordinates.csv"),
                               col_types = list(
                                 "gene" = col_character(),
                                 "chromosome" = col_character(),
                                 "transcript_ensembl" = col_character(),
                                 "transcript_refseq"	= col_character(),
                                 "gene_start" = col_double(),
                                 "gene_end" = col_double()))
  
  chromosome <- as.character(gene_coordinates[gene_coordinates$gene == gene, 2])
  
  stopifnot(chromosome != "character(0)")
  
  return(chromosome)
  
}

make_fold_change_cnv_plot <- function(df, 
                                  gene,
                                  interval = 10000, 
                                  buffer = 5000, 
                                  ymin = 0,
                                  ymax = 40,
                                  title = "") {
  
  #' Make a CNV plot with fold change on the Y axis
  #'
  #' @param df The dataframe containing CNV information to be plotted. 
  #' This should be based on the "Positive CN results" table from the PanSolid 
  #' @param gene The gene of interest
  #' @param interval The desired interval for X axis breaks
  #' @param buffer The desired buffer to add onto the width of the X axis.
  #' @param ymin The desired Y axis minimum value
  #' @param ymax The desired Y axis maximum value
  #'
  #' @return A list of values for the creation of other plots, and the rendered CNV plot.
  #' These additional values are supplied so that the output of this function can act 
  #' an input to the make_cnv_triptych_plot function.
  #' @export
  #'
  #' @examples erbb2_fold_change_plot <- make_fold_change_cnv_plot(
  #' df = validation_pos_cnv_results_collated,
  #' gene = "ERBB2")
  
  chromosome <- get_gene_chromosome(gene = {{ gene }})
  
  data_for_plot <- get_data_for_cnv_plot(df = {{ df }}, 
                                     gene = {{ gene }})
  
  plot_xmin <- get_cnv_plot_xmin(df = data_for_plot,
                             buffer = buffer)
  
  plot_xmax <- get_cnv_plot_xmax(df = data_for_plot,
                             buffer = buffer)
  
  fold_change_plot <- ggplot2::ggplot(data_for_plot, aes(x = start, y = fold_change)) +
    
    # Add theme
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    
    # Add CNV calls
    geom_segment(aes(x = start, xend = end, 
                     y = fold_change, yend = fold_change),
                 linewidth = 2,
                 colour = "#CC6677") +
    
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
      x = "",
      title = title)
  
  return(list(plot_xmin, plot_xmax, interval, fold_change_plot, chromosome))
  
}


make_labno_cnv_plot <- function(df, 
                            gene,
                            interval = 10000, 
                            buffer = 5000, 
                            yaxis = labno,
                            title = "") {
  
  #' Make a CNV plot with sample lab number on the Y axis
  #'
  #' @param df The dataframe containing CNV information to be plotted. 
  #' This should be based on the "Positive CN results" table from the PanSolid 
  #' Excel files.
  #' @param gene The gene of interest
  #' @param interval The desired interval for X axis breaks
  #' @param buffer The desired buffer to add onto the width of the X axis.
  #' @param yaxis The variable to plot on the Y axis, which defaults to labno, but can be 
  #' changed to other sample identifiers (such as NHS number) if desired.
  #'
  #' @return A list of values for the creation of other plots, and the rendered CNV plot.
  #' These additional values are supplied so that the output of this function can act 
  #' an input to the make_cnv_triptych_plot function.
  #' @export
  #'
  #' @examples erbb2_labno_plot <- make_labno_cnv_plot(
  #' df = validation_pos_cnv_results_collated,
  #' gene = "ERBB2",
  #' interval = 10000)
  
  chromosome <- get_gene_chromosome(gene = {{ gene }})
  
  data_for_plot <- get_data_for_cnv_plot(df = {{ df }}, 
                                     gene = {{ gene }})
  
  max_fold_change <- max(data_for_plot$fold_change)
  
  min_fold_change <- min(data_for_plot$fold_change)
  
  plot_xmin <- get_cnv_plot_xmin(df = data_for_plot,
                             buffer = buffer)
  
  plot_xmax <- get_cnv_plot_xmax(df = data_for_plot,
                             buffer = buffer)
  
  labno_plot <- ggplot2::ggplot(data_for_plot, aes(x = start, y = {{ yaxis }},
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
      colour = "Fold change",
      title = title)
  
  return(list(plot_xmin, plot_xmax, interval, labno_plot, chromosome))
  
}

make_primer_plot <- function(plot_xmin, plot_xmax, interval, chromosome) {
  
  #' Make a plot showing PanSolid QIAseq primer location data for inclusion in a CNV plot
  #'
  #' @param plot_xmin The desired X axis minimum value
  #' @param plot_xmax The desired X axis maximum value
  #' @param interval The desired interval for X axis breaks
  #' @param chromosome The chromosome for the gene of interest
  #'
  #' @return
  #' @export A plot showing the locations of QIAseq primers within the specified
  #' genomic region.
  #'
  #' @examples erbb2_primers <- make_primer_plot(plot_xmin = 39700064, 
  #' plot_xmax = 39728658,
  #' interval = 10000, chromosome = "17")
  
  grch38_primers <- readr::read_csv(file = paste0(data_folder,
                                           "primers/CDHS-40079Z-11284.primer3_Converted.csv"),
                             show_col_types = FALSE) |> 
    janitor::clean_names()
  
  grch38_primer_coordinates <- extract_pansolid_cnv_coordinates(df = grch38_primers,
                                                       cnv_coord_col = region)
  
  primers_filtered <- grch38_primer_coordinates |> 
    dplyr::mutate(y_value = "Primers") |> 
    dplyr::filter(chromosome == {{ chromosome }} ) |> 
    dplyr::filter(start >= {{ plot_xmin }} & end <= {{ plot_xmax }} )
  
  output <-  ggplot2::ggplot(primers_filtered, aes(x = start, y = y_value)) +
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
  
  #' Make an plot showing exon data for inclusion in a CNV plot
  #'
  #' @param plot_xmin The desired X axis minimum value
  #' @param plot_xmax The desired X axis maximum value
  #' @param interval The desired interval for X axis breaks
  #' @param chromosome The chromosome for the gene of interest
  #'
  #' @return A plot showing the exon locations for target genes within the specified
  #' genomic region.
  #' @export
  #'
  #' @examples erbb2_exons <- make_exon_plot(plot_xmin = 39700064, plot_xmax = 39728658,
  #' interval = 10000, chromosome = "17")
  
  all_transcripts <- readr::read_csv(paste0(data_folder,
                                     "transcripts/processed/",
                                     "collated_transcripts.csv"),
                              show_col_types = FALSE)
  
  exon_data_for_plot <- all_transcripts |> 
    dplyr::mutate(y_value = "Exons") |> 
    dplyr::filter(chromosome == {{ chromosome }}) |> 
    dplyr::filter(start >= {{ plot_xmin }} & end <= {{ plot_xmax }})
  
  gene_labels <- readr::read_csv(file = paste0(data_folder, 
                                        "transcripts/processed/",
                                        "gene_labels.csv"),
                          show_col_types = FALSE)
  
  labels_for_plot <- gene_labels |> 
    dplyr::filter(chromosome == {{ chromosome }} ) |> 
    dplyr::filter(start >= {{ plot_xmin }} & start <= {{ plot_xmax }})
  
  output <- ggplot2::ggplot(exon_data_for_plot, 
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
    
    labs(y = "", x = stringr::str_c("Genome coordinate (GRCh38) Chr", 
                           chromosome)) +
    
    geom_label(data = labels_for_plot, label = labels_for_plot$gene)
  
  return(output)
  
}


make_cnv_triptych_plot <- function(input_plot) {
  
  #' Create a triptych CNV plot with primer and exon tracks
  #'
  #' @param input_plot The input plot to format as a triptych. This can be combined 
  #' with the make_labno_cnv_plot and make_fold_change_cnv_plot functions above 
  #' (see example).
  #'
  #' @return A triptych plot with 3 panels: first the input plot displaying the copy number
  #' variant calls, then a plot showing the locations of PanSolid QIAseq primers, then
  #' a plot showing the locations of target genes and their individual exons. All plots
  #' have a consistently formatted X axis,
  #' @export
  #'
  #' @examples erbb2_plot <- make_cnv_triptych_plot(
  #' make_labno_cnv_plot(df = validation_pos_cnv_results_collated,
  #' gene = "ERBB2",
  #' interval = 10000))

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
