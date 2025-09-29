count_oligos_in_range <- function(chrom, coord1, coord2, df) {
  #' Count QIAseq oligos within a genomic coordinate range
  #'
  #' "Oligo" refers to an oligonucleotide, which may be a QIAseq primer or target. A QIAseq target is
  #' a region of a BED file targeted by multiple primers.
  #' 
  #' This function can be used with dplyr::rowwise() for a dataframe input.
  #'
  #' @param chrom The chromosome of the range of interest
  #' @param coord1 The first coordinate of the range of interest
  #' @param coord2 The second coordinate of the range of interest
  #' @param df The dataframe containing the oligonucleotide information.
  #' This should be in the format of one oligo (primer or target) per row with one column each for the 
  #' start and end coordinates. 
  #' 
  #' @return Returns the number of oligos within the specified region
  
  if(typeof(chrom) != "character") {
    stop("Chromosome argument must be a character")
  }
  
  if(!any(chrom %in% c(as.character(seq(1, 22, by = 1)), "X", "Y" ))){
    stop("Chromosome argument must be 1-22, X or Y")
  }
  
  if (typeof(coord1) != "double" |
      typeof(coord2) != "double") {
    stop("Coordinates must be type double")
  }
  
  if(!("start" %in% colnames(df) & "end" %in% colnames(df))) {
    stop("Oligo dataframe does not contain start and end columns")
  }
  
  region_start <- min(coord1, coord2)
  
  region_end <- max(coord1, coord2)
  
  df2 <- df |> 
    dplyr::filter(chromosome == chrom) |> 
    dplyr::mutate(in_region = case_when(
      
      (start >= region_start &
         start <= region_end) |
        (end >= region_start &
           end  <= region_end) ~"Yes",
      TRUE ~"No"
      
    ))
  
  oligos_in_region <- sum(df2$in_region == "Yes")
  
  return(oligos_in_region)
  
}