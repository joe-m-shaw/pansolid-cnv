count_primers_in_region <- function(chrom, coord1, coord2, df) {
  #' Count the QIAseq primers within a region
  #'
  #' @param chrom The chromosome of the region of interest
  #' @param coord1 The first coordinate of the region of interest
  #' @param coord2 The second coordinate of the region of interest
  #' @param df The dataframe containing the primer coordinate information.
  #' This should be in the format of one primer per row with one column each for the 
  #' primer start and end coordinates. 
  #' 
  #' @return Returns the number of primers within the specified region
  #' 
  #' @examples count_primers_in_region(chrom = "17", coord1 = 39700064,
  #' coord2 = 39728658)
  
  if(typeof(chrom) != "character") {
    stop("Chromosome argument must be a character")
  }
  
  if (typeof(coord1) != "double" |
      typeof(coord2) != "double") {
    stop("Coordinates must be type double")
  }
  
  if(!("primer_start" %in% colnames(df) & "primer_end" %in% colnames(df))) {
    stop("Primer dataframe does not contain primer_start and primer_end columns")
  }
  
  region_start <- min(coord1, coord2)
  
  region_end <- max(coord1, coord2)
  
  df2 <- df |> 
    filter(chromosome == chrom) |> 
    mutate(in_region = case_when(
      
      (primer_start >= region_start &
         primer_start <= region_end) |
        (primer_end >= region_start &
           primer_end <= region_end) ~"Yes",
      TRUE ~"No"
      
    ))
  
  primers_in_region <- sum(df2$in_region == "Yes")
  
  return(primers_in_region)
  
}