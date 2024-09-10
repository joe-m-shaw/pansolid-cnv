extract_pansolid_cnv_coordinates <- function(df, cnv_coord_col) {
  
  #' Extract coordinates from the PanSolid CNV coordinate string
  #'
  #' @param df A dataframe containing a column of CNV coordinates to parse. This function was designed
  #' mainly to be used with the "Positive CNV results" table in the PanSolid Excel output.
  #' @param cnv_coord_col The column in the dataframe containing the CNV coordinates separated by two 
  #' full stops (Example: "15889325..15946102")
  #' 
  #' @return The input dataframe with the CNV start and end coordinates as two new columns.
  #' @export
  #'
  #' @examples cnvs_with_coordinates <- extract_pansolid_cnv_coordinates(df = pos_cnv_tbl, 
  #' cnv_coord_col = cnv_co_ordinates)
  
  cnv_coord_regex <- stringr::regex(
    r"[
        (|complement\()
        (\d{1,10})     # first coordinate number (1 to 10 digits)
        \.\.           # two full stops
        (\d{1,10})     # second coordinate number (1 to 10 digits)
        ]",
    comments = TRUE
  )
  
  output <- df |> 
    dplyr::mutate(start = as.numeric(stringr::str_extract(string = {{ cnv_coord_col }}, 
                                                          pattern = cnv_coord_regex, 
                                                          group = 2)),
                  end = as.numeric(stringr::str_extract(string = {{ cnv_coord_col }}, 
                                                        pattern = cnv_coord_regex, 
                                                        group = 3))) 
  
  return(output)
  
}