extract_cnv_coordinates <- function(df, cnv_coord_col) {
  
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