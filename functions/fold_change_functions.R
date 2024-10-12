calculate_fold_change <- function(tcc_percent, sample_target_copies_per_tumour_cell) {
  
  ref_sample_copies <- 200
  
  tcc_fraction <- tcc_percent / 100
  
  sample_target_copies_in_tumour_cells <- (100 * tcc_fraction) * sample_target_copies_per_tumour_cell
  
  sample_target_copies_in_normal_cells <- (100 * (1-tcc_fraction)) * 2
  
  sample_total_target_copies <- sample_target_copies_in_tumour_cells + sample_target_copies_in_normal_cells
  
  fold_change <- sample_total_target_copies / ref_sample_copies
  
  return(fold_change)
  
}

calculate_target_copies <- function(fold_change, tcc_percent) {
  
  ref_sample_copies <- 200
  
  tcc_fraction <- tcc_percent / 100
  
  sample_total_target_copies <- fold_change * ref_sample_copies
  
  sample_target_copies_in_tumour_cells <- sample_total_target_copies - ((100 * (1-tcc_fraction)) * 2)
  
  sample_target_copies_per_tumour_cell <- sample_target_copies_in_tumour_cells / (100 * tcc_fraction)
  
  return(sample_target_copies_per_tumour_cell)
  
}
