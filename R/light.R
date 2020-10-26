#' Calculates light isotope tube shift
#'
#' @param data A dataframe containing taxon wad values per tube
#' @export

tube_shift_light = function(data = wad_long) {
  
  data = data %>%
    filter(Isotope == "12C" & wads != "NA" & wads != "NaN") %>%
    group_by(taxon) %>%
    mutate(taxon_12C_standard = median(wads),
           diff_from_12C_standard = taxon_12C_standard - wads,
           taxon_uncorrected_stress = diff_from_12C_standard^2) %>% # mean or median?
    ungroup() %>%
    group_by(SampleID) %>%
    mutate(tube_adjustment = median(diff_from_12C_standard),
           wads_tube_corrected = wads + tube_adjustment,
           taxon_corrected_stress = (taxon_12C_standard - wads_tube_corrected)^2,
           tube_rank_wads = percent_rank(wads),
           tube_rank_diff = percent_rank(-diff_from_12C_standard),
           inactive = TRUE) %>%
    ungroup()
  
  return(data)
}

