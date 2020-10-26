#' Calculates heavy isotope tube shift
#' Adds light and heavy tube shifts together
#'
#' @param data A dataframe containing taxon wad values per tube
#' @param light A dataframe containing light isotope adjusted taxon wad values averaged across all tubes
#' @export

tube_shift_heavy = function(data, light, low = 0, high = 0.1) {
  require(tidyverse)
  taxon_12C_standard = light %>%
    select(taxon, taxon_12C_standard) %>%
    unique()

  rowcount = nrow(data)

  data = data %>%
    filter(Isotope == "13C" & wads != "NA") %>%
    left_join(taxon_12C_standard, by = "taxon") %>%
    group_by(SampleID) %>%
    mutate(diff_from_12C_standard = taxon_12C_standard - wads, # if wad density is below standard, diff will be positive number
           taxon_uncorrected_stress = diff_from_12C_standard^2,
           tube_rank_wads = percent_rank(wads),
           tube_rank_diff = percent_rank(-diff_from_12C_standard),
           inactive = case_when(tube_rank_wads < low ~ FALSE,
                                tube_rank_wads > high ~ FALSE,
                                tube_rank_wads >= low & tube_rank_wads <= high ~ TRUE)) %>%
    ungroup()

  inactive = data %>%
    filter(inactive == TRUE) %>%
    group_by(SampleID) %>%
    summarize(tube_adjustment = median(diff_from_12C_standard))

  data = data %>%
    left_join(inactive, by = "SampleID") %>%
    mutate(wads_tube_corrected = wads + tube_adjustment,
           taxon_corrected_stress = (taxon_12C_standard - wads_tube_corrected)^2) %>%
    ungroup()

  data = bind_rows(light, data)
  message(rowcount - nrow(data), " ASVs were removed from original dataframe")
  return(data)
}
