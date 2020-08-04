#' Calculates and plots stress
#' Window and step size can be specified by user
#'
#' @param data A dataframe containing taxon wad values per tube
#' @param light A dataframe containing light isotope adjusted taxon wad values averaged across all tubes
#' @export
#'

stress_test = function(data = wad_long, light = df_12C, window = 0.1, step = 0.1){
  f = seq(0, 1 - window, by = step)
  df = data.frame()
  for (i in f) {
    j = i + window
    S = tube_shift_heavy(data = wad_long, light = df_12C, high = j, low = i) %>%
      filter(inactive == TRUE) %>%
      summarize(S = mean(taxon_corrected_stress)) %>%
      pull()
    df = rbind(df, data.frame("low" = i, "high" = j, "stress" = S))

  }
  return(df)
}

stress_test(window = .1, step = 0.1) %>%
  ggplot(aes(x = low, y = stress)) +
  geom_point()

options(scipen = 999)

