labeller_perc <- function(breaks) {
  breaks <- round(breaks * 100 - 100, 1)
  
  if_else(breaks <= 0,
          paste0(breaks, '%'),
          paste0('+', breaks, '%'))
}
