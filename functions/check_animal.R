source('functions/outlier_plots.R')

#' function to run `outlie_plots()` for a given animal
check_animal <- function(id, UERE = NULL, return_out = TRUE, ...) {
  tel <- filter(d, animal == id) %>%
    pull(tel) %>%
    first() %>%
    as.telemetry(mark.rm = TRUE)
  
  # add User Equivalent Range Error or a default error of 10 m
  if(! is.null(UERE)) {
    uere(tel) <- UERE
  } else {
    uere(tel) <- 10
  }
  
  outlier_plots(tel, return = return_out, ...)
}
