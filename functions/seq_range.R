seq_range <- function(x, n = 250, na.rm = TRUE) {
  MIN <- min(x, na.rm = na.rm)
  MAX <- max(x, na.rm = na.rm)
  SEQ <- seq(from = MIN, to = MAX, length.out = n)
  return(SEQ)
}
