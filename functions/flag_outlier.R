# function to flag outliers:
#' `value == 0` => point should be kept for analyses
#' `value == 1` => point to be removed: - deviation is excessive,
#'                                          - speed is excessive,
#'                                          - outside of DOP ellipses,
#'                                          - turning angle is too small...
#' `value == 2` => point *may* need to be removed (but, e.g., speed is ok)
#' **NOTE:** `as.telemetry()` removes all rows with `outlier > 0`
flag_outlier <- function(id, max_speed = 0, max_distance = 0,
                         max_angle = 0, max_dt = 0, value) {
  o <- out$speed >= max_speed &
    out$distance >= max_distance &
    (out$angle >= max_angle | is.na(out$angle)) &
    (out$dt >= max_dt | is.na(out$dt))  
  tel <- d$tel[d$animal == id][[1]] # extract telemetry
  if(is.logical(tel$outlier)) {
    tel$outlier <- if_else(tel$outlier, 1, 0) # convert logical to numeric
  }
  tel[! tel$outlier, ][o, 'outlier'] <- value # skip outliers and add flags
  d$tel[d$animal == id][[1]] <<- tel
}
