source('functions/find_angle.R')

outlier_plots <- function(telemetry, units = FALSE, ci_level = 0,
                          return = FALSE, reset_layout = TRUE,
                          cap_dt = TRUE) {
  layout(matrix(c(1, 1, 2:5), ncol = 2, byrow = TRUE))
  # plot entire telemetry
  plot(telemetry$x, telemetry$y, type = 'l', xlab = 'x (m)',
       ylab = 'y (m)')
  # plot of locations and velocities
  out <- outlie(telemetry, plot = TRUE)
  # min speed vs distance from tel center (SI units; no CIs)
  plot(out, units = units, level = ci_level)
  
  # calculate turning angles in degrees
  out$angle <- find_angle(telemetry$x, telemetry$y, radians = FALSE)
  
  # calculate time intervals
  out$dt <- out$t - lag(out$t)
  
  plot(out$angle, out$speed, pch = 19, col = '#00000050', xaxt = 'n',
       xlim = c(0, 180), xlab = 'Turning angle (degrees)',
       ylab = 'Minimum speed (m/s)')
  axis(side = 1, at = seq(0, 180, by = 45))
  
  # SLD is often not useful because dt is too often clustered
  if(cap_dt) {
    XLIM <- c(0, min(max(out$dt / (1 %#% 'hour'), na.rm = TRUE), 24))
    XLAB <- if_else(XLIM[2] == 24 & cap_dt,
                    'Time between locations (hours, capped at 24)',
                    'Time between locations (hours)')
  } else {
    XLIM <- c(0, max(out$dt / (1 %#% 'hour'), na.rm = TRUE))
    XLAB <- 'Time between locations (hours)'
  }
  plot(out$dt / (1 %#% 'hours'), out$speed, pch = 19, col = '#00000050',
       xlab = XLAB,
       ylab = 'Minimum speed (m/s)',
       xlim = XLIM)
  
  if(reset_layout) layout(1)
  if(return) return(out)
}
