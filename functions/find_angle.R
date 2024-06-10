# function to calculate angle between three points
find_angle <- function(x, y, radians = FALSE) {
  if(length(x) != length((y))) stop('x and y must be the same length.\n')
  
  # coordinates for current location
  x_1 <- x
  y_1 <- y
  
  # coordinates for previous location
  x_0 <- c(NA, x[-length(x)])
  y_0 <- c(NA, y[-length(y)])
  
  # coordinates for next location
  x_2 <- c(x[-1], NA)
  y_2 <- c(y[-1], NA)
  
  # return NA if any coordinates are missing
  l_01 <- sqrt((x_0 - x_1)^2 + (y_0 - y_1)^2) # segm. to current location
  l_12 <- sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2) # segm. from current location
  l_02 <- sqrt((x_0 - x_2)^2 + (y_0 - y_2)^2) # distance between endpoints
  
  # calculate the angle between segments (in radians) using law of cosines
  # https://stackoverflow.com/questions/1211212/how-to-calculate-an-angle-from-three-points
  # https://www.mathsisfun.com/algebra/trig-cosine-law.html
  angle <- acos((l_01^2 + l_12^2 - l_02^2) / (2 * l_01 * l_12))
  
  # calculate turning angle (in radians)
  angle <- pi - angle
  
  if(radians) {
    message('Returned angles in radians.')
    return(angle)
  } else {
    message('Returned angles in degrees.')
    return(angle / pi * 180)
  }
}

# test the function
if(FALSE) {
  x <- c(0, 0, 1, 0.5)
  y <- c(1, 0, 0, 0.5)
  plot(x, y, type = 'b', pch = as.character(1:length(x)))
  find_angle(x, y)
  
  # random angles
  x <- runif(8)
  y <- runif(8)
  plot(y ~ x, type = 'b', pch = as.character(1:length(x)))
  find_angle(x = x, y = y)
  
  # estimated distribution of all possible angles
  hist(find_angle(runif(1e6), runif(1e6)), breaks = 250)
}
