# Code written by Chris H Fleming and adapted by Stefano Mezzini
library('ctmm')       # for movement modeling
library('data.table') # for memory-efficient data.frames

# Function for de-trending empirical speed estimates from the movement model
detrend_speeds <- function(DATA, CTMM, units = FALSE) {
  
  # check if model is OUF
  if(! (grepl('OUF', toupper(summary(CTMM)$name)) |
        grepl('IOU', toupper(summary(CTMM)$name)))) {
    warning('Found a model with non-finite speed.\n')
    
    SPEEDS <- data.frame(low = NA_real_,
                         est = NA_real_,
                         high = NA_real_,
                         time = DATA$timestamp)
  } else {
    # estimate speeds conditional on the data
    EST <- speeds(DATA, CTMM, fast = TRUE, level = NULL, units = units)
    # null estimates of speed (no data, only model)
    EST.NULL <- speeds(CTMM, t = DATA$t, fast = TRUE, level = NULL,
                       units = units)
    
    # reduce DOF by the DOF of the null (i.e., mean) speed to reduce the
    # weight of speed estimates with little to no information on speed.
    DOF <- EST$DOF - EST.NULL$DOF
    DOF <- pmax(0,DOF) # haven't needed this, but just in case
    
    # assuming a chi^2 likelihood estimation, remove the background speed
    # in an amount proportional to the null DOF (since the sum of squared
    # chi random variables follows a chi^2 distribution)
    # formula used is equivalent to 
    #' `speed()^2 + d0 / ds * (speed()^2 - null^2)`
    #' where:
    #' - `speed()` is the output from `speed()`,
    #' - `null` is the null speed (mean speed not informed by nearby data), 
    #' - `d0` is the null DOF for the null speed
    #' - `ds` is the DOF for the data-informed speed
    #' this formula increases the speed if `speed() > null` and reduces it
    #' otherwise (which should keep an average of `null`). The amount by
    #' which it changes the speed depends on how informed the estimate is,
    #' based on `d0 / ds`.
    X2 <- (EST$DOF * EST$speed^2 - EST.NULL$DOF * EST.NULL$speed^2) / DOF
    X2 <- pmax(X2,0)
    
    # Calculate chi-squared CIs
    CI <- sapply(1:length(DATA$t),
                 function(i){
                   ctmm:::chisq.ci(MLE = X2[i], DOF = DOF[i],
                                   level = 0.95, robust = TRUE)
                 })

    # convert to chi distribution: sqrt(chi^2(n)) ~ chi(n)
    CI <- sqrt(CI)
    
    # transpose to a readable format
    SPEEDS <- as.data.frame(t(CI))
    
    # add timestamps
    SPEEDS$time <- DATA$timestamp
  }
  return(SPEEDS)
}

# test the function
if(FALSE) {
  library('lubridate')
  data('buffalo')
  
  TRACK <- buffalo[[1]][1:500,]
  GUESS <- ctmm.guess(TRACK, interactive = FALSE)
  FIT <- ctmm.fit(TRACK, GUESS, control = list(method = "pNewton",
                                               cores = -1))
  EST <- speeds(TRACK, FIT, fast=TRUE, level=NULL)
  TEST <- detrend_speeds(TRACK, FIT)
  
  TEST$time2 <- hour(TEST$time) + (minute(TEST$time)/60)
  EST$time2 <- hour(EST$timestamp) + (minute(EST$timestamp)/60)
  
  layout(t(1:2))
  plot(speed ~ time2, data = EST, pch = 19, cex = 0.5, ylim=c(0,1.5),
       main = 'Original')
  plot(est ~ time2, data = TEST, pch = 19, cex = 0.5, ylim=c(0,1.5),
       main = 'Detrended')
}
