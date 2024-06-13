# Code written by Chris H Fleming and adapted by Stefano Mezzini
library('ctmm')       # for movement modeling
library('data.table') # for memory-efficient data.frames

# Function for de-trending empirical speed estimates from the movement model
detrend_speeds <- function(DATA, CTMM, units = FALSE) {
  
  # check if model is OUF
  if(! grepl('OUF', toupper(summary(CTMM)$name))) {
    warning('Found a non-OUF/OUf model; cannot estimate speeds.\n')
    
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
    
    # reduce DOF by that of null distribution
    DOF <- EST$DOF - EST.NULL$DOF
    all(DOF>0) # check for validity... always has held so far
    DOF <- pmax(0,DOF) # haven't needed this, but just in case
    
    S2 <- EST$speed^2 + EST$VAR
    
    AVE <- (EST$DOF*S2)/DOF
    
    # Calculate CIs
    CI <- sapply(
      1:length(DATA$t),
      function(i){ ctmm:::chisq.ci(MLE = AVE[i], DOF = DOF[i],
                                   level = 0.95, robust = TRUE)})
    CI <- sqrt(CI)
    
    SPEEDS <- as.data.frame(t(CI))
    
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
