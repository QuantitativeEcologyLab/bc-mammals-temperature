library('dplyr')     # for data wrangling (e.g., mutate, %>%)
library('tidyr')     # for data wrangling (e.g., nest, unnest)
library('purrr')     # for functional programming (e.g., map_***())
library('furrr')     # for parallel functional programming (e.g., future_map())
library('sf')        # for working with spatial data
library('terra')     # for working with raster data
library('ncdf4')     # for nc datasets
library('ggplot2')   # for fancy plots
library('ctmm')      # for movement modeling
library('lubridate') # for working with dates
source('functions/detrend_speeds.R') # to remove baseline noise in speed

theme_set(theme_bw())

# add weather data and speeds to the telemetry data ----
# see: https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation

#' custom `which()` that returns `NA` instead of `integer(0)`
which2 <- function(logical_vector) {
  z <- which(logical_vector)
  
  if(length(z) > 0) {
    return(z)
  } else {
    return(NA_integer_)
  }
}

# estimate speeds ----
if(FALSE) {
  CORES <- availableCores(logical = FALSE) - 1
  plan(multisession, workers = CORES)
  
  d <-
    readRDS('models/movement-models-akdes-2024-06-06.rds') %>%
    mutate(
      # estimate speeds
      tel = future_imap(tel, \(.tel, i) {
        detrend_speeds(DATA = .tel, CTMM = movement_model[[i]]) %>%
          transmute(speed_lwr = low, speed_est = est, speed_upr = high) %>%
          bind_cols(data.frame(.tel), .) %>%
          # add AKDE weights
          mutate(weight = unique(akde[[i]]$DOF.H) * akde[[i]]$weights) %>%
          return()
      }, .progress = TRUE, .options = furrr_options(seed = TRUE))) %>%
    # drop unnecessary columns
    select(! c(variogram, movement_model, akde,  mm_file_name,
               akde_file_name))%>%
    unnest(tel)
  
  plan(sequential)
  saveRDS(d, paste0('data/movement-models-speed-weights-', Sys.Date(), '.rds'))
} else {
  d <- readRDS('data/movement-models-speed-weights-2024-06-10.rds')
}

# extract temperature (takes about ~1 hour on lab machine w 11 cores) ----
CORES <- future::availableCores(logical = FALSE) - 1
plan(multisession, workers = CORES)
d2 <-
  d %>%
  # round to the nearest hour
  mutate(nearest_hour = timestamp %>%
           as_datetime() %>%
           round(units = 'hours') %>%
           with_tz(tzone = 'UTC')) %>%
  # sort by time
  arrange(timestamp) %>%
  # split data by year (one nc file for each year)
  mutate(year = year(nearest_hour)) %>%
  nest(yearly_data = ! year) %>%
  mutate(yearly_data = future_map2(year, yearly_data, \(.y, .d) {
    r <- rast(paste0('data/ecmwf-era5-2m-temperature/epsg-4326/',
                     'ecmwf-era5-2m-temp-', .y, '-epsg-4326.nc'))
    
    # split the data by hour (one layer for each hour in each nc file)
    .d %>%
      nest(hourly_data = ! c(nearest_hour)) %>%
      mutate(
        layer = map_int(nearest_hour, \(nh)
                        which2(as.POSIXct(r@cpp$time, tz = 'UTC') == nh)),
        hourly_data = map2(
          layer, hourly_data,
          \(.l, hourly_d) {
            if(is.na(.l)) {
              hourly_d$temperature_K <- NA_real_
            } else {
              hourly_d$temperature_K <-
                extract(x = r[[.l]],
                        y = select(hourly_d, longitude, latitude),
                        ID = FALSE) %>%
                unlist() %>%
                unname()
            }
            return(hourly_d)
          })
      ) %>%
      select(! layer) %>%
      unnest(hourly_data) %>%
      return()
  }, .progress = TRUE)) %>%
  unnest(yearly_data) %>%
  mutate(temp_c = temperature_K - 273.15) %>% # Kelvin to Celsius
  select(! temperature_K) # to keep < 100 MB

# ensure no temperature values are NA
group_by(d2, dataset_name) %>%
  summarize(prop = mean(is.na(temp_c)))
filter(d2, is.na(temp_c))

saveRDS(d2, paste0('data/movement-models-speed-weights-temperature-',
                   Sys.Date(), '.rds'))

# check ranges in temperature
range(d2$temp_c)

# make sure times are reasonable
slice_sample(d2, n = 5e3) %>%
  mutate(hour = hour(timestamp)) %>%
  ggplot(aes(hour, temp_c)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5), color = 'red3') +
  labs(x = 'Time of day', y = 'Temperature (\u00B0C)')
