library('purrr')     # for functional programming (map_***())
library('dplyr')     # for dara wrangling (mutate, transmute, etc.)
library('tidyr')     # for data wrangling (pivot_longer, pivot_wider, nest)
library('stringi')   # for working with strings
library('lubridate') # for working with dates
library('sf')
library('terra')
source('data/bc-shapefile.R')

# to bind all climate projections from each year together into a single RDS file
d <-
  # import all files
  map_dfr(
    list.files('H:/GitHub/bc-mammals-temperature/ClimateNA_v742/bc-buffered-dem-z3/',
               full.names = TRUE, # include folder names in file name
               pattern = '@'),
    \(.fname) {
      # readr::read_csv(.fname, col_types = '?', progress = FALSE) %>%
      data.table::fread(.fname, showProgress = FALSE) %>%
        mutate(file = .fname) # add column of filename for scenario & year
    }, .progress = TRUE) %>%
  # filter to the approximate extent of the data
  filter(Longitude > -127, Longitude < -113,
         Latitude > 47, Latitude < 62) %>%
  # add scenario and year columns
  mutate(scenario = substr(
    file,
    start = stri_locate_last(file, regex = '/')[1] + 1,
    stop = stri_locate_first(file, regex = '@')[1] - 1),
    year = substr(file,
                  start = stri_locate_first(file, regex = '@')[1] + 1,
                  stop = nchar(file) - nchar('.csv'))) %>%
  # only keep necessary columns
  select(scenario, year, Latitude, Longitude, Elevation,
         Tave01:Tave12, Tmin01:Tmin12, Tmax01:Tmax12)

#' mean temperature is measured as `Tave = (Tmax + Tmin) / 2`
#' `https://www.prism.oregonstate.edu/documents/PRISM_datasets.pdf`
with(slice_sample(d, n = 5e4),
     plot(Tave01 - Tmin01, Tmax01 - Tave01, col = '#00000050', pch = '.'))
for(a in c(-0.1, 0, 0.1)) abline(a = a, b = 1, col = '#FF000040'); rm(a)

# splitting to reduce RAM costs
# pivot to long format
d <- pivot_longer(d, ! c(scenario, year, Latitude, Longitude, Elevation),
                  names_to = 'parameter', values_to = 'value')

d <- d %>%
  # extract time and parameter columns 
  mutate(month = map_chr(parameter,
                         \(.chr) substr(.chr, nchar(.chr) - 1, nchar(.chr))),
         dec_date = decimal_date(date(paste(year, month, '15', sep = '-'))),
         month = as.numeric(month),
         year = as.numeric(year),
         parameter = map_chr(parameter,
                             \(.chr) substr(.chr, 1, nchar(.chr) - 2))) %>%
  # create separate columns for temperature and precipitation
  pivot_wider(names_from = parameter, values_from = value)

d <- d %>%
  # change to names used in the models
  rename(mean_temperature = Tave,
         min_temperature = Tmin,
         max_temperature = Tmax,
         latitude = Latitude,
         longitude = Longitude,
         elevation = Elevation) %>%
  # place month and decimal date columns after the year column
  relocate(c(month, dec_date), .after = year)

# save the output for later use
saveRDS(d, paste0('data/climate-yearly-projections-', Sys.Date(), '.rds'))

# only 2025 and 2100 data ----
bc_points <- readr::read_csv('H:/GitHub/bc-mammals-temperature/ClimateNA_v742/bc-buffered-dem-z3/8GCMs_ensemble_ssp245@2025.csv') %>%
  transmute(Longitude, Latitude, z = 1) %>%
  rast(crs = 'EPSG:4326') %>%
  mask(bc_unproj) %>%
  as.data.frame(xy = TRUE) %>%
  transmute(Longitude = x, Latitude = y) %>%
  as_tibble()
bc_points
plot(bc_points)

d <-
  # import all files
  map(
    c(list.files('H:/GitHub/bc-mammals-temperature/ClimateNA_v742/bc-buffered-dem-z3/',
                 full.names = TRUE, pattern = '@2025'),
      list.files('H:/GitHub/bc-mammals-temperature/ClimateNA_v742/bc-buffered-dem-z3/',
                 full.names = TRUE, pattern = '@2100')),
    \(.fname) {
      data.table::fread(.fname, showProgress = FALSE) %>%
        mutate(file = .fname) # add column of filename for scenario & year
    }, .progress = TRUE) %>%
  bind_rows() %>%
  mutate(Longitude = round(Longitude, 3),
         Latitude = round(Latitude, 3)) %>%
  # filter to points in BC only
  inner_join(mutate(bc_points,
                    Longitude = round(Longitude, 3),
                    Latitude = round(Latitude, 3)),
             by = c('Longitude', 'Latitude')) %>%
  # add scenario and year columns
  mutate(scenario = substr(
    file,
    start = stri_locate_last(file, regex = '/')[1] + 1,
    stop = stri_locate_first(file, regex = '@')[1] - 1),
    year = substr(file,
                  start = stri_locate_first(file, regex = '@')[1] + 1,
                  stop = nchar(file) - nchar('.csv'))) %>%
  # only keep necessary columns
  select(scenario, year, Latitude, Longitude, Elevation,
         Tave01:Tave12, Tmin01:Tmin12, Tmax01:Tmax12) %>%
  pivot_longer(! c(scenario, year, Latitude, Longitude, Elevation),
               names_to = 'parameter', values_to = 'value') %>%
  # extract time and parameter columns
  mutate(month = map_chr(parameter,
                         \(.chr) substr(.chr, nchar(.chr) - 1, nchar(.chr))),
         date = date(paste(year, month, '15', sep = '-')),
         dec_date = decimal_date(date),
         month = as.numeric(month),
         year = as.numeric(year),
         parameter = map_chr(parameter,
                             \(.chr) substr(.chr, 1, nchar(.chr) - 2))) %>%
  # create separate columns for temperature for each month
  pivot_wider(names_from = parameter, values_from = value) %>%
  # change to the names used in the models
  rename(mean_temperature = Tave,
         min_temperature = Tmin,
         max_temperature = Tmax,
         latitude = Latitude,
         longitude = Longitude,
         elevation = Elevation) %>%
  # place month and decimal date columns after the year column
  relocate(c(month, dec_date), .after = year)

# save the output for later use
saveRDS(d, paste0('data/climate-yearly-projections-2025-2100-only-',
                  Sys.Date(), '.rds'))
