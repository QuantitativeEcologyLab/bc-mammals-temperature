library('purrr')     # for functional programming (map_***())
library('dplyr')     # for dara wrangling (mutate, transmute, etc.)
library('tidyr')     # for data wrangling (pivot_longer, pivot_wider, nest)
library('stringi')   # for working with strings
library('lubridate') # for working with dates

# to bind all climate projections from each year together into a single RDS file
d <-
  # import all files
  map_dfr(
    list.files('H:/GitHub/bc-mammals-speeds/ClimateNA_v742/bc-dem-z3',
      full.names = TRUE, # include folder names in file name
      pattern = '@'), # only every 10 years
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
#' *can simulate hourly temperature using a smooth sinusoidal function*
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
  # convert monthly total precip to hourly total precip
  mutate(first_day = as.Date(paste(year, month, '01', sep = '-')),
         next_month = if_else(month != '12', as.numeric(month + 1), 1),
         next_year = if_else(month != '12', year, year + 1),
         last_day = as.Date(paste(next_year, next_month, '01', sep = '-')),
         hours = as.numeric((last_day - first_day)) * 24) %>%
  # drop temporary columns
  select(-c(first_day, next_month, next_year, last_day, hours)) %>%
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
saveRDS(d, paste0('../data/climate-yearly-projections-', Sys.Date(), '.rds'))
