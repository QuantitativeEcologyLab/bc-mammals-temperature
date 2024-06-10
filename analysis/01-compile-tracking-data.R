library('dplyr')     # for data wrangling (mutate(), %>%, etc.)
library('lubridate') # for working with dates (force_tz(), etc.)
library('purrr')     # for functional programming (map_***(), etc.)
library('ctmm')      # for movement models
library('sf')        # for spatial features
library('ggplot2')   # for fancy plots
source('functions/import_rda.R') # to convert Rda files to rds files
source('data/bc-shapefile.R') # for shapefile of BC

bc_tz <- 'America/Vancouver' # to avoid typos

# import all telemetries separately ----
# mountain goats ----
goats <-
  readRDS('data/tracking-data/oreamnos-americanus-tels.rds') %>%
  mutate(animal = individual.local.identifier, # to add a column of ID
         timezone = bc_tz,
         species = 'Oreamnos americanus',
         dataset_name = 'Oreamnos_americanus',
         outlier = FALSE) # outlier column is not present

# pumas ----
pumas <- bind_rows(
  # puma 2
  import_rda('data/tracking-data/Puma_concolor_2.Rda') %>%
    as_tibble() %>%
    arrange(AnimalID, timestamp) %>% # some data out of order
    mutate(AnimalID = as.character(AnimalID),
           animal = AnimalID,
           species = 'Puma concolor',
           dataset_name = 'Puma_concolor_2',
           timestamp = ymd_hms(timestamp, tz = bc_tz),
           Time.Zone = bc_tz) %>% # unique: PST, DST, PST/PDT: assume BC tz
    rename(individual.local.identifier = AnimalID,
           location.lat = Latitude,
           location.long = Longitude,
           timezone = Time.Zone,
           outlier = Outlier) %>%
    select(species, dataset_name, animal, individual.local.identifier,
           location.lat, location.long, timestamp, timezone, outlier),
  # puma 4
  import_rda('data/tracking-data/Puma_concolor_4.Rda') %>%
    as_tibble() %>%
    mutate(animal = individual.local.identifier,
           species = 'Puma concolor',
           dataset_name = 'Puma_concolor_4',
           timezone = bc_tz,
           timestamp = ymd_hms(timestamp, tz = bc_tz),
           outlier = FALSE) %>% # outlier column is missing
    rename(temperature = external.temperature) %>%
    select(species, dataset_name, animal, individual.local.identifier,
           location.lat, location.long, timestamp, timezone,
           gps.fix.type.raw, sensor.type, outlier, temperature))

# sensor type column is not useful (only one value)
unique(pumas$sensor.type)

pumas <- select(pumas, ! sensor.type)

# southern mountain caribou in south-eastern BC ----
caribou_m <- import_rda('data/tracking-data/Rangifer_tarandus.Rda',
           object_name = 'data') %>%
  as_tibble() %>%
  # "validated" column causes "error in `[<-`:! subscript out of bounds"
  select(! 'Validated')

# two have "GMT?" as timezone, but they have similar activity patters to
# the other caribou
caribou_m %>%
  group_by(AnimalID) %>%
  summarise(GMT = any(Time.Zone == 'GMT?')) %>%
  summarize(BC = sum(! GMT), GMT = sum(GMT))

# all of them are in the same area (none in Alberta)
ggplot(caribou_m) +
  geom_sf(data = bc) +
  geom_point(aes(Longitude, Latitude, col = Time.Zone))

caribou_m <-
  caribou_m %>%
  mutate(animal = AnimalID,
         species = 'Rangifer tarandus',
         dataset_name = 'Rangifer_tarandus_southern_mountain',
         timezone = bc_tz,
         timestamp = ymd_hms(timestamp, tz = bc_tz))

# plot activity by BC time (all timestamps seem to be in BC time)
caribou_m %>%
  group_by(AnimalID) %>%
  arrange(AnimalID, timestamp) %>%
  mutate(hour = hour(timestamp),
         distance = sqrt((Longitude - lag(Longitude))^2 +
                           (Latitude - lag(Latitude))^2)) %>%
  ggplot() +
  facet_wrap(~ Time.Zone) +
  geom_smooth(aes(hour, distance), formula = y ~ s(x, k = 5, bs = 'cc'))

caribou_m <- caribou_m %>%
  rename(individual.local.identifier = AnimalID,
         outlier = Outlier,
         location.lat = Latitude,
         location.long = Longitude,
         temperature = Temperature) %>%
  select(species, dataset_name, animal, individual.local.identifier,
         location.lat, location.long, timestamp, timezone, outlier,
         DOP, PDOP, HDOP)

# PDOP column is not useful
unique(caribou_m$PDOP)
caribou_m <- select(caribou_m, ! PDOP)

# grizzly ----
grizzly <-
  import_rda('data/tracking-data/Ursus_arctos_horribilis.Rda') %>%
  as_tibble() %>%
  mutate(AnimalID = as.character(AnimalID),
         animal = AnimalID,
         species = 'Ursus arctos horribilis',
         dataset_name = 'Ursus_arctos_horribilis',
         timestamp = ymd_hms(timestamp, tz = bc_tz))

# again, timestamp seems to be corrected to Vancouver time
grizzly %>%
  mutate(hour = hour(timestamp),
         distance = sqrt((Longitude - lag(Longitude))^2 +
                           (Latitude - lag(Latitude))^2)) %>%
  ggplot() +
  facet_wrap(~ Time.Zone) +
  geom_smooth(aes(hour, distance), formula = y ~ s(x, k = 5, bs = 'cc'))

grizzly <- grizzly %>%
  rename(individual.local.identifier = AnimalID,
         location.lat = Latitude,
         location.long = Longitude,
         timezone = Time.Zone,
         outlier = Outlier) %>%
  select(species, dataset_name, animal, individual.local.identifier,
       location.lat, location.long, timestamp, timezone, outlier)

# elk in south-west AB & south-east BC ----
# https://datarepository.movebank.org/entities/datapackage/3a171441-6c9b-4bf9-b2fd-0216d962f091
#' using `read.csv()` because `read_csv()` keeps dashes in colnames
elk <- 
  read.csv('data/tracking-data/Elk in southwestern Alberta.csv') %>%
  rename(species = individual.taxon.canonical.name,
         dataset_name = study.name,
         temperature = external.temperature) %>%
  mutate(animal = as.character(individual.local.identifier),
         outlier = manually.marked.outlier == 'true',
         timezone = bc_tz,
         timestamp = ymd_hms(timestamp, tz = 'UTC') %>%
           force_tz(bc_tz)) %>% # bc_tz causes parsing issues
  # remove duplicates (too few to make a difference anyway; 0.02%)
  group_by(animal) %>%
  ungroup() %>%
  select(species, dataset_name, animal, individual.local.identifier,
         location.lat, location.long, timestamp, timezone, outlier,
         gps.dop, gps.fix.type.raw, gps.satellite.count, temperature)

# check activity levels by time of day
elk %>%
  mutate(hour = hour(timestamp),
         distance = sqrt((location.long - lag(location.long))^2 +
                           (location.lat - lag(location.lat))^2)) %>%
  ggplot() +
  geom_smooth(aes(hour, distance), formula = y ~ s(x, k = 5, bs = 'cc'))

# some duplicates in the elk dataset (~0.02 % of the data)
with(elk,
     mean(duplicated(paste(animal, timestamp)) |
            duplicated(paste(animal, timestamp), fromLast = TRUE)))

# all locations for E089 are marked as outliers
unique(filter(elk, animal == 'E089')$outlier)

elk <- elk %>%
  # drop duplicates
  filter(! (duplicated(paste(animal, timestamp)))) %>%
  # drop elk E089
  filter(animal != 'E089')

# boreal caribou ----
caribou_b <-
  readxl::read_xlsx('data/tracking-data/webexportcaribou.xlsx') %>%
  mutate(individual.local.identifier = FieldID,
         location.lat = Latitude,
         location.long = Longitude,
         timestamp = force_tz(FixDateTime, tz = bc_tz),
         timezone = bc_tz,
         animal = individual.local.identifier,
         species = 'Rangifer tarandus',
         dataset_name = 'Rangifer_tarandus_boreal',
         outlier = FALSE,
         temperature = Temperature) %>%
  select(species, dataset_name, animal, individual.local.identifier,
         location.lat, location.long, timestamp, timezone, outlier,
         NAV, HDOP, temperature)

# boreal wolves ----
wolves <-
  readxl::read_xlsx('data/tracking-data/webexportwolf.xlsx') %>%
  mutate(individual.local.identifier = FieldID,
         location.lat = Latitude,
         location.long = Longitude,
         timestamp = force_tz(FixDateTime, tz = bc_tz),
         timezone = bc_tz,
         animal = individual.local.identifier,
         species = 'Canis_lupus',
         dataset_name = 'Canis_lupus_boreal',
         outlier = FALSE,
         temperature = Temperature) %>%
  select(species, dataset_name, animal, individual.local.identifier,
         location.lat, location.long, timestamp, timezone, outlier,
         NAV, HDOP, temperature)

# bind all the telemetry objects into a single tibble ----
d <-
  bind_rows(goats, pumas, caribou_b, caribou_m, grizzly, elk, wolves) %>%
  arrange(dataset_name, animal, timestamp)

#' plot the data (`pch = '.'` is much faster than default `pch`)
plot(bc)
points(location.lat ~ location.long, goats, col = 'black', pch = '.')
points(location.lat ~ location.long, pumas, col = 'red', pch = '.')
points(location.lat ~ location.long, caribou_m, col = 'blue', pch = '.')
points(location.lat ~ location.long, grizzly, col = 'darkorange', pch = '.')
points(location.lat ~ location.long, elk, col = 'pink', pch = '.')
points(location.lat ~ location.long, caribou_b, col = 'cyan', pch = '.')
points(location.lat ~ location.long, wolves, col = 'grey', pch = '.')

#' all times have been converted to America/Vancouver
tz(d$timestamp)

# check outlier columns ----
d %>%
  group_by(dataset_name) %>%
  summarize(outliers = mean(outlier))  

# check GPS attribute columns ----
table(d$gps.fix.type.raw)
table(d$NAV)

d %>%
  group_by(dataset_name) %>%
  summarize(n_hdop = n_distinct(HDOP, na.rm = TRUE),
            n_dop = n_distinct(DOP, na.rm = TRUE),
            n_gps.dop = n_distinct(gps.dop, na.rm = TRUE),
            n_satellite = n_distinct(gps.satellite.count, na.rm = TRUE))

# check if times are correct
if(FALSE) {
  library('ggplot2')
  mov <- d %>%
    group_by(species, animal) %>%
    arrange(timestamp) %>%
    transmute(timestamp,
              time = hour(timestamp) + minute(timestamp) / 60,
              dt = as.numeric(timestamp - lag(timestamp)) / (1 %#% 'hours'),
              distance = sqrt((location.long - lag(location.long))^2 +
                                (location.lat - lag(location.lat))^2))
  
  # histograms of times between locations
  mov %>%
    ggplot() +
    facet_wrap(~ species, scales = 'free_y') +
    geom_histogram(aes(as.numeric(dt)), bins = 20) +
    scale_x_log10(expression(Time~between~locations~(log[10]~scale))) +
    ylab('Counts') +
    theme_bw()
  
  # histograms of animals' average times weighted by displacement
  summarize(mov,
            x = mean(time,
                     weights = distance / mean(distance, na.rm = TRUE),
                     na.rm = TRUE)) %>%
    ggplot() +
    facet_wrap(~ species, scales = 'free_y') +
    geom_histogram(aes(x))
  
  # hex plot of straight-line displacement by time of day
  mov %>%
    filter(distance < 100) %>% # remove locations in Germany
    ggplot() +
    facet_wrap(~ species, scales = 'free_y') +
    geom_hex(aes(time, distance), bins = 10) +
    labs(x = 'Time of day', y = 'SLD between consecutive locations') +
    theme_bw() +
    scale_fill_viridis_c('Count', limits = c(1, NA), trans = 'log10')
  
  mov %>%
    filter(distance < 100) %>% # remove locations in Germany
    ggplot() +
    facet_wrap(~ species, scales = 'free_y') +
    geom_smooth(aes(time, distance)) +
    labs(x = 'Time of day', y = 'SLD between consecutive locations') +
    scale_y_log10() +
    theme_bw() +
    scale_fill_viridis_c('Count', limits = c(1, NA), trans = 'log10')
}

# find unique time zones by species
d %>%
  group_by(species) %>%
  summarize(tz = paste(unique(timezone), collapse = ', '),
            tz2 = paste(unique(tz(timestamp)), collapse = ', '))

# convert times to UTC
d <- mutate(d, utc_timestamp = as.POSIXct(timestamp, 'UTC'))

# check times (removing NA timestamps in next script)
mean(with(d, timestamp - utc_timestamp) == 0, na.rm = TRUE)
sum(with(d, timestamp - utc_timestamp) == 0, na.rm = TRUE)

# save the full dataset ----
saveRDS(d, 'data/tracking-data/all-tracking-data-not-cleaned.rds')
