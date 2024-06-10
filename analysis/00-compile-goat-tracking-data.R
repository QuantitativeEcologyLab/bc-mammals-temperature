library('dplyr') # for data wrangling
library('purrr') # for functional programming
library('tidyr') # for data wangling
library('sf')    # for spatial features

# check layers
st_layers('data/tracking-data/oreamnos-americanus/GPS_Collar30548_20231005121804.gpx')
#' use `track_points` layer

# import all gpx files
d <- tibble(
  file = list.files('data/tracking-data/oreamnos-americanus'),
  track = map(file, function(fn) {
    st_read(paste0('data/tracking-data/oreamnos-americanus/', fn),
            layer = 'track_points', quiet = TRUE) %>%
      # add animal identifier with ctmm syntax
      mutate(individual.local.identifier = substr(fn, nchar('GPS_CollarX'),
                                                  nchar('GPS_CollarXXXXX'))) %>%
      bind_cols(., st_coordinates(.)) %>% # add coordinates
      st_drop_geometry() %>% # drop sf geometry column
      rename(location.long = X, # rename using ctmm syntax
             location.lat = Y,
             timestamp = time)
  })) %>%
  unnest(track) # turn into a single large tibble

# check which columns have more than one unique (non-NA) value
d %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(-c()) %>%
  group_by(name) %>%
  summarize(n = n_distinct(value, na.rm = TRUE)) %>%
  mutate(useful = n > 1) %>%
  arrange(desc(n)) %>%
  print(n = 28)

# drop non-NA columns with more than one unique value
d <- d %>%
  select(individual.local.identifier, location.long, location.lat,
         timestamp, fix, hdop, pdop)

saveRDS(d, 'data/tracking-data/oreamnos-americanus.rds')
