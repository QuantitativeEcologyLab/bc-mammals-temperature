library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for smother date wrangling
library('mgcv')      # for Generalized Additive Models
library('gratia')    # for useful convenience functions for GAMs
library('ggplot2')   # for fancy plots
library('khroma')    # for colorblind-friendly color palettes
library('cowplot')   # for fancy multi-panel plots
library('sf')        # for spatial data
library('ctmm')    # for movement models and utilization distributions
library('terra')   # for rasters
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids
plot_scheme(PAL, colours = TRUE)

# import data
d <- readRDS('data/hgam-speed-data.rds')

# import models
m_1 <- readRDS('models/binomial-gam.rds')
m_2 <- readRDS('models/gamma-gam.rds')
rsfs <- tibble(
  species = SPECIES,
  lab = SPECIES_LABS,
  file_name = map_chr(species, \(.sp) {
    fn <- list.files(path = 'H:/GitHub/bc-mammals-temperature/models/',
                     pattern = gsub('\\(s\\.', 'southern', .sp) %>%
                       gsub('\\(', '', .) %>%
                       gsub('\\)', '', .) %>%
                       paste0('rsf-', .),
                     full.names = TRUE)
    if(length(fn) == 0) fn <- NA_character_
    return(fn)
  })) %>%
  filter(! is.na(file_name)) %>%
  mutate(rsf = map(file_name, readRDS),
         local_newd = map(species, tibble))

# get range-resident animals
ANIMALS <- readr::read_csv('data/tracking-data/telemetry-metadata.csv',
                           show_col_types = FALSE) %>%
  filter(keep_for_RSF == 'yes') %>%
  pull(animal)

# get the extent of the range resident animals using they 99.9% UDs
uds <-
  readRDS('models/movement-models-akdes-2024-06-06.rds') %>%
  group_by(species) %>%
  # remove animals that are not range resident (invalid weights)
  filter(animal %in% ANIMALS) %>%
  mutate(species = if_else(grepl('Rangifer', dataset_name),
                           dataset_name,
                           species),
         species = stringr::str_replace_all(species, '_', ' ') %>%
           factor(),
         # get 99.9% UDs
         akde = map(akde, \(.a) {
           SpatialPolygonsDataFrame.UD(.a, level.UD = 0.999, level = 0) %>%
             st_as_sf() %>%
             filter(grepl('est', name)) %>%
             st_geometry() %>%
             st_set_crs(.a@info$projection) %>%
             st_transform('EPSG:4326') %>%
             st_as_sf()
         })) %>%
  # join within each species
  select(species, akde) %>%
  nest(akdes = ! species) %>%
  mutate(akdes = map(akdes, \(.a) {
    bind_rows(.a[[1]]) %>%
      st_make_valid() %>%
      st_union() %>%
      st_as_sf() 
  }))

# get resources for each species
if(FALSE) {
  res(rast('data/resource-rasters/forest.tif'))
  res(rast('data/resource-rasters/bc-dem-z6.tif')) # higher res than others
  res(rast('data/resource-rasters/bc-distance-from-water.tif'))
}

resources <- transmute(
  uds,
  species,
  r = map(akdes, \(.a) {
    forest <- rast('data/resource-rasters/forest.tif') %>%
      crop(.a) %>%
      mask(vect(.a))
    elev_m <- rast('data/resource-rasters/bc-buffered-dem-z6.tif') %>%
      project(forest)
    dist_w <- rast('data/resource-rasters/distance-from-water.tif') %>%
      project(forest)
    
    resources <-
      left_join(as.data.frame(forest, xy = TRUE) %>%
                  rename(forest_perc = 3),
                as.data.frame(elev_m, xy = TRUE) %>%
                  rename(elevation_m = 3),
                by = c('x', 'y')) %>%
      left_join(as.data.frame(dist_w, xy = TRUE) %>%
                  rename(dist_water_m = 3),
                by = c('x', 'y'))
  }))
resources

resources %>%
  unnest(r) %>%
  ggplot(aes(x, y, fill = forest_perc)) +
  facet_wrap(~ species, scales = 'free') +
  geom_raster() +
  labs(x = NULL, y = NULL) +
  scale_fill_gradient('Forest (%)', low = 'beige', high = 'forestgreen',
                      na.value = 'grey')

# extract temperatures for each species
temps <-
  readRDS('H:/GitHub/bc-mammals-temperature/data/weather-projections.rds')

gc() # clean up

unique_locs <-
  readRDS('data/climate-yearly-projections-2024-06-11.rds') %>%
  transmute(long = longitude, lat = latitude) %>%
  group_by(long, lat) %>%
  slice(1) %>%
  ungroup()

unique_locs_spp <- mutate(
  uds,
  subset = map(akdes, \(.ud) {
    unique_locs %>%
      transmute(
        long, lat, # drop all other columns
        inside =
          st_intersects(
            unique_locs %>%
              st_as_sf(coords = c('long', 'lat')) %>%
              st_set_crs('EPSG:4326'),
            .ud,
            sparse = TRUE) %>%
          as.numeric(),
        inside = if_else(inside > 0, TRUE, FALSE, missing = FALSE)) %>%
      filter(inside) %>%
      select(! inside)
  })) %>%
  select(! akdes) %>%
  unnest(subset)

ggplot(unique_locs_spp, aes(long, lat)) +
  facet_wrap(~ species, scales = 'free') +
  geom_point() +
  labs(x = NULL, y = NULL)

gc() # clean up

temps <- filter(temps,
                long >= min(unique_locs_spp$long),
                long <= max(unique_locs_spp$long),
                lat >= min(unique_locs_spp$lat),
                lat <= max(unique_locs_spp$lat))

map(unique(unique_locs_spp$species), \(.s) {
  .u <- unique_locs_spp %>%
    filter(species == .s)

  .t <- filter(temps,
               long >= min(.u$long),
               long <= max(.u$long),
               lat >= min(.u$lat),
               lat <= max(.u$lat))

  .u %>%
    left_join(.t, by = c('long', 'lat')) %>%
    saveRDS(paste0('data/weather-projections-', .s, '.rds'))
  
  return(as.character(.s))
})
