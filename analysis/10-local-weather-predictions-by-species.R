library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for smoother date wrangling
library('mgcv')      # for Generalized Additive Models
library('gratia')    # for useful convenience functions for GAMs
library('ggplot2')   # for fancy plots
library('khroma')    # for colorblind-friendly color palettes
library('cowplot')   # for fancy multi-panel plots
library('sf')        # for spatial data
library('terra')     # for rasters
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids
plot_scheme(PAL, colours = TRUE)

# get range-resident animals
ANIMALS <- readr::read_csv('data/tracking-data/telemetry-metadata.csv',
                           show_col_types = FALSE) %>%
  filter(keep_for_RSF == 'yes') %>%
  pull(animal)

# get the extent of the range resident animals using they 99.9% UDs
# ctmm is currently not installed on the EME linux because it needs more
# linux dependencies
if(file.exists('data/species-uds-99.9-percent.rds')) {
  uds <- readRDS('data/species-uds-99.9-percent.rds')
} else {
  library('ctmm') # for movement models and utilization distributions
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
  saveRDS(uds, 'data/species-uds-99.9-percent.rds')
}

ud_bbox <- unnest(uds, akdes) %>% st_as_sf() %>% st_bbox()

# extract temperatures for each species
temps <- readRDS('data/weather-projections.rds') %>%
  # filtering to the extent of the observed ranges in BC
  filter(long >= ud_bbox['xmin'], long <= ud_bbox['xmax'],
         lat >= ud_bbox['ymin'], lat <= ud_bbox['ymax'])

# find unique locations in BC
unique_locs <-
  temps %>%
  select(long, lat) %>%
  slice(1, .by = c(long, lat))
plot(unique_locs)

# find unique BC locations within each species' observed range
unique_locs_spp <- mutate(
  uds,
  subset = map(akdes, \(.ud) {
    unique_locs %>%
      transmute(
        long, lat, # drop all other columns
        inside =
          st_intersects(unique_locs %>%
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

plot_grid(plotlist = map(uds$species, function(.sp) {
  ggplot() +
    geom_sf(aes(geometry = x), unnest(filter(uds, species == .sp), akdes))+
    geom_point(aes(long, lat), filter(unique_locs_spp, species == .sp),
               size = 0.1) +
    geom_point(aes(long, lat), filter(unique_locs_spp, species == .sp)) +
    labs(x = NULL, y = NULL, title = .sp)
}))

map(unique(unique_locs_spp$species)[5], \(.s) {
  # find the locations for the species of interest
  .u <- filter(unique_locs_spp, species == .s)
  
  # only keep the locations within each species' set of 99.9% UDs
  # filter to the extent of the species for faster wrangling
  filter(temps,
         long >= min(.u$long) - 1, # widen the bounds by 2 degrees
         long <= max(.u$long) + 1,
         lat >= min(.u$lat) - 1,
         lat <= max(.u$lat) + 1) %>%
    right_join(.u, by = c('long', 'lat')) %>%
    saveRDS(paste0('data/weather-projections-', .s, '.rds'))
  
  return(as.character(.s))
}, .progress = TRUE)
