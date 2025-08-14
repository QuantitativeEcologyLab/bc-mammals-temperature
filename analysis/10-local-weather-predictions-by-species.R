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
library('ctmm')      # for movement models and utilization distributions
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

# find unique locations in BC
unique_locs <-
  readRDS('data/climate-yearly-projections-2025-2100-only-2025-01-21.rds') %>%
  transmute(long = longitude, lat = latitude) %>%
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
    geom_point(aes(long, lat), filter(unique_locs_spp, species == .sp)) +
    labs(x = NULL, y = NULL, title = .sp)
}))

gc() # clean up

# extract temperatures for each species
temps <- readRDS('data/weather-projections.rds') %>%
  # cut the dataset size in half by filtering to the observed ranges in BC
  filter(long >= min(unique_locs_spp$long),
         long <= max(unique_locs_spp$long),
         lat >= min(unique_locs_spp$lat),
         lat <= max(unique_locs_spp$lat))

map(unique(unique_locs_spp$species), \(.s) {
  # find the locations for the species of interest
  .u <- filter(unique_locs_spp, species == .s)
  
  # filter to the extent of the species for faster wrangling
  .t <- filter(temps,
               long >= min(.u$long),
               long <= max(.u$long),
               lat >= min(.u$lat),
               lat <= max(.u$lat))
  
  # only keep the locations within each species set of 99.9% UDs
  left_join(.u, .t, by = c('long', 'lat')) %>%
    saveRDS(paste0('data/weather-projections-', .s, '.rds'))
  
  return(as.character(.s))
})
