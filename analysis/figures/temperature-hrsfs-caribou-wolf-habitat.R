library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for smoother date wrangling
library('mgcv')      # for Generalized Additive Models
library('gratia')    # for useful convenience functions for GAMs
library('ggplot2')   # for fancy plots
library('khroma')    # for colorblind-friendly color palettes
library('cowplot')   # for fancy multi-panel plots
library('terra')   # for rasters
library('sf')      # for spatial objects
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids
source('functions/get_legend.R') # doesn't work in cowplot v.1.1.3.9000

theme_set(theme_get() + theme(text = element_text(size = 15)))

# boreal caribou RSF
rsf_bc <- readRDS('H:/GitHub/bc-mammals-temperature/models/rsf-Rangifer tarandus boreal-2025-01-21.rds')
rsf_bw <- readRDS('H:/GitHub/bc-mammals-temperature/models/rsf-Canis lupus-2025-01-20.rds')

# habitat selection strength for the current habitat ----
locs <- readRDS('data/tracking-data/all-tracking-data-cleaned-2024-02-22-13-49.rds') %>%
  filter(dataset_name == 'Rangifer_tarandus_boreal') %>%
  unnest(tel) %>%
  select(location.long, location.lat) %>%
  st_as_sf(coords = c('location.long', 'location.lat')) %>%
  st_set_crs('+proj=longlat')

bounds <- locs %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_as_sf() %>%
  st_transform('EPSG:3005')

f_rast <- rast('data/resource-rasters/forest.tif') %>%
  project('EPSG:3005') %>%
  crop(bounds)

f <- f_rast %>%
  as.data.frame(xy = TRUE) %>%
  rename(forest_perc = 3) %>%
  filter(! is.na(forest_perc))

e <- rast('data/resource-rasters/bc-buffered-dem-z6.tif') %>%
  crop(st_buffer(st_transform(bounds, crs(.)), 1e5)) %>%
  project(f_rast) %>%
  crop(bounds) %>%
  as.data.frame(xy = TRUE) %>%
  rename(elevation_m = 3) %>%
  filter(elevation_m > 0)

w <- rast('data/resource-rasters/distance-from-water.tif') %>%
  crop(st_buffer(st_transform(bounds, crs(.)), 1e5)) %>%
  project(f_rast) %>%
  crop(bounds) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE) %>%
  rename(dist_water_m = consensus_full_class_12)

resources <- left_join(f, e, by = c('x', 'y')) %>%
  left_join(w, by = c('x', 'y')) %>%
  as_tibble()

# figure of resources
p_resources <-
  cowplot::plot_grid(
    ggplot(resources, aes(x, y, fill = forest_perc)) +
      coord_sf(crs = 'EPSG:3005') +
      geom_raster() +
      scale_x_continuous(NULL, expand = c(0, 0), breaks = -c(124, 122, 120)) +
      scale_y_continuous(NULL, expand = c(0, 0)) +
      scale_fill_gradient('Forest cover (%)', low = 'white', na.value = NA,
                          high = 'darkgreen', limits = c(0, 100),
                          breaks = c(0, 100)) +
      theme(legend.position = 'top', legend.key.width = rel(0.7)),
    ggplot(resources, aes(x, y, fill = elevation_m / 1e3)) +
      coord_sf(crs = 'EPSG:3005') +
      geom_raster() +
      scale_x_continuous(NULL, expand = c(0, 0), breaks = -c(124, 122, 120)) +
      scale_y_continuous(NULL, expand = c(0, 0)) +
      scale_fill_distiller('Elevation (km) ', palette = 6, direction = 1,
                           breaks = round(range(resources$elevation_m / 1e3), 1),
                           labels = round(range(resources$elevation_m / 1e3), 1),
                           limits = round(range(resources$elevation_m / 1e3), 1)) +
      theme(legend.position = 'top', legend.key.width = rel(0.7)),
    ggplot(resources, aes(x, y, fill = dist_water_m / 1e3)) +
      coord_sf(crs = 'EPSG:3005') +
      geom_raster() +
      scale_x_continuous(NULL, expand = c(0, 0), breaks = -c(124, 122, 120)) +
      scale_y_continuous(NULL, expand = c(0, 0)) +
      scale_fill_distiller(expression(bold(atop(Distance~from,
                                                water~'(km)'~phantom(om)))),
                           na.value = NA, values = c(0, 0.05, 1),
                           limits = c(0, 18), breaks = c(0, 18)) +
      theme(legend.position = 'top', legend.key.width = rel(0.7)),
    nrow = 1)

# predict habitat selection strength at 20 degrees C ----
preds_habitat <-
  expand_grid(res = list(resources),
              temperature_C = c(-20, 0, 20)) %>%
  unnest(res) %>%
  mutate(species = 'Caribou (boreal)',
         animal = rsf_bc$model$animal[1]) %>%
  nest(dat = everything()) %>%
  bind_rows(., mutate(., dat = map(dat, \(.d) {
    mutate(.d,,
           species = c('Wolves'),
           animal = c(rsf_bw$model$animal[1]))
  }))) %>%
  mutate(dat = map2(dat, list(rsf_bc, rsf_bw), \(.d, .m) {
    .d %>%
      mutate(lambda = predict(object = .m, newdata = .,
                              se.fit = FALSE, type = 'response',
                              # exclude seasonal sampling bias
                              exclude = c('s(temperature_C)',
                                          's(temperature_C,animal)')))
  })) %>%
  unnest(dat) %>%
  mutate(temperature_C = paste0(temperature_C, '\U00B0', 'C') %>%
           factor(., levels = unique(.)))

# figure of habitat selection strength for different temperatures
# keeping elevations > 1200 m to make the point that predicting is hard
p_habitat <-
  preds_habitat %>%
  mutate(lambda = lambda * if_else(grepl('Caribou', species), 200, 1e4),
         lambda = if_else(lambda > 4, 4, lambda),
         lambda = if_else(lambda < 0.25, 0.25, lambda)) %>%
  ggplot(aes(x, y, fill = lambda)) +
  coord_sf(crs = 'EPSG:3005') +
  facet_grid(species ~ temperature_C) +
  geom_raster() +
  scale_x_continuous(NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, expand = c(0, 0)) +
  scale_fill_distiller('Relative selection strength   ', type = 'div',
                       palette = 3, direction = 1, trans = 'log2',
                       limits = c(0.25, 4),
                       breaks = c(0.25, 0.5, 1, 2, 4),
                       labels = c(0.25, 0.5, 1, 2, 4) %>%
                         as.character()) +
  theme(legend.position = 'top', legend.key.width = unit(0.7, 'in'))

p_enc <-
  preds_habitat %>%
  mutate(lambda = lambda * if_else(grepl('Caribou', species), 100, 5e3)) %>%
  select(x, y, temperature_C, species, lambda) %>%
  nest(r = ! c(species, temperature_C)) %>%
  mutate(r = map(r, \(.r) rast(.r, crs = 'EPSG:3005'))) %>%
  pivot_wider(names_from = species, values_from = r) %>%
  mutate(r = map2(`Caribou (boreal)`, Wolves, \(.c, .w) {
    as.data.frame(.c * .w, xy = TRUE) %>%
      rename(cooccupancy = lambda)
  })) %>%
  select(! c(`Caribou (boreal)`, Wolves)) %>%
  unnest(r) %>%
  mutate(cooccupancy = if_else(cooccupancy > 4, 4, cooccupancy),
         cooccupancy = if_else(cooccupancy < 0.25, 0.25, cooccupancy)) %>%
  ggplot(aes(x, y, fill = cooccupancy)) +
  coord_sf(crs = 'EPSG:3005') +
  facet_grid(. ~ temperature_C) +
  geom_raster() +
  scale_x_continuous(NULL, expand = c(0, 0), breaks = -c(124, 122, 120)) +
  scale_y_continuous(NULL, expand = c(0, 0)) +
  scale_fill_distiller('Relative co-occupancy rate   ', type = 'div',
                       palette = 2, direction = 1, trans = 'log2',
                       limits = c(0.25, 4),
                       breaks = c(0.25, 0.5, 1, 2, 4),
                       labels = c(0.25, 0.5, 1, 2, 4) %>%
                         as.character()) +
  theme(legend.position = 'top', legend.key.width = unit(0.7, 'in'))

p <- plot_grid(p_habitat,
               plot_grid(p_resources, p_enc, ncol = 1,
                         labels = c('B', 'C'), rel_heights = c(1, 1.2)),
               labels = c('A', ''), nrow = 1)

ggsave('figures/boreal-caribou-wolf-habitat.png', p,
       width = 20, height = 10, dpi = 600, bg = 'white')
