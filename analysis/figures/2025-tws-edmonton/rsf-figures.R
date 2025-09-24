library('mgcv')    # for GAMs
library('terra')   # for rasters
library('sf')      # for spatial objects
library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('ggplot2') # for fancy plots
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids
source('functions/get_legend.R') # doesn't work in cowplot v.1.1.3.9000

theme_set(theme_get() + theme(text = element_text(size = 15)))

rsf_bc <- readRDS('H:/GitHub/bc-mammals-temperature/models/rsf-Rangifer tarandus boreal-2025-01-21.rds')

# habitat selection strength for the current habitat ----
# create custom color palettes for each resource
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

# projections are unnecessary since we are not making maps
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

# ranges are quite wide for each predictor
resources %>%
  select(forest_perc, elevation_m, dist_water_m) %>%
  pivot_longer(everything()) %>%
  ggplot() +
  facet_wrap(~ name, scales = 'free') +
  geom_histogram(aes(value), bins = 10)

# predict habitat selection strength at 15 degrees C ----
preds_habitat <-
  expand_grid(res = list(resources),
              temperature_C = c(-15, 0, 15),
              animal = rsf_bc$model$animal[1]) %>%
  unnest(res) %>%
  mutate(lambda = predict(object = rsf_bc, newdata = ., se.fit = FALSE,
                          type = 'response'),
         temperature_C = paste0(temperature_C, '\U00B0', 'C') %>%
           factor(., levels = unique(.)))

# figure of resources
p_resources <-
  cowplot::plot_grid(
    ggplot(resources, aes(x, y, fill = forest_perc)) +
      coord_sf(crs = 'EPSG:3005') +
      geom_raster() +
      scale_x_continuous(NULL, expand = c(0, 0)) +
      scale_y_continuous(NULL, expand = c(0, 0)) +
      scale_fill_gradient('Tree cover (%)', low = 'white', na.value = NA,
                          high = 'darkgreen', limits = c(0, 100),
                          breaks = c(0, 100)) +
      theme(legend.position = 'top'),
    ggplot(resources, aes(x, y, fill = elevation_m / 1e3)) +
      coord_sf(crs = 'EPSG:3005') +
      geom_raster() +
      scale_x_continuous(NULL, expand = c(0, 0)) +
      scale_y_continuous(NULL, expand = c(0, 0)) +
      scale_fill_distiller('Elevation (km) ', palette = 6, direction = 1,
                           breaks = round(range(resources$elevation_m / 1e3), 1),
                           labels = round(range(resources$elevation_m / 1e3), 1),
                           limits = round(range(resources$elevation_m / 1e3), 1)) +
      theme(legend.position = 'top'),
    ggplot(resources, aes(x, y, fill = dist_water_m / 1e3)) +
      coord_sf(crs = 'EPSG:3005') +
      geom_raster() +
      scale_x_continuous(NULL, expand = c(0, 0)) +
      scale_y_continuous(NULL, expand = c(0, 0)) +
      scale_fill_distiller(expression(bold(atop(Distance~from,
                                                water~(km)~phantom(om)))),
                           na.value = NA, values = c(0, 0.05, 1),
                           limits = c(0, 18), breaks = c(0, 18)) +
      theme(legend.position = 'top'),
    nrow = 1); p_resources

ggsave('figures/2025-tws-edmonton/rsf-caribou-resources.png', p_resources,
       width = 12, height = 6, dpi = 600)

# figure of habitat selection strength for different temperatures
p_habitat <-
  preds_habitat %>%
  mutate(lambda = lambda * 100, # to make differences more visible
         lambda = if_else(lambda > 4, 4, lambda),
         lambda = if_else(lambda < 0.25, 0.25, lambda)) %>%
  ggplot(aes(x, y, fill = lambda)) +
  coord_sf(crs = 'EPSG:3005') +
  facet_wrap(~ temperature_C) +
  geom_raster() +
  scale_x_continuous(NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, expand = c(0, 0)) +
  scale_fill_distiller('Relative habitat preference   ', type = 'div',
                       palette = 3, direction = 1, trans = 'log2',
                       limits = c(0.25, 4),
                       breaks = c(0.25, 0.5, 1, 2, 4),
                       labels = c(0.25, 0.5, 1, 2, 4) %>%
                         as.character()) +
  theme(legend.position = 'top', legend.key.width = unit(0.7, 'in')); p_habitat

ggsave('figures/2025-tws-edmonton/rsf-caribou-range.png', p_habitat,
       width = 10, height = 6, dpi = 600)

# habitat preference relative to 0 degrees C
p_habitat_rel <-
  preds_habitat %>%
  # calculate values relative to 0 degrees C
  pivot_wider(names_from = temperature_C, values_from = lambda) %>%
  mutate(`-15°C` = `-15°C` / `0°C`,
         `15°C` =`15°C` / `0°C`,
         `0°C` = `0°C` /`0°C`) %>%
  pivot_longer(c(`-15°C`, `0°C`, `15°C`), names_to = 'temperature',
               values_to = 'lambda') %>%
  mutate(temperature = factor(temperature, levels = unique(temperature)),
         lambda = if_else(lambda > 4, 4, lambda),
         lambda = if_else(lambda < 0.25, 0.25, lambda)) %>%
  ggplot(aes(x, y, fill = lambda)) +
  coord_sf(crs = 'EPSG:3005') +
  facet_wrap(~ temperature) +
  geom_raster() +
  scale_x_continuous(NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, expand = c(0, 0)) +
  scale_fill_distiller(expression(atop(bold('Habitat preference   '),
                                       bold('relative to 0\U00B0', 'C'))),
                       type = 'div',
                       palette = 2, direction = 1, trans = 'log2',
                       limits = c(0.25, 4),
                       breaks = c(0.25, 0.5, 1, 2, 4),
                       labels = c(0.25, 0.5, 1, 2, 4) %>%
                         as.character()) +
  theme(legend.position = 'top', legend.key.width = unit(0.7, 'in'))

ggsave('figures/2025-tws-edmonton/rsf-caribou-range-0C.png', p_habitat_rel,
       width = 10, height = 6, dpi = 600)

readRDS('H:/GitHub/bc-mammals-temperature/data/cc-hrsf-projections-local-2100.rds') %>%
  filter(lab == 'bolditalic(Rangifer~tarandus)~bold(\"(boreal)\")') %>%
  mutate(l = if_else(l > 1.066, 1.066, l),
         scenario = gsub('"', '', scenario) %>% factor(levels = unique(.))) %>%
  ggplot(aes(x = l, fill = scenario, color = scenario)) +
  geom_density(alpha = 0.25, lwd = 1) +
  geom_vline(xintercept = 1, color = 'black', lty = 'dashed') +
  scale_x_continuous('Relative change in RSS in 2100') +
  ylab('Density') +
  scale_color_brewer('Climate change scenario', type = 'div', palette = 5,
                     direction = -1, aesthetics = c('color', 'fill')) +
  theme(legend.position = 'inside', legend.position.inside = c(3/4, 3/4))

ggsave('figures/2025-tws-edmonton/boreal-caribou-local-rss-2100-density.png',
       width = 10, height = 5, dpi = 600, bg = 'white')

# wolf figures
rsf_bw <- readRDS('H:/GitHub/bc-mammals-temperature/models/rsf-Canis lupus-2025-01-20.rds')

preds_habitat_bw <-
  expand_grid(res = list(resources),
              temperature_C = c(-15, 0, 15),
              animal = rsf_bw$model$animal[1]) %>%
  unnest(res) %>%
  mutate(lambda = predict(object = rsf_bw, newdata = ., se.fit = FALSE,
                          type = 'response'),
         temperature_C = paste0(temperature_C, '\U00B0', 'C') %>%
           factor(., levels = unique(.)))

p_habitat_bw <-
  preds_habitat_bw %>%
  mutate(lambda = lambda * 1e3, # to make differences more visible
         lambda = if_else(lambda > 4, 4, lambda),
         lambda = if_else(lambda < 0.25, 0.25, lambda)) %>%
  ggplot(aes(x, y, fill = lambda)) +
  coord_sf(crs = 'EPSG:3005') +
  facet_wrap(~ temperature_C) +
  geom_raster() +
  scale_x_continuous(NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, expand = c(0, 0)) +
  scale_fill_distiller('Relative habitat preference   ', type = 'div',
                       palette = 3, direction = 1, trans = 'log2',
                       limits = c(0.25, 4),
                       breaks = c(0.25, 0.5, 1, 2, 4),
                       labels = c(0.25, 0.5, 1, 2, 4) %>%
                         as.character()) +
  theme(legend.position = 'top', legend.key.width = unit(0.7, 'in'))

ggsave('figures/2025-tws-edmonton/rsf-wolf-range.png', p_habitat_bw,
       width = 10, height = 6, dpi = 600)

plot_grid(get_legend(p_habitat),
          plot_grid(p_habitat + theme(legend.position = 'none'),
                    p_habitat_bw + theme(legend.position = 'none'),
                    labels = c('A', 'B'), nrow = 1),
          ncol = 1, rel_heights = c(0.15, 1))

ggsave('figures/2025-tws-edmonton/rsf-both-range.png',
       width = 20, height = 5.5, dpi = 600, bg = 'white')
