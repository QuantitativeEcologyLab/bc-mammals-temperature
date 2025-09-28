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

#' created in `analysis/figures/temperature-hrsfs.R`
preds <- readRDS('models/temperature-hrsf-preds.rds') %>%
  filter(species == 'Rangifer tarandus (boreal)') %>%
  # add common species names
  mutate(lab = COMMON_NAMES[map_int(species, \(.s) {
    which(sort(SPECIES) == .s)
  })],
  variable = variable %>%
    gsub('bold\\(', '', .) %>%
    gsub('~', ' ', .) %>%
    gsub('"', '', .) %>%
    gsub('\\)\\)', '\\)', .) %>%
    factor(levels = c('Forest cover (%)',
                      'Elevation (km)',
                      'Distance from water (km)')))

slice(preds, 1)

# plot of RSS ----
LIM <- 2

p <-
  preds %>%
  filter((! too_far_2)) %>%
  select(lab, x, temperature_C, variable, lambda, too_far_1) %>%
  # make detection rates homogeneous across temperature
  group_by(lab, variable, temperature_C) %>%
  mutate(lambda = lambda / mean(lambda)) %>%
  # re-scale the center (lambda = 1) relative to the median in the data.
  # doing this because some estimated effects of distance from water and
  # elevation are extreme
  group_by(lab, variable) %>%
  mutate(lambda = lambda / median(lambda),
         lambda = if_else(variable == 'Elevation (km)', lambda * 2^(2), lambda)) %>%
  ungroup() %>%
  # cap at 2^(+/-LIM)
  mutate(log2_lambda = log2(lambda),
         log2_lambda = case_when(log2_lambda > LIM ~ LIM,
                                 log2_lambda < -LIM ~ -LIM,
                                 TRUE ~ log2_lambda)) %>%
  ggplot() +
  facet_wrap(. ~ variable, scales = 'free', strip.position = 'left') +
  geom_raster(aes(temperature_C, x, fill = log2_lambda)) +
  geom_contour(aes(temperature_C, x, z = as.numeric(too_far_1)),
               color = 'grey50', linewidth = 0.25) +
  geom_contour(aes(temperature_C, x, z = log2_lambda),
               color = 'black', breaks = c(-2, -1, 0, 1, 2)) +
  scale_x_continuous(paste0('Temperature (\U00B0', 'C)'), expand = c(0, 0),
                     breaks = c(-20, 0, 20)) +
  scale_y_continuous(NULL, expand = c(0, 0)) +
  scale_fill_PRGn(name = 'Relative selection strength   ', midpoint = 0,
                  limits = c(-LIM, LIM), breaks = -LIM:LIM,
                  labels = \(x) 2^x) +
  theme(strip.placement = 'outside', strip.background.y = element_blank(),
        strip.text.y = element_text(size = 11), legend.position = 'top',
        panel.background = element_rect(fill = 'grey90'),
        legend.key.width = rel(2))
p

ggsave('figures/2025-tws-edmonton/boreal-caribou-hrsf-surface-plots.png',
       p, width = 10, height = 4.375, units = 'in', dpi = 600, bg = 'white')

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
  scale_fill_distiller('Relative selection strength   ', type = 'div',
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

# density plot for RSS in 2100
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
  scale_fill_distiller('Relative selection strength   ', type = 'div',
                       palette = 3, direction = 1, trans = 'log2',
                       limits = c(0.25, 4),
                       breaks = c(0.25, 0.5, 1, 2, 4),
                       labels = c(0.25, 0.5, 1, 2, 4) %>%
                         as.character()) +
  theme(legend.position = 'top', legend.key.width = unit(0.7, 'in'))

ggsave('figures/2025-tws-edmonton/rsf-wolf-range.png', p_habitat_bw,
       width = 10, height = 6, dpi = 600)

plot_grid(plot_grid(p_habitat + theme(legend.position = 'none'),
                    p_habitat_bw + theme(legend.position = 'none'),
                    labels = c('A', 'B'), ncol = 1),
          get_legend(p_habitat +
                       scale_fill_distiller('RSS', type = 'div',
                                            palette = 3, direction = 1, trans = 'log2',
                                            limits = c(0.25, 4),
                                            breaks = c(0.25, 0.5, 1, 2, 4),
                                            labels = c(0.25, 0.5, 1, 2, 4) %>%
                                              as.character()) +
                       theme(legend.key.width = rel(1),
                                       legend.key.height = rel(2)),
                     position = 'right'),
          nrow = 1, rel_widths = c(1, 0.15))

ggsave('figures/2025-tws-edmonton/rsf-both-range.png',
       width = 11.5, height = 10, dpi = 600, bg = 'white')
