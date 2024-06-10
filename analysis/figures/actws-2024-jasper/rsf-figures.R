library('mgcv')    # for GAMs
library('terra')   # for rasters
library('sf')      # for spatial objects
library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('ggplot2') # for fancy plots
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids

theme_set(theme_get() + theme(text = element_text(size = 15)))

rsf <- readRDS('models/rsf-boreal-caribou-2024-03-03.rds')

# marginal effect of forest ----
# empty figure
tibble(temperature = 15,
       forest_perc = seq(0, 100, by = 0.1),
       elevation_m = 1,
       dist_water_m = 1) %>%
  bind_cols(.,
            predict(rsf, newdata = .,
                    terms = c('s(forest_perc)',
                              's(forest_perc,temperature)'),
                    type = 'link', se.fit = TRUE)) %>%
  mutate(mu = exp(fit),
         lwr = exp(fit - 1.96 * se.fit),
         upr = exp(fit + 1.96 * se.fit)) %>%
  ggplot() +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(forest_perc, ymin = lwr, ymax = upr),
              fill = 'transparent', alpha = 0.2) +
  geom_line(aes(forest_perc, mu), color = 'transparent', linewidth =  1) +
  scale_y_continuous('Relative habitat preference', trans = 'log2',
                     breaks = c(0.25, 0.5, 1, 2, 4)) +
  xlab('Forest (%)')

ggsave('figures/actws-2024-jasper/rsf-forest-empty.png',
       width = 10, height = 6, dpi = 600)

# real figure
tibble(temperature = 15,
       forest_perc = seq(0, 100, by = 0.1),
       elevation_m = 1,
       dist_water_m = 1) %>%
  bind_cols(.,
            predict(rsf, newdata = .,
                    terms = c('s(forest_perc)',
                              's(forest_perc,temperature)'),
                    type = 'link', se.fit = TRUE)) %>%
  mutate(mu = exp(fit),
         lwr = exp(fit - 1.96 * se.fit),
         upr = exp(fit + 1.96 * se.fit)) %>%
  ggplot() +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(forest_perc, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(forest_perc, mu), color = '#EE6677', linewidth =  1) +
  scale_y_continuous('Relative habitat preference', trans = 'log2',
                     breaks = c(0.25, 0.5, 1, 2, 4)) +
  xlab('Forest (%)')

ggsave('figures/actws-2024-jasper/rsf-forest.png',
       width = 10, height = 6, dpi = 600)

# marginal effect of elevation ----
elevs <- seq(200, 1100, length.out = 400) # used in later plot, too

tibble(temperature = 15,
       forest_perc = 0,
       elevation_m = elevs,
       dist_water_m = 1) %>%
  bind_cols(.,
            predict(rsf, newdata = .,
                    terms = c('s(elevation_m)',
                              'ti(elevation_m,temperature)'),
                    type = 'link', se.fit = TRUE)) %>%
  mutate(mu = exp(fit),
         lwr = exp(fit - 1.96 * se.fit),
         upr = exp(fit + 1.96 * se.fit)) %>%
  ggplot() +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(elevation_m, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(elevation_m, mu), color = '#EE6677', linewidth =  1) +
  scale_y_continuous('Relative habitat preference', trans = 'log2',
                     breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
  xlab('Elevation (m)')

ggsave('figures/actws-2024-jasper/rsf-elevation.png',
       width = 10, height = 6, dpi = 600)

# marginal effect of distance from water ----
distances <- seq(0, 12e3, length.out = 400)

tibble(temperature = 15,
       forest_perc = 0,
       elevation_m = 0,
       dist_water_m = distances) %>%
  bind_cols(.,
            predict(rsf, newdata = .,
                    terms = c('s(sqrt(dist_water_m))'),
                    type = 'link', se.fit = TRUE)) %>%
  mutate(dist_water_km = dist_water_m / 1e3,
         mu = exp(fit),
         lwr = exp(fit - 1.96 * se.fit),
         upr = exp(fit + 1.96 * se.fit)) %>%
  ggplot() +
  coord_cartesian(ylim = c(0.5, 2)) +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(dist_water_km, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(dist_water_km, mu), color = '#EE6677', linewidth =  1) +
  scale_y_continuous('Relative habitat preference', trans = 'log2',
                     breaks = c(0.25, 0.5, 1, 2, 4)) +
  scale_x_continuous('Distance from water (km)', limits = c(0, 12.5))

ggsave('figures/actws-2024-jasper/rsf-water.png',
       width = 10, height = 6, dpi = 600)

# effect of forest at different temperatures ----
expand.grid(temperature = c(-15, 0, 15),
            forest_perc = seq(0, 100, by = 0.1),
            elevation_m = 0,
            dist_water_m = 0) %>%
  bind_cols(.,
            predict(rsf, newdata = .,
                    terms = c('s(forest_perc)', 'ti(forest_perc,temperature)'),
                    type = 'link', se.fit = TRUE)) %>%
  mutate(temperature = paste0(temperature, '\U00B0', 'C'),
         mu = exp(fit),
         lwr = exp(fit - 1.96 * se.fit),
         upr = exp(fit + 1.96 * se.fit)) %>%
  ggplot() +
  facet_wrap(~ temperature) +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(forest_perc, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(forest_perc, mu), color = '#EE6677', linewidth =  1) +
  scale_y_continuous('Relative habitat preference', trans = 'log2',
                     breaks = c(0.25, 0.5, 1, 2, 4)) +
  xlab('Forest (%)')

ggsave('figures/actws-2024-jasper/rsf-forest-temperature.png',
       width = 15, height = 6, dpi = 600)

# effect of elevation at different temperatures ----
expand.grid(temperature = c(-15, 0, 15),
            forest_perc = 0,
            elevation_m = elevs,
            dist_water_m = 0) %>%
  filter(! exclude.too.far(g1 = temperature, g2 = elevation_m,
                           d1 = rsf$model$temperature,
                           d2 = rsf$model$elevation_m,
                           dist = 0.1)) %>%
  bind_cols(.,
            predict(rsf, newdata = .,
                    terms = c('s(elevation_m)', 'ti(elevation_m,temperature)'),
                    type = 'link', se.fit = TRUE)) %>%
  mutate(temperature = paste0(temperature, '\U00B0', 'C'),
         mu = exp(fit),
         lwr = exp(fit - 1.96 * se.fit),
         upr = exp(fit + 1.96 * se.fit)) %>%
  ggplot() +
  facet_wrap(~ temperature) +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(elevation_m, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(elevation_m, mu), color = '#EE6677', linewidth =  1) +
  scale_y_continuous('Relative habitat preference', trans = 'log2',
                     breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8)) +
  xlab('Elevation (m)')

ggsave('figures/actws-2024-jasper/rsf-elevation-temperature.png',
       width = 15, height = 6, dpi = 600)

# effect of distance from water at different temperatures ----
expand.grid(temperature = c(-15, 0, 15),
            forest_perc = 0,
            elevation_m = 0,
            dist_water_m = distances) %>%
  bind_cols(.,
            predict(rsf, newdata = .,
                    terms = c('s(sqrt(dist_water_m))', 'ti(sqrt(dist_water_m),temperature)'),
                    type = 'link', se.fit = TRUE)) %>%
  mutate(temperature = paste0(temperature, '\U00B0', 'C'),
         dist_water_km = dist_water_m / 1e3,
         mu = exp(fit),
         lwr = exp(fit - 1.96 * se.fit),
         upr = exp(fit + 1.96 * se.fit)) %>%
  ggplot() +
  facet_wrap(~ temperature) +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(dist_water_km, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(dist_water_km, mu), color = '#EE6677', linewidth =  1) +
  scale_y_continuous('Relative habitat preference', trans = 'log2',
                     breaks = c(0.25, 0.5, 1, 2, 4)) +
  scale_x_continuous('Distance from water (km)', limits = c(0, 12.5))

ggsave('figures/actws-2024-jasper/rsf-water-temperature.png',
       width = 15, height = 6, dpi = 600)

# habitat quality for the current habitat ----
# create custom color palettes for each resource
bounds <- readRDS('data/tracking-data/all-tracking-data-cleaned-2024-02-22-13-49.rds') %>%
  filter(dataset_name == 'Rangifer_tarandus_boreal') %>%
  unnest(tel) %>%
  select(location.long, location.lat) %>%
  st_as_sf(coords = c('location.long', 'location.lat')) %>%
  st_set_crs('+proj=longlat') %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_as_sf()

# projections are unnecessary since we are not making maps
f_rast <- rast('data/resource-rasters/bc-forest.tif') %>%
  crop(bounds)

f <- f_rast %>%
  as.data.frame(xy = TRUE) %>%
  filter(! is.na(layer)) %>%
  rename(forest_perc = layer)

e <- rast('data/resource-rasters/bc-dem-z6.tif') %>%
  crop(st_buffer(bounds, 1e5)) %>%
  project(f_rast) %>%
  crop(bounds) %>%
  as.data.frame(xy = TRUE) %>%
  rename(elevation_m = 3) %>%
  filter(elevation_m > 0)

w <- rast('data/resource-rasters/bc-distance-from-water.tif') %>%
  crop(st_buffer(bounds, 1e5)) %>%
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

# predict habitat quality at 15 degrees C ----
preds_habitat <-
  expand_grid(res = list(resources),
              temperature = c(-15, 0, 15)) %>%
  unnest(res) %>%
  mutate(lambda = predict(object = rsf, newdata = ., se.fit = FALSE,
                          type = 'response', exclude = '(Intercept)'),
         temperature = paste0(temperature, '\U00B0', 'C') %>%
           factor(., levels = unique(.)))

# figure of resources
p_resources <-
  cowplot::plot_grid(
    ggplot(resources, aes(x, y, fill = forest_perc)) +
      coord_sf(crs = '+proj=longlat') +
      geom_raster() +
      scale_x_continuous(NULL, expand = c(0, 0)) +
      scale_y_continuous(NULL, expand = c(0, 0)) +
      scale_fill_gradient('Tree cover (%)', low = 'white', na.value = NA,
                          high = 'darkgreen', limits = c(0, 100),
                          breaks = c(0, 100)) +
      theme(legend.position = 'top'),
    ggplot(resources, aes(x, y, fill = elevation_m)) +
      coord_sf(crs = '+proj=longlat') +
      geom_raster() +
      scale_x_continuous(NULL, expand = c(0, 0)) +
      scale_y_continuous(NULL, expand = c(0, 0)) +
      scale_fill_distiller('Elevation (m) ', palette = 6, direction = 1,
                           breaks = round(range(resources$elevation_m)),
                           labels = round(range(resources$elevation_m), -2),
                           limits = round(range(resources$elevation_m))) +
      theme(legend.position = 'top'),
    ggplot(resources, aes(x, y, fill = dist_water_m / 1e3)) +
      coord_sf(crs = '+proj=longlat') +
      geom_raster() +
      scale_x_continuous(NULL, expand = c(0, 0)) +
      scale_y_continuous(NULL, expand = c(0, 0)) +
      scale_fill_distiller(expression(bold(atop(Distance~from,
                                                water~(km)~phantom(om)))),
                           na.value = NA, values = c(0, 0.05, 1),
                           limits = c(0, 18), breaks = c(0, 18)) +
      theme(legend.position = 'top'),
    nrow = 1); p_resources

ggsave('figures/actws-2024-jasper/rsf-caribou-resources.png', p_resources,
       width = 12, height = 6, dpi = 600)

# figure of habitat quality for different temperatures
p_habitat <-
  preds_habitat %>%
  mutate(lambda = lambda / 4, # to make differences more visible
         lambda = if_else(lambda > 4, 4, lambda),
         lambda = if_else(lambda < 0.25, 0.25, lambda)) %>%
  ggplot(aes(x, y, fill = lambda)) +
  coord_sf(crs = '+proj=longlat') +
  facet_wrap(~ temperature) +
  geom_raster() +
  scale_x_continuous(NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, expand = c(0, 0)) +
  scale_fill_distiller('Relative habitat preference   ', type = 'div',
                       palette = 1, direction = 1, trans = 'log2',
                       limits = c(0.25, 4),
                       breaks = c(0.25, 0.5, 1, 2, 4),
                       labels = c(0.25, 0.5, 1, 2, 4) %>%
                         as.character()) +
  theme(legend.position = 'top', legend.key.width = unit(0.7, 'in'))

ggsave('figures/actws-2024-jasper/rsf-caribou-range.png', p_habitat,
       width = 10, height = 6, dpi = 600)

# habitat preference relative to 0 degrees C
p_habitat_rel <-
  preds_habitat %>%
  mutate(lambda = lambda / 4) %>% # to make differences more visible
  pivot_wider(names_from = temperature, values_from = lambda) %>%
  mutate(`-15°C` = `-15°C` / `0°C`,
         `0°C` = `0°C` /`0°C`,
         `15°C` =`15°C` / `0°C`) %>%
  pivot_longer(c(`-15°C`, `0°C`, `15°C`), names_to = 'temperature',
               values_to = 'lambda') %>%
  mutate(temperature = factor(temperature, levels = unique(temperature)),
         lambda = if_else(lambda > 4, 4, lambda),
         lambda = if_else(lambda < 0.25, 0.25, lambda)) %>%
  ggplot(aes(x, y, fill = lambda)) +
  coord_sf(crs = '+proj=longlat') +
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

ggsave('figures/actws-2024-jasper/rsf-caribou-range-0C.png', p_habitat_rel,
       width = 10, height = 6, dpi = 600)

# predict using climate change predictions ----
# projections are unnecessary since we are not making maps
f_rast <- rast('data/resource-rasters/bc-forest.tif') %>%
  crop(bounds)

e_rast <- rast('data/resource-rasters/bc-dem-z6.tif') %>%
  crop(bounds)

w_rast <- rast('data/resource-rasters/bc-distance-from-water.tif') %>%
  crop(bounds)

ecmwf_preds <-
  readRDS('data/climate-yearly-projections-2024-03-03.rds') %>%
  filter(longitude > st_bbox(bounds)['xmin'],
         longitude < st_bbox(bounds)['xmax'],
         latitude > st_bbox(bounds)['ymin'],
         latitude < st_bbox(bounds)['ymax']) %>%
  mutate(forest_perc = terra::extract(f_rast,
                                      tibble(longitude, latitude))[, 2],
         elevation_m = terra::extract(e_rast,
                                      tibble(longitude, latitude))[, 2],
         dist_water_m = terra::extract(w_rast,
                                       tibble(longitude, latitude))[, 2]) %>%
  filter(! is.na(forest_perc))

preds_cc <-
  ecmwf_preds %>%
  mutate(temperature = mean_temperature) %>%
  mutate(lambda = predict(object = rsf, newdata = ., se.fit = FALSE,
                          type = 'response', exclude = '(Intercept)')) %>%
  group_by(scenario, year, month) %>%
  summarize(lambda = mean(lambda, na.rm = TRUE), .groups = 'drop') %>%
  group_by(scenario, year) %>%
  summarize(lambda = mean(lambda, na.rm = TRUE), .groups = 'drop') %>%
  arrange(year) %>%
  mutate(lambda = lambda / first(lambda)) %>%
  mutate(scenario = case_when(
           scenario == '8GCMs_ensemble_ssp126' ~ 'Best',
           scenario == '8GCMs_ensemble_ssp245' ~ 'Good',
           scenario == '8GCMs_ensemble_ssp370' ~ 'Bad',
           scenario == '8GCMs_ensemble_ssp585' ~ 'Worst') %>%
           factor(levels = unique(.)))

p_preds_cc <-
  ggplot(preds_cc) +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_point(aes(year, lambda, color = scenario), alpha = 0.5) +
  geom_smooth(aes(year, lambda, color = scenario, fill = scenario),
              method = 'gam', formula = y ~ s(x, k = 5)) +
  scale_fill_brewer('Climate change scenario',
                    type = 'div', palette = 5, direction = -1,
                    aesthetics = c('color', 'fill')) +
  labs(x = NULL, y = 'Relative habitat quality') +
  scale_y_continuous(limits = c(0.5, 2), trans = 'log2') +
  theme(legend.position = 'top'); p_preds_cc

ggsave('figures/actws-2024-jasper/rsf-caribou-cc-preds.png', p_preds_cc,
       width = 10, height = 6, dpi = 600)
