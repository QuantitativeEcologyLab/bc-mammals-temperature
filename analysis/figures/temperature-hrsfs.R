library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for smoother date wrangling
library('mgcv')      # for Generalized Additive Models
library('gratia')    # for useful convenience functions for GAMs
library('ggplot2')   # for fancy plots
library('khroma')    # for colorblind-friendly color palettes
library('cowplot')   # for fancy multi-panel plots
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids

# import models for each species
d <- tibble(
  species = SPECIES,
  lab = SPECIES_LABS,
  file_name = map_chr(species, \(.sp) {
    fn <- list.files(path = 'H:/GitHub/bc-mammals-temperature/models/',
                     pattern = gsub('\\(s\\.', 'southern', .sp) %>%
                       gsub('\\(', '', .) %>%
                       gsub('\\)', '', .) %>%
                       paste0('rsf-', ., '-2024-'),
                     full.names = TRUE)
    if(length(fn) == 0) fn <- NA_character_
    return(fn)
  })) %>%
  filter(! is.na(file_name)) %>%
  mutate(rsf = map(file_name, readRDS),
         local_newd = map(species, tibble))
d

# get range of predictors
if(FALSE) {
  readRDS('data/tracking-data/rsf-data.rds') %>%
    filter(detected == 1) %>%
    summarize(min_elevation_m = min(elevation_m),
              max_elevation_m = max(elevation_m),
              min_dist_water_m = min(dist_water_m),
              max_dist_water_m = max(dist_water_m),
              min_temperature_C = min(temperature_C),
              max_temperature_C = max(temperature_C))
}

# surface plots of partial effects ----
newd <-
  tibble(animal = 'new animal',
         temperature_C = seq(-40, 40, length.out = 400),
         x = tibble(forest_perc = seq(0, 100, length.out = 500),
                    elevation_m = seq(0, 3000, length.out = 500),
                    dist_water_m = seq(0, 25e3, length.out = 500)) %>%
           list()) %>%
  unnest(x)

surface <- function(m, dist_1 = 0.01, dist_2 = 0.1) {
  bind_rows(
    transmute(
      newd,
      temperature_C,
      x = forest_perc,
      variable = 'bold(Forest~cover~(\'\U0025\'))',
      lambda =
        predict(object = m, newdata = newd, type = 'response',
                se.fit = FALSE, discrete = FALSE, newdata.guaranteed = TRUE,
                terms = c('s(forst_perc)',
                          'ti(forest_perc,temperature_C)')),
      too_far_1 = too_far(temperature_C, x,
                          filter(m$model, detected == 1)$temperature_C,
                          filter(m$model, detected == 1)$forest_perc,
                          dist = dist_1),
      # too_far_2 = too_far(temperature_C, x,
      #                     filter(m$model, detected == 1)$temperature_C,
      #                     filter(m$model, detected == 1)$forest_perc,
      #                     dist = dist_2)
      too_far_2 = FALSE),
    transmute(
      newd,
      temperature_C,
      x = elevation_m,
      variable = 'bold(Elevation~(km))',
      lambda =
        predict(object = m, newdata = newd, type = 'response',
                se.fit = FALSE, discrete = FALSE, newdata.guaranteed = TRUE,
                terms = c('s(elevation_m)',
                          'ti(elevation_m,temperature_C)')),
      too_far_1 = too_far(temperature_C, x,
                          filter(m$model, detected == 1)$temperature_C,
                          filter(m$model, detected == 1)$elevation_m,
                          dist = dist_1),
      # too_far_2 = too_far(temperature_C, x,
      #                     filter(m$model, detected == 1)$temperature_C,
      #                     filter(m$model, detected == 1)$elevation_m,
      #                     dist = dist_2),
      too_far_2 =
        x < min(filter(m$model, detected == 1)$elevation_m) - 500 |
        x > max(filter(m$model, detected == 1)$elevation_m) + 500,
      x = x / 1e3),
    transmute(
      newd,
      temperature_C,
      x = dist_water_m,
      variable = 'bold(Distance~from~water~(km))',
      lambda =
        predict(object = m, newdata = newd, type = 'response',
                se.fit = FALSE, discrete = FALSE, newdata.guaranteed = TRUE,
                terms = c('s(dist_water_m)',
                          'ti(dist_water_m,temperature_C)')),
      too_far_1 = too_far(temperature_C, x,
                          filter(m$model, detected == 1)$temperature_C,
                          filter(m$model, detected == 1)$dist_water_m,
                          dist = dist_1),
      # too_far_2 = too_far(temperature_C, x,
      #                     filter(m$model, detected == 1)$temperature_C,
      #                     filter(m$model, detected == 1)$dist_water_m,
      #                     dist = dist_2),
      too_far_2 =
        x < min(filter(m$model, detected == 1)$dist_water_m) - 1e3 |
        x > max(filter(m$model, detected == 1)$dist_water_m) + 1e3,
    x = x / 1e3))
}

# predict partial effects
preds <- d %>%
  transmute(species,
            lab,
            lambdas = map(rsf, surface)) %>%
  unnest(lambdas)

# need to divide by median lambda within each variable
preds %>%
  filter(! too_far_1) %>%
  group_by(species, variable) %>%
  summarize(min_lambda = min(lambda),
            median_lambda = median(lambda),
            max_lambda = max(lambda))

LIM <- 2

preds %>%
  filter(! (species == 'Canis lupus' & variable == 'bold(Elevation~(km))' &
              x > 1.2)) %>%
  # filter(species == 'Canis lupus', variable == 'bold(Elevation~(km))') %>%
  select(species, lab, x, temperature_C, variable, lambda, too_far_1,
         too_far_2) %>%
  filter((! too_far_2)) %>%
  # re-scale the center (1) relative to the median in the observed data
  # doing this because some estimated effects of distance from water and
  # elevation are extreme
  group_by(species, variable) %>%
  mutate(lambda = lambda / median(lambda[! too_far_1])) %>%
  ungroup() %>%
  # cap at 2^(+/-2)
  mutate(
    log2_lambda = log2(lambda),
    log2_lambda = case_when(log2_lambda > LIM ~ LIM,
                            log2_lambda < -LIM ~ -LIM,
                            TRUE ~ log2_lambda),
    variable = factor(variable,
                      levels = c("bold(Forest~cover~('%'))",
                                 "bold(Elevation~(km))",
                                 "bold(Distance~from~water~(km))"))) %>%
  ggplot() +
  facet_grid(variable ~ lab, scales = 'free', labeller = label_parsed,
             switch = 'y') +
  geom_raster(aes(temperature_C, x, fill = log2_lambda)) +
  geom_contour(aes(temperature_C, x, z = as.numeric(too_far_1)),
               color = 'grey50', linewidth = 0.25) +
  scale_x_continuous(paste0('Temperature (\U00B0', 'C)'), expand = c(0, 0),
                     breaks = c(-20, 0, 20)) +
  scale_y_continuous(NULL, expand = c(0, 0)) +
  scale_fill_sunset(name = 'Relative selection strength', midpoint = 0,
                    limits = c(-LIM, LIM), breaks = -LIM:LIM,
                    labels = \(x) 2^x) +
  theme(strip.placement = 'outside', strip.background.y = element_blank(),
        strip.text.y = element_text(size = 11), legend.position = 'top',
        panel.background = element_rect(fill = 'grey90'),
        legend.key.width = rel(2))

ggsave('figures/rsf-surface-plots.png', width = 17.5, height = 8,
       units = 'in', dpi = 600, bg = 'white')
