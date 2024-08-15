library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for smother date wrangling
library('mgcv')      # for Generalized Additive Models
library('gratia')    # for useful convenience functions for GAMs
library('ggplot2')   # for fancy plots
library('khroma')    # for colorblind-friendly color palettes
library('cowplot')   # for fancy multi-panel plots
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids
plot_scheme(PAL, colours = TRUE)

# import models for each species
d <- tibble(
  species = SPECIES[c(1:2, 7)],
  lab = SPECIES_LABS[c(1:2, 7)],
  file_name = map_chr(species, \(.sp) {
    list.files(path = 'H:/GitHub/bc-mammals-temperature/models/',
               pattern = gsub('\\(s\\.', 'southern', .sp) %>%
                 gsub('\\(', '', .) %>%
                 gsub('\\)', '', .) %>%
                 paste0('rsf-', .),
               full.names = TRUE)
  }),
  rsf = map(file_name, readRDS),
  local_newd = map(species, tibble))

d

# get range of predictors
d0 <- readRDS('data/tracking-data/rsf-data.rds')
range(d0$elevation_m)
range(d0$dist_water_m)
range(d0$temperature_C)

# surface plots of partial effects ----
newd <-
  tibble(animal = 'new animal',
         temperature_C = seq(-40, 40, length.out = 400),
         x = tibble(forest_perc = seq(0, 100, length.out = 400),
                    elevation_m = seq(0, 2500, length.out = 400),
                    dist_water_m = seq(0, 5e3, length.out = 400)) %>%
           list()) %>%
  unnest(x)

surface <- function(m, dist = 0.1) {
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
      too_far = too_far(temperature_C, x,
                        m$model$temperature_C, m$model$forest_perc,
                        dist = dist)),
    transmute(
      newd,
      temperature_C,
      x = elevation_m,
      variable = 'bold(Elevation~(m))',
      lambda =
        predict(object = m, newdata = newd, type = 'response',
                se.fit = FALSE, discrete = FALSE, newdata.guaranteed = TRUE,
                terms = c('s(elevation_m)',
                          'ti(elevation_m,temperature_C)'),
      too_far = too_far(temperature_C, x,
                        m$model$temperature_C, m$model$elevation_m,
                        dist = 1))),
    transmute(
      newd,
      temperature_C,
      x = dist_water_m,
      variable = 'bold(Distance~from~water~(m))',
      lambda =
        predict(object = m, newdata = newd, type = 'response',
                se.fit = FALSE, discrete = FALSE, newdata.guaranteed = TRUE,
                terms = c('s(dist_water_m)',
                          'ti(dist_water_m,temperature_C)')),
    too_far = too_far(temperature_C, x,
                      m$model$temperature_C, m$model$dist_water_m,
                      dist = dist)))
}

# predict partial effects
preds <- d %>%
  transmute(species,
            lab,
            lambdas = map(rsf, surface)) %>%
  unnest(lambdas)

preds %>%
  select(species, lab, x, temperature_C, variable, lambda, too_far) %>%
  group_by(species, variable) %>%
  mutate(lambda = lambda / median(lambda)) %>%
  ungroup() %>%
  #' **elevation gets dropped if only using too_far because it's NA**
  filter((! too_far) | variable == 'bold(Elevation~(m))') %>%
  mutate(lambda = case_when(lambda > 4 ~ 4,
                            lambda < 0.25 ~ 0.25,
                            TRUE ~ lambda)) %>%
  ggplot() +
  facet_grid(variable ~ lab, scales = 'free', labeller = label_parsed,
             switch = 'y') +
  geom_raster(aes(temperature_C, x, fill = log2(lambda))) +
  scale_x_continuous(paste0('Temperature (\U00B0', 'C)'), expand = c(0, 0)) +
  scale_y_continuous(NULL, expand = c(0, 0)) +
  scale_fill_sunset(name = 'Relative selection strength', midpoint = 0,
                    limits = c(-2, 2), breaks = -2:2, labels = 2^c(-2:2)) +
  theme(strip.placement = 'outside', strip.background.y = element_blank(),
        strip.text.y = element_text(size = 11), legend.position = 'top')

# old code -----

# predict habitat quality
hq <- mutate(hq,
             preds = map2(model, species,
                          \(.m, .sp) {
                            mutate(
                              filter(newd, species == .sp),
                              mu = predict(
                                .m, newdata = filter(newd, species == .sp),
                                type = 'response') %>%
                                # take yearly averages
                                group_by(scenario, year) %>%
                                summarize(mu = mean(mu), .groups = 'drop'))
                          }))

hq_rel <-
  hq %>%
  select(-model) %>%
  unnest(preds) %>%
  # scale to ralative change
  group_by(species, scenario) %>%
  mutate(ref = mu[2],
         mu_rel = mu / ref) %>%
  ungroup() %>%
  mutate(mu_rel = if_else(mu_rel < 0.25, 0.25, mu_rel),
         mu_rel = if_else(mu_rel > 4, 4, mu_rel))

# plot the estimated change
ggplot(hq_rel, aes(year, mu_rel, color = scenario)) +
  facet_wrap(~ species, scales = 'fixed') +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_point() +
  geom_smooth(se = FALSE, method = 'gam', formula = y ~ s(x)) +
  xlab(NULL) +
  scale_y_continuous(expression(Relative~change~'in'~habitat~quality~(log[2])),
                     trans = 'log2') +
  scale_color_brewer('Scenario', type = 'div', palette = 5, direction = -1,
                     aesthetics = c('color', 'fill')) +
  theme(legend.position = c(0.85, 0.2))

ggsave('figures/climate-change-habitat-quality-relative.png',
       width = 10, height = 5, dpi = 600, bg = 'white')

