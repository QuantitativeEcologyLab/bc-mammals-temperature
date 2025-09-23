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
         lambda = if_else(
           lab == 'Caribou (boreal)' & variable == 'bold(Elevation~"(km)")',
           lambda / exp(2), lambda)) %>%
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
  scale_fill_PRGn(name = 'Relative selection strength', midpoint = 0,
                  limits = c(-LIM, LIM), breaks = -LIM:LIM,
                  labels = \(x) 2^x) +
  theme(strip.placement = 'outside', strip.background.y = element_blank(),
        strip.text.y = element_text(size = 11), legend.position = 'top',
        panel.background = element_rect(fill = 'grey90'),
        legend.key.width = rel(2))

ggsave('figures/2025-tws-edmonton/boreal-caribou-hrsf-surface-plots.png',
       p, width = 8, height = 3.5, units = 'in', dpi = 600, bg = 'white')
