library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for working with dates
library('sf')        # for working with spatial data
library('terra')     # for working with rasters
library('mgcv')      # for generalized additive models
library('gratia')    # for useful functions for generalized additive models
library('ggplot2')   # for fancy figures
library('khroma')    # for colorblind-friendly palettes
library('stringr')   # for working with strings
source('analysis/figures/default-ggplot-theme.R') # for consistent theme
source('data/bc-shapefile.R') # import shapefile of bc

colorRampPalette(RColorBrewer::brewer.pal(11, 'PuOr'))(1e3) %>%
  plot_scheme_colorblind()

# import climate data ----
if(file.exists('data/cc-hgam-bc-projections-albers.rds')) {
  bc_preds <- readRDS('data/cc-hgam-bc-projections-albers.rds')
} else {
  # import models
  m_1 <- readRDS('models/binomial-gam.rds')
  m_2 <- readRDS('models/gamma-gam.rds')
  
  # approximate Gaussian distributions of weather
  weather <- readRDS('data/weather-projections-2025-2100-only.rds') %>%
    filter(! is.na(temp_c)) %>% # there's some NA values even after masking
    select(scenario, year, month, lat, long, temp_c, weight)
  
  # for creating a raster of predictions
  e <- rast('data/resource-rasters/bc-buffered-dem-z3.tif')
  
  # predict movement rates
  bc_preds <-
    # add columns of species and species labels
    tibble(species = SPECIES,
           lab = SPECIES_LABS,
           w = weather %>%
             filter(! is.na(temp_c)) %>% # some NA values after masking
             select(scenario, year, month, lat, long, temp_c, weight) %>%
             list()) %>%
    unnest(w) %>%
    #' ******\/ for testing *****************************************
    # filter(species == species[1], lat > 58) %>%
    #' ******/\ for testing *****************************************
    # make sure only the necessary column are kept
    transmute(scenario, year, lat, long,
              animal = m_1$model$animal[1],
              species = gsub('southern mountain', '(s. mountain)', x = species),
              lab,
              tod_pdt = 12,
              doy = yday(date_decimal(year + (month - 0.5) / 12)),
              temp_c,
              dt = 1,
              weight) %>%
    # predict P(moving), E(speed | moving) and E(speed) = displacement
    mutate(
      .,
      p = predict(m_1, newdata = ., type = 'response',se.fit = FALSE,
                  discrete = TRUE, exclude = 's(animal)'),
      s = predict(m_2, newdata = ., type = 'response', se.fit = FALSE,
                  discrete = TRUE, exclude = 's(animal)'),
      d = p * s) %>%
    # take the weighted mean of the approximated gaussian distributions 
    # after grouping all the 2025 scenarios together
    mutate(scenario = if_else(year == 2025, '2025', scenario)) %>%
    group_by(scenario, lab, long, lat) %>%
    summarize(p = weighted.mean(p, w = weight),
              s = weighted.mean(s, w = weight),
              d = weighted.mean(d, w = weight),
              .groups = 'drop') %>%
    # calculate relative change from 2025 to 2100, for each location
    arrange(lab, scenario, long, lat) %>%
    group_by(lab, long, lat) %>%
    mutate(p = p / mean(p[scenario == '2025']),
           s = s / mean(s[scenario == '2025']),
           d = d / mean(d[scenario == '2025'])) %>%
    ungroup() %>%
    filter(scenario != '2025') %>%
    mutate(scenario = case_when(grepl('126', scenario) ~ 'SSP~1-2.6',
                                grepl('245', scenario) ~ 'SSP~2-4.5',
                                grepl('370', scenario) ~ 'SSP~3-7.0',
                                grepl('585', scenario) ~ 'SSP~5-8.5') %>%
             paste0('bold(', ., ')')) %>%
    pivot_wider(values_from = c(p, s, d), names_from = c(scenario, lab)) %>%
    pivot_longer(! c(long, lat),
                 names_to = c('parameter', 'scenario', 'lab'),
                 names_sep = '_', values_to = 'value') %>%
    pivot_wider(values_from = value, names_from = parameter) %>%
    # reproject the rasters from lat-long to BC Albers
    nest(r_p = c(long, lat, p),
         r_s = c(long, lat, s),
         r_d = c(long, lat, d)) %>%
    mutate(across(c(r_p, r_s, r_d), \(a) map(a, \(.a) {
      .a %>%
        select(long, lat) %>%
        as.matrix() %>%
        rasterize(y = e, values = .a[[3]], fun = mean) %>%
        project('EPSG:3005') %>%
        crop(bc, mask = TRUE, touches = FALSE) %>%
        as.data.frame(xy = TRUE) %>%
        rename(z = 3)
    })),
    # rename back to correct columns (p, s, d) and drop repeated columns
    r_p = map(r_p, \(a) rename(a, p = z)),
    r_s = map(r_s, \(a) transmute(a, s = z)),
    r_d = map(r_d, \(a) transmute(a, d = z))) %>%
    unnest(c(r_p, r_s, r_d)) %>%
    #' prevent `label_parsed()` from removing the zero in "3-7.0"
    # improve scenario labels
    mutate(scenario = case_when(grepl('1-2.6', scenario) ~ 'bold("Best scenario (SSP 1-2.6)")',
                                grepl('2-4.5', scenario) ~ 'bold("Good scenario (SSP 2-4.5)")',
                                grepl('3-7.0', scenario) ~ 'bold("Bad scenario (SSP 3-7.0)")',
                                grepl('5-8.5', scenario) ~ 'bold("Worst scenario (SSP 5-8.5)")') %>%
             factor(levels = c('bold("Best scenario (SSP 1-2.6)")',
                               'bold("Good scenario (SSP 2-4.5)")',
                               'bold("Bad scenario (SSP 3-7.0)")',
                               'bold("Worst scenario (SSP 5-8.5)")')),
           lab = gsub('\\(s.~mountain\\)', '"\\(s. mountain\\)"', lab))
  
  gc()
  bc_preds
  
  # test the figure
  if(FALSE) {
    bc_preds %>%
      filter(lab == lab[1]) %>%
      ggplot(aes(x, y, fill = d)) +
      facet_grid(lab ~ scenario, labeller = label_parsed) +
      geom_sf(data = bc, inherit.aes = FALSE) +
      geom_raster() +
      scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = NULL) +
      labs(x = NULL, y = NULL)
  }
  
  saveRDS(bc_preds, 'data/cc-hgam-bc-projections-albers.rds')
}

# check changes with density functions
bc_preds %>%
  mutate(change = d * 100 - 100,
         change = if_else(change > 10 & grepl('boreal', lab), 10, change)) %>%
  ggplot() +
  facet_wrap(~ lab, scales = 'free', labeller = label_parsed) +
  geom_vline(xintercept = 0, color = 'grey', linetype = 'dashed') +
  geom_density(aes(change, fill = scenario, color = scenario), alpha = 0.3) +
  scale_color_manual('Scenario', values = khroma::color('sunset')(4),
                     aesthetics = c('color', 'fill'),
                     labels = scales::parse_format()) +
  labs(x = 'Change in distance travelled relative to 2025 (%)',
       y = 'Probability density') +
  theme(legend.position = 'inside', legend.position.inside = c(0.85, 0.15))

# figure of estimated speeds for each species ----
z_breaks <- seq(-log2(1.5), log2(1.5), length.out = 5)

# check color scheme
if(FALSE) {
  p_0 <- filter(bc_preds, lab == lab[1], scenario == scenario[1]) %>%
    # cap values for readability
    mutate(d = case_when(d < 2^min(z_breaks) ~ 2^min(z_breaks),
                         d > 2^max(z_breaks) ~ 2^max(z_breaks),
                         TRUE ~ d)) %>%
    ggplot() +
    geom_raster(aes(x, y, fill = log2(d))) +
    geom_sf(data = bc, fill = 'transparent') +
    scale_fill_distiller(name = 'Relative change in distance travelled',
                         palette = 'PuOr', direction = 1,
                         labels = \(x) round(2^x, 2)) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = 'top'); p_0
  
  colorblindr::cvd_grid(p_0)
}

# make figure
p <- bc_preds %>%
  # cap values for readability
  mutate(d = case_when(d < 2^min(z_breaks) ~ 2^min(z_breaks),
                       d > 2^max(z_breaks) ~ 2^max(z_breaks),
                       TRUE ~ d)) %>%
  ggplot() +
  facet_grid(scenario ~ lab, labeller = label_parsed) +
  geom_raster(aes(x, y, fill = log2(d))) +
  geom_sf(data = bc, fill = 'transparent') +
  scale_fill_distiller(name = 'Relative change in distance travelled',
                       palette = 'PuOr', limits = range(z_breaks),
                       breaks = z_breaks, labels = \(x) round(2^x, 2),
                       direction = 1) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'top', legend.key.width = rel(3))

ggsave('figures/bc-distance-2100.png', p, width = 15, height = 8,
       units = 'in', dpi = 600, bg = 'white', scale = 1.2)

# make zoomed in version
tels <- readRDS()

p <- bc_preds %>%
  filter(case_when(species == '')) %>%
  # cap values for readability
  # mutate(d = case_when(d < 2^min(z_breaks) ~ 2^min(z_breaks),
  #                      d > 2^max(z_breaks) ~ 2^max(z_breaks),
  #                      TRUE ~ d)) %>%
  ggplot() +
  facet_grid(scenario ~ lab, labeller = label_parsed) +
  geom_raster(aes(x, y, fill = log2(d))) +
  geom_sf(data = bc, fill = 'transparent') +
  scale_fill_distiller(name = 'Relative change in distance travelled',
                       palette = 'PuOr', limits = range(z_breaks),
                       breaks = z_breaks, labels = \(x) round(2^x, 2),
                       direction = 1) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'top', legend.key.width = rel(3))

ggsave('figures/ranges-dist-2100.png', p, width = 15, height = 8,
       units = 'in', dpi = 600, bg = 'white', scale = 1.2)
