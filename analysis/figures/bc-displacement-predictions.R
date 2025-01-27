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
library('rphylopic') # for animal silhouettes to plots
library('cowplot')   # to add phylopic once (not in each facet)
source('analysis/figures/default-ggplot-theme.R') # for consistent theme
source('data/bc-shapefile.R') # import shapefile of bc

d_pal <- colorRampPalette(c('#9A2600', '#EBE8DB', '#00605C'))(1e3)
plot_scheme_colorblind(d_pal)

# import climate data ----
if(file.exists('data/cc-hgam-bc-projections-albers.rds')) {
  bc_preds <- readRDS('data/cc-hgam-bc-projections-albers.rds')
} else {
  # import models
  m_1 <- readRDS('models/binomial-gam.rds')
  m_2 <- readRDS('models/gamma-gam.rds')
  
  # terms to exclude from the prediction
  SM <- unique(c(smooths(m_1), smooths(m_2)))
  EXCLUDE <- SM[grepl('tod', SM) | grepl('animal', SM)]
  
  # predict movement rates
  bc_preds <-
    readRDS('data/weather-projections-2025-2100-only.rds') %>%
    filter(! is.na(temp_c)) %>% # there's some NA values even after masking
    select(scenario, year, month, lat, long, temp_c, weight) %>%
    #' ******\/ for testing *****************************************
    # group_by(year, scenario, long) %>%
    # slice(1:3) %>%
    # ungroup() %>%
    # mutate(lat = list(tibble(lat = c(55, 56)))) %>%
    # unnest(lat) %>%
    #' ******/\ for testing *****************************************
    # add columns of species and species labels
    mutate(spp = list(tibble(species = SPECIES)),
           labs = list(tibble(lab = SPECIES_LABS))) %>%
    unnest(c(spp, labs)) %>%
    # make sure only the necessary column are kept
    transmute(scenario, year, lat, long,
              animal = m_1$model$animal[1],
              species = gsub('southern mountain', '(s. mountain)', x = species),
              lab,
              tod_pdt = 0,
              doy = yday(date_decimal(year + (month - 0.5) / 12)),
              temp_c,
              dt = 1,
              weight) %>%
    # predict P(moving), E(speed | moving) and E(speed) = displacement
    mutate(
      .,
      p = predict(m_1, newdata = ., type = 'response',se.fit = FALSE,
                  discrete = TRUE, exclude = EXCLUDE),
      s = predict(m_2, newdata = ., type = 'response', se.fit = FALSE,
                  discrete = TRUE, exclude = EXCLUDE),
      d = p * s) %>%
    # take the weighted mean of the approximated gaussian distributions 
    # after grouping all the 2025 scenarios together
    mutate(scenario = if_else(year == 2025, '2025', scenario)) %>%
    group_by(scenario, lab, long, lat) %>%
    summarize(p = weighted.mean(p, w = weight),
              s = weighted.mean(s, w = weight),
              d = weighted.mean(d, w = weight),
              .groups = 'drop') %>%
    # calculate change over space relative to the mean 2025 values 
    arrange(lab, scenario, long, lat) %>%
    group_by(lab, long, lat) %>%
    mutate(p = p / mean(p[1:4]),
           s = s / mean(s[1:4]),
           d = d / mean(d[1:4])) %>%
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
      # interpolate between points to fill the raster
      interp::interp(x = .a$long, y = .a$lat, z = .a[[3]],
                     nx = n_distinct(.a$long), ny = n_distinct(.a$lat)) %>%
        interp::interp2xyz() %>%
        as.data.frame() %>%
        setNames(c('x', 'y', 'z')) %>%
        rast(type = 'xyz', crs = 'EPSG:4326') %>%
        project('EPSG:3005') %>%
        crop(bc) %>%
        mask(bc) %>%
        as.data.frame(xy = TRUE)
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
           # fix species labs rather than having to re-run all predictions
           lab = gsub('\\(boreal\\)', '"\\(boreal\\)"', lab),
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

# check overall changes from 2025 to 2100
bc_preds %>%
  group_by(lab, scenario) %>%
  summarise(mean_change = paste0(round(mean(d) * 100 - 100, 1), '%'),
            sd_change = paste0(round(sd(d) * 100, 1), '%')) %>%
  print(n = length(SPECIES) * 4)

# check changes with density functions
bc_preds %>%
  mutate(change = d * 100 - 100) %>%
  ggplot() +
  facet_wrap(~ lab, scales = 'free_y', labeller = label_parsed) +
  geom_vline(xintercept = 0, color = 'grey', linetype = 'dashed') +
  geom_density(aes(change, fill = scenario, color = scenario), alpha = 0.3) +
  scale_color_manual('Scenario', values = khroma::color('sunset')(4),
                     aesthetics = c('color', 'fill')) +
  labs(x = 'Change in distance travelled relative to 2025 (%)',
       y = 'Probability density') +
  xlim(c(-10, 20)) +
  theme(legend.position = 'inside', legend.position.inside = c(0.85, 0.15))

# check color scheme
p_0 <- filter(bc_preds, lab == lab[1], scenario == scenario[1]) %>%
  ggplot() +
  geom_raster(aes(x, y, fill = log2(d))) +
  geom_sf(data = bc, fill = 'transparent') +
  scale_fill_gradientn(name = 'Relative change in distance travelled',
                       colors = d_pal) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'top')

colorblindr::cvd_grid(p_0)

# figure of estimated speeds for each species ----
z_breaks <- seq(-log2(1.2), log2(1.2), length.out = 5)

p <- bc_preds %>%
  # cap values at +20% for readability
  mutate(d = if_else(d > 1.2, 1.2, d)) %>%
  ggplot() +
  facet_grid(scenario ~ lab, labeller = label_parsed) +
  geom_raster(aes(x, y, fill = log2(d))) +
  geom_sf(data = bc, fill = 'transparent') +
  scale_fill_gradientn(name = 'Relative change in distance travelled',
                       colors = d_pal, limits = range(z_breaks),
                       breaks = z_breaks, labels = \(x) round(2^x, 2)) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'top', legend.key.width = rel(3))

ggsave('figures/bc-displ-2100.png', p, width = 15, height = 8,
       units = 'in', dpi = 600, bg = 'white', scale = 1.2)
