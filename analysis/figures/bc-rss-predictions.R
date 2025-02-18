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

colorRampPalette(RColorBrewer::brewer.pal(11, 'PiYG'))(1e3) %>%
  plot_scheme_colorblind()

# import climate data ----
if(file.exists('data/cc-rss-bc-projections-albers.rds')) {
  bc_preds <- readRDS('data/cc-rss-bc-projections-albers.rds')
} else {
  # get file names
  fn <- list.files(path = 'models', pattern = '^rsf-', recursive = FALSE,
                   full.names = TRUE)
  fn <- fn[! grepl('no-temperature', fn)]
  
  # terms to exclude from the prediction
  SM <- smooths(readRDS('models/rsf-Oreamnos americanus-2025-01-20.rds'))
  EXCLUDE <- SM[grepl('animal', SM)]
  
  # import resource rasters
  f <- rast('data/resource-rasters/forest.tif')
  e <- rast('data/resource-rasters/bc-buffered-dem-z3.tif')
  w <- rast('data/resource-rasters/distance-from-water.tif')
  
  layout(t(1:3))
  plot(f)
  plot(e)
  plot(w)
  layout(1)
  
  weather <- readRDS('data/weather-projections-2025-2100-only.rds') %>%
    filter(! is.na(temp_c)) %>% # there's some NA values even after masking
    select(scenario, year, month, lat, long, temp_c, weight)
  
  # predict relative selection strength
  bc_preds <-
    weather %>%
    #' ******\/ for testing *****************************************
    # filter(lat > 58) %>%
    #' ******/\ for testing *****************************************
    # extract resources for each location
    nest(temporal = ! c(long, lat)) %>% # to reduce number of extractions
    mutate(forest_perc = extract(f, tibble(long, lat))[[2]], # 1 is ID
           elevation_m = extract(e, tibble(long, lat))[[2]],
           dist_water_m = extract(w, tibble(long, lat))[[2]],
           animal = 'new animal') %>%
    filter(! is.na(forest_perc)) %>%
    unnest(temporal) %>%
    # rename to columns used in the model
    rename(temperature_C = temp_c) %>%
    # nest data for predictions for each species
    nest(newd = everything()) %>%
    # add columns of species and species labels
    mutate(spp = list(tibble(species = SPECIES)),
           labs = list(tibble(lab = SPECIES_LABS))) %>%
    unnest(c(spp, labs)) %>%
    #' ******\/ for testing *****************************************
    # filter(species == 'Oreamnos americanus') %>% # has the least data
    #' ******/\ for testing *****************************************
    # import the HRSFs for each species
    mutate(
      # predict from each HRSF
      newd = map2(species, newd, \(.sp, .d) {
        # import HRSF
        fn <- case_when(
          .sp == 'Canis lupus' ~
            'models/rsf-Canis lupus-2025-01-20.rds',
          .sp == 'Cervus canadensis' ~
            'models/rsf-Cervus canadensis-2025-01-20.rds',
          .sp == 'Oreamnos americanus' ~
            'models/rsf-Oreamnos americanus-2025-01-20.rds',
          .sp == 'Puma concolor' ~
            'models/rsf-Puma concolor-2025-01-20.rds',
          .sp == 'Rangifer tarandus (boreal)' ~
            'models/rsf-Rangifer tarandus boreal-2025-01-21.rds',
          .sp == 'Rangifer tarandus (s. mountain)' ~
            'models/rsf-Rangifer tarandus southern mountain-2025-01-20.rds',
          .sp == 'Ursus arctos horribilis' ~
            'models/rsf-Ursus arctos horribilis-2025-01-20.rds')
        
        m <- readRDS(fn)
        
        # predict RSS values
        mutate(.d,
               lambda = predict(object = m, newdata = .d,
                                type = 'response', se.fit = FALSE,
                                exclude = EXCLUDE, discrete = FALSE))
      })) %>%
    unnest(newd) %>%
    # make sure only the necessary column are kept
    select(scenario, year, lat, long, lab, weight, lambda) %>%
    # take the weighted mean of the approximated gaussian distributions
    # after grouping all the 2025 scenarios together
    mutate(scenario = if_else(year == 2025, '2025', scenario)) %>%
    group_by(scenario, lab, long, lat) %>%
    summarize(lambda = weighted.mean(lambda, w = weight),
              .groups = 'drop') %>%
    # calculate relative change from 2025 to 2100, for each location
    group_by(lab, long, lat) %>%
    mutate(lambda = lambda / median(lambda[scenario == '2025'],
                                    na.rm = TRUE)) %>%
    ungroup() %>%
    filter(scenario != '2025') %>% # drop reference
    mutate(scenario = case_when(grepl('126', scenario) ~ 'SSP~1-2.6',
                                grepl('245', scenario) ~ 'SSP~2-4.5',
                                grepl('370', scenario) ~ 'SSP~3-7.0',
                                grepl('585', scenario) ~ 'SSP~5-8.5') %>%
             paste0('bold(', ., ')')) %>%
    # reproject the rasters from lat-long to BC Albers
    nest(r_lambda = c(long, lat, lambda)) %>%
    mutate(r_lambda = map(r_lambda, \(.l) {
      # interpolate between points to fill the raster
      .l %>%
        select(long, lat) %>%
        as.matrix() %>%
        rasterize(y = e, values = .l$lambda, fun = mean) %>%
        project('EPSG:3005') %>%
        crop(bc, mask = TRUE, touches = FALSE) %>%
        as.data.frame(xy = TRUE) %>%
        rename(lambda = mean)
    })) %>%
    unnest(r_lambda) %>%
    #' to prevent `label_parsed()` from removing the zero in "3-7.0"
    # improve scenario labels
    mutate(scenario = case_when(
      grepl('1-2.6', scenario) ~ 'bold("Best scenario (SSP 1-2.6)")',
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
      ggplot(aes(x, y, fill = lambda)) +
      facet_grid(lab ~ scenario, labeller = label_parsed) +
      geom_sf(data = bc, inherit.aes = FALSE) +
      geom_raster() +
      scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = NULL) +
      labs(x = NULL, y = NULL)
  }
  
  saveRDS(bc_preds, 'data/cc-rss-bc-projections-albers.rds')
}

# check changes with density functions
bc_preds %>%
  mutate(change = lambda,
         change = case_when(change > 1.5 ~ 1.5,
                            TRUE ~ change)) %>%
  ggplot() +
  facet_wrap(~ lab, scales = 'free', labeller = label_parsed) +
  geom_vline(xintercept = 1, color = 'grey', linetype = 'dashed') +
  geom_density(aes(change, fill = scenario, color = scenario), alpha = 0.3) +
  scale_color_manual('Scenario', values = khroma::color('sunset')(4),
                     aesthetics = c('color', 'fill'), labels = scales::parse_format()) +
  labs(x = 'Relative change in RSS', y = 'Probability density') +
  theme(legend.position = 'inside', legend.position.inside = c(0.85, 0.15))

ggsave('figures/bc-rss-2100-density.png', width = 15, height = 8,
       units = 'in', dpi = 600, bg = 'white')

bc_preds %>%
  group_by(lab, scenario) %>%
  summarize(p_decrease = mean(lambda < 1)) %>%
  filter(p_decrease < 0.5)

bc_preds %>%
  group_by(lab) %>%
  summarise(q_0.05 = round(quantile(lambda, 0.05) * 100) - 100,
            q_0.95 = round(quantile(lambda, 0.95) * 100) - 100)

# check color scheme
p_0 <- filter(bc_preds, lab == lab[1], scenario == scenario[1]) %>%
  mutate(log2_l = log2(lambda),
         log2_l = if_else(abs(log2_l) > 2, 2 * sign(log2_l), log2_l)) %>%
  ggplot() +
  geom_raster(aes(x, y, fill = log2_l)) +
  geom_sf(data = bc, fill = 'transparent') +
  scale_fill_distiller(name = 'Relative selection strength',
                       type = 'div', palette = 2, direction = 1) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'top'); p_0

colorblindr::cvd_grid(p_0)

# figure of estimated speeds for each species ----
z_breaks <- seq(-log2(1.5), log2(1.5), length.out = 5)

p <- bc_preds %>%
  # cap log2(values) for readability
  mutate(log2_l = log2(lambda),
         log2_l = if_else(abs(log2_l) > max(z_breaks),
                          max(z_breaks) * sign(log2_l), log2_l)) %>%
  ggplot() +
  facet_grid(scenario ~ lab, labeller = label_parsed) +
  geom_raster(aes(x, y, fill = log2_l)) +
  geom_sf(data = bc, fill = 'transparent') +
  scale_fill_distiller(name = 'Relative change in selection strength',
                       type = 'div', palette = 2, direction = 1,
                       labels = \(x) round(2^x, 2),  limits = range(z_breaks),
                       breaks = z_breaks) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'top', legend.key.width = rel(3))

ggsave('figures/bc-rss-2100.png', p, width = 15, height = 8,
       units = 'in', dpi = 600, bg = 'white', scale = 1.2)
