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

# import climate data ----
if(file.exists('data/cc-hrsf-bc-projections-albers.rds')) {
  bc_preds <- readRDS('data/cc-hrsf-bc-projections-albers.rds')
} else {
  # get file names
  fn <- list.files(path = 'models', pattern = '^rsf-', recursive = FALSE,
                   full.names = TRUE)
  fn <- fn[! grepl('no-temperature', fn)]
  
  # terms to exclude from the prediction
  SM <- smooths(readRDS('models/rsf-Oreamnos americanus-2025-01-20.rds'))
  EXCLUDE <- SM[grepl('animal', SM)]
  
  # import resource rasters
  f <- rast('data/resource-rasters/bc-forest.tif')
  e <- rast('data/resource-rasters/bc-dem-z3.tif')
  w <- rast('data/resource-rasters/bc-distance-from-water.tif')
  
  layout(t(1:3))
  plot(f)
  plot(e)
  plot(w)
  layout(1)
  
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
    # extract resources for each location
    mutate(forest_perc = extract(f, tibble(long, lat))[[2]], # 1 is ID
           elevation_m = extract(e, tibble(long, lat))[[2]],
           dist_water_m = extract(w, tibble(long, lat))[[2]],
           animal = 'new animal') %>%
    filter(! is.na(forest_perc)) %>%
    # rename to columns used in the model
    rename(temperature_C = temp_c) %>%
    # nest data for predictions for each species
    nest(newd = everything()) %>%
    # add columns of species and species labels
    mutate(spp = list(tibble(species = SPECIES)),
           labs = list(tibble(lab = SPECIES_LABS))) %>%
    unnest(c(spp, labs)) %>%
    arrange(species) %>% # to match file names 
    # import the HRSFs for each species
    mutate(file = fn) %>%
    # predict from each HRSF
    mutate(newd = map2(file, newd, \(.fn, .d) {
      mutate(.d,
             lambda = predict(object = readRDS(.fn), newdata = .d,
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
    # calculate change over space relative to the mean 2025 values
    arrange(lab, scenario, long, lat) %>%
    group_by(lab, long, lat) %>%
    mutate(lambda = lambda / mean(lambda[1:4])) %>%
    ungroup() %>%
    filter(scenario != '2025') %>%
    mutate(scenario = case_when(grepl('126', scenario) ~ 'SSP~1-2.6',
                                grepl('245', scenario) ~ 'SSP~2-4.5',
                                grepl('370', scenario) ~ 'SSP~3-7.0',
                                grepl('585', scenario) ~ 'SSP~5-8.5') %>%
             paste0('bold(', ., ')')) %>%
    # drop any NAs
    filter(! is.na(lambda)) %>%
    # reproject the rasters from lat-long to BC Albers
    nest(r_lambda = c(long, lat, lambda)) %>%
    mutate(r_lambda = map(r_lambda, \(.l) {
      # interpolate between points to fill the raster
      interp::interp(x = .l$long, y = .l$lat, z = .l$lambda,
                     nx = n_distinct(.l$long), ny = n_distinct(.l$lat)) %>%
        interp::interp2xyz() %>%
        as.data.frame() %>%
        setNames(c('x', 'y', 'z')) %>%
        rast(type = 'xyz', crs = 'EPSG:4326') %>%
        project('EPSG:3005') %>%
        crop(bc) %>%
        mask(bc) %>%
        as.data.frame(xy = TRUE) %>%
        rename(lambda = z)
    })) %>%
    unnest(r_lambda) %>%
    #' prevent `label_parsed()` from removing the zero in "3-7.0"
    # improve scenario labels
    mutate(scenario = case_when(grepl('1-2.6', scenario) ~ 'bold("Best scenario (SSP 1-2.6)")',
                                grepl('2-4.5', scenario) ~ 'bold("Good scenario (SSP 2-4.5)")',
                                grepl('3-7.0', scenario) ~ 'bold("Bad scenario (SSP 3-7.0)")',
                                grepl('5-8.5', scenario) ~ 'bold("Worst scenario (SSP 5-8.5)")') %>%
             factor(levels = c('bold("Best scenario (SSP 1-2.6)")',
                               'bold("Good scenario (SSP 2-4.5)")',
                               'bold("Bad scenario (SSP 3-7.0)")',
                               'bold("Worst scenario (SSP 5-8.5)")')))
  
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
  
  saveRDS(bc_preds, 'data/cc-hrsf-bc-projections-albers.rds')
}

# check overall changes from 2025 to 2100
bc_preds %>%
  group_by(lab, scenario) %>%
  summarise(mean_change = paste0(round(mean(lambda) * 100 - 100, 1), '%'),
            sd_change = paste0(round(sd(lambda) * 100, 1), '%')) %>%
  print(n = length(SPECIES) * 4)

# check changes with density functions
bc_preds %>%
  mutate(change = lambda * 100 - 100) %>%
  ggplot() +
  facet_wrap(~ lab, scales = 'free_y', labeller = label_parsed) +
  geom_vline(xintercept = 0, color = 'grey', linetype = 'dashed') +
  geom_density(aes(change, fill = scenario, color = scenario), alpha = 0.3) +
  scale_color_manual('Scenario', values = khroma::color('sunset')(4),
                     aesthetics = c('color', 'fill'), labels = scales::parse_format()) +
  labs(x = 'Change in RSS relative to 2025 (%)',
       y = 'Probability density') +
  theme(legend.position = 'inside', legend.position.inside = c(0.85, 0.15))

# check color scheme
p_0 <- filter(bc_preds, lab == lab[1], scenario == scenario[1]) %>%
  ggplot() +
  geom_raster(aes(x, y, fill = log2(lambda))) +
  geom_sf(data = bc, fill = 'transparent') +
  scale_fill_PRGn(name = 'Relative selection strength', midpoint = 0,
                 labels = \(x) round(2^x, 2), limits = c(-1, 1)) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'top')

colorblindr::cvd_grid(p_0)

# figure of estimated speeds for each species ----
z_breaks <- seq(-2, 2, length.out = 5)

p <- bc_preds %>%
  # cap values at +20% for readability
  mutate(log2_l = log2(lambda),
         log2_l = if_else(abs(log2_l) > 2, 2 * sign(log2_l), log2_l)) %>%
  ggplot() +
  facet_grid(scenario ~ lab, labeller = label_parsed) +
  geom_raster(aes(x, y, fill = log2_l)) +
  geom_sf(data = bc, fill = 'transparent') +
  scale_fill_PRGn(name = 'Relative selection strength', midpoint = 0,
                  labels = \(x) round(2^x, 2),  limits = range(z_breaks),
                  breaks = z_breaks) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'top', legend.key.width = rel(3))

ggsave('figures/bc-rss-2100.png', p, width = 15, height = 8,
       units = 'in', dpi = 600, bg = 'white', scale = 1.2)
