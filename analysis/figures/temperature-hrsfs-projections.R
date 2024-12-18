library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for smoother date wrangling
library('mgcv')      # for Generalized Additive Models
library('terra')     # for rasters
library('ggplot2')   # for fancy plots
library('khroma')    # for colorblind-friendly color palettes
library('cowplot')   # for fancy multi-panel plots
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids

# import resource rasters ----
f <- rast('data/resource-rasters/forest.tif')
e <- rast('data/resource-rasters/bc-buffered-dem-z6.tif')
w <- rast('data/resource-rasters/distance-from-water.tif')

if(FALSE) {
  layout(t(1:3))
  plot(f)
  plot(e)
  plot(w)
  layout(1)
}

# import models and predict ----
if(file.exists('data/cc-rsf-projections.rds')) {
  cc_proj <- readRDS('data/cc-rsf-projections.rds')
} else {
  # terms to exclude from the predictions
  EXCLUDE <- c('s(animal)', 's(forest_perc,animal)',
               's(elevation_m,animal)', 's(dist_water_m,animal)')
  MODEL_FILES <- list.files('models', '^rsf-', full.names = TRUE)
  MODEL_FILES <- MODEL_FILES[! grepl('no-temp', MODEL_FILES)]
  
  map_chr(MODEL_FILES, function(fn) {
    rsf <- readRDS(fn)
    sp <- tolower(fn) %>%
      gsub('models/rsf-', '', .) %>%
      gsub('-2.*', '', .) %>%
      gsub(' ', '-', .) %>%
      gsub('\\(', '', .) %>%
      gsub('\\)', '', .) %>%
      gsub('\\.', '', .)
    
    rm(x) # remove previous object, if it exists
    
    x <<- # to keep in global environment in case saving fails
      tibble(
        species = fn %>%
          gsub(pattern = 'models/rsf-', replacement = '', x = .) %>%
          gsub(pattern = '(-).*', replacement = '', x = ., perl = FALSE),
        lab = gsub('southern mountain', '(s. mountain)', x = species) %>%
          gsub('boreal', '(boreal)', x = .)) %>%
      mutate(
        lab = SPECIES_LABS[map_int(lab, \(x) which(SPECIES == x))],
        preds = map(
          paste0('data/weather-projections-', species, '.rds'),
          function(p_fn) {
            .nd <-
              readRDS(p_fn) %>%
              filter(year >= 2020) %>%
              # make sure only the necessary column are kept
              transmute(
                scenario,
                year,
                long,
                lat,
                animal = rsf$model$animal[1],
                temperature_C = temp_c,
                weight, # for weather quantiles
                forest_perc = terra::extract(f, tibble(long, lat))[, 2],
                elevation_m = terra::extract(e, tibble(long, lat))[, 2],
                dist_water_m = terra::extract(w, tibble(long, lat))[, 2]) %>%
              ungroup() %>%
              select(! species)
            
            #' use an `animal` from the model to use `discrete = TRUE`
            .nd %>%
              mutate(l = predict.bam(rsf, newdata = .nd, type = 'response',
                                     se.fit = FALSE, exclude = EXCLUDE,
                                     discrete = TRUE)) %>%
              return()
          })) %>%
      select(lab, preds) %>%
      unnest(preds) %>%
      # calculate habitat suitability for each pixel in each year 
      group_by(lab, scenario, year, long, lat) %>%
      summarize(l = weighted.mean(l, w = weight),
                .groups = 'drop') %>%
      # find median and 90% percentile interval of the predicted means
      # not including CIs because averaging them is not straightforward
      group_by(lab, scenario, year) %>%
      summarize(l_lwr_05 = quantile(l, 0.05),
                l_median = quantile(l, 0.50),
                l_upr_95 = quantile(l, 0.95),
                .groups = 'drop') %>%
      # divide by mean of 2025 to find relative change since 2025
      mutate(l_ref = mean(l_median[year == 2025])) %>%
      mutate(scenario = case_when(grepl('126', scenario) ~ 'Best scenario (SSP 1-2.6)',
                                  grepl('245', scenario) ~ 'Good scenario (SSP 2-4.5)',
                                  grepl('370', scenario) ~ 'Bad scenario (SSP 3-7.0)',
                                  grepl('585', scenario) ~ 'Worst scenario (SSP 5-8.5)') %>%
               factor(., levels = unique(.)))
    saveRDS(x, paste0('data/cc-rsf-projections-', sp, '.rds'))
    
    return(sp)
  })
  
  cc_proj <- map_dfr(list.files('data', 'cc-rsf-projections-', full.names = TRUE), readRDS)
  saveRDS(cc_proj, 'data/cc-rsf-projections.rds')
}

# figures of habitat quality ----
cc_proj %>%
  ggplot() +
  facet_wrap(~ lab, scales = 'free_y', labeller = label_parsed) +
  geom_ribbon(aes(year, ymin = l_lwr_05 / l_ref, ymax = l_upr_95 / l_ref,
                  fill = scenario), alpha = 0.2) +
  geom_line(aes(year, l_upr_95 / l_ref, color = scenario), lwd = 0.5) +
  geom_line(aes(year, l_lwr_05 / l_ref, color = scenario), lwd = 0.5) +
  geom_hline(yintercept = 1, color = 'black', lty = 'dashed') +
  geom_line(aes(year, l_median / l_ref, color = scenario), lwd = 1) +
  xlab(NULL) +
  scale_y_continuous('Change in relative selection strength') +
  scale_color_brewer('Scenario', type = 'div', palette = 5, direction = -1,
                     aesthetics = c('color', 'fill')) +
  theme(legend.position = 'inside', legend.position.inside = c(5/6, 1/6))

ggsave('figures/rss-local-cc-predictions.png',
       width = 10, height = 5, dpi = 600, bg = 'white')

ggsave('figures/2024-ubco-grad-symposium/rss-local-cc-predictions.png',
       width = 17.5, height = 9.5, dpi = 300, bg = 'white', scale = 0.75)
