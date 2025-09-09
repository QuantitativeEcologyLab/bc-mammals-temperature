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
library('cowplot')   # for plots in grids (for free scales)
source('analysis/figures/default-ggplot-theme.R') # for consistent theme
source('functions/get_legend.R') # to extract legens from ggplots
source('data/bc-shapefile.R') # import shapefile of bc

colorRampPalette(RColorBrewer::brewer.pal(11, 'PiYG'))(1e3) %>%
  plot_scheme_colorblind()

if(file.exists('data/cc-hrsf-projections-local-2100.rds')) {
  cc_proj <- readRDS('data/cc-hrsf-projections-local-2100.rds')
} else {
  # import resource rasters
  f <- rast('data/resource-rasters/forest.tif')
  e <- rast('data/resource-rasters/bc-buffered-dem-z6.tif')
  w <- rast('data/resource-rasters/distance-from-water.tif')
  
  # terms to exclude from the predictions
  SM <- gratia::smooths(readRDS('models/rsf-Oreamnos americanus-2025-01-20.rds'))
  EXCLUDE <- SM[grepl('animal', SM)] # drop all individual-specific smooths
  
  # find model files and drop those with no temperature smooths
  MODEL_FILES <- list.files('models', '^rsf-', full.names = TRUE)
  MODEL_FILES <- MODEL_FILES[! grepl('no-temp', MODEL_FILES)]
  
  cc_proj <-
    map(MODEL_FILES, function(fn) {
      rsf <- readRDS(fn)
      sp <- tolower(fn) %>%
        gsub('models/rsf-', '', .) %>%
        gsub('-2.*', '', .) %>%
        gsub(' ', '-', .) %>%
        gsub('\\(', '', .) %>%
        gsub('\\)', '', .) %>%
        gsub('\\.', '', .)
      
      tibble(
        species = fn %>%
          gsub(pattern = 'models/rsf-', replacement = '', x = .) %>%
          gsub(pattern = '(-).*', replacement = '', x = ., perl = FALSE),
        lab = gsub('southern mountain', '(s. mountain)', x = species) %>%
          gsub('boreal', '(boreal)', x = .)) %>%
        mutate(
          # add bold and italic to labels
          lab = SPECIES_LABS[map_int(lab, \(x) which(SPECIES == x))],
          # get predictions based on the file name of the cc projections
          preds = map(
            paste0('data/weather-projections-', species, '.rds'),
            function(p_fn) {
              .nd <-
                readRDS(p_fn) %>%
                filter(year == 2025 | year == 2100) %>%
                # make sure only the necessary column are kept
                transmute(
                  scenario,
                  year,
                  long,
                  lat,
                  #' use an `animal` from the model to use `discrete = TRUE`
                  animal = rsf$model$animal[1],
                  temperature_C = temp_c,
                  weight, # for weather quantiles
                  forest_perc = terra::extract(f, tibble(long, lat))[, 2],
                  elevation_m = terra::extract(e, tibble(long, lat))[, 2],
                  dist_water_m = terra::extract(w, tibble(long, lat))[, 2]) %>%
                ungroup()
              
              .nd %>%
                mutate(l = predict.bam(rsf, newdata = .nd, type = 'response',
                                       se.fit = FALSE, exclude = EXCLUDE,
                                       discrete = TRUE)) %>%
                return()
            })) %>%
        select(lab, preds)
    }) %>%
    bind_rows() %>%
    unnest(preds) %>%
    # average across day of year
    group_by(scenario, year, lab, long, lat) %>%
    summarize(forest_perc = unique(forest_perc),
              elevation_m = unique(elevation_m),
              dist_water_m = unique(dist_water_m),
              l = weighted.mean(l, w = weight), .groups = 'drop') %>%
    # divide by mean of 2025 to find relative change since 2025
    group_by(lab, long, lat) %>%
    mutate(l_ref = mean(l[year == 2025])) %>%
    filter(year != 2025) %>%
    mutate(l = l / l_ref) %>%
    ungroup() %>%
    mutate(scenario = case_when(
      grepl('126', scenario) ~ '"Best scenario (SSP 1-2.6)"',
      grepl('245', scenario) ~ '"Good scenario (SSP 2-4.5)"',
      grepl('370', scenario) ~ '"Bad scenario (SSP 3-7.0)"',
      grepl('585', scenario) ~ '"Worst scenario (SSP 5-8.5)"') %>%
        factor(., levels = unique(.)))
  cc_proj
  
  saveRDS(cc_proj, 'data/cc-hrsf-projections-local-2100.rds')
  beepr::beep(2)
}

cc_proj <- cc_proj %>%
  mutate(name = COMMON_NAMES[map_int(lab, \(.lab) {
    which(sort(SPECIES_LABS) == .lab)
  })],
  scenario = gsub('"', '', scenario) %>%
    factor(., levels = unique(.)))

slice(cc_proj, 1, .by = lab) %>%
  select(lab, name)

# figure of estimated speeds for each species ----
z_breaks <- seq(-log2(1.25), log2(1.25), length.out = 5)

make_plot <- function(sn, y_facets = FALSE, get_legend = FALSE,
                      reproject = TRUE) {
  sp_data <- cc_proj %>%
    filter(name == sn) %>%
    select(scenario, name, long, lat, l)
  
  if(reproject) {
    sp_data <- sp_data %>%
      nest(rast = ! c(scenario, name)) %>%
      mutate(rast = map(rast, function(.r) {
        rast(.r) %>%
          `crs<-`('EPSG:4326') %>%
          project('EPSG:3005') %>%
          as.data.frame(xy = TRUE) %>%
          rename(long = x, lat = y) %>%
          return()
      })) %>%
      unnest(rast)
  }
  
  p <- sp_data %>%
    # cap values for readability
    mutate(l = case_when(l < 2^min(z_breaks) ~ 2^min(z_breaks) + 1e-10,
                         l > 2^max(z_breaks) ~ 2^max(z_breaks) - 1e-10,
                         TRUE ~ l)) %>%
    ggplot() +
    coord_sf(crs = 'EPSG:3005') +
    facet_grid(scenario ~ name) +
    geom_raster(aes(long, lat, fill = log2(l))) +
    scale_fill_distiller(name = 'Pixel-level relative change in selection strength',
                         palette = 'PiYG', limits = range(z_breaks),
                         breaks = z_breaks,
                         labels = \(x) round(2^x, 2),
                         direction = 1) +
    labs(x = NULL, y = NULL) +
    ggspatial::annotation_scale(style = 'ticks', text_cex = 0.6,
                                location = 'bl', text_face = 'bold') +
    theme(legend.position = 'none', axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill = 'grey'))
  
  if(! y_facets) {
    p <- p + theme(strip.background.y = element_blank(),
                   strip.text.y = element_blank())
  }
  
  if(get_legend) {
    p <- p + theme(legend.position = 'top', legend.key.width = rel(3))
    p <- get_legend(p)
  }
  
  return(p)
}

plot_list <- map(sort(COMMON_NAMES), make_plot)
plot_list[[7]] <- make_plot('Wolves', y_facets = TRUE) # add vertical facet labels

rel_w <- c(1.015, 1.23, 1.36, 0.95, 1.35, 1.42, 1.1)

plot_grid(make_plot(COMMON_NAMES[1], get_legend = TRUE), 
          plot_grid(plotlist = plot_list, nrow = 1, rel_widths = rel_w),
          ncol = 1, rel_heights = c(0.05, 1))

ggsave('figures/local-rss-2100.png',
       width = 17.8, height = 12.4, units = 'in', dpi = 600, bg = 'white')

# make density plots of the current range
cc_proj %>%
  mutate(l = if_else(l > 1.3, 1.3, l)) %>%
  ggplot(aes(x = l, fill = scenario, color = scenario)) +
  facet_wrap(~ name, scales = 'free_y', nrow = 2) +
  geom_density(alpha = 0.25) +
  geom_vline(xintercept = 1, color = 'black', lty = 'dashed') +
  scale_x_continuous('Relative change in RSS in 2100') +
  ylab('Density') +
  scale_color_brewer('Climate change scenario', type = 'div', palette = 5,
                     direction = -1, aesthetics = c('color', 'fill')) +
  theme(legend.position = 'inside', legend.position.inside = c(7/8, 1/4))

ggsave('figures/local-rss-2100-density.png',
       width = 12, height = 6.67, dpi = 600, bg = 'white')

# maps of resources ----
plot_resources <- function(common_name, resource) {
  .d <-
    cc_proj %>%
    filter(scenario == scenario[1], name == common_name)
  
  .d <- rast(.d[, c('long', 'lat', resource)], crs = 'EPSG:4326') %>%
    project('EPSG:3005') %>%
    as.data.frame(xy = TRUE) %>%
    mutate(name = unique(.d$name)) # add species name back in
  
  colnames(.d)[which(colnames(.d) == resource)] <- 'z'
  
  if(resource == 'elevation_m') .d$z <- .d$z / 1e3 # convert to km
  if(resource == 'dist_water_m') .d$z <- .d$z / 1e3 # convert to km
  
  p <-
    ggplot(.d) +
    coord_sf(crs = 'EPSG:3005') +
    facet_grid(. ~ name) +
    geom_raster(aes(x, y, fill = z)) +
    labs(x = NULL, y = NULL) +
    ggspatial::annotation_scale(style = 'ticks', text_cex = 0.6,
                                location = 'bl', text_face = 'bold') +
    theme(legend.position = 'none', axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill = 'grey'))
  
  if(resource == 'forest_perc') {
    p <- p +
      scale_fill_gradient('Tree cover (%)', low = 'white', na.value = NA,
                          high = 'darkgreen', limits = c(0, 100))
  } else if(resource == 'elevation_m') {
    p <- p +
      scale_fill_distiller('Elevation (km)', palette = 6, direction = 1,
                           limits = range(cc_proj$elevation_m) / 1e3)
  } else if(resource == 'dist_water_m') {
    p <- p +
      scale_fill_distiller('Distance from water (km)')
  }
  
  return(p)
}

make_row <- function(resource) {
  plot_list <- map(sort(COMMON_NAMES),
                   \(.z) plot_resources(.z, resource = resource))
  
  return(plot_list)
}

trees <- make_row('forest_perc')
elevs <- make_row('elevation_m')
water <- make_row('dist_water_m')

plot_grid(
  get_legend(trees[[7]]),
  plot_grid(plotlist = trees, nrow = 1, rel_widths = rel_w),
  get_legend(elevs[[7]]),
  plot_grid(plotlist = elevs, nrow = 1, rel_widths = rel_w),
  get_legend(water[[7]]),
  plot_grid(plotlist = water, nrow = 1, rel_widths = rel_w),
  ncol = 1, rel_heights = rep(c(0.2, 1), 3))

ggsave('figures/local-resources.png',
       width = 20, height = 12, units = 'in', dpi = 600, bg = 'white')
