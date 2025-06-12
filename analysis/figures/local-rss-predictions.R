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
source('functions/get_legend.R') # 
source('data/bc-shapefile.R') # import shapefile of bc

colorRampPalette(RColorBrewer::brewer.pal(11, 'PiYG'))(1e3) %>%
  plot_scheme_colorblind()

if(file.exists('data/cc-hrsf-projections-local-2100.rds')) {
  local_2100_proj <- readRDS('data/cc-hrsf-projections-local-2100.rds')
} else {
  # import resource rasters
  f <- rast('data/resource-rasters/forest.tif')
  e <- rast('data/resource-rasters/bc-buffered-dem-z6.tif')
  w <- rast('data/resource-rasters/distance-from-water.tif')
  
  # terms to exclude from the predictions
  SM <- gratia::smooths(readRDS('models/rsf-Oreamnos americanus-2025-01-20.rds'))
  EXCLUDE <- SM[grepl('animal', SM)]
  
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

# figure of estimated speeds for each species ----
z_breaks <- seq(-log2(1.5), log2(1.5), length.out = 5)

make_plot <- function(sp, y_facets = FALSE, get_legend = FALSE,
                      reproject = TRUE) {
  sp_data <- cc_proj %>%
    filter(lab == sp) %>%
    select(scenario, lab, long, lat, l)
  
  if(reproject) {
    sp_data <- sp_data %>%
      nest(rast = ! c(scenario, lab)) %>%
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
    facet_grid(scenario ~ lab, labeller = label_parsed) +
    geom_raster(aes(long, lat, fill = log2(l))) +
    scale_fill_distiller(name = 'Relative change in selection strength',
                         palette = 'PiYG', limits = range(z_breaks),
                         breaks = z_breaks,
                         labels = \(x) round(2^x, 2),
                         direction = 1) +
    scale_x_continuous(NULL) +
    scale_y_continuous(NULL) +
    ggspatial::annotation_scale(style = 'ticks', text_cex = 0.6,
                                location = 'tr', text_face = 'bold') +
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

plot_list <- map(sort(as.character(SPECIES_LABS)), make_plot)
plot_list[[7]] <- make_plot('bolditalic(Ursus~arctos~horribilis)',
                            y_facets = TRUE)

# get approximate relative widths
cc_proj %>%
  arrange(lab) %>%
  summarize(ratio = diff(range(long)) / diff(range(lat)),
            .by = lab) %>%
  pull(ratio) %>%
  round(2) %>%
  cat(sep = ', ')

plot_grid(
  make_plot(SPECIES_LABS[2], get_legend = TRUE), 
  plot_grid(plotlist = plot_list, nrow = 1,
            rel_widths = c(1.32, 1.1, 1.67, 1.64, 1.57, 1.58, 1.64)),
  ncol = 1, rel_heights = c(0.05, 1))

ggsave('figures/local-rss-2100.png', width = 15, height = 9.65,
       units = 'in', dpi = 600, bg = 'white', scale = 1.2)

cc_proj %>%
  summarize(mean = quantile(l, 0.95) > 1,
            .by = c(lab, scenario)) %>%
  filter(mean) %>%
  arrange(lab)

# maps of resources ----
plot_resources <- function(sp, resource = 'forest_perc') {
  .d <-
    cc_proj %>%
    filter(scenario == scenario[1], lab == sp)
  
  .d <- rast(.d[, c('long', 'lat', resource)], crs = 'EPSG:4326') %>%
    project('EPSG:3005') %>%
    as.data.frame(xy = TRUE) %>%
    mutate(lab = unique(.d$lab))
  
  colnames(.d)[which(colnames(.d) == resource)] <- 'z'
  
  if(resource == 'dist_water_m') .d$z <- .d$z / 1e3
  
  p <-
    ggplot(.d) +
    coord_sf(crs = 'EPSG:3005') +
    facet_grid(. ~ lab, labeller = label_parsed) +
    geom_raster(aes(x, y, fill = z)) +
    scale_x_continuous(NULL) +
    scale_y_continuous(NULL) +
    ggspatial::annotation_scale(style = 'ticks', text_cex = 0.6,
                                location = 'tr', text_face = 'bold') +
    theme(legend.position = 'none', axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill = 'grey'))
  
  if(resource == 'forest_perc') {
    p <- p +
      scale_fill_gradient('Tree cover (%)', low = 'white', na.value = NA,
                          high = 'darkgreen', limits = c(0, 100))
  } else if(resource == 'elevation_m') {
    p <- p +
      scale_fill_distiller('Elevation (m)', palette = 6, direction = 1,
                           limits = range(cc_proj$elevation_m))
  } else if(resource == 'dist_water_m') {
    p <- p +
      scale_fill_distiller('Distance from water (km)')
  }
  
  return(p)
}

trees <- map(sort(as.character(SPECIES_LABS)), plot_resources)
elevs <- map(sort(as.character(SPECIES_LABS)),
             \(.z) plot_resources(.z, resource = 'elevation_m'))
water <- map(sort(as.character(SPECIES_LABS)),
             \(.z) plot_resources(.z, resource = 'dist_water_m'))

plot_grid(
  get_legend(trees[[7]]),
  plot_grid(plotlist = trees, nrow = 1,
            rel_widths = c(1.32, 1.1, 1.67, 1.64, 1.57, 1.58, 1.64)),
  get_legend(elevs[[7]]),
  plot_grid(plotlist = elevs, nrow = 1,
            rel_widths = c(1.32, 1.1, 1.67, 1.64, 1.57, 1.58, 1.64)),
  get_legend(water[[7]]),
  plot_grid(plotlist = water, nrow = 1,
            rel_widths = c(1.32, 1.1, 1.67, 1.64, 1.57, 1.58, 1.64)),
  ncol = 1, rel_heights = rep(c(0.2, 1), 3))

ggsave('figures/local-resources.png', width = 15, height = 8.5,
       units = 'in', dpi = 600, bg = 'white', scale = 1.5)
