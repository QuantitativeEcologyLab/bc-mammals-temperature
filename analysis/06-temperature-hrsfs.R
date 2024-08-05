library('mgcv')    # for GAMs
library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('sf')      # for spatial data
library('ctmm')    # for movement models and utilization distributions
library('terra')   # for rasters
library('purrr')   # for functional programming
library('ggplot2') # for fancy plots
library('gratia')  # for graceful GAM plots
source('analysis/figures/default-ggplot-theme.R')

# K <- 1e6 # scaling constant for weights

# import data ----
if(file.exists('data/tracking-data/rsf-data.rds')) {
  d <- readRDS('data/tracking-data/rsf-data.rds')
} else {
  ANIMALS <- readr::read_csv('data/tracking-data/telemetry-metadata.csv',
                             show_col_types = FALSE) %>%
    filter(keep_for_RSF == 'yes') %>%
    pull(animal)
  
  # import movement models
  mm <- readRDS('models/movement-models-akdes-2024-06-06.rds') %>%
    # remove animals that are not range resident (invalid weights)
    filter(animal %in% ANIMALS) %>%
    mutate(species = if_else(grepl('Rangifer', dataset_name),
                             dataset_name,
                             species),
           species = stringr::str_replace_all(species, '_', ' ') %>%
             factor())
  
  # import tracking data with annotated temperatures
  d_1 <-
    readRDS('data/movement-models-speed-weights-temperature-2024-06-10.rds') %>%
    # remove animals that are not range resident (invalid weights)
    filter(animal %in% ANIMALS) %>%
    mutate(species = if_else(grepl('Rangifer', dataset_name),
                             dataset_name,
                             species),
           species = stringr::str_replace_all(species, '_', ' ') %>%
             factor()) %>%
    select(species, animal, longitude, latitude, temperature_C, weight)
  
  # minimum convex polygon to compare data extent to AKDEs
  locs_mcp <-
    map_dfr(mm$tel, \(.t) {
      .t %>%
        SpatialPoints.telemetry() %>%
        adehabitatHR::mcp() %>%
        st_as_sf() %>%
        st_geometry() %>%
        #' st_make_valid() %>% #' use if `Error: Loop x edge x has duplicate...`
        st_union() %>% # unite estimate with lower and upper 95% CIs
        st_transform('EPSG:4326') %>%
        st_as_sf()
    }) %>%
    st_make_valid() %>%
    st_union() %>%
    st_as_sf()
  
  ud_bounds <-
    map_dfr(mm$akde, \(.a) {
      SpatialPolygonsDataFrame.UD(.a, level.UD = 0.95, level = 0) %>%
        st_as_sf() %>%
        st_geometry() %>%
        #' st_make_valid() %>% #' use if `Error: Loop x edge x has duplicate...`
        st_union() %>% # unite estimate with lower and upper 95% CIs
        st_set_crs(.a@info$projection) %>%
        st_transform('EPSG:4326') %>%
        st_as_sf()
    }) %>%
    st_make_valid() %>%
    st_union() %>%
    st_as_sf()
  
  wide_bounds <- map_dfr(mm$akde, \(.a) {
    SpatialPolygonsDataFrame.UD(.a, level.UD = 0.999, level = 0) %>%
      st_as_sf() %>%
      st_geometry() %>%
      #' st_make_valid() %>% #' use if `Error: Loop x edge x has duplicate...`
      st_union() %>% # unite estimate with lower and upper 95% CIs
      st_set_crs(.a@info$projection) %>%
      st_transform('EPSG:4326') %>%
      st_as_sf()
  }) %>%
    st_make_valid() %>%
    st_union() %>%
    st_as_sf()
  
  # import rasters of resources
  f <- rast('data/resource-rasters/uncropped-rasters/consensus_full_class_1.tif') %>%
    crop(wide_bounds) %>%
    mask(wide_bounds)
  e <- rast('data/resource-rasters/bc-buffered-dem-z6.tif') %>%
    crop(wide_bounds) %>%
    mask(wide_bounds)
  w <- rast('data/resource-rasters/distance-from-water.tif') %>%
    crop(wide_bounds) %>%
    mask(wide_bounds)
  
  if(FALSE) {
    layout(t(1:3))
    plot(f)
    # plot(locs, add = TRUE)
    plot(wide_bounds, add = TRUE)
    plot(ud_bounds, add = TRUE)
    plot(locs_mcp, add = TRUE, col = 'red3')
    plot(e)
    plot(wide_bounds, add = TRUE)
    plot(ud_bounds, add = TRUE)
    plot(locs_mcp, add = TRUE, col = 'red3')
    plot(w)
    plot(wide_bounds, add = TRUE)
    plot(ud_bounds, add = TRUE)
    plot(locs_mcp, add = TRUE, col = 'red3')
    layout(1)
  }
  
  # add resources, responses, and weights observed locations ----
  d_1 <- d_1 %>%
    mutate(forest_perc = extract(f, tibble(longitude, latitude))[, 2],
           elevation_m = extract(e, tibble(longitude, latitude))[, 2],
           dist_water_m = extract(w, tibble(longitude, latitude))[, 2],
           detected = 1,
           weight = weight)
  
  # data frame of null locations ----
  # check range of temperatures
  range(d_1$temperature_C, na.rm = TRUE)
  
  if(! file.exists('data/quadrature-data-2024-08-02.rds')) {
    d_0 <- mm %>%
      transmute(
        animal, # for each animal
        zeros = map2(akde, species, function(.a, .species) {
          # extract the extent of each AKDE
          .a <- SpatialPolygonsDataFrame.UD(.a) %>%
            st_as_sf() %>%
            st_transform('EPSG:4326')
          
          # uniform temperatures to compare to no selection
          null_temperatures <- filter(d_1, species == .species) %>%
            pull(temperature_C) %>%
            quantile(probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE) %>%
            unname()
          
          q <-
            # crop and mask forest raster with the AKDE
            f %>%
            crop(.a) %>%
            mask(.a) %>%
            as.data.frame(xy = TRUE) %>% # convert to a data.frame
            rename(forest_perc = consensus_full_class_1) %>%
            mutate( # extract corresponding elevation and distance values
              elevation_m = st_as_sf(tibble(x, y), coords = c('x', 'y')) %>%
                extract(x = e, y = .) %>%
                pull(2),
              dist_water_m = st_as_sf(tibble(x, y),
                                      coords = c('x', 'y')) %>%
                extract(x = w, y = .) %>%
                pull(2)) %>%
            rename(long = x, lat = y) %>%
            # drop NAs (masked values)
            filter(! is.na(forest_perc + elevation_m + dist_water_m)) %>%
            mutate(detected = 0, # because it's a quadrature point
                   weight = 1) # all quadrature points are independent
          
          # make a final tibble of the species, temperatures, and q points
          tibble(species = .species,
                 temperature_C = null_temperatures,
                 quadrature_points = list(q)) %>%
            unnest(quadrature_points)
        })) %>%
      unnest(zeros)
    saveRDS(d_0, paste0('data/quadrature-data-', Sys.Date(), '.rds'))
  } else {
    d_0 <- readRDS('data/quadrature-data-2024-08-02.rds')
  }
  
  # bind observed and quadrature locations together
  d <- bind_rows(d_1, d_0) %>%
    mutate(species = factor(species),
           animal = factor(animal))
  
  d %>%
    group_by(species) %>%
    summarise(ratio = sum(weight * (1 - detected)))
  
  saveRDS(d, 'data/tracking-data/rsf-data.rds')
  
  # check coverage ----
  if(FALSE) {
    layout(matrix(1:6, ncol = 2))
    hist(d$forest_perc, breaks = 10)
    hist(d$elevation_m, breaks = 10)
    hist(d$dist_water_m, breaks = 10)
    hist(d_1$forest_perc, xlim = range(d$forest_perc), breaks = 10)
    hist(d_1$elevation_m, xlim = range(d$elevation_m), breaks = 10)
    hist(d_1$dist_water_m, xlim = range(d$dist_water_m), breaks = 10)
    layout(1)
  }
  
  rm(d_0, d_1, e, f, ud_bounds, mm, w, wide_bounds)
}

# fit the RSFs ----
#' - adding `log()` gives too much leverage to low elevations and
#' distances from water.
#' - dividing `detected` by `K` and adding `K` in the weights does not
#' improve the model fit
#' - adding the AKDE weights makes everything flat
#' took 18 minutes for goats and bears only (1.5% of the data)
SPECIES <- unique(d$species)

# find number of quadrature points per each detection
d %>%
  group_by(species) %>%
  summarize(unweighted = 1 / mean(d$detected),
            weighted = 1 / mean(detected * weight))

for(sp in SPECIES) {
  d %>%
    filter(species == sp) %>%
    pivot_longer(c(forest_perc, elevation_m, dist_water_m)) %>%
    ggplot() +
    facet_wrap(~ name, scales = 'free') +
    geom_histogram(aes(value, fill = factor(detected)),
                   position = 'identity', alpha = 0.5, bins = 10) +
    scale_fill_brewer('Detected', type = 'qual', palette = 6)
  
  rsf <- bam(
    detected ~
      # species-level average resource preference
      s(forest_perc, k = 6, bs = 'cr') +
      s(elevation_m, k = 6, bs = 'cr') +
      s(dist_water_m, k = 6, bs = 'cr') +
      # animal-level deviations from the species-level average
      s(animal, bs = 're') +
      s(forest_perc, animal, k = 6, bs = 'fs', xt = list(bs = 'cr')) +
      s(elevation_m, animal, k = 6, bs = 'fs', xt = list(bs = 'cr')) +
      s(dist_water_m, animal, k = 6, bs = 'fs', xt = list(bs = 'cr')) +
      # changes in preference with temperature
      ti(forest_perc, temperature_C, k = 6, bs = 'cr') +
      ti(elevation_m, temperature_C, k = 6, bs = 'cr') +
      ti(dist_water_m, temperature_C, k = 6, bs = 'cr'),
    family = poisson(link = 'log'),
    data = d,
    weights = weight,
    subset = species == sp,
    method = 'fREML',
    discrete = TRUE,
    control = gam.control(trace = TRUE))
  
  saveRDS(rsf, paste0('models/rsf-', sp, '-', Sys.Date(), '.rds'))
  draw(rsf, scales = 'free', rug = FALSE, ci_alpha = 0.05,
       discrete_colour = scale_color_manual(values = rep('#00000080', 300)),
       discrete_fill = scale_fill_manual(values = rep('black', 300)))
  
  ggsave(filename = paste0('figures/rsf-', sp, '.png'), width = 8, height = 8,
         units = 'in', dpi = 600, bg = 'white')
  
  summary(rsf, re.test = FALSE)
}

if(FALSE) {
  appraise(rsf, type = 'pearson')
}
