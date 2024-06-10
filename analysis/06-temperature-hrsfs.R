library('mgcv')  # for GAMs
library('dplyr') # for data wrangling
library('tidyr') # for data wrangling
library('sf')    # for spatial data
library('ctmm')  # for movement models and utilization distributions
library('terra') # for rasters
library('purrr') # for functional programming

# K <- 1e6 # scaling constant for weights

# import data ----
if(file.exists('data/tracking-data/rsf-data.rds')) {
  d <- readRDS('data/tracking-data/rsf-data.rds')
} else {
  # import movement models
  mm <- readRDS('models/movement-models-akdes-2024-03-18.rds') %>%
    # remove animals that are not range resident (invalid weights)
    filter(! animal %in% c('SCEK014', 'SCEK014b', 'BW028')) %>%
    mutate(species = if_else(grepl('Rangifer', dataset_name),
                             dataset_name,
                             species),
           species = stringr::str_replace_all(species, '_', ' ') %>%
             factor())
  
  # import tracking data with annonated temperatures
  d_1 <-
    readRDS('data/movement-models-speed-weights-temperature-2024-04-05.rds') %>%
    # remove animals that are not range resident (invalid weights)
    filter(! animal %in% c('SCEK014', 'SCEK014b', 'BW028')) %>%
    mutate(species = if_else(grepl('Rangifer', dataset_name),
                             dataset_name,
                             species),
           species = stringr::str_replace_all(species, '_', ' ') %>%
             factor()) %>%
    select(species, animal, longitude, latitude, temperature_C, weight)
  
  locs_mcp <- map_dfr(mm$tel, \(.t) {
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
  
  hr_bounds <- map_dfr(mm$akde, \(.a) {
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
  plot(wide_bounds)
  
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
    plot(hr_bounds, add = TRUE)
    plot(locs_mcp, add = TRUE, col = 'red3')
    plot(e)
    plot(wide_bounds, add = TRUE)
    plot(hr_bounds, add = TRUE)
    plot(locs_mcp, add = TRUE, col = 'red3')
    plot(w)
    plot(wide_bounds, add = TRUE)
    plot(hr_bounds, add = TRUE)
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
  
  if(! file.exists('data/quadrature-data-2024-05-07.rds')) {
    d_0 <- mm %>%
      transmute(
        animal,
        zeros = map2(akde, species, function(.a, .species) {
          .a <- SpatialPolygonsDataFrame.UD(.a) %>%
            st_as_sf() %>%
            st_transform('EPSG:4326')
          
          null_temperatures <- filter(d_1, species == .species) %>%
            pull(temperature_C) %>%
            quantile(probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE) %>%
            unname()
          
          q <- f %>%
            crop(.a) %>%
            mask(.a) %>%
            as.data.frame(xy = TRUE) %>%
            rename(forest_perc = consensus_full_class_1) %>%
            mutate(
              elevation_m = st_as_sf(tibble(x, y), coords = c('x', 'y')) %>%
                extract(x = e, y = .) %>%
                pull(2),
              dist_water_m = st_as_sf(tibble(x, y),
                                      coords = c('x', 'y')) %>%
                extract(x = w, y = .) %>%
                pull(2)) %>%
            rename(long = x, lat = y) %>%
            filter(! is.na(forest_perc + elevation_m + dist_water_m)) %>%
            mutate(detected = 0,
                   weight = 1)
          
          tibble(species = .species,
                 temperature_C = null_temperatures,
                 quadrature_points = list(q)) %>%
            unnest(quadrature_points)
        })) %>%
      unnest(zeros)
    saveRDS(d_0, paste0('data/quadrature-data-', Sys.Date(), '.rds'))
  } else {
    d_0 <- readRDS('data/quadrature-data-2024-05-07.rds')
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
  
  rm(d_0, d_1, e, f, hr_bounds, mm, w, wide_bounds)
}

# fit the RSF ----
#' - adding `log()` gives too much leverage to low elevations and
#' distances from water.
#' - dividing `detected` by `K` and adding `K` in the weights does not
#' improve the model fit
#' - adding the AKDE weights makes everything flat
#' took 18 minutes for goats and bears only (1.5% of the data)
#' **need to fit separate RSFs**
rsf <- bam(
  detected ~
    # average resource preference
    s(forest_perc, by = species, k = 6, bs = 'cr') +
    s(elevation_m, by = species, k = 6, bs = 'cr') +
    s(dist_water_m, by = species, k = 6, bs = 'cr') +
    # animal-level deviations from the species average
    s(forest_perc, animal, by = species, k = 6, bs = 'sz', xt = list(bs = 'cr')) +
    s(elevation_m, animal, by = species, k = 6, bs = 'sz', xt = list(bs = 'cr')) +
    s(dist_water_m, animal, by = species, k = 6, bs = 'sz', xt = list(bs = 'cr')) +
    # changes in preference with temperature
    ti(forest_perc, by = species, temperature_C, k = 6, bs = 'cr') +
    ti(elevation_m, by = species, temperature_C, k = 6, bs = 'cr') +
    ti(dist_water_m, by = species, temperature_C, k = 6, bs = 'cr'),
  family = poisson(link = 'log'),
  data = d,
  method = 'fREML',
  discrete = TRUE,
  control = gam.control(trace = TRUE))

saveRDS(rsf, paste0('models/rsf-', Sys.Date(), '.rds'))
plot(rsf, pages = 1, scheme = 3, scale = 0)
summary(rsf, re.test = FALSE)

if(FALSE) {
  layout(matrix(1:4, ncol = 2))
  gam.check(rsf, type = 'pearson')
  layout(1)
}
