library('mgcv')    # for GAMs
library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('sf')      # for spatial data
library('ctmm')    # for movement models and utilization distributions
library('terra')   # for rasters
library('purrr')   # for functional programming
library('ggplot2') # for fancy plots
library('gratia')  # for graceful GAM plots
library('khroma')  # for colorblind-friendly figures
source('analysis/figures/default-ggplot-theme.R')

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
  f <- rast('data/resource-rasters/forest.tif') %>%
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
    plot(e)
    plot(wide_bounds, add = TRUE)
    plot(ud_bounds, add = TRUE)
    plot(w)
    plot(wide_bounds, add = TRUE)
    plot(ud_bounds, add = TRUE)
    layout(1)
  }
  
  # add resources, responses, and weights observed locations ----
  d_1 <- d_1 %>%
    mutate(forest_perc = extract(f, tibble(longitude, latitude))[, 2],
           elevation_m = extract(e, tibble(longitude, latitude))[, 2],
           dist_water_m = extract(w, tibble(longitude, latitude))[, 2],
           detected = 1)
  
  # make a data frame of null locations ----
  # check range of temperatures
  range(d_1$temperature_C, na.rm = TRUE)
  
  if(file.exists('data/quadrature-data-2025-01-17.rds')) {
    d_0 <- readRDS('data/quadrature-data-2025-01-17.rds')
  } else {
    d_0 <- mm %>%
      transmute(
        animal, # for each animal
        zeros = map2(akde, species, function(.a, .species) {
          # extract the extent of each AKDE
          .a <- SpatialPolygonsDataFrame.UD(.a, level.UD = 0.999, level = 0) %>%
            st_as_sf() %>%
            st_transform('EPSG:4326')
          
          q <-
            # crop and mask forest raster with the AKDE
            f %>%
            crop(.a) %>%
            mask(.a) %>%
            as.data.frame(xy = TRUE) %>% # convert to a data.frame
            rename(forest_perc = consensus_full_class_1) %>%
            mutate(
              # extract corresponding elevation and distance values
              elevation_m = st_as_sf(tibble(x, y), coords = c('x', 'y')) %>%
                extract(x = e, y = .) %>%
                pull(2),
              dist_water_m = st_as_sf(tibble(x, y),
                                      coords = c('x', 'y')) %>%
                extract(x = w, y = .) %>%
                pull(2)) %>%
            rename(longitude = x, latitude = y) %>%
            # drop NAs (masked values)
            filter(! is.na(forest_perc + elevation_m + dist_water_m)) %>%
            mutate(detected = 0, # because it's a quadrature point
                   weight = 1) # all quadrature points are independent
          
          # make a final tibble of the species, temperatures, and q points
          tibble(species = .species,
                 quadrature_points = list(q)) %>%
            unnest(quadrature_points) %>%
            # uniform temperatures to compare to no selection
            mutate(temperature_C = runif(n(), min = -40, max = 40))
        })) %>%
      unnest(zeros)
    saveRDS(d_0, paste0('data/quadrature-data-', Sys.Date(), '.rds'))
  }
  
  # bind observed and quadrature locations together
  d <- bind_rows(d_1, d_0) %>%
    mutate(species = factor(species),
           animal = factor(animal))
  
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

# check goat UDs as an example
if(FALSE) {
  goat_akdes <- readRDS('models/movement-models-akdes-goats-2024-06-06.rds')
  goat_akdes_095 <- map_dfr(goat_akdes$akde, \(.a) {
    SpatialPolygonsDataFrame.UD(.a, level.UD = 0.95) %>%
      st_as_sf() %>%
      slice(2) %>%
      st_transform('EPSG:4326') %>%
      st_geometry() %>%
      st_union() %>%
      st_as_sf()
  })
  
  goat_akdes_0999 <- map_dfr(goat_akdes$akde, \(.a) {
    SpatialPolygonsDataFrame.UD(.a, level.UD = 0.999) %>%
      st_as_sf() %>%
      slice(2) %>%
      st_transform('EPSG:4326') %>%
      st_geometry() %>%
      st_union() %>%
      st_as_sf()
  })
  
  ggplot(goat_akdes) +
    geom_raster(aes(longitude, latitude, fill = dist_water_m),
                readRDS('data/quadrature-data-2025-01-17.rds') %>%
                  filter(species == 'Oreamnos americanus')) +
    geom_sf(data = goat_akdes_0999, fill = 'transparent', lwd = 2,
            color = 'red') +
    geom_sf(data = goat_akdes_095, fill = 'transparent', lwd = 2,
            color = 'black') +
    scale_fill_viridis_c()
  
  rm(goat_akdes, goat_akdes_095, goat_akdes_0999)
}

# all sets of quadrature points have forest values 0-100%
range(d$forest_perc)
d %>%
  filter(detected == 0) %>%
  group_by(animal) %>%
  summarize(min = min(forest_perc),
            max = max(forest_perc),
            .groups = 'drop') %>%
  # summarize across animals
  summarize(min = min(min),
            max = max(max))

# find n(locations) with resource values outside the quadrature values
problematic <-
  d %>%
  group_by(animal) %>%
  mutate(min_ele = min(elevation_m[detected == 0]),
         max_ele = max(elevation_m[detected == 0]),
         min_dis = min(dist_water_m[detected == 0]),
         max_dis = max(dist_water_m[detected == 0]),
  ) %>%
  filter(elevation_m < min_ele |
           elevation_m > max_ele |
           dist_water_m < min_dis |
           dist_water_m > max_dis) %>%
  filter(detected == 1) %>%
  select(c(species, animal, weight,
           min_ele, elevation_m, max_ele,
           min_dis, dist_water_m, max_dis)) %>%
  mutate(problematic = paste(if_else(elevation_m < min_ele |
                                       elevation_m > max_ele,
                                     'elevation', ''),
                             if_else(dist_water_m < min_dis |
                                       dist_water_m > max_dis,
                                     'water', '')))

# only three problematic points because the elevation is slightly too low
problematic

# drop two wolf locations with unusually high elevations
filter(d, species == 'Canis lupus', detected == 1) %>%
  pull(elevation_m) %>%
  quantile(c(0.9, 0.95, 0.99, 0.999, 1))

filter(d, species == 'Canis lupus', detected == 1) %>%
  pull(elevation_m) %>%
  hist(ylim = c(0, 100))
abline(v = 1200, col = 'red', lwd = 2)

d <- filter(d, ! (species == 'Canis lupus' & elevation_m > 1200))

filter(d, species == 'Canis lupus', detected == 1) %>%
  pull(elevation_m) %>%
  hist(ylim = c(0, 200), xlim = c(200, 1200))

# the only caribou (likely tracked twice) that moves up in elevation
filter(d, species == 'Rangifer tarandus boreal',
       elevation_m > 1500) %>% pull(animal) %>% unique()

d <- filter(d, animal != 'SCEK014' & animal != 'SCEK014b')

# fit the RSFs ----
#' *NOTES:*
#' - adding `log()` gives too much leverage to low elevations and
#' distances from water.
#' - dividing `detected` by a large number like `K = 1e6` and multiplying
#'   the weights by `K` does not improve the model fit

# find number of quadrature points per each detection
d %>%
  group_by(species) %>%
  summarize(unweighted = 1 / mean(detected),
            weighted = 1 / mean(detected * weight))

# arrange by number of rows
SPECIES <- d %>%
  group_by(species) %>%
  summarise(n = n()) %>%
  arrange(n) %>%
  pull(species)

COMPLETED <- character(0)

# fit HRSFs with temperature
for(sp in as.character(SPECIES)) {
  cat('Working on ', as.character(sp), '...\n', sep = '')
  print(
    d %>%
      filter(species == sp) %>%
      pivot_longer(c(forest_perc, elevation_m, dist_water_m)) %>%
      mutate(name = factor(name, levels = unique(name))) %>%
      ggplot() +
      facet_wrap(~ name, scales = 'free') +
      geom_histogram(aes(value, fill = factor(detected)),
                     position = 'identity', alpha = 0.5, bins = 10) +
      scale_fill_brewer('Detected', type = 'qual', palette = 6))
  
  #' not using `log1p(dist_water_m)` because many times it jumps from 0 m
  #' to 450 m, which gives a gap from log(0 + 1) = 0 to about
  #' log(451) = 6.1. Using log(x + 100) doesn't change smooths much.
  #' not using `sqrt(dist_water_m)` because trends stay about the same,
  #' but `ti()` terms get too extreme.
  #' using smooth terms gets around the need for nonlinear transformations.
  rsf <- bam(
    detected ~
      # species-level average resource preference
      s(forest_perc, k = 4, bs = 'cr') +
      s(elevation_m, k = 4, bs = 'cr') +
      s(dist_water_m, k = 4, bs = 'cr') +
      # animal-level deviations from the species-level average
      s(animal, bs = 're') +
      s(forest_perc, animal, k = 4, bs = 'fs', xt = list(bc = 'cr')) +
      s(elevation_m, animal, k = 4, bs = 'fs', xt = list(bc = 'cr')) +
      s(dist_water_m, animal, k = 4, bs = 'fs', xt = list(bc = 'cr')) +
      # changes in preference with temperature
      ti(forest_perc, temperature_C, k = 4, bs = 'cr') +
      ti(elevation_m, temperature_C, k = 4, bs = 'cr') +
      ti(dist_water_m, temperature_C, k = 4, bs = 'cr') +
      # include marginals of temperature to post-stratify over afterwards
      s(temperature_C, k = 4, bs = 'cr') +
      s(temperature_C, animal, k = 4, bs = 'fs', xt = list(bc = 'cr')),
    family = poisson(link = 'log'),
    data = d,
    weights = weight,
    subset = species == sp,
    method = 'fREML',
    discrete = TRUE,
    control = gam.control(trace = TRUE))
  
  unique(warnings())
  
  saveRDS(rsf, paste0('models/rsf-', sp, '-', Sys.Date(), '.rds'))
  
  draw(rsf, scales = 'free', rug = FALSE, ci_alpha = 0.05,
       discrete_colour = scale_color_manual(values = rep('#00000080', 300)),
       discrete_fill = scale_fill_manual(values = rep('black', 300)),
       overall_uncertainty = FALSE) %>%
    print()
  
  ggsave(paste0('figures/hrsf-partial-effects/rsf-', sp, '.png'),
         width = 16, height = 8, units = 'in', dpi = 300, bg = 'white',
         scale = 1.5)
  
  summary(rsf, re.test = FALSE)
  
  COMPLETED <- c(COMPLETED, sp)
  cat(paste0('Completed: ', paste(COMPLETED, collapse = ', ')), '\n')
  
  if(FALSE) {
    appraise(rsf, type = 'pearson')
  }
  rm(rsf)
}

# fit HRSFs without temperature
for(sp in as.character(SPECIES)) {
  cat('Working on ', as.character(sp), '...\n', sep = '')
  
  rsf <- bam(
    detected ~
      # species-level average resource preference
      s(forest_perc, k = 4, bs = 'cr') +
      s(elevation_m, k = 4, bs = 'cr') +
      s(dist_water_m, k = 4, bs = 'cr') +
      # animal-level deviations from the species-level average
      s(animal, bs = 're') +
      s(forest_perc, animal, k = 4, bs = 'fs', xt = list(bc = 'cr')) +
      s(elevation_m, animal, k = 4, bs = 'fs', xt = list(bc = 'cr')) +
      s(dist_water_m, animal, k = 4, bs = 'fs', xt = list(bc = 'cr')),
    family = poisson(link = 'log'),
    data = d,
    weights = weight,
    subset = species == sp,
    method = 'fREML',
    discrete = TRUE,
    control = gam.control(trace = TRUE))
  
  saveRDS(rsf,
          paste0('models/rsf-', sp, '-no-temperature-', Sys.Date(), '.rds'))
  
  draw(rsf, scales = 'free', rug = FALSE, ci_alpha = 0.05,
       discrete_colour = scale_color_manual(values = rep('#00000080', 300)),
       discrete_fill = scale_fill_manual(values = rep('black', 300)),
       overall_uncertainty = FALSE) %>%
    print()
  
  ggsave(paste0('figures/hrsf-partial-effects/rsf-', sp,
                '-no-temperature.png'),
         width = 16, height = 8, units = 'in', dpi = 300, bg = 'white',
         scale = 1.5)
  
  summary(rsf, re.test = FALSE)
  
  COMPLETED <- c(COMPLETED, sp)
  cat(paste0('Completed: ', paste(COMPLETED, collapse = ', ')), '\n')
  
  if(FALSE) {
    appraise(rsf, type = 'pearson')
  }
  rm(rsf)
}

# find change in AIC values and deviance explained
change <- tibble(
  species = SPECIES %>%
    gsub('\\(s\\.', 'southern', .) %>%
    gsub('\\(', '', .) %>%
    gsub('\\)', '', .),
  w = map_chr(species, \(.sp) {
    list.files(path = 'models',
               pattern = paste0('rsf-', .sp, '-2025-'),
               full.names = TRUE)
  }),
  wo = map_chr(species, \(.sp) {
    fn <- list.files(path = 'models',
                     pattern = paste0('rsf-', .sp, '-no-temperature-2025-'),
                     full.names = TRUE)
    return(fn)
  })) %>%
  filter(! is.na(wo)) %>%
  mutate(delta_AIC = map2_dbl(w, wo, \(.w, .wo) {
    AIC(readRDS(.w)) - AIC(readRDS(.wo))
  }),
  delta_de = map2_dbl(w, wo, \(.w, .wo) {
    (summary(readRDS(.w), re.test = FALSE)$dev.expl -
       summary(readRDS(.wo), re.test = FALSE)$dev.expl) * 100
  })) %>%
  select(! c(w, wo)) %>%
  arrange(species)

change

# create figure of agreement between pairs of models w and w/o temperature
fits <- tibble(
  species = SPECIES %>%
    gsub('\\(s\\.', 'southern', .) %>%
    gsub('\\(', '', .) %>%
    gsub('\\)', '', .),
  lab = SPECIES_LABS,
  w = map_chr(species, \(.sp) {
    list.files(path = 'models',
               pattern = paste0('rsf-', .sp, '-2025-'),
               full.names = TRUE)
  }),
  wo = map_chr(species, \(.sp) {
    fn <- list.files(path = 'models',
                     pattern = paste0('rsf-', .sp, '-no-temperature-2025-'),
                     full.names = TRUE)
    return(fn)
  })) %>%
  mutate(
    w = map(w, \(.fn) {
      .m <- readRDS(.fn)
      tibble(fits_w = predict(.m, type = 'link'))
    }),
    wo = map(wo, \(.fn) {
      .m <- readRDS(.fn)
      tibble(fits_wo = predict(.m, type = 'link'))
    })) %>%
  unnest(c(w, wo))

p_fits <-
  fits %>%
  ggplot(aes(fits_wo, fits_w)) +
  facet_wrap(~ lab, labeller = label_parsed, scales = 'free') +
  geom_hex(aes(fill = log10(after_stat(count)))) +
  geom_abline(slope = 1, intercept = 0, color = 'black') +
  labs(x = 'log(RSS) without temperature',
       y = 'log(RSS) with temperature') +
  scale_fill_iridescent(name = expression(bold(Count~(log['10']~scale))),
                        range = c(0.3, 1), breaks = (0:3) * 2,
                        labels = round(10^(c(0:3) * 2))) +
  theme(legend.position = 'inside', legend.position.inside = c(5/6, 1/6))

ggsave('figures/hrsf-with-without-temp-prediction-agreement.png',
       plot = p_fits, width = 8, height = 8, units = 'in', dpi = 600,
       bg = 'white')
beepr::beep()
