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

# import climate data ----
if(file.exists('data/cc-hgam-bc-projections-albers.rds')) {
  cc_proj <- readRDS('data/cc-hgam-bc-projections-albers.rds')
} else {
  # import models
  m_1 <- readRDS('models/binomial-gam.rds')
  m_2 <- readRDS('models/gamma-gam.rds')
  
  weather <- readRDS('data/weather-projections.rds') %>%
    filter(year == 2020 | year == 2100)
  
  # terms to exclude from the prediction
  SM <- unique(c(smooths(m_1), smooths(m_2)))
  EXCLUDE <- SM[grepl('tod', SM) | grepl('dt', SM) | grepl('animal', SM)]
  
  # predict movement rates
  bc_preds <- weather %>%
    filter(! is.na(temp_c)) %>%
    select(scenario, year, month, lat, long, temp_c, weight) %>%
    #' #' ******\/ for testing *****************************************
    #' group_by(year, scenario, long) %>%
    #' slice(1:3) %>%
    #' ungroup() %>%
    #' mutate(lat = list(tibble(lat = 45:46))) %>%
    #' unnest(lat) %>%
    #' #' ******/\ for testing *****************************************
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
    # after grouping all the 2020 scenarios together
    mutate(scenario = if_else(year == 2020, '2020', scenario)) %>%
    group_by(scenario, lab, long, lat) %>%
    summarize(p = weighted.mean(p, w = weight),
              s = weighted.mean(s, w = weight),
              d = weighted.mean(d, w = weight),
              .groups = 'drop') %>%
    # calculate change over space relative to the mean 2020 values 
    arrange(lab, scenario, long, lat) %>%
    group_by(lab, long, lat) %>%
    mutate(p = p / mean(p[1:4]),
           s = s / mean(s[1:4]),
           d = d / mean(d[1:4])) %>%
    ungroup() %>%
    filter(scenario != '2020') %>%
    mutate(scenario = case_when(grepl('126', scenario) ~ 'SSP~1-2.6',
                                grepl('245', scenario) ~ 'SSP~2-4.5',
                                grepl('370', scenario) ~ 'SSP~3-7.0',
                                grepl('585', scenario) ~ 'SSP~5-8.5') %>%
             paste0('bold(', ., ')')) %>%
    pivot_wider(values_from = c(p, s, d), names_from = c(scenario, lab)) %>%
    rast(crs = 'EPSG:4326') %>%
    project('EPSG:3005') %>%
    as.data.frame(xy = TRUE) %>%
    pivot_longer(! c(x, y),
                 names_to = c('parameter', 'scenario', 'lab'),
                 names_sep = '_', values_to = 'value') %>%
    pivot_wider(values_from = value, names_from = parameter)
  
  bc_preds
  
  # test the figure
  ggplot(bc_preds, aes(x, y, fill = d)) +
    facet_grid(lab ~ scenario, labeller = label_parsed) +
    geom_raster() +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL)
  
  saveRDS(bc_preds, 'data/cc-hgam-bc-projections-albers.rds')
}

#' **left to edit**
#' # check overall changes from 2020 to 2100
#' preds %>%
#'   select(species, scenario, long, lat, displ) %>%
#'   pivot_wider(names_from = scenario, values_from = displ) %>%
#'   pivot_longer(cols = -c(species, long, lat, `2020`),
#'                names_to = 'scenario', values_to = 'displ') %>%
#'   mutate(change = displ / `2020`) %>%
#'   group_by(species, scenario) %>%
#'   summarise(mean_change = paste0(round(mean(change) * 100 - 100, 1), '%'),
#'             sd_change = paste0(round(sd(change) * 100, 1), '%')) %>%
#'   print(n = length(SPECIES) * 4)
#' 
#' # check changes with density functions
#' preds %>%
#'   select(species, scenario, longitude, latitude, displ) %>%
#'   pivot_wider(names_from = scenario, values_from = displ) %>%
#'   pivot_longer(cols = -c(species, longitude, latitude, `2020`),
#'                names_to = 'scenario', values_to = 'displ') %>%
#'   mutate(change = displ / `2020` * 100 - 100) %>%
#'   ggplot() +
#'   facet_wrap(~ species, scales = 'free_y') +
#'   geom_vline(xintercept = 0, color = 'grey') +
#'   geom_density(aes(change, fill = scenario, color = scenario), alpha = 0.3) +
#'   scale_color_manual('Scenario', values = khroma::color('sunset')(4),
#'                      aesthetics = c('color', 'fill')) +
#'   labs(x = 'Change in distance travelled relative to 2020 (%)',
#'        y = 'Density') +
#'   theme(legend.position = 'inside', legend.position.inside = c(0.85, 0.15))
#' 
#' # figure of estimated speeds for each species ----
#' for(SPECIES in unique(preds$species)) {
#'   cat('Making figure for ', SPECIES, '...\n', sep = '')
#'   
#'   # get phylopic for the species
#'   pic <- case_when(
#'     SPECIES == 'Elk' ~ 'cc03f5c2-933f-4c40-9c64-7f8727556fdb',
#'     SPECIES == 'Mountain goat' ~ 'd6b7ba72-a51b-47ce-a1c9-d8489483ea4c',
#'     SPECIES == 'Cougar' ~ 'ba8012a3-cfc9-4403-b4bb-1959988766ca',
#'     SPECIES == 'Caribou' ~ 'e6e864fd-8e3d-435f-9db3-dc6869c589f1',
#'     SPECIES == 'Grizzly bear' ~ '0cd82109-bb1c-4e08-ab11-c845d8a82eba') %>%
#'     get_phylopic()
#'   
#'   pic <- ggplot() +
#'     add_phylopic(pic) +
#'     theme(panel.border = element_blank())
#'   
#'   preds_sp <-
#'     preds %>%
#'     filter(species == SPECIES)
#'   # # cap displacements at 97.5% quantile
#'   #   mutate(displ = if_else(displ > round(quantile(displ, 0.975), 4),
#'   #                          round(quantile(displ, 0.975), 4), displ)) %>%
#'   #   # cap displacements at 2.5% quantile
#'   #   mutate(displ = if_else(displ < round(quantile(displ, 0.025), 4),
#'   #                          round(quantile(displ, 0.025), 4), displ))
#'   
#'   plt_0 <-
#'     ggplot() +
#'     facet_wrap(~ scenario) +
#'     geom_sf(data = bc_unproj) +
#'     geom_raster(aes(longitude, latitude, fill = displ), preds_sp) +
#'     scale_fill_iridescent(name = 'Distance travelled (km/day)') +
#'     labs(x = NULL, y = NULL) +
#'     coord_sf(xlim = c(-138.5, NA)) +
#'     theme(legend.position = 'none')
#'   
#'   if(SPECIES == 'Caribou') {
#'     bounds <-
#'       read_sf('data/species-ranges-shapefiles/caribou/Herd_Boundaries.shp') %>%
#'       st_geometry() %>%
#'       st_transform(crs = '+proj=longlat') %>%
#'       st_as_sf() %>%
#'       st_make_valid()
#'     
#'     plt_0 <- plt_0 + geom_sf(data = bounds, color = 'black', fill = NA) +
#'       coord_sf(xlim = c(-138.5, NA))
#'   }
#'   
#'   plt <-
#'     ggdraw(plt_0) +
#'     # add the phylopic
#'     draw_plot(pic,
#'               x = case_when(SPECIES == 'Cougar' ~ 0.8,
#'                             SPECIES == 'Grizzly bear' ~ 0.8,
#'                             SPECIES == 'Mountain goat' ~ 0.8,
#'                             TRUE ~ 0.75),
#'               y = 0.05,
#'               width = if_else(SPECIES %in% c('Elk', 'Caribou'), 0.3, 0.2),
#'               height = if_else(SPECIES %in% c('Elk', 'Caribou'), 0.3, 0.2)) +
#'     # add the legend over the phylopic
#'     draw_grob(get_legend(plt_0 + theme(legend.position = 'right',
#'                                        legend.background = element_blank())),
#'               x = 0.35, y = -0.15)
#'   
#'   #' edit `SPECIES` appropriately for a file name
#'   SPECIES <- tolower(SPECIES) %>%
#'     str_replace_all(pattern = ' ', replacement = '-')
#'   
#'   ggsave(paste0('figures/bc-speeds-', SPECIES, '.png'),
#'          plot = plt, scale = 0.75, width = 13, height = 8, dpi = 300,
#'          bg = 'white', device = 'png', type = 'cairo')
#' }
