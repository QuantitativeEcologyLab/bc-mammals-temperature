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
source('functions/get_legend.R') #' cowplot's `get_legend()` fails (v1.1.3)
source('data/bc-shapefile.R') # import shapefile of bc

colorRampPalette(RColorBrewer::brewer.pal(11, 'PuOr'))(1e3) %>%
  plot_scheme_colorblind()

if(file.exists('data/cc-hgam-projections-local-2100.rds')) {
  cc_proj <- readRDS('data/cc-hgam-projections-local-2100.rds')
} else {
  # import models
  m_1 <- readRDS('models/binomial-gam.rds')
  m_2 <- readRDS('models/gamma-gam.rds')
  
  # when predicting: only excluding animal-level random effect and choosing
  # specific values for the rest
  
  # import prediction data for each species
  #' from `analysis/figures/temperature-hgams-projections.R`
  cc_newd <- readRDS('data/hgam-cc_newd.rds')
  
  cc_proj <- cc_newd %>%
    filter(year == 2025 | year == 2100) %>%
    mutate(.,
           p = predict(m_1, newdata = ., type = 'response', se.fit = FALSE,
                       discrete = TRUE, exclude = 's(animal)'),
           s = predict(m_2, newdata = ., type = 'response', se.fit = FALSE,
                       discrete = TRUE, exclude = 's(animal)'),
           d = p * s) %>%
    # average across day of year
    group_by(scenario, year, species, long, lat) %>%
    summarize(p = weighted.mean(p, w = weight),
              s = weighted.mean(s, w = weight),
              d = weighted.mean(d, w = weight),
              .groups = 'drop') %>%
    # divide by mean in 2025 at each pixel to find rel change since 2025
    # grouping by pixel so that we look at the pixel-level change
    # if we don't then values that are higher than the average will show as
    # higher despite not being an increase
    group_by(species, long, lat) %>%
    mutate(p_ref = mean(p[year == 2025]),
           s_ref = mean(s[year == 2025]),
           d_ref = mean(d[year == 2025])) %>%
    filter(year != 2025) %>%
    mutate(p = p / p_ref,
           s = s / s_ref,
           d = d / d_ref) %>%
    ungroup() %>%
    mutate(scenario = case_when(
      grepl('126', scenario) ~ '"Best scenario (SSP 1-2.6)"',
      grepl('245', scenario) ~ '"Good scenario (SSP 2-4.5)"',
      grepl('370', scenario) ~ '"Bad scenario (SSP 3-7.0)"',
      grepl('585', scenario) ~ '"Worst scenario (SSP 5-8.5)"') %>%
        factor(., levels = unique(.)),
      species = gsub(' ', '~', species) %>%
        gsub('~\\(', '\\)~bold\\((', .) %>%
        paste0('bolditalic(', ., ')') %>%
        gsub('\\(boreal\\)', '"\\(boreal\\)"', .) %>%
        gsub('\\(s.~mountain\\)', '"\\(s. mountain\\)"', .) %>%
        factor())
  cc_proj
  
  saveRDS(cc_proj, 'data/cc-hgam-projections-local-2100.rds')
  beepr::beep(2)
}

# add common species names
cc_proj <- cc_proj %>%
  mutate(lab = COMMON_NAMES[map_int(species, \(.s) {
    which(sort(SPECIES_LABS) == .s)
  })],
  scenario = gsub('"', '', scenario) %>%
    factor(., levels = unique(.)))

# figure of estimated speeds for each species ----
range(c(range(cc_proj$p), range(cc_proj$s), range(cc_proj$d)))

z_breaks <- seq(log2(0.75), -log2(0.75), length.out = 5)

make_plot <- function(.lab, y_facets = FALSE, get_legend = FALSE,
                      reproject = TRUE, variable) {
  sp_data <- filter(cc_proj, lab == .lab)
  
  if(variable == 'p') {
    sp_data <- sp_data %>% select(scenario, lab, long, lat, p)
    z_lab <- 'probability of moving'
  } else if(variable == 's') {
    sp_data <- sp_data %>% select(scenario, lab, long, lat, s)
    z_lab <- 'speed when moving'
  } else if(variable == 'd') {
    sp_data <- sp_data %>% select(scenario, lab, long, lat, d)
    z_lab <- 'distance travelled'
  } else stop('Please choose a variable among p, s, or d.')
  
  sp_data <- rename(sp_data, z = ncol(sp_data))
  z_lab <- paste('Pixel-level relative change in', z_lab)
  
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
    mutate(z = case_when(z < 2^min(z_breaks) ~ 2^min(z_breaks) + 1e-10,
                         z > 2^max(z_breaks) ~ 2^max(z_breaks) - 1e-10,
                         TRUE ~ z)) %>%
    ggplot() +
    coord_sf(crs = 'EPSG:3005') +
    facet_grid(scenario ~ lab) +
    geom_raster(aes(long, lat, fill = log2(z))) +
    scale_fill_distiller(name = z_lab,
                         palette = 'PuOr', limits = range(z_breaks),
                         breaks = z_breaks,
                         labels = \(x) round(2^x, 2),
                         direction = 1) +
    ggspatial::annotation_scale(style = 'ticks', text_cex = 0.6,
                                location = 'bl', text_face = 'bold') +
    labs(x = NULL, y = NULL) +
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

make_full_plot <- function(variable) {
  plot_list <- map(sort(COMMON_NAMES),
                   \(.cn) make_plot(.lab = .cn, variable = variable))
  
  # add y facets to rightmost plot
  plot_list[[7]] <- make_plot(.lab = levels(COMMON_NAMES)[7],
                              y_facets = TRUE, variable = variable)
  
  plot_grid(
    make_plot(COMMON_NAMES[1], get_legend = TRUE, variable = variable),
    plot_grid(plotlist = plot_list, nrow = 1,
              rel_widths = c(1.015, 1.23, 1.36, 0.95, 1.35, 1.42, 1.1)),
  ncol = 1, rel_heights = c(0.05, 1))
}

ggsave('figures/local-p-moving-2100.png', make_full_plot('p'),
       width = 17.8, height = 12.4, units = 'in', dpi = 600, bg = 'white')

ggsave('figures/local-speed-2100.png', make_full_plot('s'),
       width = 17.8, height = 12.4, units = 'in', dpi = 600, bg = 'white')

ggsave('figures/local-distance-2100.png', make_full_plot('d'),
       width = 17.8, height = 12.4, units = 'in', dpi = 600, bg = 'white')
