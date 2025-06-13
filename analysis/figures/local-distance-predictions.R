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
    # divide by mean of 2025 to find relative change since 2025
    group_by(species) %>%
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

# figure of estimated speeds for each species ----
z_breaks <- seq(-log2(1.5), log2(1.5), length.out = 5)

make_plot <- function(sp, y_facets = FALSE, get_legend = FALSE,
                      reproject = TRUE) {
  sp_data <- cc_proj %>%
    filter(species == sp) %>%
    select(scenario, species, long, lat, d)
  
  if(reproject) {
    sp_data <- sp_data %>%
      nest(rast = ! c(scenario, species)) %>%
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
    mutate(d = case_when(d < 2^min(z_breaks) ~ 2^min(z_breaks) + 1e-10,
                         d > 2^max(z_breaks) ~ 2^max(z_breaks) - 1e-10,
                         TRUE ~ d)) %>%
    ggplot() +
    coord_sf(crs = 'EPSG:3005') +
    facet_grid(scenario ~ species, labeller = label_parsed) +
    geom_raster(aes(long, lat, fill = log2(d))) +
    scale_fill_distiller(name = 'Relative change in distance travelled',
                         palette = 'PuOr', limits = range(z_breaks),
                         breaks = z_breaks,
                         labels = \(x) round(2^x, 2),
                         direction = 1) +
    scale_x_continuous(NULL) +
    scale_y_continuous(NULL, expand = c(0.2, 0)) +
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
  arrange(species) %>%
  summarize(ratio = diff(range(long)) / diff(range(lat)),
            .by = species) %>%
  pull(ratio) %>%
  round(2) %>%
  cat(sep = ', ')

plot_grid(
  make_plot(SPECIES_LABS[1], get_legend = TRUE), 
  plot_grid(plotlist = plot_list, nrow = 1,
          rel_widths = c(1.32, 1.1, 1.67, 1.64, 1.57, 1.58, 1.64)),
  ncol = 1, rel_heights = c(0.05, 1))

ggsave('figures/local-distance-2100.png', width = 15, height = 10,
       units = 'in', dpi = 600, bg = 'white', scale = 1.5)
