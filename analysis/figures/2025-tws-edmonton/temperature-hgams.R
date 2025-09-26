library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for smoother date wrangling
library('mgcv')      # for Generalized Additive Models
library('gratia')    # for useful convenience functions for GAMs
library('ggplot2')   # for fancy plots
library('khroma')    # for colorblind-friendly color palettes
library('cowplot')   # for fancy multi-panel plots
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids
source('functions/get_legend.R') # to extract legends from ggplot plots
plot_scheme(PAL, colours = TRUE)

# import data
d <- readRDS('data/hgam-speed-data.rds') %>%
  filter(species %in% c('Rangifer tarandus (boreal)', 'Canis lupus'))

SPECIES <- unique(d$species)
COMMON_NAMES <- COMMON_NAMES[c(5, 1)]

# import models
m_1 <- readRDS('models/binomial-gam.rds')
m_2 <- readRDS('models/gamma-gam.rds')

# backgrounds for plots ----
tod <- tibble(x = c(seq(0, 6, by = 0.01),
                    seq(18, 24, by = 0.01)),
              g = c(rep(1, 601), rep(2, 601)))

seasons <- tibble(x = seq(1, 366, by = 0.1),
                  y = 0.49,
                  season = case_when(x < yday('2025-03-21') ~ 'w2',
                                     x < yday('2025-06-21') ~ 'sp',
                                     x < yday('2025-09-21') ~ 'su',
                                     x < yday('2025-12-21') ~ 'f',
                                     TRUE ~ 'w1') %>%
                    factor(levels = c('sp', 'su', 'f', 'w1', 'w2')))

PAL_SEASONS <- c('forestgreen', 'goldenrod', 'darkorange4', 'grey90', 'grey90')

season_breaks <-
  as.Date(paste0('2024-', c('03-20', '06-20', '09-22', '12-21'))) + 45

# marginal term plots ----
#' not averaging across uniform `tod_pdt` and `doy` because I want to
#' preserve the uneaven sampling of `doy` for `E(speed | moving)`
marginal <- function(newd, term) {
  newd <- newd %>%
    mutate(lab = case_when(species == 'Rangifer tarandus (s. mountain)' ~ 'bolditalic(Rangifer~tarandus)~bold("(s. mountain)")',
                           species == 'Ursus arctos horribilis' ~ 'bolditalic(Ursus~arctos~horribilis)',
                           species == 'Puma concolor' ~ 'bolditalic(Puma~concolor)',
                           species == 'Cervus canadensis' ~ 'bolditalic(Cervus~canadensis)',
                           species == 'Rangifer tarandus (boreal)' ~ 'bolditalic(Rangifer~tarandus)~bold("(boreal)")',
                           species == 'Canis lupus' ~ 'bolditalic(Canis~lupus)',
                           species == 'Oreamnos americanus' ~ 'bolditalic(Oreamnos~americanus)'))
  
  preds_1 <-
    predict(object = m_1, newdata = newd, type = 'link', se.fit = TRUE,
            exclude = 's(animal)', discrete = TRUE) %>%
    as.data.frame() %>%
    transmute(p_mu = m_1$family$linkinv(fit),
              p_lwr = m_1$family$linkinv(fit - 1.96 * se.fit),
              p_upr = m_1$family$linkinv(fit + 1.96 * se.fit))
  
  preds_2 <-
    predict(object = m_2, newdata = newd, type = 'link', se.fit = TRUE,
            exclude = 's(animal)', discrete = TRUE) %>%
    as.data.frame() %>%
    transmute(s_mu = exp(fit),
              s_lwr = exp(fit - 1.96 * se.fit),
              s_upr = exp(fit + 1.96 * se.fit))
  
  # add columns of distance travelled
  bind_cols(newd, preds_1, preds_2) %>%
    mutate(d_mu = p_mu * s_mu,
           d_lwr = p_lwr * s_lwr,
           d_upr = p_upr * s_upr) %>%
    group_by(species) %>%
    mutate(across(c(s_mu, s_lwr, s_upr), \(.col) .col / mean(s_mu)),
           across(c(d_mu, d_lwr, d_upr), \(.col) .col / mean(d_mu))) %>%
    ungroup() %>%
    mutate(lab = COMMON_NAMES[map_int(species,
                                      \(.s) which(SPECIES == .s))]) %>%
    return()
}

# interaction term plots ----
# color palettes
d_pal <- colorRampPalette(c('#9A2600', '#EBE8DB', '#00605C'))(1e3)

#' exclude values further from an observation than `DIST * 100%` of the
#' range of the observed data
DIST <- 0.1
DOY <- yday('2025-10-06')

# make more appropriate scale breaks
z_breaks <- seq(-1, 1, length.out = 5)
tod_breaks <- seq(0, 24, by = 6)
tod_labs <- paste0(c('0', '0', '', '', ''), tod_breaks, ':00')
doys <- as.Date(paste0('2025-', c('03', '06', '09', '12'), '-01'))
doy_breaks <- yday(doys)
doy_labs <- format(doys, format = '%b 1'); doy_labs

doy_breaks
# function to calculate predictions
surface <- function(newd, term) {
  newd <- newd %>%
    mutate(animal = m_2$model$animal[1],
           lab = case_when(species == 'Rangifer tarandus (s. mountain)' ~ 'bolditalic(Rangifer~tarandus)~bold("(s. mountain)")',
                           species == 'Ursus arctos horribilis' ~ 'bolditalic(Ursus~arctos~horribilis)',
                           species == 'Puma concolor' ~ 'bolditalic(Puma~concolor)',
                           species == 'Cervus canadensis' ~ 'bolditalic(Cervus~canadensis)',
                           species == 'Rangifer tarandus (boreal)' ~ 'bolditalic(Rangifer~tarandus)~bold("(boreal)")',
                           species == 'Canis lupus' ~ 'bolditalic(Canis~lupus)',
                           species == 'Oreamnos americanus' ~ 'bolditalic(Oreamnos~americanus)'))
  
  preds_1 <-
    predict(object = m_1, newdata = newd, type = 'link', se.fit = FALSE,
            exclude = 's(animal)', discrete = TRUE) %>%
    as.data.frame() %>%
    transmute(p_mu = m_1$family$linkinv(`.`))
  
  preds_2 <-
    predict(object = m_2, newdata = newd, type = 'link', se.fit = FALSE,
            exclude = 's(animal)', discrete = TRUE) %>%
    as.data.frame() %>%
    rename(fit = '.') %>%
    transmute(s_mu = exp(fit))
  
  # add columns of distance travelled
  bind_cols(newd, preds_1, preds_2) %>%
    mutate(d_mu = p_mu * s_mu) %>%
    group_by(species) %>%
    mutate(s_mu = s_mu / mean(s_mu),
           d_mu = d_mu / mean(d_mu),
           lab = COMMON_NAMES[map_int(species,
                                      \(.s) which(SPECIES == .s))]) %>%
    ungroup() %>%
    return()
}

# create datasets for prediction ----
newd_tod <- expand_grid(animal = m_2$model$animal[1],
                        species = SPECIES,
                        tod_pdt = seq(0, 24, length.out = 200),
                        doy = DOY,
                        temp_c = seq(-40, 40, length.out = 200),
                        dt = 1)

newd_doy <- expand_grid(animal = m_2$model$animal[1],
                        species = SPECIES,
                        tod_pdt = 12,
                        doy = seq(1, 365, length.out = 200),
                        temp_c = seq(-40, 40, length.out = 200),
                        dt = 1) %>%
  nest(dat = -species) %>%
  mutate(dat = map2(dat, species, \(.d, .s) {
    ref <- filter(m_1$model, species == .s)
    
    filter(.d, ! too_far(x = doy, y = temp_c, ref_1 = ref$doy,
                         ref_2 = ref$temp_c, dist = DIST)) %>%
      return()
  })) %>%
  unnest(dat)

# create predictions ----
ti_tod <- surface(newd_tod, term = 'tod_pdt')
ti_doy <- surface(newd_doy, term = 'doy')

# make caribou figures ----
p_tod_c <-
  ti_tod %>%
  filter(grepl('Caribou', lab)) %>%
  mutate(d_mu = case_when(d_mu < 0.25 ~ 0.25,
                          d_mu > 4 ~ 4,
                          TRUE ~ d_mu)) %>%
  ggplot(aes(temp_c, tod_pdt, fill = log2(d_mu))) +
  geom_raster() +
  geom_contour(aes(temp_c, tod_pdt, z = log2(d_mu)), color = 'black',
               inherit.aes = FALSE, bins = 5) +
  scale_x_continuous(paste0('Temperature (\U00B0', 'C)'),
                     breaks = c(-20, 0, 20)) +
  scale_y_continuous('Time of day (PDT)', expand = c(0, 0),
                     breaks = tod_breaks, labels = tod_labs) +
  scale_fill_gradientn(name = 'Relative change in distance travelled',
                       colors = d_pal, limits = range(z_breaks) * 2,
                       breaks = z_breaks * 2,
                       labels = \(x) 2^x) +
  theme(panel.background = element_rect(fill = 'grey90'),
        legend.position = 'none', legend.key.width = rel(1.5),
        legend.justification = 'center', legend.direction = 'horizontal')

# distance travelled
p_doy_c <-
  ti_doy %>%
  filter(grepl('Caribou', lab)) %>%
  mutate(d_mu = case_when(d_mu < 0.25 ~ 0.25,
                          d_mu > 4 ~ 4,
                          TRUE ~ d_mu)) %>%
  ggplot(aes(temp_c, doy, fill = log2(d_mu))) +
  geom_raster() +
  geom_contour(aes(temp_c, doy, z = log2(d_mu)), color = 'black',
               inherit.aes = FALSE, bins = 5) +
  scale_x_continuous(paste0('Temperature (\U00B0', 'C)'),
                     breaks = c(-20, 0, 20)) +
  scale_y_continuous('Day of year', expand = c(0, 0),
                     breaks = doy_breaks, labels = doy_labs) +
  scale_fill_gradientn(name = 'Relative change in distance travelled',
                       colors = d_pal, limits = range(z_breaks) * 2,
                       breaks = z_breaks * 2,
                       labels = \(x) round(2^x, 2)) +
  theme(panel.background = element_rect(fill = 'grey90'),
        legend.position = 'none', legend.key.width = rel(1.5),
        legend.justification = 'center', legend.direction = 'horizontal')

plot_grid(
  get_legend(p_doy_c),
  plot_grid(p_tod_c, p_doy_c, labels = c('A', 'B'), nrow = 1),
  ncol = 1, rel_heights = c(0.2, 1))

ggsave('figures/2025-tws-edmonton/distance-caribou.png',
       width = 8, height = 4, dpi = 600, bg = 'white')

# caribou and wolves ----
p_tod_cw <-
  ti_tod %>%
  mutate(d_mu = case_when(d_mu < 0.25 ~ 0.25,
                          d_mu > 4 ~ 4,
                          TRUE ~ d_mu)) %>%
  ggplot(aes(temp_c, tod_pdt, fill = log2(d_mu))) +
  facet_wrap(~ lab, nrow = 1) +
  geom_raster() +
  geom_contour(aes(temp_c, tod_pdt, z = log2(d_mu)), color = 'black',
               inherit.aes = FALSE, bins = 5) +
  scale_x_continuous(paste0('Temperature (\U00B0', 'C)'),
                     breaks = c(-20, 0, 20)) +
  scale_y_continuous('Time of day (PDT)', expand = c(0, 0),
                     breaks = tod_breaks, labels = tod_labs) +
  scale_fill_gradientn(name = 'Relative change in distance travelled',
                       colors = d_pal, limits = range(z_breaks) * 2,
                       breaks = z_breaks * 2,
                       labels = \(x) 2^x) +
  theme(panel.background = element_rect(fill = 'grey90'),
        legend.position = 'none', legend.key.width = rel(1.5),
        legend.justification = 'center', legend.direction = 'horizontal')

p_doy_cw <-
  ti_doy %>%
  mutate(d_mu = case_when(d_mu < 0.25 ~ 0.25,
                          d_mu > 4 ~ 4,
                          TRUE ~ d_mu)) %>%
  ggplot(aes(temp_c, doy, fill = log2(d_mu))) +
  facet_wrap(~ lab, nrow = 1) +
  geom_raster() +
  geom_contour(aes(temp_c, doy, z = log2(d_mu)), color = 'black',
               inherit.aes = FALSE, bins = 5) +
  scale_x_continuous(paste0('Temperature (\U00B0', 'C)'),
                     breaks = c(-20, 0, 20)) +
  scale_y_continuous('Day of year', expand = c(0, 0),
                     breaks = doy_breaks, labels = doy_labs) +
  scale_fill_gradientn(name = 'Relative change in distance travelled',
                       colors = d_pal, limits = range(z_breaks) * 2,
                       breaks = z_breaks * 2,
                       labels = \(x) round(2^x, 2)) +
  theme(panel.background = element_rect(fill = 'grey90'),
        legend.position = 'none', legend.key.width = rel(1.5),
        legend.justification = 'center', legend.direction = 'horizontal')

plot_grid(
  get_legend(p_doy_cw),
  plot_grid(p_tod_cw, p_doy_cw, labels = c('A', 'B'), nrow = 1),
  ncol = 1, rel_heights = c(0.2, 1))

ggsave('figures/2025-tws-edmonton/distance-caribou-wolves.png',
       width = 10, height = 4, dpi = 600, bg = 'white')
