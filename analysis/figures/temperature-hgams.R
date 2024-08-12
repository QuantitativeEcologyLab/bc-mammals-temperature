library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for smother date wrangling
library('mgcv')      # for Generalized Additive Models
library('gratia')    # for useful convenience functions for GAMs
library('ggplot2')   # for fancy plots
library('khroma')    # for colorblind-friendly color palettes
library('cowplot')   # for fancy multi-panel plots
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids
plot_scheme(PAL, colours = TRUE)

# import data
d <- readRDS('data/hgam-speed-data.rds')

# import models
m_1 <- readRDS('models/binomial-gam.rds')
m_2 <- readRDS('models/gamma-gam.rds')

# backgrounds for plots ----
tod <- tibble(x = c(seq(0, 6, by = 0.01),
                    seq(18, 24, by = 0.01)),
              g = c(rep(1, 601), rep(2, 601)))

seasons <- tibble(x = seq(1, 366, by = 0.1),
                  y = 0.49,
                  season = case_when(x < 80 ~ 'w2',
                                     x < 172 ~ 'sp',
                                     x < 266 ~ 'su',
                                     x < 356 ~ 'f',
                                     TRUE ~ 'w1') %>%
                    factor(levels = c('sp', 'su', 'f', 'w1', 'w2')))

PAL_SEASONS <- c('forestgreen', 'goldenrod', 'darkorange4', 'grey90', 'grey90')

season_breaks <-
  as.Date(paste0('2024-', c('03-20', '06-20', '09-22', '12-21'))) + 45

# marginal term plots ----
#' not averaging across uniform `tod_pdt` and `doy` because I want to
#' preserve the uneaven sampling of `doy` for `E(speed | moving)`
marginal <- function(newd, term) {
  preds_1 <-
    predict(object = m_1, newdata = newd, type = 'link', se.fit = TRUE,
            terms = c('(Intercept)', 'species',
                      paste0('s(', term, '):species', SPECIES)),
            discrete = FALSE) %>%
    as.data.frame() %>%
    transmute(p_mu = m_1$family$linkinv(fit),
              p_lwr = m_1$family$linkinv(fit - 1.96 * se.fit),
              p_upr = m_1$family$linkinv(fit + 1.96 * se.fit))
  
  preds_2 <-
    predict(object = m_2, newdata = newd, type = 'link', se.fit = TRUE,
            terms = paste0('s(', term, '):species', SPECIES),
            discrete = FALSE) %>%
    as.data.frame() %>%
    transmute(s_mu = exp(fit),
              s_lwr = exp(fit - 1.96 * se.fit),
              s_upr = exp(fit + 1.96 * se.fit))
  
  # add columns of distance travelled
  bind_cols(newd, preds_1, preds_2) %>%
    mutate(species = gsub(' ', '~', species),
           species = gsub('~\\(', '\\)~bold\\((', species),
           species = paste0('bolditalic(', species, ')'),
           d_mu = p_mu * s_mu,
           d_lwr = p_lwr * s_lwr,
           d_upr = p_upr * s_upr) %>%
    group_by(species) %>%
    mutate(across(c(s_mu, s_lwr, s_upr), \(.col) .col / mean(s_mu)),
           across(c(d_mu, d_lwr, d_upr), \(.col) .col / mean(d_mu))) %>%
    ungroup() %>%
    return()
}

# s(tod_pdt) ----
d %>%
  group_by(species) %>%
  summarize(min_tod = round(min(tod_pdt), 2),
            max_tod = round(max(tod_pdt), 2))

newd_tod <-
  expand_grid(animal = 'new animal',
              species = SPECIES,
              tod_pdt = seq(0, 24, length.out = 1e3),
              doy = 0,
              temp_c = 0,
              dt = 0)

preds_tod <- marginal(newd = newd_tod, term = 'tod_pdt')

p_mov_tod <-
  preds_tod %>%
  mutate(p_upr = if_else(p_upr > 0.75, 0.75, p_upr)) %>%
  ggplot() +
  coord_polar() +
  facet_wrap(~ species, nrow = 1, labeller = label_parsed) +
  geom_area(aes(x, 0.75, group = g), tod, fill = 'black', alpha = 0.2) +
  geom_line(aes(tod_pdt, p_mu, color = species), linewidth = 1) +
  geom_ribbon(aes(tod_pdt, ymin = p_lwr, ymax = p_upr, fill = species),
              alpha = 0.2) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 24),
                     breaks = c(0, 6, 12, 18),
                     labels = c('00:00', '06:00', '12:00', '18:00')) +
  scale_y_continuous(limits = c(0, 0.75), expand = c(0, 0),
                     breaks = c(0, 0.25, 0.5, 0.75)) +
  labs(x = 'Time of day (PDT)', y = 'P(moving)') +
  theme(legend.position = 'none',
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = c(rep('grey', 4), NA)))

speed_tod <-
  ggplot(preds_tod) +
  coord_polar() +
  facet_wrap(~ species, nrow = 1, labeller = label_parsed) +
  geom_area(aes(x, 1.5, group = g), tod, fill = 'black', alpha = 0.2) +
  geom_line(aes(tod_pdt, s_mu, color = species), linewidth = 1) +
  geom_ribbon(aes(tod_pdt, ymin = s_lwr, ymax = s_upr, fill = species),
              alpha = 0.2) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 24),
                     breaks = c(0, 6, 12, 18),
                     labels = c('00:00', '06:00', '12:00', '18:00')) +
  scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
  labs(x = 'Time of day (PDT)', y = 'Relative change in speed') +
  theme(legend.position = 'none',
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = c(rep('grey', 4), NA)))

distance_tod <-
  preds_tod %>%
  mutate(d_upr = if_else(d_upr > 4, 4, d_upr)) %>%
  ggplot() +
  coord_polar() +
  facet_wrap(~ species, nrow = 1, labeller = label_parsed) +
  geom_area(aes(x, 4, group = g), tod, fill = 'black', alpha = 0.2) +
  geom_line(aes(tod_pdt, d_mu, color = species), linewidth = 1) +
  geom_ribbon(aes(tod_pdt, ymin = d_lwr, ymax = d_upr, fill = species),
              alpha = 0.2) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 24),
                     breaks = c(0, 6, 12, 18),
                     labels = c('00:00', '06:00', '12:00', '18:00')) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  labs(x = 'Time of day (PDT)', y = 'Relative change in distance travelled') +
  theme(legend.position = 'none',
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = c(rep('grey', 5), NA)))

# s(doy) ----
newd_doy <- expand_grid(animal = 'new animal',
                        species = SPECIES,
                        tod_pdt = 12,
                        doy = seq(1, 366, by = 0.1),
                        temp_c = 0,
                        dt = 0)

preds_doy <- marginal(newd = newd_doy, term = 'doy')

p_mov_doy <-
  ggplot(preds_doy) +
  coord_polar(start = 11 / 366 * 2 * pi) + # offset to make axes vertical
  facet_wrap(~ species, nrow = 1, labeller = label_parsed) +
  geom_area(aes(x, 1, fill = season), seasons, alpha = 0.3,
            show.legend = FALSE) +
  scale_fill_manual(values = PAL_SEASONS) +
  ggnewscale::new_scale('fill') + # to remove the season scale for fill
  geom_ribbon(aes(doy, ymin = p_lwr, ymax = p_upr, fill = species),
              alpha = 0.2) +
  geom_line(aes(doy, p_mu, color = species), linewidth = 1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(breaks = yday(season_breaks),
                     labels = c('Spring', 'Summer', 'Fall', 'Winter')) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = 'Day of year', y = 'P(moving)') +
  theme(legend.position = 'none',
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = c(rep('grey', 5), NA)))

speed_doy <-
  ggplot(preds_doy) +
  coord_polar(start = 11 / 366 * 2 * pi) + # offset to make axes vertical
  facet_wrap(~ species, nrow = 1, labeller = label_parsed) +
  geom_area(aes(x, 2, fill = season), seasons, alpha = 0.3,
            show.legend = FALSE) +
  scale_fill_manual(values = PAL_SEASONS) +
  ggnewscale::new_scale('fill') + # to remove the season scale for fill
  geom_ribbon(aes(doy, ymin = s_lwr, ymax = s_upr, fill = species),
              alpha = 0.2) +
  geom_line(aes(doy, s_mu, color = species), linewidth = 1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(breaks = yday(season_breaks),
                     labels = c('Spring', 'Summer', 'Fall', 'Winter')) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
  labs(x = 'Day of year', y = 'Relative change in speed') +
  theme(legend.position = 'none',
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = c(rep('grey', 5), NA)))

distance_doy <-
  preds_doy %>%
  mutate(d_upr = if_else(d_upr > 3, 3, d_upr)) %>%
  ggplot() +
  coord_polar(start = 11 / 366 * 2 * pi) + # offset to make axes vertical
  facet_wrap(~ species, nrow = 1, labeller = label_parsed) +
  geom_area(aes(x, 3, fill = season), seasons, alpha = 0.3,
            show.legend = FALSE) +
  scale_fill_manual(values = PAL_SEASONS) +
  ggnewscale::new_scale('fill') + # to remove the season scale for fill
  geom_ribbon(aes(doy, ymin = d_lwr, ymax = d_upr, fill = species),
              alpha = 0.2) +
  geom_line(aes(doy, d_mu, color = species), linewidth = 1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(breaks = yday(season_breaks),
                     labels = c('Spring', 'Summer', 'Fall', 'Winter')) +
  scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
  labs(x = 'Day of year', y = 'Relative change in distance travelled') +
  theme(legend.position = 'none',
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = c(rep('grey', 4), NA)))

# s(temp_c) ----
newd_temp_c <- tibble(animal = 'new animal',
                      species = SPECIES,
                      tod_pdt = 12,
                      doy = 0,
                      dt = 0) %>%
  mutate(temp_data = purrr::map(species, \(.s) {
    .temp <- filter(d, species == .s)$temp_c
    
    tibble(temp_c = seq(quantile(.temp, 0, na.rm = TRUE),
                        quantile(.temp, 1, na.rm = TRUE),
                        length.out = 400))
  })) %>%
  unnest(temp_data)

preds_temp_c <- marginal(newd = newd_temp_c, term = 'temp_c')

# change species names for rug plots
d <- mutate(d,
            species = gsub(' ', '~', species),
            species = gsub('~\\(', '\\)~bold\\((', species),
            species = paste0('bolditalic(', species, ')'))

p_mov_temp <-
  ggplot(preds_temp_c) +
  facet_wrap(~ species, nrow = 1, labeller = label_parsed,
             scales = 'free') +
  geom_rug(aes(temp_c), d, alpha = 0.1) +
  geom_ribbon(aes(temp_c, ymin = p_lwr, ymax = p_upr, fill = species),
              alpha = .2) +
  geom_line(aes(temp_c, p_mu, color = species), linewidth =  1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)'))+
  scale_y_continuous('P(moving)') +
  theme(legend.position = 'none')

speed_temp <-
  ggplot(preds_temp_c) +
  facet_wrap(~ species, nrow = 1, labeller = label_parsed) +
  geom_rug(aes(temp_c), filter(d, moving), alpha = 0.1) +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(temp_c, ymin = s_lwr, ymax = s_upr, fill = species),
              alpha = .2) +
  geom_line(aes(temp_c, s_mu, color = species), linewidth =  1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)'))+
  scale_y_continuous('Relative change in speed') +
  theme(legend.position = 'none')

distance_temp <-
  ggplot(preds_temp_c) +
  facet_wrap(~ species, nrow = 1, labeller = label_parsed,
             scales = 'free_y') +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(temp_c, ymin = d_lwr, ymax = d_upr, fill = species),
              alpha = .2) +
  geom_line(aes(temp_c, d_mu, color = species), linewidth =  1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)'))+
  scale_y_continuous('Relative change in distance travelled') +
  theme(legend.position = 'none')

# full figures ----
p_mov_full <- plot_grid(p_mov_temp, p_mov_tod, p_mov_doy,
                        labels = 'AUTO', ncol = 1)
ggsave('figures/p-moving-marginals.png', p_mov_full,
       width = 25, height = 13, dpi = 600, bg = 'white')

speed_full <- plot_grid(speed_temp, speed_tod, speed_doy,
                        labels = 'AUTO', ncol = 1)
ggsave('figures/speed-marginals.png', speed_full,
       width = 25, height = 13, dpi = 600, bg = 'white')

distance_full <- plot_grid(distance_temp, distance_tod, distance_doy,
                           labels = 'AUTO', ncol = 1)
ggsave('figures/distance-travelled-marginals.png', distance_full,
       width = 25, height = 13, dpi = 600, bg = 'white')

temp_full <- plot_grid(p_mov_temp, speed_temp, distance_temp,
                       labels = 'AUTO', ncol = 1)
ggsave('figures/temp_c-all.png', temp_full,
       width = 25, height = 13, dpi = 600, bg = 'white')

# interaction term plots ----
# color palettes
s_pal <- colorRampPalette(c('#BD4301', '#EBE8DB', '#560A63'))(1e3)
plot_scheme_colorblind(s_pal)
d_pal <- colorRampPalette(c('#9A2600', '#EBE8DB', '#00605C'))(1e3)
plot_scheme_colorblind(d_pal)
plot_scheme_colorblind(c(color('acton')(1e3), s_pal, d_pal))

#' exclude values further from an observation than `DIST * 100%` of the
#' range of the observed data
DIST <- 0.1

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
  preds_1 <-
    predict(object = m_1, newdata = newd, type = 'link',
            se.fit = FALSE,
            terms = c('(Intercept)', 's(species)',
                      paste0('s(', term, '):species', SPECIES),
                      paste0('ti(temp_c,', term, '):species', SPECIES)),
            discrete = FALSE) %>%
    as.data.frame() %>%
    transmute(p_mu = m_1$family$linkinv(`.`))
  
  preds_2 <-
    predict(object = m_2, newdata = newd, type = 'link', se.fit = FALSE,
            terms = c(paste0('s(', term, '):species', SPECIES),
                      paste0('ti(temp_c,', term, '):species', SPECIES)),
            discrete = FALSE) %>%
    as.data.frame() %>%
    rename(fit = '.') %>%
    transmute(s_mu = exp(fit))
  
  # add columns of distance travelled
  bind_cols(newd, preds_1, preds_2) %>%
    mutate(species = gsub(' ', '~', species),
           species = gsub('~\\(', '\\)~bold\\((', species),
           species = paste0('bolditalic(', species, ')'),
           d_mu = p_mu * s_mu) %>%
    group_by(species) %>%
    mutate(s_mu = s_mu / mean(s_mu),
           d_mu = d_mu / mean(d_mu)) %>%
    ungroup() %>%
    return()
}

# time of day ----
# P(moving)
newd_tod <- expand_grid(animal = 'new animal',
                        species = SPECIES,
                        tod_pdt = seq(0, 24, length.out = 200),
                        doy = 0,
                        temp_c = seq(-40, 40, length.out = 200),
                        dt = 0) %>%
  nest(dat = -species) %>%
  mutate(dat = map2(dat, species, \(.d, .s) {
    ref <- filter(m_1$model, species == .s)
    
    filter(.d, ! too_far(x = tod_pdt, y = temp_c, ref_1 = ref$tod_pdt,
                         ref_2 = ref$temp_c, dist = DIST)) %>%
      return()
  })) %>%
  unnest(dat)

ti_tod <- surface(newd_tod, term = 'tod_pdt')

p_mov_tod_int <-
  mutate(ti_tod,
         p_mu = if_else(p_mu > 0.4, 0.4, p_mu)) %>%
  ggplot(aes(temp_c, tod_pdt, fill = p_mu)) +
  facet_wrap(~ species, labeller = label_parsed, nrow = 3) +
  geom_raster() +
  geom_contour(aes(temp_c, tod_pdt, z = log2(p_mu)), color = 'grey',
               inherit.aes = FALSE) +
  scale_x_continuous(paste0('Temperature (\U00B0', 'C)'),
                     breaks = c(-20, 0, 20)) +
  scale_y_continuous('Time of day (PDT)', expand = c(0, 0),
                     breaks = tod_breaks, labels = tod_labs) +
  scale_fill_acton(name = 'P(moving)', limits = c(0, 0.3),
                   breaks = seq(0, 0.3, by = 0.1)) +
  theme(panel.background = element_rect(fill = 'grey'),
        legend.position = 'inside', legend.position.inside = c(0.66, 0.15),
        legend.key.width = rel(1.5),  legend.justification = 'center',
        legend.direction = 'horizontal')

ggsave('figures/p-moving-tod-interaction.png', p_mov_tod_int,
       width = 8, height = 8, units = 'in', dpi = 600, bg = 'white')

# E(speed | moving)
s_tod_int <-
  ggplot(ti_tod, aes(temp_c, tod_pdt, fill = log2(s_mu))) +
  facet_wrap(~ species, labeller = label_parsed, nrow = 3) +
  geom_raster() +
  geom_contour(aes(temp_c, tod_pdt, z = log2(s_mu)), color = 'black',
               inherit.aes = FALSE, bins = 5) +
  scale_x_continuous(paste0('Temperature (\U00B0', 'C)'),
                     breaks = c(-20, 0, 20)) +
  scale_y_continuous('Time of day (PDT)', expand = c(0, 0),
                     breaks = tod_breaks, labels = tod_labs) +
  scale_fill_gradientn(name = 'Relative change in speed', colors = s_pal,
                       limits = range(z_breaks), breaks = z_breaks,
                       labels = round(2^z_breaks, 2)) +
  theme(panel.background = element_rect(fill = 'grey'),
        legend.position = 'inside', legend.position.inside = c(0.66, 0.15),
        legend.key.width = rel(1.5),  legend.justification = 'center',
        legend.direction = 'horizontal')

ggsave('figures/speed-tod-interaction.png', s_tod_int,
       width = 8, height = 8, units = 'in', dpi = 600, bg = 'white')

# E(distance)
d_tod_int <-
  ti_tod %>%
  mutate(d_mu = case_when(d_mu < 0.25 ~ 0.25,
                          d_mu > 4 ~ 4,
                          TRUE ~ d_mu)) %>%
  ggplot(aes(temp_c, tod_pdt, fill = log2(d_mu))) +
  facet_wrap(~ species, labeller = label_parsed, nrow = 3) +
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
                       labels = round(2^(z_breaks * 2), 2)) +
  theme(panel.background = element_rect(fill = 'grey'),
        legend.position = 'inside', legend.position.inside = c(0.66, 0.15),
        legend.key.width = rel(1.5),  legend.justification = 'center',
        legend.direction = 'horizontal')

ggsave('figures/distance-travelled-tod-interaction.png', d_tod_int,
       width = 8, height = 8, units = 'in', dpi = 600, bg = 'white')

# doy ----
newd_doy <- expand_grid(animal = 'new animal',
                        species = SPECIES,
                        tod_pdt = 0,
                        doy = seq(1, 365, length.out = 200),
                        temp_c = seq(-40, 40, length.out = 200),
                        dt = 0) %>%
  nest(dat = -species) %>%
  mutate(dat = map2(dat, species, \(.d, .s) {
    ref <- filter(m_1$model, species == .s)
    
    filter(.d, ! too_far(x = doy, y = temp_c, ref_1 = ref$doy,
                         ref_2 = ref$temp_c, dist = DIST)) %>%
      return()
  })) %>%
  unnest(dat)

ti_doy <- surface(newd_doy, term = 'doy')

p_mov_doy_int <-
  mutate(ti_doy,
         p_mu = if_else(p_mu > 0.4, 0.4, p_mu)) %>%
  ggplot(aes(temp_c, doy, fill = p_mu)) +
  facet_wrap(~ species, labeller = label_parsed, nrow = 3) +
  geom_raster() +
  geom_contour(aes(temp_c, doy, z = log2(p_mu)), color = 'grey',
               inherit.aes = FALSE, bins = 5) +
  scale_x_continuous(paste0('Temperature (\U00B0', 'C)'),
                     breaks = c(-20, 0, 20)) +
  scale_y_continuous('Day of year', expand = c(0, 0),
                     breaks = doy_breaks, labels = doy_labs) +
  scale_fill_acton(name = 'P(moving)', limits = c(0, NA)) +
  theme(panel.background = element_rect(fill = 'grey'),
        legend.position = 'inside', legend.position.inside = c(0.66, 0.15),
        legend.key.width = rel(1.5),  legend.justification = 'center',
        legend.direction = 'horizontal')

ggsave('figures/p-moving-doy-interaction.png', p_mov_doy_int,
       width = 8, height = 8, units = 'in', dpi = 600, bg = 'white')

# E(speed | moving)
s_doy_int <-
  ggplot(ti_doy, aes(temp_c, doy, fill = log2(s_mu))) +
  facet_wrap(~ species, labeller = label_parsed, nrow = 3) +
  geom_raster() +
  geom_contour(aes(temp_c, doy, z = log2(s_mu)), color = 'black',
               inherit.aes = FALSE, bins = 5) +
  scale_x_continuous(paste0('Temperature (\U00B0', 'C)'),
                     breaks = c(-20, 0, 20)) +
  scale_y_continuous('Day of year', expand = c(0, 0),
                     breaks = doy_breaks, labels = doy_labs) +
  scale_fill_gradientn(name = 'Relative change in speed', colors = s_pal,
                       limits = range(z_breaks), breaks = z_breaks,
                       labels = round(2^z_breaks, 2)) +
  theme(panel.background = element_rect(fill = 'grey'),
        legend.position = 'inside', legend.position.inside = c(0.66, 0.15),
        legend.key.width = rel(1.5),  legend.justification = 'center',
        legend.direction = 'horizontal')

ggsave('figures/speed-doy-interaction.png', s_doy_int,
       width = 8, height = 8, units = 'in', dpi = 600, bg = 'white')

# distance travelled
d_doy_int <-
  ti_doy %>%
  mutate(d_mu = case_when(d_mu < 0.25 ~ 0.25,
                          d_mu > 4 ~ 4,
                          TRUE ~ d_mu)) %>%
  ggplot(aes(temp_c, doy, fill = log2(d_mu))) +
  facet_wrap(~ species, labeller = label_parsed, nrow = 3) +
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
                       labels = round(2^(z_breaks * 2), 2)) +
  theme(panel.background = element_rect(fill = 'grey'),
        legend.position = 'inside', legend.position.inside = c(0.66, 0.15),
        legend.key.width = rel(1.5),  legend.justification = 'center',
        legend.direction = 'horizontal')

ggsave('figures/distance-travelled-doy-interaction.png', d_doy_int,
       width = 8, height = 8, units = 'in', dpi = 600, bg = 'white')

# figures by response parameter ----
#' `get_legend()` for `{cowplot}` version 1.1.3 returns an empty plot 
get_legend <- function(.plot) {
  get_plot_component(.plot + theme(legend.position = 'top'),
                     pattern = 'guide-box-top', return_all = TRUE)
}

p_mov <- plot_grid(
  get_legend(p_mov_tod_int),
  p_mov_tod_int +
    facet_wrap(~ species, labeller = label_parsed, nrow = 1) +
    theme(legend.position = 'none'),
  p_mov_doy_int +
    facet_wrap(~ species, labeller = label_parsed, nrow = 1) +
    theme(legend.position = 'none'),
  labels = c('', 'A', 'B'), ncol = 1, rel_heights = c(0.2, 1, 1))
ggsave('figures/p-moving.png', p_mov,
       width = 20, height = 10, dpi = 600, bg = 'white')

speed <- plot_grid(
  get_legend(s_tod_int),
  s_tod_int +
    facet_wrap(~ species, labeller = label_parsed, nrow = 1) +
    theme(legend.position = 'none'),
  s_doy_int +
    facet_wrap(~ species, labeller = label_parsed, nrow = 1) +
    theme(legend.position = 'none'),
  labels = c('', 'A', 'B'), ncol = 1, rel_heights = c(0.2, 1, 1))
ggsave('figures/speed.png', speed,
       width = 20, height = 10, dpi = 600, bg = 'white')

distance <- plot_grid(
  get_legend(d_tod_int),
  d_tod_int +
    facet_wrap(~ species, labeller = label_parsed, nrow = 1) +
    theme(legend.position = 'none'),
  d_doy_int +
    facet_wrap(~ species, labeller = label_parsed, nrow = 1) +
    theme(legend.position = 'none'),
  labels = c('', 'A', 'B'), ncol = 1, rel_heights = c(0.2, 1, 1))
ggsave('figures/distance.png', distance,
       width = 20, height = 10, dpi = 600, bg = 'white')
