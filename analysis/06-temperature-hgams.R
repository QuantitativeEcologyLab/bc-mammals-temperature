library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for smother date wrangling
library('mgcv')      # for Generalized Additive Models
library('gratia')    # for ggplot-based figures for GAMs
library('ggplot2')   # for fancy plots
library('khroma')    # for colorblind-friendly color palettes
library('cowplot')   # for fancy multi-panel plots
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids
plot_scheme(PAL, colours = TRUE)

# using Pacific Time for time of day
with_tz(as.POSIXct('2024-01-01 12:00'), 'America/Vancouver')
with_tz(as.POSIXct('2024-01-01 12:00'), 'UTC')

# import data
d <-
  readRDS('data/movement-models-speed-weights-temperature-2024-06-10.rds') %>%
  filter(! is.na(speed_est), is.finite(speed_est)) %>%
  mutate(animal = factor(animal),
         species = if_else(grepl('Rangifer', dataset_name),
                           dataset_name,
                           species),
         species = stringr::str_replace_all(species, '_', ' ') %>%
           stringr::str_replace_all('boreal', '(boreal)') %>%
           stringr::str_replace_all('southern mountain',
                                    '(southern mountain)') %>%
           factor()) %>%
  #' cluster speeds using `kmeans()`
  nest(speeds = ! species) %>%
  mutate(speeds = map(speeds, \(.d) {
    km <- kmeans(.d$speed_est, centers = 2)
    
    moving_id <- which.max(km$centers)
    
    .d <- mutate(.d,
                 moving = km$cluster == moving_id)
    
    return(.d)
  })) %>%
  unnest(speeds) %>%
  mutate(tod_pdt = (hour(timestamp) +
                      minute(timestamp) / 60 +
                      second(timestamp) / 60 / 60),
         # UTC to PDT
         tod_pdt = if_else(tod_pdt < 8, tod_pdt + 16, tod_pdt - 8), 
         doy = yday(timestamp)) %>% # using UTC doy
  rename(temp_c = temperature_C)

SPECIES <- unique(d$species)
N_SPECIES <- length(SPECIES)

# check speed classification for each species
d %>%
  group_by(species) %>%
  summarise(prop_moving = mean(moving),
            min_speed = min(speed_est[moving]))

ggplot(d) +
  facet_wrap(~ species + moving, scales = 'free', ncol = 2) +
  geom_density(aes(speed_est), fill = 'grey')

# drop single data point in February for bears
plot(as.numeric(moving) ~ doy,
     filter(d, species == 'Ursus arctos horribilis'))
d <- filter(d, ! (species == 'Ursus arctos horribilis' & doy < 50))
points(as.numeric(moving) ~ doy,
       filter(d, species == 'Ursus arctos horribilis'), col = 'red')

# some species have data with extreme temperatures but little data ----
# but clipping temperature data to the 0.8 quantile doesn't change the high
# P(moving) at cold temperatures for goats
ggplot(d, aes(temp_c)) +
  facet_wrap(~ species, scales = 'free') +
  geom_histogram(color = 'black', binwidth = 1, fill = 'grey') +
  geom_vline(aes(xintercept = value), color = 'red',
             d %>%
               group_by(species) %>%
               summarize(lwr = quantile(temp_c, 0.025),
                         upr = quantile(temp_c, 0.975)) %>%
               pivot_longer(-species)) +
  ylim(c(0, NA)) +
  xlab(paste0('Temperature (\U00B0', 'C)'))

# model P(movement) ----
#' using `bs = 'fs'` rather than `by` smooths because to keep the models
#' reasonably constrained. Using `by` gives extremely high `P(moving)`
#' in spring for grizzly bears
#' fits in 6 minutes
if(file.exists('models/binomial-gam-2024-06-11rds')) {
  m_1 <- readRDS('models/binomial-gam-2024-06-11.rds')
} else {
  #' `ti(doy, temp_c)` does not increase BIC or AIC
  m_1 <-
    bam(moving ~
          # random effect for each animal
          s(animal, bs = 're') +
          # to account for changes in behavior within days
          s(tod_pdt, species, k = 10, bs = 'fs', xt = list(bs = 'cc')) +
          # to account for changes in behavior within years
          s(doy, species, k = 10, bs = 'fs', xt = list(bs = 'cc')) +
          # species-level effect of temperature
          s(temp_c, species, k = 5, bs = 'fs', xt = list(bs = 'tp')) +
          # to account for seasonal changes in day length
          ti(doy, tod_pdt, species, k = 5, bs = c('cc', 'cc', 're')) +
          # to account for changes in day nocturnality with temperature
          ti(temp_c, tod_pdt, species, k = 5, bs = c('cc', 'cc', 're')),
        family = binomial(link = 'logit'),
        data = d,
        method = 'fREML', # fast REML
        discrete = TRUE, # discretize the posterior for faster computation
        knots = list(tod_pdt = c(0, 1), doy = c(0.5, 366.5)),# for bs = 'cc'
        control = gam.control(trace = TRUE))
  saveRDS(m_1, paste0('models/binomial-gam-', Sys.Date(), '.rds'))
  
  # southern mountain and boreal caribou move differently
  draw(m_1, contour = FALSE, rug = FALSE,
       discrete_colour = scale_color_manual(values = PAL))
  
  summary(m_1)
  
  qq.gam(m_1, type = 'pearson') # good q-q plot
}

# check predictions
transmute(d,
          species,
          moving,
          est = predict(m_1, type = 'response') %>%
            round(2)) %>%
  group_by(species, est) %>%
  summarise(empirical_mean = mean(moving),
            n = n(),
            .groups = 'drop') %>%
  ggplot(aes(est, empirical_mean, color = log2(n))) +
  facet_wrap(~ species) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_point() +
  lims(x = c(0, 1), y = c(0, 1)) +
  scale_color_viridis_c(breaks = 4 * (0:4),
                        labels = 2^(2^(0:4)))

# s(tod_pdt) ----
d %>%
  group_by(species) %>%
  summarize(min_tod = round(min(tod_pdt), 2),
            max_tod = round(max(tod_pdt), 2))

newd_tod <- expand_grid(animal = 'new animal',
                        species = SPECIES,
                        tod_pdt = seq(0, 24, length.out = 1e3),
                        doy = 0,
                        temp_c = 0)

p_mov_tod <-
  bind_cols(newd_tod,
            predict(m_1, newdata = newd_tod,
                    terms = 's(tod_pdt,species)',
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_1$family$linkinv(fit),
         lwr = m_1$family$linkinv(fit - 1.96 * se.fit),
         upr = m_1$family$linkinv(fit + 1.96 * se.fit)) %>%
  ggplot() +
  coord_polar() +
  facet_wrap(~ species, nrow = 1) +
  geom_area(aes(x, y), tibble(x = seq(0, 6, by = 0.01), y = 1),
            fill = 'black', alpha = 0.3) +
  geom_area(aes(x, y), tibble(x = seq(18, 24, by = 0.01), y = 1),
            fill = 'black', alpha = 0.3) +
  geom_area(aes(tod_pdt, mu, color = species, fill = species), linewidth = 1,
            alpha = 0.5) +
  geom_ribbon(aes(tod_pdt, ymin = lwr, ymax = upr, fill = species),
              alpha = 0.2) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 24),
                     breaks = c(0, 6, 12, 18),
                     labels = c('00:00', '06:00', '12:00', '18:00')) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = 'Time of day', y = 'P(moving)') +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold.italic'),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = c(rep('grey', 5), NA)))

ggsave('figures/p-moving-tod.png', p_mov_tod,
       width = 25, height = 4, dpi = 600, bg = 'white')

# s(doy) ----
seasons <- tibble(x = seq(1, 366, by = 0.1),
                  y = 1,
                  season = case_when(x < 80 ~ 'w2',
                                     x < 172 ~ 'sp',
                                     x < 266 ~ 'su',
                                     x < 356 ~ 'f',
                                     TRUE ~ 'w1') %>%
                    factor(levels = c('sp', 'su', 'f', 'w1', 'w2'))) %>%
  list() %>%
  tibble(data = ., species = SPECIES) %>%
  unnest(data)

PAL_SEASONS <- c('forestgreen', 'goldenrod', 'darkorange4', 'white', 'white')

# ensure smoooth for grizzly bears is reasonable
d %>%
  group_by(species) %>%
  summarize(min_doy = round(min(doy), 2),
            max_doy = round(max(doy), 2))

newd_doy <- expand_grid(animal = 'new animal',
                        species = SPECIES,
                        tod_pdt = 12,
                        doy = seq(1, 366, by = 0.1),
                        temp_c = 0)

p_mov_doy <-
  bind_cols(newd_doy,
            predict(m_1, newdata = newd_doy,
                    terms = 's(doy,species)',
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_1$family$linkinv(fit),
         lwr = m_1$family$linkinv(fit - 1.96 * se.fit),
         upr = m_1$family$linkinv(fit + 1.96 * se.fit)) %>%
  ggplot() +
  coord_polar(start = 11 / 366 * 2 * pi) + # offset to make axes vertical
  facet_wrap(~ species, nrow = 1) +
  geom_area(aes(x, y, fill = season), seasons, alpha = 0.3,
            show.legend = FALSE) +
  scale_fill_manual(values = PAL_SEASONS) +
  ggnewscale::new_scale('fill') + # to remove the season scale for fill
  geom_ribbon(aes(doy, ymin = lwr, ymax = upr, fill = species),
              alpha = 0.2) +
  geom_area(aes(doy, mu, color = species, fill = species), linewidth = 1,
            alpha = 0.5) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(breaks = yday(paste0('2024-', c('03-20', '06-20',
                                                     '09-22', '12-21')))) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = 'Day of year', y = 'P(moving)') +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold.italic'),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_line(color = 'grey'),
        panel.grid.major.y = element_line(colour = c(rep('grey', 5), NA)))

ggsave('figures/p-moving-doy.png', p_mov_doy,
       width = 25, height = 4, dpi = 600, bg = 'white')

# s(temp_c) ----
newd_temp_c <- tibble(animal = 'new animal',
                      species = SPECIES,
                      tod_pdt = 12,
                      doy = 0) %>%
  mutate(temp_data = purrr::map(species, \(.s) {
    .temp <- filter(d, species == .s)$temp_c
    
    tibble(temp_c = seq(quantile(.temp, 0, na.rm = TRUE),
                        quantile(.temp, 1, na.rm = TRUE),
                        length.out = 400))
  })) %>%
  unnest(temp_data)

p_mov_temp <-
  bind_cols(newd_temp_c,
            predict(m_1, newdata = newd_temp_c,
                    terms = 's(temp_c,species)',
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_1$family$linkinv(fit),
         lwr = m_1$family$linkinv(fit - 1.96 * se.fit),
         upr = m_1$family$linkinv(fit + 1.96 * se.fit)) %>%
  ggplot() +
  facet_wrap(~ species, nrow = 1) +
  geom_ribbon(aes(temp_c, ymin = lwr, ymax = upr, fill = species),
              alpha = .2) +
  geom_line(aes(temp_c, mu, color = species), linewidth =  1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)'))+
  scale_y_continuous('P(moving)', expand = c(0, 0), limits = c(0, 1)) +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold.italic'))

ggsave('figures/p-moving-temperature.png', p_mov_temp,
       width = 25, height = 4, dpi = 600, bg = 'white')

# full figure for P(moving) ----
p_mov <- plot_grid(p_mov_temp, p_mov_tod, p_mov_doy, labels = 'AUTO',
                   ncol = 1)

ggsave('figures/p-moving-all.png', p_mov,
       width = 25, height = 12, dpi = 600, bg = 'white')

# model mean speed given that animals are moving ----
d_2 <- filter(d, moving)

if(file.exists('models/gamma-gam-2024-05-03.rds')) {
  m_2 <- readRDS('models/gamma-gam-2024-05-03.rds')
} else {
  # regular HGAM has non-gaussian residuals each species
  m_2 <-
    bam(speed_est ~
          # random effect for each animal
          s(animal, bs = 're') +
          # to account for changes in behavior within days
          s(tod_pdt, species, k = 10, bs = 'fs', xt = list(bs = 'cc')) +
          # to account for changes in behavior within years
          s(doy, species, k = 10, bs = 'fs', xt = list(bs = 'cc')) +
          # species-level effect of temperature
          s(temp_c, species, k = 5, bs = 'fs', xt = list(bs = 'tp')) +
          # to account for seasonal changes in day length
          ti(doy, tod_pdt, species, k = 5, bs = c('cc', 'cc', 're')) +
          # to account for changes in day nocturnality with temperature
          ti(temp_c, tod_pdt, species, k = 5, bs = c('cc', 'cc', 're')),
        family = Gamma(link = 'log'), # can use Gamma because no zeros
        data = d_2,
        method = 'fREML', # fast REML
        discrete = TRUE, # discretize the posterior for faster computation
        knots = list(tod_pdt = c(0, 1), doy = c(0.5, 366.5)),
        control = gam.control(trace = TRUE))
  
  draw(m_2, contour = FALSE, rug = FALSE,
       discrete_colour = scale_color_manual(values = PAL))
  
  summary(m_2)
  
  mutate(d_2, e = resid(m_2)) %>%
    ggplot(aes(e, fill = species)) +
    facet_wrap(~ species, scales = 'free') +
    geom_histogram(show.legend = FALSE) +
    scale_fill_manual('Species', values = PAL)
  
  m_2 <-
    gam(list(
      speed_est ~
          # random effect for each animal
          s(animal, bs = 're') +
          # to account for changes in behavior within days
          s(tod_pdt, species, k = 10, bs = 'fs', xt = list(bs = 'cc')) +
          # to account for changes in behavior within years
          s(doy, species, k = 10, bs = 'fs', xt = list(bs = 'cc')) +
          # species-level effect of temperature
          s(temp_c, species, k = 5, bs = 'fs', xt = list(bs = 'tp')) +
          # to account for seasonal changes in day length
          ti(doy, tod_pdt, species, k = 5, bs = c('cc', 'cc', 're')) +
          # to account for changes in day nocturnality with temperature
          ti(temp_c, tod_pdt, species, k = 5, bs = c('cc', 'cc', 're')),
      ~
        # to account for changes in behavior within days
        s(tod_pdt, species, k = 10, bs = 'fs', xt = list(bs = 'cc')) +
        # to account for changes in behavior within years
        s(doy, species, k = 10, bs = 'fs', xt = list(bs = 'cc')) +
        # species-level effect of temperature
        s(temp_c, species, k = 5, bs = 'fs', xt = list(bs = 'tp')) +
        # to account for seasonal changes in day length
        ti(doy, tod_pdt, species, k = 5, bs = c('cc', 'cc', 're')) +
        # to account for changes in day nocturnality with temperature
        ti(temp_c, tod_pdt, species, k = 5, bs = c('cc', 'cc', 're'))),
        family = gammals(),
        data = d_2,
        method = 'REML', # cannot use fREML
        knots = list(tod_pdt = c(0, 1), doy = c(0.5, 366.5)),
        control = gam.control(trace = TRUE))
  
  saveRDS(m_2, paste0('models/gammals-gam-', Sys.Date(), '.rds'))
  
  draw(m_2, contour = FALSE, rug = FALSE,
       discrete_colour = scale_color_manual(values = PAL))
  
  summary(m_2)
  
  mutate(d_2, e = resid(m_2)) %>%
    ggplot(aes(e, fill = species)) +
    facet_wrap(~ species, scales = 'free') +
    geom_histogram(show.legend = FALSE) +
    scale_fill_manual('Species', values = PAL)
  
  layout(matrix(1:4, ncol = 2))
  gam.check(m_2)
  layout(1)
}

# s(tod_pdt) ----
p_s_tod <-
  bind_cols(newd_tod,
            predict(m_2, newdata = newd_tod,
                    terms = c('species',
                              paste0('s(tod_pdt):species', SPECIES)),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_2$family$linkinv(fit),
         lwr = m_2$family$linkinv(fit - 1.96 * se.fit),
         upr = m_2$family$linkinv(fit + 1.96 * se.fit)) %>%
  ggplot() +
  facet_wrap(~ species, nrow = 1, scales = 'free_y') +
  # geom_point(aes(tod_pdt, speed_est), d_2, alpha = 0.01) +
  geom_ribbon(aes(tod_pdt, ymin = lwr, ymax = upr, fill = species),
              alpha = 0.2) +
  geom_line(aes(tod_pdt, mu, color = species), linewidth = 1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 24),
                     breaks = c(0, 6, 12, 18),
                     labels = c('00:00', '06:00', '12:00', '18:00')) +
  labs(x = 'Time of day', y = 'Speed (m/s)') +
  theme(legend.position = 'none')

ggsave('figures/speed-tod.png', p_s_tod,
       width = 25, height = 4, dpi = 600, bg = 'white')

# s(doy) ----
p_s_doy <-
  bind_cols(newd_doy,
            predict(m_2, newdata = newd_doy,
                    terms = c('species',
                              paste0('s(doy):species', SPECIES)),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_2$family$linkinv(fit),
         lwr = m_2$family$linkinv(fit - 1.96 * se.fit),
         upr = m_2$family$linkinv(fit + 1.96 * se.fit)) %>%
  ggplot() +
  facet_wrap(~ species, nrow = 1, scales = 'free_y') +
  # geom_point(aes(doy, speed_est), d_2, alpha = 0.01) +
  geom_ribbon(aes(doy, ymin = lwr, ymax = upr, fill = species),
              alpha = 0.2) +
  geom_line(aes(doy, mu, color = species), linewidth = 1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = yday(paste0('2024-', c('03-20', '06-20',
                                                     '09-22', '12-21')))) +
  labs(x = 'Day of year', y = 'Speed (m/s)') +
  theme(legend.position = 'none')

ggsave('figures/speed-doy.png', p_s_doy,
       width = 25, height = 4, dpi = 600, bg = 'white')

# s(temp_c) ----
p_s_temp <-
  bind_cols(newd_temp_c,
            predict(m_2, newdata = newd_temp_c,
                    terms = c('species',
                              paste0('s(temp_c):species', SPECIES)),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_2$family$linkinv(fit),
         lwr = m_2$family$linkinv(fit - 1.96 * se.fit),
         upr = m_2$family$linkinv(fit + 1.96 * se.fit)) %>%
  ggplot() +
  facet_wrap(~ species, nrow = 1, scales = 'free_y') +
  # geom_point(aes(temp_c, speed_est), d_2, alpha = 0.01) +
  geom_ribbon(aes(temp_c, ymin = lwr, ymax = upr, fill = species),
              alpha = .2) +
  geom_line(aes(temp_c, mu, color = species), linewidth =  1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)'))+
  scale_y_continuous('Speed (m/s)') +
  theme(legend.position = 'none')

ggsave('figures/speed-temperature.png', p_s_temp,
       width = 25, height = 4, dpi = 600, bg = 'white')

# full figure for speed | moving ----
p_s <- plot_grid(p_s_temp, p_s_tod, p_s_doy,
                 labels = 'AUTO', ncol = 1)

ggsave('figures/speed-all.png',
       width = 25, height = 12, dpi = 600, bg = 'white')

# caribou plots only
filter(newd_tod, species == 'Rangifer tarandus (boreal)') %>%
  bind_cols(.,
            predict(m_1, newdata = .,
                    terms = c('species',
                              paste0('s(tod_pdt):species', SPECIES),
                              paste0('ti(doy:tod_pdt):species', SPECIES)),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_1$family$linkinv(fit),
         lwr = m_1$family$linkinv(fit - 1.96 * se.fit),
         upr = m_1$family$linkinv(fit + 1.96 * se.fit)) %>%
  ggplot() +
  coord_polar() +
  # geom_area(aes(x, y), tibble(x = seq(0, 6, by = 0.01), y = 1),
  #           fill = 'black', alpha = 0.3) +
  # geom_area(aes(x, y), tibble(x = seq(18, 24, by = 0.01), y = 1),
  #           fill = 'black', alpha = 0.3) +
  geom_area(aes(tod_pdt, mu, color = species, fill = species), linewidth = 1,
            alpha = 0.5) +
  geom_ribbon(aes(tod_pdt, ymin = lwr, ymax = upr, fill = species),
              alpha = 0.2) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 24),
                     breaks = c(0, 6, 12, 18),
                     labels = c('00:00', '06:00', '12:00', '18:00')) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Time of day', y = 'P(moving)') +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold.italic'),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = c(rep('grey', 5), NA)))
