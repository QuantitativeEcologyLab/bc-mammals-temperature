library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('lubridate') # for smother date wrangling
library('mgcv')      # for Generalized Additive Models
library('ggplot2')   # for fancy plots
library('khroma')    # for colorblind-friendly color palettes
library('cowplot')   # for fancy multi-panel plots
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids

#' clipping temperature data to the 0.8 quantile doesn't fix the high
#' probabilities for goats
d <-
  readRDS('data/movement-models-speed-weights-temperature-2024-04-05.rds') %>%
  filter(! is.na(speed_est), is.finite(speed_est)) %>%
  mutate(speed_est = round(speed_est, 2), # some values are ~1e-10
         animal = factor(animal),
         species = if_else(grepl('Rangifer', dataset_name),
                           dataset_name,
                           species),
         species = stringr::str_replace_all(species, '_', ' ') %>%
           factor(),
         # using a higher limit doesn't change much
         moving = factor(if_else(speed_est == 0, 'no', 'yes')),
         timestamp_bc = with_tz(timestamp, tz = 'America/Vancouver'),
         tod = (hour(timestamp) / 24 +
                  minute(timestamp) / 60 / 24 +
                  second(timestamp) / 60 / 60 / 24),
         doy = yday(timestamp)) %>%
  rename(temp_c = temperature_C)

SPECIES <- unique(d$species)
N_SPECIES <- length(SPECIES)

# drop single data point in February for bears
plot(as.numeric(moving) ~ doy,
     filter(d, species == 'Ursus arctos horribilis'))
d <- filter(d, ! (species == 'Ursus arctos horribilis' & doy < 50))

# drop extreme temperatures ----
# some species have data with extreme temperatures but little data
ggplot(d, aes(temp_c)) +
  facet_wrap(~ species, scales = 'free') +
  geom_histogram(color = 'black', binwidth = 1) +
  geom_vline(aes(xintercept = value), color = 'red',
             d %>%
               group_by(species) %>%
               summarize(lwr = quantile(temp_c, 0.025),
                         upr = quantile(temp_c, 0.975)) %>%
               pivot_longer(-species)) +
  ylim(c(0, NA))

if(FALSE) {
  # add some observations where bears are not moving in winter to inform
  # the model better
  quantile(filter(d, species == 'Ursus arctos horribilis')$doy,
           c(0, 0.01, 0.02, 0.05, 0.95, 0.98, 0.99, 1))
  
  # fraction of added data points
  N_HYB <- 100
  N_HYB / sum(d$species == 'Ursus arctos horribilis')
  
  set.seed(1)
  d <- d %>%
    bind_rows(tibble(moving = FALSE,
                     speed_est = 0,
                     species = factor('Ursus arctos horribilis',
                                      levels = levels(d$species)),
                     animal = 'winter_prior',
                     temp_c = -20,
                     tod = runif(N_HYB),
                     doy = sample(x = c(1:127, 290:365), size = N_HYB)))
}

range(d$speed_est)
layout(t(1:3))
hist(d$speed_est, breaks = 100)
hist(log10(d$speed_est), breaks = 5000)
hist(log10(d$speed_est), breaks = 5000, xlim = c(-3, 1))
layout(1)
mean(d$moving == 'yes')

# model P(movement) ----
#' using `bs = 'fs'` rather than `by` smooths because to keep the models
#' reasonably constrained. Using `by` gives extremely high `P(moving)`
#' in spring for grizzly bears
#' fits in ~90 s
if(FALSE) {
  m_1 <-
    bam(moving ~
          # random effect for each animal
          s(animal, bs = 're') +
          # fixed intercept of species without one species as a default
          species +
          # to account for changes in behavior within days
          s(tod, by = species, k = 10, bs = 'cc') +
          # to account for changes in behavior within years
          s(doy, by = species, k = 10, bs = 'cc') +
          # species-level effect of temperature
          s(temp_c, by = species, k = 6, bs = 'tp'),
        family = binomial(link = 'logit'),
        data = d,
        method = 'fREML', # fast REML
        discrete = TRUE,  # discretize the posterior for faster computation
        select = TRUE, # perform model shrinkage
        knots = list(tod = c(0, 1), doy = c(0.5, 366.5)),# for bs = 'cc'
        control = gam.control(trace = TRUE))
  saveRDS(m_1, paste0('models/binomial-hgam-', Sys.Date(), '.rds'))
  
  # southern mountain and boreal caribou move differently
  layout(matrix(c(rep(0, 6), 1:29), ncol = 7, byrow = TRUE))
  plot(m_1, scheme = 3, scale = 0)
  layout(1)
  
  summary(m_1)
  
  layout(matrix(1:4, ncol = 2))
  gam.check(m_1, type = 'pearson') # good q-q plot
  layout(1)
} else {
  m_1 <- readRDS('models/binomial-hgam-2024-05-03.rds')
}

# s(tod) ----
d %>%
  group_by(species) %>%
  summarize(min_tod = round(min(tod), 2),
            max_tod = round(max(tod), 2))

newd_tod <- expand_grid(animal = 'new animal',
                        species = SPECIES,
                        tod = seq(0, 1, length.out = 1e3),
                        doy = 0,
                        temp_c = 0) %>%
  mutate(bc_tod = tod * 24 - 7,
         bc_tod = if_else(bc_tod < 0, bc_tod + 24, bc_tod))

p_mov_tod <-
  bind_cols(newd_tod,
            predict(m_1, newdata = newd_tod,
                    terms = c('species',
                              paste0('s(tod):species', SPECIES)),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_1$family$linkinv(fit),
         lwr = m_1$family$linkinv(fit - 1.96 * se.fit),
         upr = m_1$family$linkinv(fit + 1.96 * se.fit)) %>%
  ggplot() +
  coord_polar() +
  facet_wrap(~ species, nrow = 2) +
  geom_ribbon(aes(bc_tod, ymin = lwr, ymax = upr, fill = species),
              alpha = 0.2) +
  geom_area(aes(x, y), tibble(x = seq(0, 6, by = 0.001), y = 1),
            fill = 'black', alpha = 0.3) +
  geom_area(aes(x, y), tibble(x = seq(18, 24, by = 0.001), y = 1),
            fill = 'black', alpha = 0.3) +
  geom_area(aes(bc_tod, mu, color = species, fill = species), linewidth = 1,
            alpha = 0.5) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 6, 12, 18),
                     labels = c('00:00', '06:00', '12:00', '18:00')) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = 'Time of day', y = 'P(moving)') +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold.italic'),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = c(rep('grey', 5), NA)))

ggsave('figures/p-moving-tod.png', p_mov_tod,
       width = 15, height = 7.5, dpi = 600, bg = 'white')

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

d %>%
  group_by(species) %>%
  summarize(min_doy = round(min(doy), 2),
            max_doy = round(max(doy), 2))

newd_doy <- expand_grid(animal = 'new animal',
                        species = SPECIES,
                        tod = 12,
                        doy = seq(1, 366, by = 0.1),
                        temp_c = 0)

p_mov_doy <-
  bind_cols(newd_doy,
            predict(m_1, newdata = newd_doy,
                    terms = c('species',
                              paste0('s(doy):species', SPECIES)),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_1$family$linkinv(fit),
         lwr = m_1$family$linkinv(fit - 1.96 * se.fit),
         upr = m_1$family$linkinv(fit + 1.96 * se.fit)) %>%
  ggplot() +
  coord_polar(start = 11 / 366 * 2 * pi) + # offset to make axes vertical
  facet_wrap(~ species, nrow = 2) +
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
       width = 15, height = 7.5, dpi = 600, bg = 'white')

# s(temp_c) ----
newd_temp_c <- tibble(animal = 'new animal',
                      species = SPECIES,
                      tod = 12,
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
                    terms = c('species',
                              paste0('s(temp_c):species', SPECIES)),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_1$family$linkinv(fit),
         lwr = m_1$family$linkinv(fit - 1.96 * se.fit),
         upr = m_1$family$linkinv(fit + 1.96 * se.fit)) %>%
  ggplot() +
  facet_wrap(~ species, nrow = 2) +
  geom_ribbon(aes(temp_c, ymin = lwr, ymax = upr, fill = species),
              alpha = .2) +
  geom_line(aes(temp_c, mu, color = species), linewidth =  1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)'))+
  scale_y_continuous('P(moving)', expand = c(0, 0), limits = c(0, 1)) +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'italic'))

ggsave('figures/p-moving-temperature.png', p_mov_temp,
       width = 15, height = 7.5, dpi = 600, bg = 'white')

# full figure ----
p_mov <- plot_grid(p_mov_temp,
                   plot_grid(p_mov_tod, p_mov_doy, labels = c('B', 'C'),
                             nrow = 1),
                   labels = c('A', ''), nrow = 2)

ggsave('figures/p-moving-all.png', p_mov,
       width = 30, height = 15, dpi = 600, bg = 'white')

# model mean speed given that animals are moving ----
d_2 <- filter(d, moving == 'yes')

if(FALSE) {
  # takes <1 minute to fit
  # non-gaussian residuals for single species, too
  m_2 <-
    bam(speed_est ~
          # random effect for each animal
          s(animal, bs = 're') +
          # fixed effect of species without one species as the default
          species +
          # to account for changes in behavior within days
          s(tod, by = species, k = 10, bs = 'cc') +
          # to account for changes in behavior within years
          s(doy, by = species, k = 10, bs = 'cc') +
          # species-level effect of temperature
          s(temp_c, by = species, k = 6, bs = 'tp'),
        family = Gamma(link = 'log'), # can use Gamma because no zeros
        data = d_2,
        method = 'fREML', # fast REML
        discrete = TRUE, # discretize the posterior for faster computation
        knots = list(tod = c(0, 1), doy = c(0.5, 366.5)),
        control = gam.control(trace = TRUE))
  
  saveRDS(m_2, paste0('models/gamma-hgam-', Sys.Date(), '.rds'))
  
  layout(matrix(c(rep(0, 6), 1:22), ncol = 7, byrow = TRUE))
  plot(m_2, scheme = 1, scale = 0)
  
  summary(m_2)
  
  layout(matrix(1:4, ncol = 2))
  gam.check(m_2)
  layout(1)
} else {
  m_2 <- readRDS('models/gamma-hgam-2024-04-26.rds')
}

# s(tod) ----
p_s_tod <-
  bind_cols(newd_tod,
            predict(m_2, newdata = newd_tod,
                    terms = c('species',
                              paste0('s(tod):species', SPECIES)),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_2$family$linkinv(fit),
         lwr = m_2$family$linkinv(fit - 1.96 * se.fit),
         upr = m_2$family$linkinv(fit + 1.96 * se.fit)) %>%
  ggplot() +
  facet_wrap(~ species, nrow = 2, scales = 'free_y') +
  geom_ribbon(aes(bc_tod, ymin = lwr, ymax = upr, fill = species),
              alpha = 0.2) +
  geom_line(aes(bc_tod, mu, color = species), linewidth = 1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 24),
                     breaks = c(0, 6, 12, 18),
                     labels = c('00:00', '06:00', '12:00', '18:00')) +
  labs(x = 'Time of day', y = 'Speed (m/s)') +
  theme(legend.position = 'none')

ggsave('figures/speed-tod.png', p_s_tod,
       width = 15, height = 7.5, dpi = 600, bg = 'white')

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
  facet_wrap(~ species, nrow = 2, scales = 'free_y') +
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
       width = 15, height = 7.5, dpi = 600, bg = 'white')

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
  facet_wrap(~ species, nrow = 2, scales = 'free_y') +
  geom_ribbon(aes(temp_c, ymin = lwr, ymax = upr, fill = species),
              alpha = .2) +
  geom_line(aes(temp_c, mu, color = species), linewidth =  1) +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)'))+
  scale_y_continuous('Speed (m/s)') +
  theme(legend.position = 'none')

ggsave('figures/speed-temperature.png', p_s_temp,
       width = 15, height = 7.5, dpi = 600, bg = 'white')

# full figure ----
p_s <- plot_grid(p_s_temp,
                 plot_grid(p_s_tod, p_s_doy, labels = c('B', 'C'),
                           nrow = 1),
                 labels = c('A', ''), nrow = 2)

ggsave('figures/speed-all.png',
       p_s, width = 30, height = 15, dpi = 600, bg = 'white')
