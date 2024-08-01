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
if(file.exists('data/hgam-speed-data.rds')) {
  d <- readRDS('data/hgam-speed-data.rds')
} else {
  d <-
    readRDS('data/movement-models-speed-weights-temperature-2024-06-10.rds') %>%
    rename(temp_c = temperature_C) %>%
    filter(! is.na(speed_est), is.finite(speed_est)) %>%
    mutate(animal = factor(animal),
           species = if_else(grepl('Rangifer', dataset_name),
                             dataset_name,
                             species),
           species = stringr::str_replace_all(species, '_', ' ') %>%
             stringr::str_replace_all('boreal', '(boreal)') %>%
             stringr::str_replace_all('southern mountain',
                                      '(s. mountain)') %>%
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
    group_by(animal) %>%
    mutate(dt = as.numeric(timestamp - lag(timestamp), units = 'hours')) %>%
    ungroup() %>%
    filter(! is.na(dt))
  
  # drop single data point in February for bears
  plot(as.numeric(moving) ~ doy,
       filter(d, species == 'Ursus arctos horribilis'))
  d <- filter(d, ! (species == 'Ursus arctos horribilis' & doy < 50))
  points(as.numeric(moving) ~ doy,
         filter(d, species == 'Ursus arctos horribilis'), pch = 'x', col = 'red')
  
  saveRDS(d, 'data/hgam-speed-data.rds')
}

SPECIES <- unique(d$species)
N_SPECIES <- length(SPECIES)

#' `kmeans()` splits data poorly
manual_splits <- c(0.06, 0.05, 0.08, 0.07, 0.1, 0.2, 0.05)
km_splits <- rep(NA, N_SPECIES)

layout(matrix(c(1:7, 0), ncol = 4, byrow = TRUE))
for(i in 1:N_SPECIES) {
  km <- kmeans(filter(d, species == SPECIES[i])$speed_est, centers = 2)
  km_splits[i] <- max(min(filter(d, species == SPECIES[i]) %>%
                            filter(km$cluster == 1) %>%
                            pull(speed_est)),
                      min(filter(d, species == SPECIES[i]) %>%
                            filter(km$cluster == 2) %>%
                            pull(speed_est)))
  hist(filter(d, species == SPECIES[i])$speed_est, breaks = 100,
       xlim = c(0,
                max(filter(d, species == SPECIES[i])$speed_est)) * 0.3,
       ylim = c(0, nrow(filter(d, species == SPECIES[i])) / 40),
       main = SPECIES[i], xlab = 'Speed (m/s)')
  abline(v = km_splits[i], col = 'red', lwd = 2)
  abline(v = manual_splits[i], col = 'blue', lwd= 2)
}
layout(1)

d <- mutate(d,
            moving =
              (species == 'Rangifer tarandus (s. mountain)' & speed_est > 0.06) |
              (species == 'Ursus arctos horribilis' & speed_est > 0.05) |
              (species == 'Puma concolor' & speed_est > 0.08) |
              (species == 'Cervus canadensis' & speed_est > 0.07) |
              (species == 'Rangifer tarandus (boreal)' & speed_est > 0.10) |
              (species == 'Canis lupus' & speed_est > 0.20) |
              (species == 'Oreamnos americanus' & speed_est > 0.05))

# check speed classification for each species
d %>%
  group_by(species) %>%
  summarise(prop_moving = mean(moving),
            min_speed = min(speed_est[moving]))

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

# dataset of non-zero speeds only
d_2 <- filter(d, moving)

ggplot(d_2) +
  facet_wrap(~ species, scales = 'free', ncol = 2) +
  geom_density(aes(speed_est), fill = 'grey')

# check spread for tod, doy, and temp for each species except for grizzly
# should be ok to keep cs = 'cc'
plot_grid(
  d %>%
    select(species, tod_pdt, doy, temp_c, moving) %>%
    pivot_longer(-c(species, moving)) %>%
    ggplot(aes(value)) +
    facet_grid(species ~ name, scales = 'free') +
    geom_histogram(aes(fill = moving), position = 'stack') +
    ggtitle('Full dataset') +
    scale_fill_brewer(type = 'qual', palette = 6) +
    theme(legend.position = 'top'),
  d_2 %>%
    select(species, tod_pdt, doy, temp_c) %>%
    pivot_longer(-species) %>%
    ggplot(aes(value)) +
    facet_grid(species ~ name, scales = 'free') +
    geom_histogram() +
    ggtitle('Moving only'),
  labels = 'AUTO', nrow = 1)

# model P(movement) ----
#' using `bs = 'fs'` rather than `by` smooths because to keep the models
#' reasonably constrained. Using `by` gives extremely high `P(moving)`
#' in spring for grizzly bears
#' fits in 3-4 minutes
if(file.exists('models/binomial-gam.rds')) {
  m_1 <- readRDS('models/binomial-gam.rds')
} else {
m_1 <-
    bam(moving ~
          # random intercept for each animal
          s(animal, bs = 're') +
          # fixed intercept for each species
          species +
          # to account for changes in behavior within days
          s(tod_pdt, by = species, k = 5, bs = 'cc') +
          # to account for changes in behavior within years
          s(doy, by = species, k = 5, bs = 'cc') +
          # species-level effect of temperature
          s(temp_c, by = species, k = 5, bs = 'tp') +
          # to account for seasonal changes in day length
          ti(doy, tod_pdt, by = species, k = 5, bs = c('cc', 'cc')) +
          # to account for changes in day nocturnality with temperature
          ti(temp_c, tod_pdt, by = species, k = 5, bs = c('tp', 'cc')) +
          # to account for changes in fur coats seasonally
          ti(temp_c, doy, by = species, k = 5, bs = c('tp', 'cc')) +
          # larger sampling intervals underestimate movement speed
          s(log(dt), k = 3) +
          s(log(dt), species, k = 3, bs = 'fs'),
        family = binomial(link = 'logit'),
        data = d,
        method = 'fREML', # fast REML
        discrete = TRUE, # discretize the posterior for faster computation
        knots = list(tod_pdt = c(0, 1), doy = c(0.5, 366.5)),# for bs = 'cc'
        control = gam.control(trace = TRUE))
  saveRDS(m_1, 'models/binomial-gam.rds')
  
  qq.gam(m_1, type = 'pearson') # good q-q plot
  
  summary(m_1)
  
  # southern mountain and boreal caribou move differently
  draw(m_1, contour = FALSE, rug = FALSE,
       discrete_colour = scale_color_manual(values = PAL))
}

# check predictions ----
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
  scale_color_viridis_c(breaks = 4 * (0:5), labels = 2^(2 * 0:5),
                        limits = c(0, 16)) +
  labs(x = 'Predicted P(moving)', y = 'Observed P(moving)')

ggsave('figures/p-moving-observed-vs-predicted.png',
       width = 8, height = 8, dpi = 600, bg = 'white')

# model mean speed given that animals are moving ----
# sampling is relatively consistent over DOY with some additional bursts
ggplot(d_2, aes(doy)) +
  coord_cartesian(ylim = c(0, 0.025)) +
  facet_wrap(~ species, scales = 'free') +
  geom_density(aes(group = animal, color = species), show.legend = FALSE) +
  scale_color_manual('Species', values = paste0(PAL, '40'))

ggplot(d_2, aes(doy)) +
  facet_wrap(~ species, scales = 'free') +
  geom_density(aes(group = animal, color = species), show.legend = FALSE) +
  scale_color_manual('Species', values = paste0(PAL, '40'))

# each species has stochastic speed but only one state
ggplot(d_2, aes(speed_est)) +
  facet_wrap(~ species, scales = 'free') +
  geom_density(aes(group = animal, color = species), show.legend = FALSE) +
  scale_x_continuous(expression(Estiamted~'speed,'~log[2]~scale),
                     transform = 'log2') +
  scale_color_manual('Species', values = paste0(PAL, '40'))

if(file.exists('models/gamma-gam.rds')) {
  m_2 <- readRDS('models/gamma-gam.rds')
} else {
  m_2 <-
    bam(
      speed_est ~
        # random intercept for each animal
        s(animal, bs = 're') +
        # fixed intercept for each species
        species +
        # to account for changes in behavior within days
        s(tod_pdt, by = species, k = 5, bs = 'cc') +
        # to account for changes in behavior within years
        s(doy, by = species, k = 5, bs = 'cc') +
        # species-level effect of temperature
        s(temp_c, by = species, k = 5, bs = 'tp') +
        # to account for seasonal changes in day length
        ti(doy, tod_pdt, by = species, k = 5, bs = c('cc', 'cc')) +
        # to account for changes in day nocturnality with temperature
        ti(temp_c, tod_pdt, by = species, k = 5, bs = c('tp', 'cc')) +
        # to account for changes in fur coats seasonally
        ti(temp_c, doy, by = species, k = 5, bs = c('tp', 'cc')) +
        # larger sampling intervals underestimate movement speed
        s(log(dt), k = 3) +
        s(log(dt), species, k = 3, bs = 'fs'),
      family = Gamma(link = 'log'), # can use Gamma because no zeros
      data = d_2,
      method = 'fREML', # fast REML
      discrete = TRUE, # discretize the posterior for faster computation
      knots = list(tod_pdt = c(0, 1), doy = c(0.5, 366.5)),
      control = gam.control(trace = TRUE))
  hist(resid(m_2))
  qqnorm(resid(m_2))
  qqline(resid(m_2), col = 'red')
  beepr::beep()
  saveRDS(m_2, 'models/gamma-gam.rds')
  
  draw(m_2, contour = FALSE, rug = FALSE,
       discrete_colour = scale_color_manual(values = PAL))
  
  appraise(m_2, point_alpha = 0.05)
  
  summary(m_2)
  
  # do not need different scale parameters for each species
  mutate(d_2, e = resid(m_2)) %>%
    group_by(species) %>%
    summarise(mean = mean(e), median = median(e),
              sd = sd(e), iqr = IQR(e)) %>%
    arrange(sd)
  
  mutate(d_2, e = resid(m_2)) %>%
    ggplot(aes(e, fill = species)) +
    facet_wrap(~ species, scales = 'free_y') +
    geom_histogram(show.legend = FALSE, color = 'black') +
    scale_fill_manual('Species', values = PAL) +
    labs(x = 'Deviance residuals', y = 'Count')
  
  mutate(d_2, e = resid(m_2)) %>%
    ggplot(aes(sample = e)) +
    facet_wrap(~ species) +
    geom_qq_line(color = 'red') +
    geom_qq(aes(color = species), show.legend = FALSE, alpha = 0.1) +
    scale_color_manual('Species', values = PAL,
                       aesthetics = c('color', 'fill')) +
    labs(x = 'Deviance residuals', y = 'Density')
  
  mutate(d_2, mu = predict(m_2, type = 'response')) %>%
    ggplot(aes(mu, speed_est, color = species)) +
    facet_wrap(~ species, scales = 'free') +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(alpha = 0.1) +
    scale_color_manual('Species', values = PAL,
                       aesthetics = c('color', 'fill')) +
    theme(legend.position = 'none')
}

mutate(d_2, mu = predict(m_2, type = 'response')) %>%
  ggplot(aes(mu, speed_est, color = species)) +
  facet_wrap(~ species, scales = 'free') +
  geom_abline(intercept = 0, slope = 1, color = 'grey') +
  geom_point(alpha = 0.3, shape = '.', size = 2) +
  theme(legend.position = 'none') +
  scale_color_manual('Species', values = PAL,
                     aesthetics = c('color', 'fill')) +
  labs(x = 'Predicted', y = 'Observed')

ggsave('figures/speed-observed-vs-predicted.png',
       width = 8, height = 8, dpi = 600, bg = 'white')

# sampling interval has a small effect on estimated state and speed ----
dt_breaks <- c(-5, 0, 5)
dt_labs <- round(exp(dt_breaks), c(3, 0, 0))

p_dt <-
  d %>%
  mutate(species = gsub(' ', '~', species),
         species = gsub('~\\(', '\\)~bold\\((', species),
         species = paste0('bolditalic(', species, ')')) %>%
  ggplot() +
  facet_wrap(~ species, labeller = label_parsed, nrow = 1) +
  geom_histogram(aes(log(dt), color = species, fill = species), bins = 15,
                 alpha = 0.2) +
  scale_x_continuous(name = expression(bold(paste('\U0394', 't (hours, log scale)'))),
                     breaks = dt_breaks, labels = dt_labs) +
  scale_y_log10(expression(bold(Count~(log[10]~scale)))) +
  scale_color_manual(values = PAL, aesthetics = c('color', 'fill')) +
  theme(legend.position = 'none')

preds_dt <-
  d %>%
  group_by(species) %>%
  summarize(min_log_dt = min(log(dt)),
            max_log_dt = max(log(dt)),
            animal = 'new animal',
            tod_pdt = 0,
            doy = 0,
            temp_c = 0) %>%
  mutate(dt = map2(min_log_dt, max_log_dt, \(.a, .b) {
    tibble(dt = exp(seq(.a, .b, by = 0.01)))
    })) %>%
  unnest(dt) %>%
  select(! min_log_dt, ! max_log_dt) %>%
  bind_cols(
    predict(m_1, type = 'link', newdata = ., se.fit = TRUE,
            newdata.guaranteed = TRUE, discrete = FALSE,
            terms = c('s(log(dt))', 's(log(dt),species)',
                      '(Intercept)', 's(species)')) %>%
      bind_cols() %>%
      transmute(p = m_1$family$linkinv(fit),
                p_lwr = m_1$family$linkinv(fit - 1.96 * se.fit),
                p_upr = m_1$family$linkinv(fit + 1.96 * se.fit)),
    mu_2 = predict(m_2, type = 'link', newdata = ., se.fit = TRUE,
                   newdata.guaranteed = TRUE, discrete = FALSE,
                   terms = c('s(log(dt))', 's(log(dt),species)')) %>%
      bind_cols() %>%
      transmute(s = m_2$family$linkinv(fit),
                s_lwr = m_2$family$linkinv(fit - 1.96 * se.fit),
                s_upr = m_2$family$linkinv(fit + 1.96 * se.fit))) %>%
  mutate(species = gsub(' ', '~', species),
         species = gsub('~\\(', '\\)~bold\\((', species),
         species = paste0('bolditalic(', species, ')'))

p_1_dt <-
  ggplot(preds_dt) +
  facet_wrap(~ species, labeller = label_parsed, nrow = 1) +
  geom_ribbon(aes(log(dt), ymin = p_lwr, ymax = p_upr, fill = species),
              alpha = 0.2) +
  geom_line(aes(log(dt), p, color = species)) +
  scale_x_continuous(name = expression(bold(paste('\U0394', 't (hours, log scale)'))),
                     breaks = dt_breaks, labels = dt_labs) +
  ylab('P(moving)') +
  scale_color_manual(values = PAL, aesthetics = c('color', 'fill')) +
  theme(legend.position = 'none')

p_2_dt <-
  ggplot(preds_dt) +
  facet_wrap(~ species, labeller = label_parsed, nrow = 1) +
  geom_hline(yintercept = 1, color = 'grey', linetype = 'dashed') +
  geom_ribbon(aes(log(dt), ymin = s_lwr, ymax = s_upr, fill = species),
              alpha = 0.2) +
  geom_line(aes(log(dt), s, color = species)) +
  scale_x_continuous(name = expression(bold(paste('\U0394', 't (hours, log scale)'))),
                     breaks = dt_breaks, labels = dt_labs) +
  scale_y_continuous('Relative change in speed',
                     breaks = c(0.5, 0.75, 1, 1.25)) +
  scale_color_manual(values = PAL, aesthetics = c('color', 'fill')) +
  theme(legend.position = 'none')

plot_grid(p_dt, p_1_dt, p_2_dt, labels = 'AUTO', ncol = 1)
ggsave('figures/dt-smooths.png', width = 16, height = 8, units = 'in',
       dpi = 600, bg = 'white')
