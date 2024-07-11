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
    readRDS('data/movement-models-speed-weights-temperature-2024-06-21.rds') %>%
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
           doy = yday(timestamp)) # using UTC doy
  
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

# check speed classification for each species
d %>%
  group_by(species) %>%
  summarise(prop_moving = mean(moving),
            min_speed = min(speed_est[moving]))

ggplot(d) +
  facet_wrap(~ species + moving, scales = 'free', ncol = 2) +
  geom_density(aes(speed_est), fill = 'grey')

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

# check spread for tod, doy, and temp for each species except for grizzly
# should be ok to keep cs = 'cc'
plot_grid(
  d %>%
    select(species, tod_pdt, doy, temp_c) %>%
    pivot_longer(-species) %>%
    ggplot(aes(value)) +
    facet_grid(species ~ name, scales = 'free') +
    geom_histogram() +
    ggtitle('Full dataset'),
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
  #' `ti(doy, temp_c)` does not increase BIC or AIC
  m_1 <-
    bam(moving ~
          # random intercept for each animal
          s(animal, bs = 're') +
          # random effect for each species
          s(species, bs = 're') +
          # to account for changes in behavior within days
          s(tod_pdt, species, k = 5, bs = 'fs', xt = list(bs = 'cc')) +
          # to account for changes in behavior within years
          s(doy, species, k = 5, bs = 'fs', xt = list(bs = 'cc')) +
          # species-level effect of temperature
          s(temp_c, species, k = 5, bs = 'fs', xt = list(bs = 'tp')) +
          # to account for seasonal changes in day length
          ti(doy, tod_pdt, species, k = 5, bs = c('cc', 'cc', 're')) +
          # to account for changes in day nocturnality with temperature
          ti(temp_c, tod_pdt, species, k = 5, bs = c('tp', 'cc', 're')) +
          # to account for changes in fur coats seasonally
          ti(temp_c, doy, species, k = 5, bs = c('tp', 'cc', 're')),
        family = binomial(link = 'logit'),
        data = d,
        method = 'fREML', # fast REML
        discrete = TRUE, # discretize the posterior for faster computation
        knots = list(tod_pdt = c(0, 1), doy = c(0.5, 366.5)),# for bs = 'cc'
        control = gam.control(trace = TRUE))
  saveRDS(m_1, paste0('models/binomial-gam.rds'))
  
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

if(file.exists('models/gamma-gamls.rds')) {
  m_2 <- readRDS('models/gamma-gamls.rds')
} else {
  # fir an HGAM with constant and common scale parameter to show issues
  if(file.exists('models/gamma-gam.rds')) {
    m_2 <- readRDS('models/gamma-gam.rds')
  } else {
    m_2 <-
      bam(
        speed_est ~
          # random effect for each animal
          s(animal, bs = 're') +
          # random effect for each species
          s(species, bs = 're') +
          # to account for changes in behavior within days
          s(tod_pdt, species, k = 5, bs = 'fs', xt = list(bs = 'cc')) +
          # to account for changes in behavior within years
          s(doy, species, k = 5, bs = 'fs', xt = list(bs = 'cc')) +
          # species-level effect of temperature
          s(temp_c, species, k = 5, bs = 'fs', xt = list(bs = 'tp')) +
          # to account for seasonal changes in day length
          ti(doy, tod_pdt, species, k = 5, bs = c('cc', 'cc', 're')) +
          # to account for changes in day nocturnality with temperature
          ti(temp_c, tod_pdt, species, k = 5, bs = c('tp', 'cc', 're')) +
          # to account for changes in fur coats seasonally
          ti(temp_c, doy, species, k = 5, bs = c('tp', 'cc', 're')),
        family = Gamma(link = 'log'), # can use Gamma because no zeros
        data = d_2,
        method = 'fREML', # fast REML
        discrete = TRUE, # discretize the posterior for faster computation
        knots = list(tod_pdt = c(0, 1), doy = c(0.5, 366.5)),
        control = gam.control(trace = TRUE))
    saveRDS(m_2, paste0('models/gamma-gam.rds'))
    
    draw(m_2, contour = FALSE, rug = FALSE,
         discrete_colour = scale_color_manual(values = PAL))
    
    appraise(m_2, point_alpha = 0.05)
    
    summary(m_2)
    
    # need different scale parameters for each species
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
      theme(legend.position = 'none')
  }
  
  m_2 <- gam(
    list(
      # mean parameter; E(Y)
      speed_est ~
        # random effect for each animal
        s(animal, bs = 're') +
        # random effect for each species
        s(species, bs = 're') +
        # to account for changes in behavior within days
        s(tod_pdt, species, k = 5, bs = 'fs', xt = list(bs = 'cc')) +
        # to account for changes in behavior within years
        s(doy, species, k = 5, bs = 'fs', xt = list(bs = 'cc')) +
        # species-level effect of temperature
        s(temp_c, species, k = 5, bs = 'fs', xt = list(bs = 'tp')) +
        # to account for seasonal changes in day length
        ti(doy, tod_pdt, species, k = 5, bs = c('cc', 'cc', 're')) +
        # to account for changes in day nocturnality with temperature
        ti(temp_c, tod_pdt, species, k = 5, bs = c('tp', 'cc', 're')) +
        # to account for changes in fur coats seasonally
        ti(temp_c, doy, species, k = 5, bs = c('tp', 'cc', 're')),
      # scale parameter; Var(Y) = E(Y)^2 * scale
      ~
        # random effect for each species
        s(species, bs = 're') +
        # to account for changes in behavior within days
        s(tod_pdt, species, k = 5, bs = 'fs', xt = list(bs = 'cc')) +
        # to account for changes in behavior within years
        s(doy, species, k = 5, bs = 'fs', xt = list(bs = 'cc')) +
        # species-level effect of temperature
        s(temp_c, species, k = 5, bs = 'fs', xt = list(bs = 'tp'))),
    family = gammals(), # can use Gamma because no zeros
    data = d_2,
    method = 'REML',
    knots = list(tod_pdt = c(0, 1), doy = c(0.5, 366.5)),
    control = gam.control(trace = TRUE))
  
  saveRDS(m_2, 'models/gammals-gam.rds')
  
  draw(m_2, contour = FALSE, rug = FALSE, parametric = TRUE,
       discrete_colour = scale_color_manual(values = PAL))
  
  appraise(m_2, point_alpha = 0.05)
  
  summary(m_2)
}

mutate(d_2, mu = predict(m_2, type = 'response')[ , 1]) %>%
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
