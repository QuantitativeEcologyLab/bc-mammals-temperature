library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for smoother date wrangling
library('mgcv')      # for Generalized Additive Models
library('gratia')    # for ggplot-based figures for GAMs
library('ggplot2')   # for fancy plots
library('khroma')    # for colorblind-friendly color palettes
library('cowplot')   # for fancy multi-panel plots
library('ctmm')      #'for unit conversions with `%#%`
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
    #' cluster speeds using `kmeans()` (does a poor job, see below)
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

#' `kmeans()` splits data poorly
splits <-
  tibble(species = SPECIES,
         lab = SPECIES_LABS,
         manual = c(0.06, 0.05, 0.08, 0.07, 0.10, 0.20, 0.05),
         `k means` = map_dbl(species, \(sp) {
           km <- kmeans(filter(d, species == sp)$speed_est, centers = 2)
           max(min(filter(d, species == sp) %>%
                     filter(km$cluster == 1) %>%
                     pull(speed_est)),
               min(filter(d, species == sp) %>%
                     filter(km$cluster == 2) %>%
                     pull(speed_est)))
         }),
  ) %>%
  pivot_longer(c(manual, `k means`), values_to = 'split_m_s',
               names_to = 'method') %>%
  mutate(split_km_h = 'km/h' %#% split_m_s)

hists <- list()
for(i in (1:N_SPECIES)[order(SPECIES)]) {
  sp <- as.character(SPECIES[i])
  XLIM <- c(0, quantile(filter(d, species == sp)$speed_est, 0.99))
  XLIM <- 'km/h' %#% XLIM # convert to km/h
  YLIM <- c(0, case_when(
    sp == 'Rangifer tarandus (s. mountain)' ~ 400,
    sp == 'Ursus arctos horribilis' ~ 1e3,
    sp == 'Puma concolor' ~ 1e3,
    sp == 'Cervus canadensis' ~ 25e3,
    sp == 'Rangifer tarandus (boreal)' ~ 5e3,
    sp == 'Canis lupus' ~ 4e3,
    sp == 'Oreamnos americanus' ~ 500))
  
  hists[[i]] <-
    ggplot() +
    geom_histogram(aes('km/h' %#% speed_est), filter(d, species == sp),
                   bins = 200, fill = 'grey', color = 'black') +
    geom_vline(aes(xintercept = split_km_h, color = method,
                   lty = method), filter(splits, species == sp), lwd = 0.75) +
    coord_cartesian(xlim = XLIM, ylim = YLIM) +
    labs(x = 'Speed (km/h)', y = 'Count',
         title = scales::parse_format()(SPECIES_LABS[i])) +
    scale_color_bright(name = 'Method') +
    scale_linetype_manual(name = 'Method', values = c(2, 1)) +
    theme(legend.position = 'none')
  rm(sp, XLIM, YLIM)
}

hists[1:7] <- hists[order(SPECIES)]
hists[[8]] <-
  get_legend(
    ggplot() +
      geom_vline(aes(xintercept = split_km_h, color = method,
                     lty = method), splits, lwd = 0.75) +
      scale_color_bright(name = 'Method') +
      scale_linetype_manual(name = 'Method', values = c(2, 1)) +
      theme(legend.position = 'right'))
plot_grid(plotlist = hists, nrow = 2)
ggsave(filename = 'figures/moving-0-1-split.png',
       width = 16, height = 5.5, units = 'in', dpi = 600, bg = 'white')

# assign moving/not moving status
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

# dataset of non-zero speeds only
d_2 <- filter(d, moving)

ggplot(d_2) +
  facet_wrap(~ species, scales = 'free', ncol = 2) +
  geom_density(aes(speed_est), fill = 'grey')

d <- mutate(d,
            lab = case_when(species == SPECIES[1] ~ SPECIES_LABS[1],
                            species == SPECIES[2] ~ SPECIES_LABS[2],
                            species == SPECIES[3] ~ SPECIES_LABS[3],
                            species == SPECIES[4] ~ SPECIES_LABS[4],
                            species == SPECIES[5] ~ SPECIES_LABS[5],
                            species == SPECIES[6] ~ SPECIES_LABS[6],
                            species == SPECIES[7] ~ SPECIES_LABS[7]))

d_2 <- mutate(d_2,
              lab = case_when(species == SPECIES[1] ~ SPECIES_LABS[1],
                              species == SPECIES[2] ~ SPECIES_LABS[2],
                              species == SPECIES[3] ~ SPECIES_LABS[3],
                              species == SPECIES[4] ~ SPECIES_LABS[4],
                              species == SPECIES[5] ~ SPECIES_LABS[5],
                              species == SPECIES[6] ~ SPECIES_LABS[6],
                              species == SPECIES[7] ~ SPECIES_LABS[7]))


# check spread for tod, doy, and temp for each species except for grizzly
# should be ok to keep cs = 'cc'
hist_1 <- d %>%
  select(lab, tod_pdt, doy, temp_c, moving) %>%
  pivot_longer(-c(lab, moving), names_to = 'variable') %>%
  mutate(variable = case_when(
    variable == 'doy' ~ 'Day~of~year',
    variable == 'temp_c' ~ paste0('Temperature~(degree*C)'),
    variable == 'tod_pdt' ~ 'Time~of~day~(PDT)') %>%
      paste0('bold(', ., ')') %>%
      factor(., levels = unique(.))) %>%
  ggplot(aes(value)) +
  facet_grid(lab ~ variable, scales = 'free', switch = 'x',
             labeller = label_parsed) +
  geom_histogram(aes(fill = moving), position = 'stack', bins = 24) +
  scale_x_continuous(NULL, expand = c(0, 0)) +
  ylab('Count') +
  scale_fill_brewer('Moving', type = 'qual', palette = 6,
                    breaks = c(FALSE, TRUE), labels = c('No', 'Yes')) +
  theme(legend.position = 'none', strip.background.x = element_blank(),
        strip.text.x = element_text(size = 11), strip.placement = 'outside')

hist_2 <- d_2 %>%
  select(lab, tod_pdt, doy, temp_c, moving) %>%
  pivot_longer(-c(lab, moving), names_to = 'variable') %>%
  mutate(variable = case_when(
    variable == 'doy' ~ 'Day~of~year',
    variable == 'temp_c' ~ paste0('Temperature~(degree*C)'),
    variable == 'tod_pdt' ~ 'Time~of~day~(PDT)') %>%
      paste0('bold(', ., ')') %>%
      factor(., levels = unique(.))) %>%
  ggplot(aes(value)) +
  facet_grid(lab ~ variable, scales = 'free', switch = 'x',
             labeller = label_parsed) +
  geom_histogram(fill = '#377EB8', position = 'stack', bins = 24) +
  scale_x_continuous(NULL, expand = c(0, 0)) +
  ylab('Count') +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_text(size = 11), strip.placement = 'outside')

plot_grid(get_plot_component(hist_1 + theme(legend.position = 'top'),
                             pattern = 'guide-box-top', return_all = TRUE),
          plot_grid(hist_1, hist_2, labels = 'AUTO', nrow = 1),
          ncol = 1, rel_heights = c(1, 15))

ggsave('figures/temperature-movement-rates-hist.png',
       width = 16, height = 17, units = 'in', dpi = 600, bg = 'white')

# model P(movement) ----
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
        discrete = TRUE, # discretize the covariates for faster computation
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
p_op <- 
  transmute(d,
            lab,
            moving,
            est = predict(m_1, type = 'response') %>%
              round(2)) %>%
  group_by(lab, est) %>%
  summarise(empirical_mean = mean(moving),
            n = n(),
            se = sd(moving) / sqrt(n),
            lwr = empirical_mean - se,
            upr = empirical_mean + se,
            .groups = 'drop') %>%
  ggplot(aes(est, empirical_mean, color = log2(n))) +
  facet_wrap(~ lab, nrow = 2, labeller = label_parsed) +
  geom_abline(slope = 1, intercept = 0, color = 'grey') +
  geom_errorbar(aes(est, ymin = lwr, ymax = upr), width = 0, alpha = 0.3) +
  geom_point() +
  lims(x = c(0, 1), y = c(0, 1)) +
  scale_color_viridis_c(breaks = 4 * (0:5), labels = 2^(2 * 0:5),
                        limits = c(0, 16)) +
  labs(x = 'Predicted P(moving)', y = 'Observed P(moving)') +
  theme(legend.position = 'inside', legend.position.inside = c(0.875, 0.25))

ggsave('figures/p-moving-observed-vs-predicted.png', p_op,
       width = 12, height = 6, dpi = 600, bg = 'white')

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
      discrete = TRUE, # discretize the covariates for faster computation
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

s_op <-
  mutate(d_2, mu = predict(m_2, type = 'response')) %>%
  ggplot(aes(mu, speed_est)) +
  facet_wrap(~ lab, scales = 'free', nrow = 2, labeller = label_parsed) +
  geom_hex(aes(fill = log10(after_stat(count)))) +
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  theme(legend.position = 'none') +
  scale_fill_iridescent(name = expression(bold(Count~(log['10']~scale))),
                        range = c(0.3, 1), breaks = 0:3,
                        labels = round(10^(c(0:3)))) +
  ylim(c(0, NA)) +
  labs(x = 'Predicted speed (m/s)', y = 'Observed speed (m/s)') +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  theme(legend.position = 'top')

ggsave('figures/speed-observed-vs-predicted.png', s_op,
       width = 12, height = 6, dpi = 600, bg = 'white')

op <- plot_grid(p_op, s_op, labels = 'AUTO', ncol = 1)

ggsave('figures/p-and-speed-observed-vs-predicted.png', op,
       width = 12, height = 12, dpi = 600, bg = 'white')

# sampling interval has a small effect on estimated state and speed ----
dt_breaks <- c(-5, 0, 5)
dt_labs <- round(exp(dt_breaks), c(3, 0, 0))

p_dt <-
  ggplot(d) +
  facet_wrap(~ lab, labeller = label_parsed, nrow = 1) +
  geom_histogram(aes(log(dt), color = species, fill = species), bins = 15,
                 alpha = 0.2) +
  scale_x_continuous(name = expression(bold(paste('\U0394', 't (hours, log scale)'))),
                     breaks = dt_breaks, labels = dt_labs) +
  scale_y_log10(expression(bold(Count~(log["10"]~scale)))) +
  scale_color_manual(values = PAL, aesthetics = c('color', 'fill')) +
  theme(legend.position = 'none')

preds_dt <-
  d %>%
  group_by(species, lab) %>%
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
                s_upr = m_2$family$linkinv(fit + 1.96 * se.fit)))

p_1_dt <-
  ggplot(preds_dt) +
  facet_wrap(~ lab, labeller = label_parsed, nrow = 1) +
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
  facet_wrap(~ lab, labeller = label_parsed, nrow = 1) +
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

# fit HGAMs without temperature ----
# P(moving)
if(file.exists('models/binomial-gam-without-temperature.rds')) {
  m_1_no_t <- readRDS('models/binomial-gam-without-temperature.rds')
} else {
  m_1_no_t <-
    bam(moving ~
          # random intercept for each animal
          s(animal, bs = 're') +
          # fixed intercept for each species
          species +
          # to account for changes in behavior within days
          s(tod_pdt, by = species, k = 5, bs = 'cc') +
          # to account for changes in behavior within years
          s(doy, by = species, k = 5, bs = 'cc') +
          # to account for seasonal changes in day length
          ti(doy, tod_pdt, by = species, k = 5, bs = c('cc', 'cc')) +
          # larger sampling intervals underestimate movement speed
          s(log(dt), k = 3) +
          s(log(dt), species, k = 3, bs = 'fs'),
        family = binomial(link = 'logit'),
        data = d,
        method = 'fREML', # fast REML
        discrete = TRUE, # discretize the covariates for faster computation
        knots = list(tod_pdt = c(0, 1), doy = c(0.5, 366.5)),# for bs = 'cc'
        control = gam.control(trace = TRUE))
  
  saveRDS(m_1_no_t, 'models/binomial-gam-without-temperature.rds')
  qq.gam(m_1_no_t, type = 'pearson') # some extreme outliers
}

# speed
if(file.exists('models/gamma-gam-without-temperature.rds')) {
  m_2_no_t <- readRDS('models/gamma-gam-without-temperature.rds')
} else {
  m_2_no_t <-
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
        # to account for seasonal changes in day length
        ti(doy, tod_pdt, by = species, k = 5, bs = c('cc', 'cc')) +
        # larger sampling intervals underestimate movement speed
        s(log(dt), k = 3) +
        s(log(dt), species, k = 3, bs = 'fs'),
      family = Gamma(link = 'log'), # can use Gamma because no zeros
      data = d_2,
      method = 'fREML', # fast REML
      discrete = TRUE, # discretize the covariates for faster computation
      knots = list(tod_pdt = c(0, 1), doy = c(0.5, 366.5)),
      control = gam.control(trace = TRUE))
  
  saveRDS(m_2_no_t, 'models/gamma-gam-without-temperature.rds')
}

# unsurprisingly, deviance explained does not change much since predictors
# are correlated
de_1 <- round((1 - m_1$deviance / m_1$null.deviance) * 100, 2)
de_1_nt <- round((1 - m_1_no_t$deviance / m_1_no_t$null.deviance) * 100, 2)
de_1
de_1_nt
de_1 / de_1_nt

de_2 <- round((1 - m_2$deviance / m_2$null.deviance) * 100, 2)
de_2_nt <- round((1 - m_2_no_t$deviance / m_2_no_t$null.deviance) * 100, 2)
de_2
de_2_nt
de_2 / de_2_nt

# RMSE does not change much
rmse <- function(.m) {
  sqrt(mean(resid(.m)^2))
}

tibble(model = c('P(moving)', 'speed'),
       rmse_without_temp = c(rmse(m_1_no_t), rmse(m_2_no_t)),
       rmse_with_temp = c(rmse(m_1), rmse(m_2)),
       perc_change = (((rmse_with_temp / rmse_without_temp) - 1) * 100) %>%
         round(2) %>%
         paste0('%')) %>%
  knitr::kable()

# P(moving)
AIC(m_1, m_1_no_t) %>%
  mutate(delta = AIC - min(AIC))

# speed
AIC(m_2, m_2_no_t) %>%
  mutate(delta = AIC - min(AIC))

# figure of comparison between predictions
p_preds <-
  plot_grid(ggplot() +
              geom_hex(aes(fitted(m_1), fitted(m_1_no_t),
                           fill = log10(after_stat(count)))) +
              geom_abline(intercept = 0, slope = 1, color = 'black') +
              labs(x = 'Predictions without temperature',
                   y = 'Predictions with temperature') +
              scale_fill_iridescent(name = expression(bold(Count~(log['10']~scale))),
                                    range = c(0.3, 1), breaks = (0:3) * 2,
                                    labels = round(10^(c(0:3) * 2))) +
              theme(legend.position = 'top'),
            ggplot() +
              geom_hex(aes(fitted(m_2), fitted(m_2_no_t),
                           fill = log10(after_stat(count)))) +
              geom_abline(intercept = 0, slope = 1, color = 'black') +
              labs(x = 'Predictions without temperature',
                   y = 'Predictions with temperature') +
              scale_fill_iridescent(name = expression(bold(Count~(log['10']~scale))),
                                    range = c(0.3, 1), breaks = (0:3) * 2,
                                    labels = round(10^(c(0:3) * 2))) +
              theme(legend.position = 'top'),
            labels = 'AUTO')

ggsave('figures/hgam-with-without-temp-prediction-agreement.png',
       plot = p_preds, width = 12, height = 6, units = 'in', dpi = 600,
       bg = 'white')
