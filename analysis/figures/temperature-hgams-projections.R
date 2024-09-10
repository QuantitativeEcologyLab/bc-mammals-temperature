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
plot_scheme(PAL, colours = TRUE)

# functions for calculating odds and back-transforming
odds <- function(p) p / (1 - p)
inv_odds <- function(o) 1 / (1 + 1/o)
if(FALSE) inv_odds(odds(0.2)) # check the functions

# import models
m_1 <- readRDS('models/binomial-gam.rds')
m_2 <- readRDS('models/gamma-gam.rds')

# terms to exclude from the prediction
EXCLUDE <- c('s(animal)',
             paste0('s(tod_pdt):species', SPECIES),
             paste0('ti(doy,tod_pdt):species', SPECIES),
             paste0('ti(temp_c,tod_pdt):species', SPECIES),
             's(log(dt))', 's(log(dt), species)')

# import prediction data for each species in the data's extent ----
cc_newd <- tibble(
  wp = list.files('data', 'weather-projections-', full.names = TRUE) %>%
    map(readRDS)) %>%
  unnest(wp) %>%
  filter(year >= 2020) %>%
  # make sure only the necessary column are kept
  transmute(
    scenario,
    year,
    animal = m_1$model$animal[1],
    species = gsub('boreal', '(boreal)', species) %>%
      gsub('southern mountain', '(s. mountain)', x = .),
    tod_pdt = 0,
    doy = yday(date_decimal(year + (month - 0.5) / 12)),
    temp_c,
    dt = 1,
    weight)

if(file.exists('data/cc-hgam-projections.rds')) {
  cc_proj <- readRDS('data/cc-hgam-projections.rds')
} else {
  cc_proj <- bind_cols(
    cc_newd,
    predict(m_1, newdata = cc_newd, type = 'link', se.fit = TRUE,
            discrete = TRUE, exclude = EXCLUDE) %>%
      as.data.frame() %>%
      transmute(p_lwr = inv_link(m_1)(fit - se.fit * 1.96),
                p = inv_link(m_1)(fit),
                p_upr = inv_link(m_1)(fit + se.fit * 1.96)),
    predict(m_2, newdata = cc_newd, type = 'link', se.fit = TRUE,
            discrete = TRUE, exclude = EXCLUDE) %>%
      as.data.frame() %>%
      transmute(s_lwr = inv_link(m_2)(fit - se.fit * 1.96),
                s = inv_link(m_2)(fit),
                s_upr = inv_link(m_2)(fit + se.fit * 1.96))) %>%
    mutate(d_lwr = p_lwr * s_lwr,
           d = p * s,
           d_upr = p_upr * s_upr) %>%
    group_by(scenario, year, species) %>%
    summarize(p_lwr = mean(p_lwr, weight = weight),
              s_lwr = mean(s_lwr, weight = weight),
              d_lwr = mean(d_lwr, weight = weight),
              p = mean(p, weight = weight),
              s = mean(s, weight = weight),
              d = mean(d, weight = weight),
              p_upr = mean(p_upr, weight = weight),
              s_upr = mean(s_upr, weight = weight),
              d_upr = mean(d_upr, weight = weight),
              .groups = 'drop') %>%
    # center by the mean of the 2020 estimates for each species
    arrange(year, species, scenario) %>%
    group_by(species) %>%
    # convert probabilities to odds for figures
    mutate(o_lwr = odds(p_lwr) / mean(odds(p_lwr[1:4])),
           s_lwr = s_lwr / mean(s_lwr[1:4]),
           d_lwr = d_lwr / mean(d_lwr[1:4]),
           o = odds(p) / mean(odds(p[1:4])),
           s = s / mean(s[1:4]),
           d = d / mean(d[1:4]),
           o_upr = odds(p_upr) / mean(odds(p_upr[1:4])),
           s_upr = s_upr / mean(s_upr[1:4]),
           d_upr = d_upr / mean(d_upr[1:4])) %>%
    ungroup() %>%
    mutate(scenario = case_when(grepl('126', scenario) ~ 'SSP 1-2.6',
                                grepl('245', scenario) ~ 'SSP 2-4.5',
                                grepl('370', scenario) ~ 'SSP 3-7.0',
                                grepl('585', scenario) ~ 'SSP 5-8.5'),
           species = gsub(' ', '~', species) %>%
             gsub('~\\(', '\\)~bold\\((', .) %>%
             paste0('bolditalic(', ., ')') %>%
             factor(levels = SPECIES_LABS))
  saveRDS(cc_proj, 'data/cc-hgam-projections.rds')
}
beepr::beep(2)

# make figures ----
labeller_perc <- function(breaks) {
  breaks <- round(breaks * 100 - 100, 1)
  
  if_else(breaks <= 0,
          paste0(breaks, '%'),
          paste0('+', breaks, '%'))
}

p_o_mov <-
  ggplot(cc_proj) +
  facet_wrap(~ species, scales = 'free_y', drop = FALSE,
             labeller = label_parsed) +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(year, ymin = o_lwr, ymax = o_upr, fill = scenario),
              alpha = 0.2) +
  geom_line(aes(year, o, color = scenario)) +
  scale_color_brewer('Climate change scenario', type = 'div',
                     palette = 5, direction = -1,
                     aesthetics = c('color', 'fill')) +
  labs(x = NULL, y = 'Relative change in odds of moving') +
  scale_y_continuous(labels = labeller_perc) +
  theme(legend.position = 'inside',
        legend.position.inside = c(5/6, 1/6)); p_o_mov
ggsave('figures/odds-moving-local-cc-predictions.png', p_o_mov,
       width = 12, height = 8, dpi = 600, bg = 'white')

p_s <-
  ggplot(cc_proj) +
  facet_wrap(~ species, scales = 'free_y', drop = FALSE,
             labeller = label_parsed) +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(year, ymin = s_lwr, ymax = s_upr, fill = scenario),
              alpha = 0.2) +
  geom_line(aes(year, s, color = scenario)) +
  scale_color_brewer('Climate change scenario', type = 'div',
                     palette = 5, direction = -1,
                     aesthetics = c('color', 'fill')) +
  labs(x = NULL, y = 'Relative change in speed when moving') +
  scale_y_continuous(labels = labeller_perc) +
  theme(legend.position = 'inside',
        legend.position.inside = c(5/6, 1/6)); p_s
ggsave('figures/speed-local-cc-predictions.png', p_s,
       width = 12, height = 8, dpi = 600, bg = 'white')

p_d <-
  ggplot(cc_proj) +
  facet_wrap(~ species, scales = 'free_y', drop = FALSE,
             labeller = label_parsed) +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(year, ymin = d_lwr, ymax = d_upr, fill = scenario),
              alpha = 0.2) +
  geom_line(aes(year, d, color = scenario)) +
  scale_color_brewer('Climate change scenario', type = 'div',
                     palette = 5, direction = -1,
                     aesthetics = c('color', 'fill')) +
  labs(x = NULL, y = 'Relative change in distance travelled') +
  scale_y_continuous(labels = labeller_perc) +
  theme(legend.position = 'inside',
        legend.position.inside = c(5/6, 1/6)); p_d
ggsave('figures/distance-travelled-local-cc-predictions.png', p_d,
       width = 12, height = 8, dpi = 600, bg = 'white')
