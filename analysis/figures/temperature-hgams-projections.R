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

if(file.exists('data/cc-hgam-projections.rds')) {
  cc_proj <- readRDS('data/cc-hgam-projections.rds')
} else {
  # import models
  m_1 <- readRDS('models/binomial-gam.rds')
  m_2 <- readRDS('models/gamma-gam.rds')
  
  # terms to exclude from the prediction
  SM <- smooths(m_1)
  TERMS <- c('(Intercept)', 'species',
             SM[!(grepl('tod', SM) | grepl('dt', SM) | grepl('animal', SM))])
  
  # import prediction data for each species in the data's extent ----
  if(file.exists('data/hgam-cc_newd.rds')) {
    cc_newd <- readRDS('data/hgam-cc_newd.rds')
  } else {
  cc_newd <- tibble(
    wp = list.files('data', 'weather-projections-', full.names = TRUE) %>%
      map(readRDS)) %>%
    unnest(wp) %>%
    filter(year >= 2025) %>%
    # make sure only the necessary columns are kept
    transmute(
      scenario,
      year,
      animal = m_1$model$animal[1], #' to use `discrete = TRUE`
      species = gsub('boreal', '(boreal)', species) %>%
        gsub('southern mountain', '(s. mountain)', x = .),
      tod_pdt = 0,
      doy = yday(date_decimal(year + (month - 0.5) / 12)),
      temp_c,
      dt = 1,
      weight,
      long, lat)
  saveRDS(cc_newd, 'data/hgam-cc_newd.rds')
  }
  
  cc_proj <- cc_newd %>%
    mutate(.,
           p = predict(m_1, newdata = ., type = 'response', se.fit = FALSE,
                       discrete = TRUE, terms = TERMS),
           s = predict(m_2, newdata = ., type = 'response', se.fit = FALSE,
                       discrete = TRUE, terms = TERMS),
           d = p * s) %>%
    # average across day of year
    group_by(scenario, year, species, long, lat) %>%
    summarize(p = weighted.mean(p, w = weight),
              s = weighted.mean(s, w = weight),
              d = weighted.mean(d, w = weight),
              .groups = 'drop') %>%
    # find median and 90% percentile interval of the predicted means
    # not including CIs because averaging them is not straightforward
    group_by(scenario, year, species) %>%
    summarize(p_lwr_05 = quantile(p, 0.05),
              s_lwr_05 = quantile(s, 0.05),
              d_lwr_05 = quantile(d, 0.05),
              p_median = quantile(p, 0.50),
              s_median = quantile(s, 0.50),
              d_median = quantile(d, 0.50),
              p_upr_95 = quantile(p, 0.95),
              s_upr_95 = quantile(s, 0.95),
              d_upr_95 = quantile(d, 0.95),
              .groups = 'drop') %>%
    # divide by mean of 2025 to find relative change since 2025
    group_by(species) %>%
    mutate(p_ref = mean(p_median[year == 2025]),
           s_ref = mean(s_median[year == 2025]),
           d_ref = mean(d_median[year == 2025])) %>%
    ungroup() %>%
    mutate(scenario = case_when(grepl('126', scenario) ~ 'Best scenario (SSP 1-2.6)',
                                grepl('245', scenario) ~ 'Good scenario (SSP 2-4.5)',
                                grepl('370', scenario) ~ 'Bad scenario (SSP 3-7.0)',
                                grepl('585', scenario) ~ 'Worst scenario (SSP 5-8.5)') %>%
             factor(., levels = unique(.)),
           species = gsub(' ', '~', species) %>%
             gsub('~\\(', '\\)~bold\\((', .) %>%
             paste0('bolditalic(', ., ')') %>%
             factor())
  cc_proj
  
  # sanity check
  if(FALSE) {
    ggplot(cc_proj, aes(year, p_median, group = scenario)) +
      facet_wrap(~ species, scales = 'free') +
      geom_ribbon(aes(year, ymin = p_lwr_05, ymax = p_upr_95,
                      fill = scenario), alpha = 0.2) +
      geom_hline(aes(yintercept = p_ref), lty = 'dashed') +
      geom_line(color = 'black', lwd = 1.5) +
      geom_line(aes(color = scenario), lwd = 1) +
      scale_color_brewer('Climate change scenario', type = 'div',
                         palette = 5, direction = -1,
                         aesthetics = c('color', 'fill'))
    ggplot(cc_proj) +
      facet_wrap(~ species, scales = 'free') +
      geom_ribbon(aes(year, ymin = s_lwr_05, ymax = s_upr_95,
                      fill = scenario), alpha = 0.2) +
      geom_line(aes(year, s_median, color = scenario)) +
      scale_color_brewer('Climate change scenario', type = 'div',
                         palette = 5, direction = -1,
                         aesthetics = c('color', 'fill'))
    ggplot(cc_proj) +
      facet_wrap(~ species, scales = 'free') +
      geom_ribbon(aes(year, ymin = d_lwr_05, ymax = d_upr_95,
                      fill = scenario), alpha = 0.2) +
      geom_line(aes(year, d_median, color = scenario)) +
      scale_color_brewer('Climate change scenario', type = 'div',
                         palette = 5, direction = -1,
                         aesthetics = c('color', 'fill'))
  }
  
  saveRDS(cc_proj, 'data/cc-hgam-projections.rds')
  beepr::beep(2)
}

# make figures ----
p_p_mov <-
  ggplot(cc_proj, aes(year, p_median / p_ref, group = scenario)) +
  facet_wrap(~ species, scales = 'free_y', drop = FALSE,
             labeller = label_parsed) +
  geom_ribbon(aes(ymin = p_lwr_05 / p_ref, ymax = p_upr_95 / p_ref,
                  fill = scenario), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 'dashed') +
  geom_line(color = 'black', lwd = 1.5) +
  geom_line(aes(color = scenario), lwd = 1) +
  scale_color_brewer('Climate change scenario', type = 'div',
                     palette = 5, direction = -1,
                     aesthetics = c('color', 'fill')) +
  labs(x = NULL, y = 'Relative change in P(moving)') +
  theme(legend.position = 'inside',
        legend.position.inside = c(5/6, 1/6)); p_p_mov
ggsave('figures/odds-moving-local-cc-predictions.png', p_p_mov,
       width = 12, height = 8, dpi = 600, bg = 'white')

p_s <-
  ggplot(cc_proj, aes(year, s_median / s_ref, group = scenario)) +
  facet_wrap(~ species, scales = 'free_y', drop = FALSE,
             labeller = label_parsed) +
  geom_ribbon(aes(ymin = s_lwr_05 / s_ref, ymax = s_upr_95 / s_ref,
                  fill = scenario), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 'dashed') +
  geom_line(color = 'black', lwd = 1.5) +
  geom_line(aes(color = scenario), lwd = 1) +
  scale_color_brewer('Climate change scenario', type = 'div',
                     palette = 5, direction = -1,
                     aesthetics = c('color', 'fill')) +
  labs(x = NULL, y = 'Relative change in speed when moving') +
  theme(legend.position = 'inside',
        legend.position.inside = c(5/6, 1/6)); p_s
ggsave('figures/speed-local-cc-predictions.png', p_s,
       width = 12, height = 8, dpi = 600, bg = 'white')

p_d <-
  ggplot(cc_proj, aes(year, d_median / d_ref, group = scenario)) +
  facet_wrap(~ species, scales = 'free_y', drop = FALSE,
             labeller = label_parsed) +
  geom_ribbon(aes(ymin = d_lwr_05 / d_ref, ymax = d_upr_95 / d_ref,
                  fill = scenario), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 'dashed') +
  geom_line(color = 'black', lwd = 1.5) +
  geom_line(aes(color = scenario), lwd = 1) +
  scale_color_brewer('Climate change scenario', type = 'div',
                     palette = 5, direction = -1,
                     aesthetics = c('color', 'fill')) +
  labs(x = NULL, y = 'Relative change in distance travelled') +
  theme(legend.position = 'inside',
        legend.position.inside = c(5/6, 1/6)); p_d
ggsave('figures/distance-travelled-local-cc-predictions.png', p_d,
       width = 12, height = 8, dpi = 600, bg = 'white')

# for poster
ggsave('figures/2024-ubco-grad-symposium/distance-travelled-local-cc-predictions.png',
       p_d, width = 17.5, height = 9.5, dpi = 300, bg = 'white', scale = 0.75)
