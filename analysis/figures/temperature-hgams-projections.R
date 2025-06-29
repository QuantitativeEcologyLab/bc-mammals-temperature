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
  
  # when predicting: only excluding animal-level random effect and choosing
  # specific values for the rest
  
  # import prediction data for each species in the data's extent ----
  if(file.exists('data/hgam-cc_newd.rds')) {
    cc_newd <- readRDS('data/hgam-cc_newd.rds')
  } else {
    cc_newd <- tibble(
      # list species files (not starting with numbers)
      wp = list.files('data', 'weather-projections-[^0123456789]',
                      full.names = TRUE) %>%
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
        tod_pdt = 12,
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
    # divide by mean of 2025 to find pixel-level relative change since 2025
    group_by(species, long, lat) %>%
    mutate(p = p / mean(p[year == 2025]),
           s = s / mean(s[year == 2025]),
           d = d / mean(d[year == 2025])) %>%
    ungroup() %>%
    # find median and 90% percentile interval of the predicted means
    # not including CIs because averaging them is not straightforward
    group_by(scenario, year, species) %>%
    summarize(p_lwr_05 = quantile(p, 0.05), # bottom 5%
              s_lwr_05 = quantile(s, 0.05),
              d_lwr_05 = quantile(d, 0.05),
              p_median = quantile(p, 0.50), # median
              s_median = quantile(s, 0.50),
              d_median = quantile(d, 0.50),
              p_upr_95 = quantile(p, 0.95), # top 5%
              s_upr_95 = quantile(s, 0.95),
              d_upr_95 = quantile(d, 0.95),
              .groups = 'drop') %>%
    mutate(scenario = case_when(grepl('126', scenario) ~ 'Best scenario (SSP 1-2.6)',
                                grepl('245', scenario) ~ 'Good scenario (SSP 2-4.5)',
                                grepl('370', scenario) ~ 'Bad scenario (SSP 3-7.0)',
                                grepl('585', scenario) ~ 'Worst scenario (SSP 5-8.5)') %>%
             factor(., levels = unique(.)),
           species = gsub(' ', '~', species) %>%
             gsub('~\\(', '\\)~bold\\((', .) %>%
             paste0('bolditalic(', ., ')') %>%
             gsub('\\(boreal\\)', '"\\(boreal\\)"', .) %>%
             gsub('\\(s.~mountain\\)', '"\\(s. mountain\\)"', .) %>%
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
  ggplot(cc_proj, aes(year, p_median, group = scenario)) +
  coord_cartesian(ylim = c(0.75, 1.25)) +
  facet_wrap(~ species, scales = 'fixed', drop = FALSE,
             labeller = label_parsed) +
  geom_hline(aes(yintercept = 1), lty = 'dashed') +
  geom_ribbon(aes(ymin = p_lwr_05, ymax = p_upr_95,
                  fill = factor(scenario, levels = rev(levels(scenario))),
                  color = scenario), linewidth = 0.2, alpha = 0.25) +
  geom_line(color = 'black', lwd = 1.5) +
  geom_line(aes(color = scenario), lwd = 1) +
  scale_color_brewer('Climate change scenario', type = 'div',
                     palette = 5, direction = -1,
                     aesthetics = c('color', 'fill')) +
  labs(x = NULL, y = 'Relative change in annual P(moving)') +
  scale_x_continuous(breaks = c(2025, 2050, 2075, 2100)) +
  theme(legend.position = 'inside',
        legend.position.inside = c(5/6, 1/6)); p_p_mov
ggsave('figures/p-moving-local-cc-predictions.png', p_p_mov,
       width = 10, height = 6.67, dpi = 600, bg = 'white')

p_s <-
  ggplot(cc_proj, aes(year, s_median, group = scenario)) +
  coord_cartesian(ylim = c(NA, 1.021)) +
  facet_wrap(~ species, scales = 'fixed', drop = FALSE,
             labeller = label_parsed) +
  geom_ribbon(aes(ymin = s_lwr_05, ymax = s_upr_95,
                  fill = factor(scenario, levels = rev(levels(scenario))),
                  color = scenario), linewidth = 0.2, alpha = 0.25) +
  geom_hline(aes(yintercept = 1), lty = 'dashed') +
  geom_line(color = 'black', lwd = 1.5) +
  geom_line(aes(color = scenario), lwd = 1) +
  scale_color_brewer('Climate change scenario', type = 'div',
                     palette = 5, direction = -1,
                     aesthetics = c('color', 'fill')) +
  labs(x = NULL, y = 'Relative change in annual speed when moving') +
  scale_x_continuous(breaks = c(2025, 2050, 2075, 2100)) +
  theme(legend.position = 'inside',
        legend.position.inside = c(5/6, 1/6)); p_s
ggsave('figures/speed-local-cc-predictions.png', p_s,
       width = 10, height = 6.67, dpi = 600, bg = 'white')

p_d <-
  ggplot(cc_proj, aes(year, d_median, group = scenario)) +
  coord_cartesian(ylim = c(0.75, 1.255)) +
  facet_wrap(~ species, scales = 'fixed', drop = FALSE,
             labeller = label_parsed) +
  geom_ribbon(aes(ymin = d_lwr_05, ymax = d_upr_95,
                  fill = factor(scenario, levels = rev(levels(scenario))),
                  color = scenario), linewidth = 0.2, alpha = 0.25) +
  geom_hline(aes(yintercept = 1), lty = 'dashed') +
  geom_line(color = 'black', lwd = 1.5) +
  geom_line(aes(color = scenario), lwd = 1) +
  scale_color_brewer('Climate change scenario', type = 'div',
                     palette = 5, direction = -1,
                     aesthetics = c('color', 'fill')) +
  labs(x = NULL, y = 'Relative change in annual distance travelled') +
  scale_x_continuous(breaks = c(2025, 2050, 2075, 2100)) +
  theme(legend.position = 'inside',
        legend.position.inside = c(5/6, 1/6)); p_d
ggsave('figures/distance-travelled-local-cc-predictions.png', p_d,
       width = 10, height = 6.67, dpi = 600, bg = 'white')
