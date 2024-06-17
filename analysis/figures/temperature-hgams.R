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

# import models
m_1 <- readRDS('models/binomial-gam-2024-06-13.rds')
m_2 <- readRDS('models/gammals-gam-2024-06-15.rds')

# find unique species
SPECIES <- unique(m_1$model$species)
N_SPECIES <- length(SPECIES)

# crete new data for predictions for each species
newd_sp <- tibble(animal = 'new animal',
                  species = c(SPECIES),
                  tod_pdt = 12,
                  doy = 0) %>%
  mutate(temp_data = purrr::map(species, \(.s) {
    .temp <- filter(m_1$model, species == .s)$temp_c
    
    tibble(temp_c = seq(quantile(.temp, 0, na.rm = TRUE),
                        quantile(.temp, 1, na.rm = TRUE),
                        length.out = 400))
  })) %>%
  unnest(temp_data)

pred_1_sp <- bind_cols(
  newd_sp,
  predict(m_1, newdata = newd_sp,
          terms = c('s(temp_c,species)', 's(temp_c)'),
          type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_1$family$linkinv(fit),
         lwr = m_1$family$linkinv(fit - 1.96 * se.fit),
         upr = m_1$family$linkinv(fit + 1.96 * se.fit)) %>%
  mutate(species = if_else(species == 'Rangifer tarandus (southern mountain)',
                           'Rangifer tarandus (s. m.)', species))

pred_2_sp <- bind_cols(
  newd_sp,
  predict(m_2, newdata = newd_sp,
          terms = c('s(temp_c,species)', 's(temp_c)'),
          type = 'link', se.fit = TRUE, discrete = FALSE) %>%
    as.data.frame()) %>%
  mutate(mu = exp(fit.1),
         lwr = exp(fit.1 - 1.96 * se.fit.1),
         upr = exp(fit.1 + 1.96 * se.fit.1)) %>%
  mutate(species = if_else(species == 'Rangifer tarandus (southern mountain)',
                           'Rangifer tarandus (s. m.)', species))

# total distance travelled
pred_3_sp <- mutate(pred_2_sp,
                    mu = pred_1_sp$mu * pred_2_sp$mu,
                    lwr = pred_1_sp$lwr * pred_2_sp$lwr,
                    upr = pred_1_sp$upr * pred_2_sp$upr)

# predictions for the average trend across species
newd_ave <- expand_grid(
  animal = 'new animal',
  species = 'Average',
  tod_pdt = 12,
  doy = 0,
  temp_c = seq(quantile(m_1$model$temp_c, 0, na.rm = TRUE),
               quantile(m_1$model$temp_c, 1, na.rm = TRUE),
               length.out = 400))

pred_1_ave <- bind_cols(
  newd_ave,
  predict(m_1, newdata = newd_ave,
          terms = 's(temp_c)',
          type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(mu = m_1$family$linkinv(fit),
         lwr = m_1$family$linkinv(fit - 1.96 * se.fit),
         upr = m_1$family$linkinv(fit + 1.96 * se.fit))

pred_2_ave <- bind_cols(
  newd_ave,
  predict(m_2, newdata = newd_ave,
          terms = 's(temp_c)',
          type = 'link', se.fit = TRUE, discrete = FALSE) %>%
    as.data.frame()) %>%
  mutate(mu = exp(fit.1),
         lwr = exp(fit.1 - 1.96 * se.fit.1),
         upr = exp(fit.1 + 1.96 * se.fit.1))

pred_3_ave <- mutate(pred_2_ave,
                     mu = pred_1_ave$mu * pred_2_ave$mu,
                     lwr = pred_1_ave$lwr * pred_2_ave$lwr,
                     upr = pred_1_ave$upr * pred_2_ave$upr)

# figures ----
p_1 <-
  ggplot() +
  facet_wrap( ~ species, ncol = 4) +
  geom_ribbon(aes(temp_c, ymin = lwr, ymax = upr, fill = species),
              pred_1_sp, alpha = .1) +
  geom_ribbon(aes(temp_c, ymin = lwr, ymax = upr),
              pred_1_ave, alpha = .1) +
  geom_line(aes(temp_c, mu, color = species), pred_1_sp, linewidth =  1) +
  geom_line(aes(temp_c, mu), pred_1_ave, linewidth =  1) +
  scale_color_manual('Species', values = c(PAL),
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)'))+
  scale_y_continuous('P(moving)', expand = c(0, 0), limits = c(0, 1)) +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold.italic')); p_1

p_2 <-
  ggplot() +
  facet_wrap( ~ species, ncol = 4) +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_ribbon(aes(temp_c, ymin = lwr, ymax = upr, fill = species),
              pred_2_sp, alpha = .1) +
  geom_ribbon(aes(temp_c, ymin = lwr, ymax = upr),
              pred_2_ave, alpha = .1) +
  geom_line(aes(temp_c, mu, color = species), pred_2_sp, linewidth =  1) +
  geom_line(aes(temp_c, mu), pred_2_ave, linewidth =  1) +
  scale_color_manual('Species', values = c(PAL),
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)'))+
  scale_y_continuous('Relative change in movement speed') +
  ggtitle('Need to account for uncertainty in scale parameter') +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold.italic')); p_2

p_3 <- ggplot() +
  facet_wrap( ~ species, ncol = 4) +
  geom_ribbon(aes(temp_c, ymin = lwr, ymax = upr, fill = species),
              pred_3_sp, alpha = .1) +
  geom_ribbon(aes(temp_c, ymin = lwr, ymax = upr),
              pred_3_ave, alpha = .1) +
  geom_line(aes(temp_c, mu, color = species), pred_3_sp, linewidth =  1) +
  geom_line(aes(temp_c, mu), pred_3_ave, linewidth =  1) +
  scale_color_manual('Species', values = c(PAL),
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)'))+
  scale_y_continuous('Relative change in distance travelled') +
  ggtitle('Need to account for uncertainty in scale parameter') +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold.italic')); p_3

# save the figures ----
p <- plot_grid(p_1, p_2, p_3, labels = 'AUTO', ncol = 1)

ggsave('figures/hgams-temp_c.png', plot = p, width = 5, height = 7,
       units = 'in', dpi = 600, bg = 'white', scale = 1.5)
