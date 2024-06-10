library('dplyr')     # for data wrangling
library('lubridate') # for smother date wrangling
library('mgcv')      # for Generalized Additive Models
library('ggplot2')   # for fancy plots
library('khroma')    # for colorblind-friendly color palettes
library('sf')        # for spatial objects
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids

theme_set(theme_get() +
            theme(text = element_text(size = 15)))

# using boreal caribou
readRDS('data/tracking-data/all-tracking-data-cleaned-2024-02-22-13-49.rds') %>%
  group_by(dataset_name) %>%
  tidyr::unnest(tel) %>%
  summarise(prop = mean(! is.na(temperature)),
            total = sum(! is.na(temperature))) %>%
  arrange(desc(total))

# import GAMs
m_1 <- readRDS('models/binomial-gam-2024-03-06-boreal-caribou.rds')
m_2 <- readRDS('models/gamma-gam-2024-03-06-boreal-caribou.rds')

# P(moving) ----
#' check ranges of `temperature`
m_1$model %>%
  summarise(min = min(temperature),
            max = max(temperature))

p_change <- function(.data) {
  inv_logit <- m_1$family$linkinv
  
  mutate(.data,
         mu = (inv_logit(fit) - 0.5) * 100,
         lwr = (inv_logit(fit - 1.96 * se.fit) - 0.5) * 100,
         upr = (inv_logit(fit + 1.96 * se.fit) - 0.5) * 100) %>%
    return()
}

#' marginal effect of `temperature`
tibble(temperature = seq(-30, 40, by = 0.1),
       time_of_day = 12,
       doy = 100,
       animal = 'new animal') %>%
  bind_cols(.,
            predict(m_1, newdata = ., terms = c('s(temperature)'),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  p_change() %>%
  ggplot() +
  geom_hline(yintercept = 0, col = 'grey') +
  geom_ribbon(aes(temperature, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(temperature, mu), color = '#EE6677', linewidth =  1) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)')) +
  scale_y_continuous('Change in P(moving) (%)', limits = c(-5, 5))

ggsave('figures/actws-2024-jasper/p-moving-temperature-caribou.png',
       width = 10, height = 6, dpi = 600)

#' effect of `time_of_day`
tibble(time_of_day = seq(0, 24, by = 0.001),
       doy = 0,
       temperature = 0,
       animal = 'new animal') %>%
  bind_cols(.,
            predict(m_1, newdata = ., terms = c('s(time_of_day)'),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  p_change() %>%
  ggplot() +
  geom_hline(yintercept = 0, col = 'grey') +
  geom_ribbon(aes(time_of_day, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(time_of_day, mu), color = '#EE6677', linewidth =  1) +
  scale_color_bright(name = 'Species') +
  scale_fill_bright(name = 'Species') +
  scale_x_continuous('Time of day', expand = c(0, 0),
                     breaks = c(6, 12, 18),
                     labels = paste0(c(6, 12, 18), ':00')) +
  scale_y_continuous('Change in P(moving) (%)', limits = c(-5, 5))

ggsave('figures/actws-2024-jasper/p-moving-tod-caribou.png',
       width = 10, height = 6, dpi = 600)

#' effect of `doy`
tibble(time_of_day = 12,
       date = seq(as.Date('2024-01-01'), as.Date('2024-12-31'), by = 1),
       doy = yday(date),
       temperature = 0,
       animal = 'new animal') %>%
  bind_cols(.,
            predict(m_1, newdata = ., terms = c('s(doy)'),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  p_change() %>%
  ggplot() +
  geom_hline(yintercept = 0, col = 'grey') +
  geom_ribbon(aes(date, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(date, mu), color = '#EE6677', linewidth =  1) +
  scale_color_bright(name = 'Species') +
  scale_fill_bright(name = 'Species') +
  scale_x_date(expand = c(0, 0)) +
  scale_y_continuous() +
  labs(x = 'Day of year', y = 'Change in P(moving) (%)')

ggsave('figures/actws-2024-jasper/p-moving-doy-caribou.png',
       width = 10, height = 6, dpi = 600)

#' effects of `temperature` on different `doy`s
expand.grid(temperature = seq(-30, 40, by = 0.1),
            time_of_day = 12,
            date = as.Date(c('2024-03-19', '2024-06-20',
                             '2024-09-22', '2024-12-21')),
            animal = 'new animal') %>%
  mutate(doy = yday(date)) %>%
  filter(! exclude.too.far(g1 = temperature, g2 = doy,
                           d1 = m_1$model$temperature, d2 = m_1$model$doy,
                           dist = 0.01)) %>%
  bind_cols(.,
            predict(m_1, newdata = .,
                    terms = c('s(temperature)', 's(doy)',
                              'ti(doy,temperature)'),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(day = format(date, '%B %d') %>%
           factor(., levels = unique(.))) %>%
  p_change() %>%
  ggplot() +
  facet_wrap(~ day) +
  geom_hline(yintercept = 0, col = 'grey') +
  geom_ribbon(aes(temperature, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(temperature, mu), color = '#EE6677', linewidth =  1) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)')) +
  scale_y_continuous('Change in P(moving) (%)')

ggsave('figures/actws-2024-jasper/p-moving-temperature-doy-caribou.png',
       width = 10, height = 6, dpi = 600)

#' effect of `time_of_day` on different `doy`s
expand.grid(temperature = 0,
            time_of_day = seq(0, 24, length.out = 400),
            date = as.Date(c('2024-03-19', '2024-06-20',
                             '2024-09-22', '2024-12-21')),
            animal = 'new animal') %>%
  mutate(doy = yday(date)) %>%
  bind_cols(.,
            predict(m_1, newdata = ., terms = c('s(doy)', 's(time_of_day)',
                                                'ti(doy,time_of_day)'),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(day = format(date, '%B %d') %>%
           factor(., levels = unique(.))) %>%
  p_change() %>%
  ggplot() +
  facet_wrap(~ day) +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_ribbon(aes(time_of_day, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(time_of_day, mu), color = '#EE6677', linewidth =  1) +
  labs(x = 'Time of day', y = 'Change in P(moving) (%)') +
  scale_x_continuous(expand = c(0, 0))

ggsave('figures/actws-2024-jasper/p-moving-tod-doy-caribou.png',
       width = 10, height = 6, dpi = 600)

# speed given animal is moving ----
if(FALSE) {
  # check predicted vs fitted
  tibble(speed = m_2$model$sld,#speed_est,
         mu = m_2$fitted.values,#[, 1],
         species = m_2$model$species) %>%
    # group_by(species) %>%
    # summarise(speed = mean(speed), mu = mean(mu)) %>%
    ggplot() +
    facet_wrap(~ species) +
    # geom_point(aes(speed, mu, col = species), d_2, alpha = 0.5) +
    geom_point(aes(speed, mu)) +
    geom_point(aes(speed, mu, col = species)) +
    scale_color_bright(name = 'Species') +
    geom_abline(slope = 1, intercept = 0) +
    labs(x = 'Observed speeds', y = 'Estimated speeds')
}

#' check ranges of `temperature`
m_2$model %>%
  summarise(min = min(temperature),
            max = max(temperature))

speed_change <- function(.data) {
  .data %>%
    mutate(mu = m_2$family$linkinv(fit),
           lwr = m_2$family$linkinv(fit - 1.96 * se.fit),
           upr = m_2$family$linkinv(fit + 1.96 * se.fit)) %>%
    mutate(mu = (mu - 1) * 100,
           lwr = (lwr - 1) * 100,
           upr = (upr - 1) * 100) %>%
    return()
}

#' marginal effect of `temperature`
tibble(temperature = seq(-30, 40, by = 0.1),
       time_of_day = 12,
       doy = 0,
       animal = 'new animal') %>%
  bind_cols(.,
            predict(m_2, newdata = .,
                    terms = c('s(temperature)'),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(fit = fit - fit[which(temperature == 0)]) %>%
  speed_change() %>%
  ggplot() +
  geom_hline(yintercept = 0, color = 'grey') +
  # geom_vline(xintercept = 0, color = 'grey', lty = 'dashed') +
  geom_ribbon(aes(temperature, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(temperature, mu), color = '#EE6677', linewidth =  1) +
  scale_x_continuous(paste0('Temperature (', '\U00B0', 'C)'))+
  ylab('Change in speed (%)')

ggsave('figures/actws-2024-jasper/speed-temperature-caribou.png',
       width = 10, height = 6, dpi = 600)

#' effect of `time_of_day`
tibble(time_of_day = seq(0, 24, by = 0.001),
       doy = 0,
       temperature = 0,
       animal = 'new animal') %>%
  bind_cols(.,
            predict(m_2, newdata = .,
                    terms = c('s(time_of_day)'),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  speed_change() %>%
  ggplot() +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_ribbon(aes(time_of_day, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(time_of_day, mu), color = '#EE6677', linewidth =  1) +
  scale_x_continuous('Time of day', expand = c(0, 0),
                     breaks = c(0, 6, 12, 18),
                     labels = paste0(c(0, 6, 12, 18), ':00')) +
  ylab('Change in speed (%)')

ggsave('figures/actws-2024-jasper/speed-tod-caribou.png',
       width = 10, height = 6, dpi = 600)

#' effect of `doy`
tibble(time_of_day = 12,
       date = seq(as.Date('2024-01-01'), as.Date('2024-12-31'), by = 1),
       doy = yday(date),
       temperature = 0,
       animal = 'new animal') %>%
  bind_cols(.,
            predict(m_2, newdata = ., terms = c('s(doy)'),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  speed_change() %>%
  ggplot() +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_ribbon(aes(date, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(date, mu), color = '#EE6677', linewidth =  1) +
  scale_color_bright(name = 'Species') +
  scale_fill_bright(name = 'Species') +
  scale_x_date(expand = c(0, 0)) +
  labs(x = 'Day of year', y = 'Change in speed (%)')

ggsave('figures/actws-2024-jasper/speed-doy-caribou.png',
       width = 10, height = 6, dpi = 600)

#' effects of `time_of_day` on different `doy`s
expand.grid(temperature = 0,
            time_of_day = seq(0, 24, length.out = 400),
            date = as.Date(c('2024-03-19', '2024-06-20',
                             '2024-09-22', '2024-12-21')),
            animal = 'new animal') %>%
  mutate(doy = yday(date)) %>%
  bind_cols(.,
            predict(m_2, newdata = .,
                    terms = c('s(doy)', 's(time_of_day)',
                              'ti(doy,time_of_day)'),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(day = format(date, '%B %d') %>%
           factor(., levels = unique(.))) %>%
  speed_change() %>%
  ggplot() +
  facet_wrap(~ date) +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_ribbon(aes(time_of_day, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(time_of_day, mu), color = '#EE6677', linewidth =  1) +
  labs(x = 'Time of day', y = 'Change in speed (%)')

ggsave('figures/actws-2024-jasper/speed-tod-doy-caribou.png',
       width = 10, height = 6, dpi = 600)

#' effects of `temperature` on different DOYs
expand.grid(temperature = seq(-30, 40, by = 0.1),
            time_of_day = 12,
            date = as.Date(c('2024-03-19', '2024-06-20',
                             '2024-09-22', '2024-12-21')),
            animal = 'new animal') %>%
  mutate(doy = yday(date)) %>%
  filter(! exclude.too.far(g1 = temperature, g2 = doy,
                           d1 = m_1$model$temperature, d2 = m_1$model$doy,
                           dist = 0.01)) %>%
  bind_cols(.,
            predict(m_2, newdata = .,
                    terms = c('s(temperature)', 's(doy)',
                              'ti(doy,temperature)'),
                    type = 'link', se.fit = TRUE, discrete = FALSE)) %>%
  mutate(day = format(date, '%B %d') %>%
           factor(., levels = unique(.))) %>%
  speed_change() %>%
  ggplot() +
  facet_wrap(~ day) +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_ribbon(aes(temperature, ymin = lwr, ymax = upr), fill = '#EE6677',
              alpha = 0.2) +
  geom_line(aes(temperature, mu), color = '#EE6677', linewidth =  1) +
  labs(x = paste0('Temperature (', '\U00B0', 'C)'),
       y = 'Change in speed (%)')

ggsave('figures/actws-2024-jasper/speed-temperature-doy-caribou.png',
       width = 10, height = 6, dpi = 600)

# future predictions ----
locs <- readRDS('data/tracking-data/all-tracking-data-cleaned-2024-02-22-13-49.rds') %>%
  filter(dataset_name == 'Rangifer_tarandus_boreal') %>%
  filter(! grepl('SCEK014', animal)) %>% # drop outlier caribou
  tidyr::unnest(tel) %>%
  select(animal, location.long, location.lat) %>%
  st_as_sf(coords = c('location.long', 'location.lat'))

bounds <- locs %>%
  group_by(animal) %>%
  st_set_crs('+proj=longlat') %>%
  as_Spatial() %>%
  adehabitatHR::mcp(percent = 100) %>%
  st_as_sfc() %>%
  st_as_sf() %>%
  st_union() %>%
  st_as_sf()

if(FALSE) { # to check coordinates
  plot(st_geometry(locs))
  plot(bounds, add = TRUE, col = '#FF000030')
}

# monthly predictions
preds_cc_month <-
  readRDS('data/climate-yearly-projections-2024-03-03.rds') %>%
  filter(longitude > st_bbox(bounds)['xmin'],
         longitude < st_bbox(bounds)['xmax'],
         latitude > st_bbox(bounds)['ymin'],
         latitude < st_bbox(bounds)['ymax']) %>%
  filter(year >= 2024) %>%
  mutate(temperature = mean_temperature,
         doy = 0,
         time_of_day = 0,
         animal = 'new animal') %>%
  mutate(p_mov = predict(object = m_1, newdata = ., se.fit = FALSE,
                         type = 'response', discrete = FALSE,
                         # can't include time of day
                         terms = c('s(temperature)', 's(doy)',
                                   'ti(doy,temperature)')),
         speed = predict(object = m_2, newdata = ., se.fit = FALSE,
                         type = 'response', discrete = FALSE,
                         terms = c('s(temperature)', 's(doy)',
                                   'ti(doy,temperature)')),
         distance = p_mov * speed) %>%
  # take the average over the whole area
  group_by(scenario, year, month) %>%
  summarize(p_mov = mean(p_mov, na.rm = TRUE),
            speed = mean(speed, na.rm = TRUE),
            distance = mean(distance, na.rm = TRUE),
            .groups = 'drop') %>%
  # make scenarios more understandable
  mutate(scenario = case_when(
    scenario == '8GCMs_ensemble_ssp126' ~ 'Best',
    scenario == '8GCMs_ensemble_ssp245' ~ 'Good',
    scenario == '8GCMs_ensemble_ssp370' ~ 'Bad',
    scenario == '8GCMs_ensemble_ssp585' ~ 'Worst'),
    scenario = factor(scenario,
                      levels = c('Best', 'Good', 'Bad', 'Worst')))

# interesting but too weak and complex for the presentation
ggplot(preds_cc_month,
       aes(year, month, fill = (p_mov - 0.5) * 100)) +
  geom_raster() +
  scale_fill_distiller('Change in P(moving) (%)', type = 'div',
                       palette = 4, limits = c(-8, 8)) +
  scale_x_continuous(NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, expand = c(0, 0), breaks = c(1, 4, 7, 10),
                     labels = c('Jan', 'Apr', 'Jul', 'Oct')) +
  theme(legend.position = 'top')

# average within years
preds_cc <-
  preds_cc_month %>%
  group_by(scenario, year) %>%
  summarize(p_mov = mean(p_mov, na.rm = TRUE),
            speed = mean(speed, na.rm = TRUE),
            distance = sum(distance, na.rm = TRUE),
            .groups = 'drop') %>%
  group_by(scenario) %>% # convert to relative to make figure look cleaner
  arrange(year) %>%
  mutate(p_mov = p_mov / first(p_mov),
         speed = speed / first(speed),
         distance = distance / first(distance),
         # convert to percentage
         p_mov = (p_mov - 1) * 100,
         speed = (speed - 1) * 100,
         distance = (distance - 1) * 100)

# P(moving)
p_cc_mov <-
  ggplot(preds_cc) +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_point(aes(year, p_mov, color = scenario), alpha = 0.5) +
  geom_smooth(aes(year, p_mov, color = scenario, fill = scenario),
              method = 'gam', formula = y ~ s(x, k = 5)) +
  scale_fill_brewer('Climate change scenario',
                    type = 'div', palette = 5, direction = -1,
                    aesthetics = c('color', 'fill')) +
  scale_y_continuous('Change in P(moving) (%)', limits = c(-5, 5)) +
  scale_x_continuous(NULL, breaks = c(2024, 2040, 2060, 2080, 2100)) +
  theme(legend.position = 'top'); p_cc_mov

ggsave('figures/actws-2024-jasper/p-moving-caribou-cc-preds.png', p_cc_mov,
       width = 10, height = 6, dpi = 600)

# average speed trends by year
p_cc_speed <- 
  ggplot(preds_cc) +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_point(aes(year, speed, color = scenario), alpha = 0.5) +
  geom_smooth(aes(year, speed, color = scenario, fill = scenario),
              method = 'gam', formula = y ~ s(x, k = 5)) +
  scale_fill_brewer('Climate change scenario',
                    type = 'div', palette = 5, direction = -1,
                    aesthetics = c('color', 'fill')) +
  scale_x_continuous(breaks = c(2024, 2040, 2060, 2080, 2100)) +
  labs(x = NULL, y = 'Change in speed (%)') +
  ylim(c(-5, NA)) +
  theme(legend.position = 'top'); p_cc_speed

ggsave('figures/actws-2024-jasper/speed-caribou-cc-preds.png',
       p_cc_speed, width = 10, height = 6, dpi = 600)

# average speed trends by month for four years
preds_cc_month %>%
  filter(year %in% c(2040, 2060, 2080, 2100)) %>%
  mutate(p_mov = p_mov / first(p_mov),
         speed = speed / first(speed),
         distance = distance / first(distance),
         # convert to percentage
         speed = (speed - 1) * 100) %>%
  ggplot() +
  facet_wrap(~ year) +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_point(aes(month, speed, color = scenario), alpha = 0.5) +
  geom_smooth(aes(month, speed, color = scenario, fill = scenario),
              method = 'gam', formula = y ~ s(x, k = 5)) +
  scale_fill_brewer('Climate change scenario',
                    type = 'div', palette = 5, direction = -1,
                    aesthetics = c('color', 'fill')) +
  scale_x_continuous(breaks = c(2024, 2040, 2060, 2080, 2100)) +
  labs(x = NULL, y = 'Change in speed (%)') +
  theme(legend.position = 'top')

ggsave('figures/actws-2024-jasper/speed-caribou-cc-preds-doy.png',
       p_cc_speed, width = 10, height = 6, dpi = 600)

# total movement
p_cc_distance <- 
  ggplot(preds_cc) +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_point(aes(year, distance, color = scenario), alpha = 0.5) +
  geom_smooth(aes(year, distance, color = scenario, fill = scenario),
              method = 'gam', formula = y ~ s(x, k = 5)) +
  scale_fill_brewer('Climate change scenario',
                    type = 'div', palette = 5, direction = -1,
                    aesthetics = c('color', 'fill')) +
  scale_x_continuous(breaks = c(2024, 2040, 2060, 2080, 2100)) +
  labs(x = NULL, y = 'Change in yearly movement (%)') +
  ylim(c(-5, NA)) +
  theme(legend.position = 'top'); p_cc_distance

ggsave('figures/actws-2024-jasper/distance-caribou-cc-preds.png',
       p_cc_distance, width = 10, height = 6, dpi = 600)
