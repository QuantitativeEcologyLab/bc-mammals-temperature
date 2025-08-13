library('dplyr')     #' for data wrangling
library('tidyr')     #' for data wrangling
library('lubridate') #' for smoother date wrangling
library('mgcv')      #' for Generalized Additive Models
library('ggplot2')   #' for fancy plots
library('gratia')    #' for ggplot GAM graphics
library('terra')     #' for working with rasters
library('ggridges')  #' for ridgeline plots
library('khroma')    # for colorblind-friendly color palettes
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids
source('data/bc-shapefile.R')

# import temperature data ----
if(! file.exists('data/ecmwf-era5-2m-temperature/aggregated-data.rds')) {
  NCHAR <- nchar('data/ecmwf-era5-2m-temperature/epsg-4326/ecmwf-era5-2m-temp-X')
  
  dem <- rast('data/resource-rasters/bc-dem-z6.tif')
  
  d <- rast(list.files('data/ecmwf-era5-2m-temperature/epsg-4326',
                       full.names = TRUE)) %>%
    crop(bc_unproj) %>%
    mask(bc_unproj) %>%
    as.data.frame(xy = TRUE) %>%
    pivot_longer(- c(x, y), names_to = 'name', values_to = 'temp_k') %>%
    mutate(year = as.numeric(substr(name, nchar('ecmwf-era5-2m-temp-X'),
                                    nchar('ecmwf-era5-2m-temp-XXXX'))),
           # get raster number hour (same as hour number using UTC)
           h = as.numeric(substr(name,
                                 nchar('ecmwf-era5-2m-temp-1998-epsg-4326_X'),
                                 nchar(name))),
           date = date_decimal(year + h / 365 / 24),
           temp_c = temp_k - 273.15) %>%
    select(x, y, date, temp_c) %>%
    mutate(year = year(date),
           month = month(date)) %>%
    group_by(x, y, year, month) %>%
    summarize(monthly_mean = mean(temp_c, na.rm = TRUE),
              monthly_var = var(temp_c, na.rm = TRUE),
              n = n(),
              .groups = 'drop') %>%
    rename(long = x, lat = y) %>%
    mutate(elev_m = extract(dem, tibble(long, lat))[, 2]) %>%
    filter(! is.na(elev_m)) %>%
    filter(! is.na(monthly_var)) %>%
    mutate(sqrt_prec = sqrt(1 / monthly_var))
  
  saveRDS(d, 'data/ecmwf-era5-2m-temperature/aggregated-data.rds')
} else {
  d <- readRDS('data/ecmwf-era5-2m-temperature/aggregated-data.rds')
}

range(d$monthly_mean)

d

#' no `NA` variances; all values have at least 672 samples
mean(is.na(d$monthly_mean))
mean(is.na(d$monthly_var))
quantile(d$n, probs = c(0, 0.003, 0.005, 0.25, 0.5, 1))

# Gamma model for the variance around the mean, given the mean ----
#' Assuming `temp_c ~ N(mu, sigma2)` implies sigma^2 is independent of mu,
#' but extreme temperatures (`abs(E(T)) > 15`) are more variable than
#' temperatures near 0, so we should include `s(monthly_mean)`
#' not including `s(year)` because it only contributes to ~2% of dev.expl,
#' and the predictions for 2100 would likely not be realistic
#' fits in ~5 minutes
#' inverse gaussian distribution fits much worse than Gamma
if(! file.exists('models/bc-temperature-var-gamma-2024-06-14.rds')) {
  m_gamma <-
    bam(monthly_var ~
          s(monthly_mean, bs = 'tp', k = 8) +
          s(elev_m, bs = 'tp', k = 8) +
          s(month, bs = 'cc', k = 10) +
          s(long, lat, bs = 'ds', k = 500) +
          ti(long, lat, month, bs = c('ds', 'cc'), d = c(2, 1),
             k = c(200, 10)),
        family = Gamma(link = 'log'),
        data = d,
        weights = sqrt(n),
        method = 'fREML',
        discrete = TRUE,
        knots = list(month = c(0.5, 12.5)),
        control = gam.control(trace = TRUE))
  
  saveRDS(m_gamma,
          paste0('models/bc-temperature-var-gamma-', Sys.Date(), '.rds'))
  
  # great fit; ignore q-q line
  layout(matrix(1:4, ncol = 2))
  gam.check(m_gamma, pch = 19, col = '#00000020')
  layout(1)
  
  summary(m_gamma) # dev.expl is quite high
  plot(m_gamma, pages = 1, scheme = c(1, 1, 1, 2, 2), too.far = 0.02)
  plot(m_gamma, pages = 1, scheme = c(1, 1, 1, 2, 2), too.far = 0.02,
       scale = 0)
  plot(m_gamma, select = 1, scale = 0, n = 300)
  plot(m_gamma, select = 3, scale = 0, n = 300, xlim = c(0.5, 12.5))
  
  # but variance in residuals varies between months
  d %>%
    mutate(month = factor(month, levels = 12:1),
           e = resid(m_gamma)) %>%
    ggplot() +
    geom_density_ridges(aes(x = e, y = month, fill = month),
                        alpha = 0.75, scale = 2) +
    geom_vline(xintercept = 0, color = 'black') +
    xlab('Deviance residuals') +
    scale_y_discrete(expand = c(0, 1.2)) +
    scale_fill_manual('Month', values = color('BuRd')(6)[c(1:6, 6:1)])
  
  plot(factor(m_gamma$model$month), resid(m_gamma, type = 'deviance'),
       xlab = NULL, ylab = 'Residuals')
  abline(h = 0, lwd = 1.5, col = 'red3')
  
  # checking spatial distributions
  gamma_resids <- mutate(d, e_gamma = resid(m_gamma)) %>%
    group_by(long, lat) %>%
    summarize(mean_e_gamma = mean(e_gamma),
              var_e_gamma = var(e_gamma),
              .groups = 'drop')
  
  # very little spatial bias
  # spatial bias along the coast should average to zero  
  ggplot(gamma_resids, aes(long, lat, fill = mean_e_gamma)) +
    geom_raster() +
    geom_sf(data = bc_unproj, inherit.aes = FALSE, fill = 'transparent') +
    scale_fill_viridis_c(expression(paste('E(e) (\U00B0', C^2, ')')),
                         option = 'D', limits = c(-0.5, 0.5), direction = -1) +
    theme(legend.position = c(0.9, 0.8))
  
  # some spatial structure, but coastal effects are ok
  ggplot(gamma_resids, aes(long, lat, fill = var_e_gamma)) +
    geom_raster() +
    geom_sf(data = bc_unproj, inherit.aes = FALSE, fill = 'transparent') +
    scale_fill_viridis_c(expression(paste('Var(e) (\U00B0', C^4, ')')),
                         option = 'A', limits = c(0, NA), direction = -1) +
    theme(legend.position = c(0.9, 0.8))
} else {
  m_gamma <- readRDS('models/bc-temperature-var-gamma-2024-06-14.rds')
}

# location-scale Gamma model for Var(T|E(T|X),X) ----
#' account for seasonal and spatial heteroskedasticity with considerations
#' from above
#' fits in < 4-5 days
if(! file.exists('models/bc-temperature-var-gammals-2024-06-16.rds')) {
  tictoc::tic()
  m_gammals <- gam(list(
    monthly_var ~
      s(monthly_mean, bs = 'cr', k = 10) + # using cr to reduce fitting time
      s(month, bs = 'cc', k = 10) +
      s(long, lat, bs = 'ds', k = 500) +
      ti(long, lat, month, bs = c('ds', 'cc'), d = c(2, 1), k = c(100, 10)),
    ~ s(monthly_mean, bs = 'cr', k = 10) +
      s(month, bs = 'cc', k = 10) +
      s(long, lat, bs = 'ds', k = 100)),
    family = gammals(),
    data = d,
    method = 'REML',
    knots = list(month = c(0.5, 12.5)),
    control = gam.control(trace = TRUE))
  tictoc::toc(); beepr::beep()
  
  saveRDS(m_gammals,
          paste0('models/bc-temperature-var-gammals-', Sys.Date(),
                 '.rds'))
  
  # good fit, but but residuals have heavy tails
  layout(matrix(1:4, ncol = 2))
  gam.check(m_gammals, pch = 19, col = '#00000020')
  layout(1)
  
} else {
  m_gammals <- readRDS('models/bc-temperature-var-gammals-2024-06-16.rds')
}

summary(m_gammals) # dev.expl is quite high
plot(m_gammals, pages = 1, scheme = c(1, 1, 2, 2, 1, 1, 2), too.far = 0.02)

# but variance in residuals varies between months
d %>%
  mutate(month = factor(month, levels = 12:1),
         e = resid(m_gammals)) %>%
  ggplot() +
  geom_density_ridges(aes(x = e, y = month, fill = month),
                      alpha = 0.75, scale = 2) +
  geom_vline(xintercept = 0, color = 'black') +
  xlab('Deviance residuals') +
  scale_fill_manual('Month', values = color('BuRd')(6)[c(1:6, 6:1)])

plot(factor(m_gammals$model$month), resid(m_gammals, type = 'deviance'),
     xlab = NULL, ylab = 'Residuals')
abline(h = 0, lwd = 1.5, col = 'red3')

# still, AIC is substantially better for the gammals model
AIC(m_gamma, m_gammals) %>%
  mutate(delta_AIC = AIC - min(AIC))

# get SD(temperature) based on SSP monthly averages ----
preds <-
  # import monthly averagaes based on climate change projections
  readRDS('data/climate-yearly-projections-2024-06-11.rds') %>%
  rename(elev_m = elevation,
         monthly_mean = mean_temperature,
         lat = latitude,
         long = longitude) %>%
  mutate(monthly_mean = if_else(monthly_mean == -9999, NA_real_, monthly_mean)) %>%
  select(! c(min_temperature, max_temperature)) %>%
  mutate(sd = sqrt(predict(m_gammals, newdata = ., type = 'response')[, 1]))

range(preds$monthly_mean, na.rm = TRUE)
range(preds$sd, na.rm = TRUE)
hist(rnorm(nrow(preds), mean = preds$monthly_mean, sd = preds$sd),
     main = expression('Gaussian temperatures (n = 1 for each {'~
                         hat(mu)~','~hat(sigma)~'})'))

# save the file with added mean and standard deviation
saveRDS(preds,
        paste0('data/climate-yearly-projections-mean-sd-', Sys.Date(),
               '.rds'))

# simulate weather for each month ----
if(! exists('preds')) {
  preds <- readRDS('data/climate-yearly-projections-mean-sd-2024-06-26.rds')
}

weather_proj <-
  preds %>%
  transmute(scenario, year, month, monthly_mean, long, lat,
            monthly_sd = sd) %>%
  # add quantiles for standard normal distribution
  mutate(., qs = list(tibble(q = qnorm(seq(0.1, 0.9, by = 0.1))))) %>%
  unnest(qs) %>%
  mutate(
    # convert Z quantiles to the monthly quantiles
    temp_c = monthly_mean + monthly_sd * q,
    # weigh each quantile by the probability density
    weight = dnorm(x = temp_c, mean = monthly_mean, sd = monthly_sd)) %>%
  # ensure each month has a weight of 1
  group_by(scenario, long, lat, year, month) %>%
  mutate(weight = weight / sum(weight)) %>%
  ungroup()

# check weights for a few temperatures
weather_proj %>%
  slice(1:(9 * 3)) %>%
  ggplot(aes(temp_c, weight)) +
  facet_wrap(~ monthly_mean, ncol = 1) +
  geom_smooth() +
  geom_point()

saveRDS(weather_proj, 'data/weather-projections.rds')

# simulated weather for only 2025 and 2100 ----
# get SD(temperature) based on SSP monthly averages
preds_2025_2100 <-
  # import monthly averages based on CC projections for 2025 and 2100
  readRDS('data/climate-yearly-projections-2025-2100-only-2025-08-13.rds') %>%
  filter(mean_temperature != -9999) %>%
  rename(elev_m = elevation,
         monthly_mean = mean_temperature,
         lat = latitude,
         long = longitude) %>%
  select(! c(min_temperature, max_temperature)) %>%
  mutate(sd = sqrt(predict(m_gammals, newdata = ., type = 'response')[, 1]))

range(preds_2025_2100$monthly_mean, na.rm = TRUE)
range(preds_2025_2100$sd, na.rm = TRUE)
hist(rnorm(nrow(preds_2025_2100),
           mean = preds_2025_2100$monthly_mean, sd = preds_2025_2100$sd),
     main = expression('Gaussian temperatures (n = 1 for each {'~
                         hat(mu)~','~hat(sigma)~'})'))

saveRDS(preds_2025_2100,
        paste0('data/climate-yearly-projections-mean-sd-2025-2100-only-',
               Sys.Date(), '.rds'))

# simulate weather for each month ----
if(! exists('preds_2025_2100')) {
  preds_2025_2100 <- readRDS('data/climate-yearly-projections-mean-sd-2025-2100-only-2025-01-21.rds')
}

weather_proj_2025_2100 <-
  preds_2025_2100 %>%
  transmute(scenario, year, month, dec_date, long, lat, monthly_mean,
            monthly_sd = sd) %>%
    # add quantiles for standard normal distribution
  mutate(., qs = list(tibble(q = qnorm(seq(0.1, 0.9, by = 0.1))))) %>%
  unnest(qs) %>%
  mutate(
    # convert Z quantiles to the monthly quantiles
    temp_c = monthly_mean + monthly_sd * q,
    # weigh each quantile by the probability density
    weight = dnorm(x = temp_c, mean = monthly_mean, sd = monthly_sd)) %>%
  # ensure each month has a weight of 1
  group_by(scenario, long, lat, year, month) %>%
  mutate(weight = weight / sum(weight)) %>%
  ungroup()

# check weights for a few temperatures
weather_proj_2025_2100 %>%
  slice(1:(9 * 3)) %>%
  ggplot(aes(temp_c, weight)) +
  facet_wrap(~ monthly_mean, ncol = 1) +
  geom_smooth() +
  geom_point()

saveRDS(weather_proj_2025_2100,
        'data/weather-projections-2025-2100-only.rds')
