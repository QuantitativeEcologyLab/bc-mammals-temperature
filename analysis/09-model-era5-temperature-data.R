library('dplyr')     #' for data wrangling
library('tidyr')     #' for data wrangling
library('lubridate') #' for smother date wrangling
library('mgcv')      #' for Generalized Additive Models
library('ggplot2')   #' for fancy plots
library('gratia')    #' for ggplot GAM graphics
library('terra')     #' for working with rasters
library('ggridges')  #' for ridgeline plots
library('khroma')    # for colorblind-friendly color palettes
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids
source('data/bc-shapefile.R')

# import temperature data ----
SAMPLE <- seq(1, 365 * 24 * 2, by = 50) # sample of rasters for test

if(FALSE) {
  NCHAR <- nchar('data/ecmwf-era5-2m-temperature/epsg-4326/ecmwf-era5-2m-temp-X')
  
  d <- rast(list.files('data/ecmwf-era5-2m-temperature/epsg-4326',
                       full.names = TRUE)) %>%
    crop(bc_unproj) %>%
    mask(bc_unproj) %>%
    as.data.frame(xy = TRUE) %>%
    # mutate(elev_m = terra::extract(dem, tibble(x, y))[, 2]) %>%
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
    rename(long = x, lat = y)
  
  saveRDS(d, 'data/ecmwf-era5-2m-temperature/aggregated-data.rds')
} else {
  d <- readRDS('data/ecmwf-era5-2m-temperature/aggregated-data.rds')
}

d

#' some variances are `NA` because there is only 1 value for some (x, y)
mean(is.na(d$monthly_mean))
mean(is.na(d$monthly_var))
quantile(d$n, probs = c(0, 0.003, 0.005, 0.25, 0.5, 1))

d <- filter(d, ! is.na(monthly_var))

# Gamma model for the variance around the mean, given the mean ----
# Assuming `temp_c ~ N(mu, sigma2) implies sigma^2 is independent of mu,
# but extreme temperatures (E(T) < -15, E(T) > 15) are more variable than
# temperatures near 0
# fit in < 90 seconds
m_gamma <-
  bam(monthly_var ~
        s(monthly_mean, bs = 'tp', k = 10) +
        s(month, bs = 'cc', k = 10) +
        s(year, bs = 'tp', k = 8) +
        s(long, lat, bs = 'ds', k = 500) +
        ti(long, lat, month, bs = c('ds', 'cc'), d = c(2, 1),
           k = c(100, 10)),
      family = Gamma(link = 'log'),
      data = d,
      method = 'fREML',
      discrete = TRUE,
      knots = list(month = c(0.5, 12.5)),
      control = gam.control(trace = TRUE))

# good fit, but but residuals have heavy tails with only 2 years
layout(matrix(1:4, ncol = 2))
gam.check(m_gamma, pch = 19, col = '#00000020')
layout(1)

summary(m_gamma) # dev.expl is quite high
plot(m_gamma, pages = 1, scheme = c(1, 1, 1, 2, 2))
plot(m_gamma, select = 1, scheme = 1, xlim = c(0.5, 12.5), n = 300)

# but variance in residuals varies between months
d %>%
  mutate(month = factor(month, levels = 12:1),
         e = resid(m_gamma)) %>%
  ggplot() +
  geom_vline(xintercept = 0, color = 'black') +
  geom_density_ridges(aes(x = e, y = month, fill = month),
                      alpha = 0.75, scale = 5) +
  labs(x = 'Deviance residuals') +
  scale_fill_manual('Month', values = color('BuRd')(6)[c(1:6, 6:1)])

plot(factor(m_gamma$model$month), resid(m_gamma, type = 'deviance'),
     xlab = NULL, ylab = 'Residuals')
abline(h = 0, lwd = 1.5, col = 'red3')

# bias is also non-zero
plot(factor(m_gamma$model$month), abs(resid(m_gamma)),
     xlab = NULL, ylab = 'Absolute residuals')

# checking spatial distributions
d <- mutate(d, e_gamma = resid(m_gamma)) %>%
  group_by(long, lat) %>%
  summarize(mean_e_gamma = mean(e_gamma),
            var_e_gamma = var(e_gamma),
            .groups = 'drop')

# very little spatial bias
ggplot(d, aes(long, lat, fill = mean_e_gamma)) +
  geom_raster() +
  scale_fill_distiller(paste0('E(e) (\U00B0', 'C)'), type = 'div',
                       palette = 5, limits = c(-0.5, 0.5))

# some spatial structure
ggplot(d, aes(long, lat, fill = var_e_gamma)) +
  geom_raster() +
  scale_fill_distiller(paste0('E(e) (\U00B0', 'C)'), type = 'seq',
                       palette = 4, direction = 1, limits = c(0, 1))

# fitting a gammals model to account for seasonal heteroskedasticity ----
m_gammals <-
  gam(list(
    # linear predictor for the mean
    monthly_var ~
      s(month, bs = 'cc', k = 10) +
      s(long, lat, bs = 'ds', k = 500),
    # linear predictor for the scale (variance = mean^2 * scale)
    ~
      s(month, bs = 'cc', k = 10) +
      s(long, lat, bs = 'ds', k = 100)),
    family = gammals(),
    data = d,
    method = 'REML',
    knots = list(month = c(0.5, 12.5)),
    control = gam.control(trace = TRUE))

# check model diagnostics
layout(matrix(1:4, ncol = 2))
gam.check(m_gammals, pch = 19, col = '#00000020')
layout(1)

summary(m_gammals) # 
plot(m_gammals, pages = 1, scheme = c(1, 1, 2, 2))

# but variance in residuals varies between months
plot(factor(m_gammals$model$month), resid(m_gammals))
