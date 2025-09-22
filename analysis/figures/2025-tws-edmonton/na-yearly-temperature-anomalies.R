library('dplyr')     # for data wrangling
library('ggplot2')   # for fancy plots
library('lubridate') # for working with dates
source('analysis/figures/2024-actws-jasper/2024-actws-jasper-theme.R')

# study area temperature anomalies relative to average from 
d <- read.csv('data/presentations/climate-anomalies/data-2025.csv',
              skip = 4) %>%
  mutate(Date = as.Date(paste0(Date, '15'), format = '%Y%m%d'))

ggplot(d, aes(Date, Anomaly, color = Anomaly > 0)) +
  geom_bar(stat = 'identity', fill = 'transparent') +
  geom_hline(yintercept = 0) +
  labs(x = NULL, y = paste0('Global temperature anomaly (\U00B0', 'C)')) +
  scale_x_date(expand = c(0, 365 * 4),
               breaks = c(1850.5, 1900.5, 1950.5, 2000.5, 2025.5) %>%
                 date_decimal() %>% as.Date(), labels = year) +
  scale_color_manual(values = c('#0571B0', '#CA0020'))

ggsave('figures/2025-tws-edmonton/temperature-anomalies.png',
       scale = 1, width = 10, height = 6, dpi = 600)
