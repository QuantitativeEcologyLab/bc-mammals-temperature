library('dplyr')   # for data wrangling
library('ggplot2') # for fancy plots
source('analysis/figures/actws-2024-jasper/actws-2024-jasper-theme.R')

# temperature anomalies relative to average from 
d <- read.csv('data/climate-anomalies/data.csv', skip = 4) %>%
  mutate(Date = as.Date(paste0(Date, '15'), format = '%Y%m%d'))

ggplot(d, aes(Date, Anomaly, color = Anomaly > 0)) +
  geom_bar(stat = 'identity', fill = 'transparent') +
  geom_hline(yintercept = 0) +
  labs(x = NULL, y = paste0('Temperature Anomaly (\U00B0', 'C)')) +
  ylim(c(-2.25, 2.25)) +
  scale_x_date(expand = c(0, 365 * 2)) +
  scale_color_manual(values = c('#0571B0', '#CA0020'))

ggsave('figures/actws-2024-jasper/temperature-anomalies.png',
       scale = 1, width = 10, height = 6, dpi = 600)
