library('dplyr')   # for data wrangling
library('ggplot2') # for fancy plots
source('analysis/figures/2024-actws-jasper/2024-actws-jasper-theme.R')

# north america temperature anomalies relative to average from 1910-2000
d <- read.csv('data/presentations/climate-anomalies/data-2024.csv',
              skip = 4) %>%
  mutate(Date = as.Date(paste0(Date, '15'), format = '%Y%m%d'))

ggplot(d, aes(Date, Anomaly, color = Anomaly > 0)) +
  geom_bar(stat = 'identity', fill = 'transparent') +
  geom_hline(yintercept = 0) +
  labs(x = NULL, y = paste0('Temperature Anomaly (\U00B0', 'C)')) +
  ylim(c(-2.25, 2.25)) +
  scale_x_date(expand = c(0, 365 * 2)) +
  scale_color_manual(values = c('#0571B0', '#CA0020'))

ggsave('figures/2024-actws-jasper/temperature-anomalies.png',
       scale = 1, width = 10, height = 6, dpi = 600)
