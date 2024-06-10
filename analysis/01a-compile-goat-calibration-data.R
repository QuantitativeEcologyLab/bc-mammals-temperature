library('dplyr')     # for data wrangling (mutate(), %>%, etc.)
library('tidyr')     # for data wrangling (nest(), unnest(), pivot_*, etc.)
library('purrr')     # for functional programming (map_***(), etc.)
library('ctmm')      # for movement models
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
theme_set(theme_bw())

# source custom functions
source('functions/outlier_plots.R') # to plot outlier diagnostic plots
source('functions/check_animal.R') # to run diagnostic plots
source('functions/plot_adj.R') # to plot 20 adjacent locations
source('functions/flag_outlier.R') # to mark outliers

# import goat data
goats <- readRDS('data/tracking-data/all-tracking-data-not-cleaned.rds') %>%
  filter(species == 'Oreamnos americanus') %>%
  select(! c(hdop, individual.local.identifier, species, dataset_name,
             gps.fix.type.raw:gps.satellite.count))

# only keep calibration data
d <- goats %>%
  filter(location.long > -119.7) %>%
  mutate(animal = paste(year(timestamp), month(timestamp), animal,
                        sep = '_')) %>%
  group_by(animal) %>%
  mutate(individual.local.identifier = animal) %>% # for nesting ID 
  nest(tel = ! animal) %>%
  arrange(animal)

# look at number of rows for each tag
hist(map_int(d$tel, nrow), breaks = 20)

# look at variograms for each set
variograms <- function() {
  N <- ceiling(sqrt(nrow(d)))^2
  layout(matrix(1:N, ncol = sqrt(N), byrow = TRUE))
  for(i in 1:nrow(d)) {
    ctmm.guess(as.telemetry(d$tel[[i]]), interactive = FALSE)
    title(d$animal[i])
  }
  layout(1)
}
variograms()

# collars turned on for use in 2019-06; dropping all >= 2019 to be safe
d <- filter(d, map_int(tel, ~ unique(year(.x$timestamp))) < 2019)
variograms()

d <- d %>%
  mutate(tel = map(tel,
                   ~ suppressMessages(as.telemetry(.x) %>%
                                        data.frame()))) %>%
  unnest(tel) %>%
  group_by(animal) %>%
  mutate(hours = 'hours' %#% (t - first(t)))

# plot projected telemetries
ggplot(d, aes(x, y, color = hours)) +
  facet_wrap(~ animal, scales = 'free') +
  geom_path(linewidth = 1) +
  geom_point(aes(size = HDOP), alpha = 0.3) +
  scale_size(limits = range(d$HDOP)) +
  scale_color_gradientn(colors = khroma::color('bright')(6))

#' drop the last 3 points of `2018_8_30561`
nrow(d)
tail(filter(d, animal == '2018_8_30561')$timestamp)
d <- filter(d, ! (animal == '2018_8_30561' &
                    timestamp > as.POSIXct('2018-08-07 02:00:20 PDT'))) %>%
  mutate(distance = sqrt(x^2 + y^2))
nrow(d)

ggplot(d, aes(x, y, color = hours)) +
  facet_wrap(~ animal, scales = 'free') +
  geom_path(size = 1) +
  geom_point(aes(size = HDOP), alpha = 0.3) +
  scale_size(limits = range(d$HDOP)) +
  scale_color_gradientn(colors = khroma::color('bright')(6))

ggplot(d) +
  facet_wrap(~ animal, scales = 'free_y') +
  geom_density(aes(distance), fill = 'grey', bw = 10)

# points look independent
d %>%
  mutate(distance = sqrt(x^2 + y^2)) %>%
  filter(distance < 50) %>% # drop larger errors to focus on small patterns
  ggplot(aes(x, y, color = hours)) +
  facet_wrap(~ animal, scales = 'free') +
  geom_path(size = 1) +
  geom_point(aes(size = HDOP), alpha = 0.3) +
  scale_size(limits = range(d$HDOP)) +
  scale_color_gradientn(colors = khroma::color('bright')(6))

d <- select(d, ! c(x, y, hours, distance))

# all locations are 3D locations
table(d$class)

# test the UERE object
goat_uere <- uere.fit(as.telemetry(d))

# adding the UERE does not affect the plot
goats %>%
  mutate(outlier = location.long > 0) %>%
  as.telemetry(mark.rm = TRUE) %>%
  plot(ylim = c(-30e3, 20e3))
goats %>%
  mutate(outlier = location.long > 0) %>%
  as.telemetry(mark.rm = TRUE) %>%
  `uere<-`(goat_uere) %>%
  plot(ylim = c(-30e3, 20e3), add = TRUE, col = 'black')

write.csv(x = d, 'data/tracking-data/goat-calibration-data.csv',
          row.names = FALSE)
