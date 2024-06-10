library('dplyr')        # for data wrangling
library('tidyr')        # for data wrangling
library('ggplot2')      # for fancy plots
library('khroma')       # for color blind-friendly palettes
source('data/bc-shapefile.R')
source('analysis/figures/actws-2024-jasper/actws-2024-jasper-theme.R')

locs <- tibble(city = c('Kelowna', 'Jasper'),
               long = c(-119.4966, -118.0814),
               lat = c(49.8863, 52.8737)) %>%
  st_as_sf(coords = c('long', 'lat')) %>%
  st_set_crs('+proj=longlat') %>%
  st_transform('EPSG:3005') %>%
  bind_cols(.,
            st_coordinates(.) %>%
              as.data.frame() %>%
              rename(alb_long = X, alb_lat = Y))

# boreal datasets
d <- readRDS('data/tracking-data/all-tracking-data-cleaned-2023-11-17-16-04.rds') %>%
  filter(dataset_name == 'Rangifer_tarandus_boreal') %>%
  unnest(tel) %>%
  filter(! outlier) %>%
  select(location.long, location.lat, dataset_name, species) %>%
  st_as_sf(coords = c('location.long', 'location.lat')) %>%
  st_set_crs('+proj=longlat') %>%
  st_transform('EPSG:3005')

# phylopics added manually using hex codes below
pal <- c(color('bright')(7), brown = '#654321')

# caribou data only
p_caribou <-
  ggplot() +
  theme_void() +
  theme(legend.position = 'none') +
  geom_sf(data = prov) +
  geom_sf(data = d, color = '#EE6677') +
  geom_sf(data = locs, size = 2) +
  geom_label(aes(x = alb_long, y = alb_lat, label = city), locs,
             nudge_x = 5e4, nudge_y = 5e4, label.size = 0, size = 6,
             bg = 'transparent') +
  coord_sf(xlim = c(2.5e5, 26e5), ylim = c(4.2e5, 18e5)) +
  scale_color_bright(name = 'Species') +
  labs(x = NULL, y = NULL)

ggsave('figures/actws-2024-jasper/caribou-map.png', p_caribou,
       width = 33.87, height = 19.05, units = 'cm', dpi = 600,
       bg = 'transparent')

# all species
d <- readRDS('data/tracking-data/all-tracking-data-cleaned-2023-11-17-16-04.rds') %>%
  unnest(tel) %>%
  filter(! outlier) %>%
  select(location.long, location.lat, dataset_name, species) %>%
  st_as_sf(coords = c('location.long', 'location.lat')) %>%
  st_set_crs('+proj=longlat') %>%
  st_transform('EPSG:3005') %>%
  mutate(dataset_name = if_else(grepl('Puma_concolor', dataset_name),
                                'Puma_concolor', dataset_name),
         dataset_name =
           factor(dataset_name,
                  levels = c('Canis_lupus_boreal',
                             'Rangifer_tarandus_boreal',
                             'Rangifer_tarandus_southern_mountain',
                             'Puma_concolor',
                             'Elk in southwestern Alberta',
                             'Oreamnos_americanus',
                             'Ursus_arctos_horribilis')))

# plot all datasets
p_all <-
  ggplot() +
  theme_void() +
  theme(legend.position = 'none') +
  geom_sf(data = prov) +
  geom_sf(aes(color = dataset_name), d) +
  geom_sf(data = locs, size = 2) +
  geom_label(aes(x = alb_long, y = alb_lat, label = city), locs,
             nudge_x = 5e4, nudge_y = 5e4, label.size = 0, size = 6,
             bg = 'transparent') +
  coord_sf(xlim = c(2.5e5, 26e5), ylim = c(3.8e5, 18e5)) +
  scale_color_manual('Species', values = pal) +
  labs(x = NULL, y = NULL)

ggsave('figures/actws-2024-jasper/telemetries-all-map.png', plot = p_all,
       width = 33.87, height = 19.05, units = 'cm', dpi = 600,
       bg = 'transparent')
