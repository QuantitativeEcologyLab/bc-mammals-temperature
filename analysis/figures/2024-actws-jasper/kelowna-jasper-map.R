library('sf')      # for spatial objects
library('ggplot2') # for fancy plots
library('dplyr')   # for data wrangling

provs <- canadianmaps::PROV %>%
  st_geometry() %>%
  st_transform('EPSG:3005')

plot(provs)

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

plot(locs, add = TRUE)

ggplot() +
  geom_sf(data = provs) +
  geom_sf(data = locs, size = 2) +
  geom_label(aes(x = alb_long, y = alb_lat, label = city), locs,
             nudge_x = 5e4, nudge_y = 5e4, label.size = 0, size = 6,
             bg = 'transparent') +
  coord_sf(xlim = c(25e4, 26e5), ylim = c(37e4, 175e4)) +
  theme(legend.position = 'none')

ggsave('figures/actws-2024-jasper/kelowna-jasper-map.png',
       width = 33.87, height = 19.05, units = 'cm', dpi = 600,
       bg = 'transparent')
