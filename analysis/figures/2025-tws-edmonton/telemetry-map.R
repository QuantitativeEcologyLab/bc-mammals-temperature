library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('sf')      # for working with simple features
library('purrr')   # for functional programming
library('ctmm')    # for working with akdes
library('ggplot2') # for fancy plots
library('scales')  #' for `parse_format()`
library('terra')   # for rasters
source('analysis/figures/default-ggplot-theme.R')
source('data/bc-shapefile.R')

tels <- readRDS('models/movement-models-akdes-2024-06-06.rds') %>%
  mutate(dataset_name = case_when(
    dataset_name == 'Canis_lupus_boreal' ~ 'Canis lupus',
    dataset_name == 'Rangifer_tarandus_boreal' ~ 'Rangifer tarandus (boreal)',
    dataset_name == 'Rangifer_tarandus_southern_mountain' ~ 'Rangifer tarandus (southern mountain)',
    dataset_name == 'Puma_concolor_2' ~ 'Puma concolor',
    dataset_name == 'Puma_concolor_4' ~ 'Puma concolor',
    dataset_name == 'Elk in southwestern Alberta' ~ 'Cervus canadensis',
    dataset_name == 'Oreamnos_americanus' ~ 'Oreamnos americanus',
    dataset_name == 'Ursus_arctos_horribilis' ~ 'Ursus arctos horribilis')) %>%
  mutate(dataset_name = gsub(' ', '~', dataset_name) %>%
           gsub('southern~', 'southern ', .) %>%
           gsub('~\\(', '\\)~bold\\("(', .),
         dataset_name = if_else(grepl('bold\\(', dataset_name),
                                paste0('bolditalic(', dataset_name, '")'),
                                paste0('bolditalic(', dataset_name, ')'))) %>%
  transmute(dataset_name,
            tel = map(tel, \(.t) {
              .t %>%
                data.frame() %>%
                select(longitude, latitude) %>%
                st_as_sf(coords = c('longitude', 'latitude')) %>%
                st_set_crs('+proj=longlat') %>%
                st_transform('EPSG:3005') %>%
                as.data.frame()
            })) %>%
  unnest(tel) %>%
  st_as_sf()

# create the figures ----
us <- tigris::states(cb = TRUE, resolution = '5m') %>%
  st_geometry() %>%
  st_transform(crs(bc)) %>%
  st_as_sf()
plot(us)

na <- rbind(prov, us)

plot(na)

na <- rbind(prov, us) %>%
  slice(c(1, 5, 6, 10, 11, 21, 32, 44, 47, 60))

dem <- readRDS('data/presentations/study-area-dem.rds') %>%
  mutate(elev_m = if_else(elev_m < 0, 0, elev_m))

# coords for edmonton conference center
ecc <- tibble(x = -113.4864, y = 53.5417) %>%
  st_as_sf(coords = c('x', 'y')) %>%
  st_set_crs('EPSG:4326') %>%
  st_transform('EPSG:3005')

# create and save the figures
p_dem <-
  ggplot() +
  geom_raster(aes(x, y, fill = elev_m), dem) +
  geom_sf(aes(geometry = geometry), color = 'transparent',
          filter(tels, dataset_name == SPECIES_LABS[5])) +
  geom_sf(data = ecc, color = 'white') +
  coord_sf(xlim = c(95e4, 110e4), clip = 'off',
           ylim = c(33e4, 185e4), crs = 'EPSG:3005') +
  scale_color_manual(name = ' ', values = PAL, labels = parse_format()) +
  scale_fill_distiller(name = '\nElevation (m)', palette = 6) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = 'inside', legend.position.inside = c(-3, 0.05),
        legend.justification = c(1, 0)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1,
                                                  order = 2)),
         fill = guide_colorbar(order = 1))

ggsave('figures/2025-tws-edmonton/map-dem-edmonton.png', plot = p_dem,
       width = 15.11, height = 8.5, units = 'in', dpi = 250, bg = 'white')

p_tels <-
  ggplot() +
  geom_raster(aes(x, y, fill = elev_m), dem) +
  geom_sf(aes(geometry = geometry, color = dataset_name), tels,
          size = 0.1) +
  geom_sf(data = ecc, color = 'white') +
  coord_sf(xlim = c(95e4, 110e4), clip = 'off',
           ylim = c(33e4, 185e4), crs = 'EPSG:3005') +
  scale_color_manual(name = ' ', values = PAL, labels = parse_format()) +
  scale_fill_distiller(name = '\nElevation (m)', palette = 6) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = 'inside', legend.position.inside = c(-3, 0.05),
        legend.justification = c(1, 0)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1,
                                                  order = 2)),
         fill = guide_colorbar(order = 1))

ggsave('figures/2025-tws-edmonton/tels-map-dem-edmonton.png', plot = p_dem,
       width = 15.11, height = 8.5, units = 'in', dpi = 250, bg = 'white')

# boreal caribou only
p_dem_bc <-
  ggplot() +
  geom_raster(aes(x, y, fill = elev_m), dem) +
  geom_sf(aes(geometry = geometry), color = PAL[5],
          filter(tels, dataset_name == SPECIES_LABS[5]), size = 0.1) +
  geom_sf(data = ecc, color = 'white') +
  coord_sf(xlim = c(95e4, 110e4), clip = 'off',
           ylim = c(33e4, 185e4), crs = 'EPSG:3005') +
  scale_color_manual(name = ' ', values = PAL, labels = parse_format()) +
  scale_fill_distiller(name = '\nElevation (m)', palette = 6) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = 'inside', legend.position.inside = c(-3, 0.05),
        legend.justification = c(1, 0)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1,
                                                  order = 2)),
         fill = guide_colorbar(order = 1))

ggsave('figures/2025-tws-edmonton/tels-map-dem-edmonton-boreal-caribou.png',
       plot = p_dem_bc, width = 15.11, height = 8.5, units = 'in',
       dpi = 250, bg = 'white')
