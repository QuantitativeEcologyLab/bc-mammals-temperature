library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('sf')      # for working with simple features
library('purrr')   # for functional programming
library('ctmm')    # for working with akdes
library('ggplot2') # for fancy plots
library('scales')  #' for `parse_format()`
library('terra')   # for rasters
library('cowplot') # for inset
source('analysis/figures/default-ggplot-theme.R')
source('data/bc-shapefile.R')

d <- readRDS('models/movement-models-akdes-2024-06-06.rds') %>%
  mutate(dataset_name = case_when(
    dataset_name == 'Canis_lupus_boreal' ~ 'Wolves',
    dataset_name == 'Rangifer_tarandus_boreal' ~ 'Caribou (boreal)',
    dataset_name == 'Rangifer_tarandus_southern_mountain' ~ 'Caribou (s. mountain)',
    dataset_name == 'Puma_concolor_2' ~ 'Cougars',
    dataset_name == 'Puma_concolor_4' ~ 'Cougars',
    dataset_name == 'Elk in southwestern Alberta' ~ 'Elk',
    dataset_name == 'Oreamnos_americanus' ~ 'Mountain goats',
    dataset_name == 'Ursus_arctos_horribilis' ~ 'Grizzly bears'))

tels <- transmute(d,
                  dataset_name,
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

na <- rbind(prov, us)

dem <- rast('data/resource-rasters/fig-1-dem-z-6.tif') %>%
  project(crs(bc)) %>%
  mask(na) %>%
  as.data.frame(dem, xy = TRUE) %>%
  rename(elev_m = 3)

# some values < 0, but not enough to be useful
mean(dem$elev_m < 0)
mean(dem$elev_m < -10)
mean(dem$elev_m < -100)

dem <- mutate(dem, elev_m = if_else(elev_m < 0, 0, elev_m))

# tels with dem
# Global Change Biology discourages the use of boundaries in figures:
# https://onlinelibrary.wiley.com/page/journal/13652486/homepage/forauthors.html#4
p_dem <-
  ggplot() +
  geom_raster(aes(x, y, fill = elev_m), dem) +
  # geom_sf(data = na, color = 'white', fill = 'transparent', lwd = 0.1) +
  geom_sf(aes(geometry = geometry, color = dataset_name), tels,
          size = 0.1) +
  coord_sf(xlim = c(100e4, 110e4), clip = 'off',
           ylim = c(33e4, 185e4), crs = 'EPSG:3005') +
  scale_color_manual(name = ' ', values = PAL) +
  scale_fill_distiller(name = '\nElevation (m)', palette = 6) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = 'inside', legend.position.inside = c(-4.5, 0.2),
        legend.justification = c(1, 0), text = element_text(face = 'bold')) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2,
                                                  order = 2)),
         fill = guide_colorbar(order = 1))

p_inset <-
  spData::world %>%
  filter(continent == 'North America') %>%
  st_transform('EPSG:3005') %>%
  ggplot() +
  geom_sf(fill = 'black', color = 'black') +
  geom_sf(aes(geometry = geometry, color = dataset_name), tels,
          color = 'darkorange', size = 0.1) +
  theme_void() +
  theme(panel.border = element_rect(colour = 'black', fill = 'transparent'))

p_full <- ggdraw() +
  draw_plot(p_dem) +
  draw_plot(p_inset, height = 0.2, width = 0.2, x = 0.01, y = 0)

ggsave('figures/tels-map-dem.png', plot = p_full,
       width = 8, height = 6.8, units = 'in', dpi = 600, bg = 'white')
