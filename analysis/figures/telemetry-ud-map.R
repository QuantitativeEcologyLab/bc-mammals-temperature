library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('sf')      # for working with simple features
library('purrr')   # for functional programming
library('ctmm')    # for working with akdes
library('ggplot2') # for fancy plots
library('scales')  #' for `parse_format()`
source('analysis/figures/default-ggplot-theme.R')
source('data/bc-shapefile.R')

d <- readRDS('models/movement-models-akdes-2024-06-06.rds') %>%
  mutate(dataset_name = case_when(
    dataset_name == 'Canis_lupus_boreal' ~ 'Canis lupus',
    dataset_name == 'Rangifer_tarandus_boreal' ~ 'Rangifer tarandus (boreal)',
    dataset_name == 'Rangifer_tarandus_southern_mountain' ~ 'Rangifer tarandus (southern mountain)',
    dataset_name == 'Puma_concolor_2' ~ 'Puma concolor',
    dataset_name == 'Puma_concolor_4' ~ 'Puma concolor',
    dataset_name == 'Elk in southwestern Alberta' ~ 'Cervus canadensis',
    dataset_name == 'Oreamnos_americanus' ~ 'Oreamnos americanus',
    dataset_name == 'Ursus_arctos_horribilis' ~ 'Ursus arctos horribilis')) %>%
  mutate(dataset_name = gsub(' ', '~', dataset_name),
         dataset_name = gsub('~\\(', '\\)~bold\\((', dataset_name),
         dataset_name = paste0('bolditalic(', dataset_name, ')'))

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
  unnest(tel)

uds <- d %>%
  transmute(
    dataset_name,
    akde =  map(akde, \(.a) {
      SpatialPolygonsDataFrame.UD(.a, level = 0, level.UD = 0.95) %>%
        st_as_sf() %>%
        st_set_crs(.a@info$projection) %>%
        filter(grepl(pattern = 'est', x = name)) %>%
        st_transform('EPSG:3005')
    }) %>%
      bind_rows() %>%
      pull(geometry))

# create the figures ----
# tels only
p_tels <-
  ggplot() +
  geom_sf(data = bc) +
  geom_sf(aes(geometry = geometry, color = dataset_name), tels,
          size = 0.1) +
  scale_fill_manual(name = NULL, values = PAL, labels = parse_format(),
                    aesthetics = c('color', 'fill')) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = 'inside', legend.position.inside = c(-0.025, 0.5),
        legend.justification = c(0, 0)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1)))

ggsave('figures/tels-map.png', plot = p_tels,
       width = 10, height = 8.5, units = 'in', dpi = 600, bg = 'white')

# UDs only
p_uds <-
  ggplot() +
  geom_sf(data = bc) +
  geom_sf(aes(geometry = akde, fill = dataset_name, color = dataset_name),
          uds, alpha = 0.3) +
  scale_fill_manual(name = NULL, values = PAL, labels = parse_format(),
                    aesthetics = c('color', 'fill')) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = 'inside', legend.position.inside = c(0, 0.45),
        legend.justification = c(0, 0))

ggsave('figures/uds-map.png', plot = p_uds,
       width = 10, height = 10.5, units = 'in', dpi = 600, bg = 'white')

# telemetries and UDs
p_tels <- p_uds +
  geom_sf(aes(geometry = geometry, color = dataset_name), tels, pch = '.')

ggsave('figures/tels-and-uds-map.png', plot = p_tels,
       width = 10, height = 10.5, units = 'in', dpi = 600, bg = 'white')
