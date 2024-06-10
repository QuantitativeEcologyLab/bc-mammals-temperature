library('sf')      # for working with spatial data
library('terra')  # for working with raster data (tifs)
library('dplyr')   # for data wrangling (e.g., mutate, %>%, ...)
library('ggplot2') # for fancy plots
source('analysis/figures/actws-2024-jasper/actws-2024-jasper-theme.R')
source('data/bc-shapefile.R')

theme_set(theme_void() +
            theme(axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  legend.position = c(0.85, 0.7),
                  legend.text = element_text(face = 'bold'),
                  legend.title = element_text(face = 'bold')))

# create custom color palettes for each predictor
forest_pal <- colorRampPalette(c('#C2B280', 'darkgreen'))(100)
elevation_pal <- colorRampPalette(c('white', 'brown4'))(100)
water_pal <- c('transparent', 'blue')

f <- rast('data/resource-rasters/bc-forest.tif') %>%
  project(crs(bc)) %>%
  as.data.frame(xy = TRUE) %>%
  filter(! is.na(layer))

e <- rast('data/resource-rasters/bc-dem.tif') %>%
  project(crs(bc)) %>%
  as.data.frame(xy = TRUE) %>%
  filter(layer > 0)

w <- rast('data/resource-rasters/bc-distance-from-water.tif') %>%
  project(crs(bc)) %>%
  crop(bc) %>%
  mask(bc) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE) %>%
  rename(layer = consensus_full_class_12)

# check ranges of each variable
range(f$layer)
range(e$layer)
range(w$layer)

if(any(f$layer < 0) | any(f$layer > 100) | any(w$layer < 0)) {
  stop('Check ranges.\n')
} else {
  cat('Ranges are ok.\n')
}

# plot the rasters
p_f <-
  ggplot() +
  geom_raster(aes(x, y, fill = layer), f) +
  geom_sf(data = bc, fill = NA, color = 'black') +
  scale_x_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_fill_gradient('Tree cover (%)', low = 'white', na.value = NA,
                      high = 'darkgreen', limits = c(0, 100),
                      breaks = c(0, 100))
ggsave('figures/actws-2024-jasper/bc-forest.png', plot = p_f,
       width = 6, height = 6, dpi = 600, bg = 'transparent')

p_e <-
  ggplot() +
  geom_raster(aes(x, y, fill = layer), e) +
  geom_sf(data = bc, fill = NA, color = 'black') +
  scale_x_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_fill_distiller('Elevation (m)', palette = 6, direction = 1,
                       breaks = round(range(e$layer, na.rm = TRUE)),
                       labels = round(range(e$layer, na.rm = TRUE), -2),
                       limits = round(range(e$layer, na.rm = TRUE)))
ggsave('figures/actws-2024-jasper/bc-elevation.png', plot = p_e,
       width = 6, height = 6, dpi = 600, bg = 'transparent')

p_w <-
  ggplot() +
  geom_raster(aes(x, y, fill = layer / 1e3), w) +
  geom_sf(data = bc, fill = NA, color = 'black') +
  scale_x_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_fill_distiller(expression(bold(atop(Distance~from,
                                            water~(km)~phantom(om)))),
                       na.value = NA, values = c(0, 0.05, 1),
                       breaks = range(w$layer) / 1e3,
                       labels = round(range(w$layer) / 1e3))
ggsave('figures/actws-2024-jasper/bc-water.png', plot = p_w,
       width = 6, height = 6, dpi = 600, bg = 'transparent')
