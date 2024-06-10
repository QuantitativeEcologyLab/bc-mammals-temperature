library('dplyr') # for data wrangling
library('tidyr') # for data wrangling
library('purrr') # for functional programming
library('sf')    # for working with spatial data
library('sp')    # for working with spatial data
library('terra') # for raster data (faster than raster)
library('canadianmaps') # to download a shapefile of BC

# import a shapefile of the locations
locs <-
  readRDS('data/tracking-data/all-tracking-data-cleaned-2024-02-22-13-49.rds') %>%
  unnest(tel) %>%
  filter(! outlier) %>%
  select(location.long, location.lat, animal) %>%
  st_as_sf(coords = c('location.long', 'location.lat'))

locs_ext <- locs %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_as_sf()

# original rasters from http://www.earthenv.org/landcover
dem <- rast('data/resource-rasters/bc-buffered-dem-z6.tif')
plot(dem)
plot(locs_ext, add = TRUE, lwd = 3)

# forest raster
rast('data/resource-rasters/uncropped-rasters/consensus_full_class_1.tif') %>%
  crop(dem) %>%
  writeRaster('data/resource-rasters/forest.tif')
plot(rast('data/resource-rasters/forest.tif'))
plot(locs_ext, add = TRUE, lwd = 3)

# raster of open water
layout(matrix(1:8, ncol = 2, byrow = TRUE))
w <- rast('data/resource-rasters/uncropped-rasters/consensus_full_class_12.tif') %>%
  crop(dem)
plot(w, main = 'Percent water')
plot(locs_ext, add = TRUE, lwd = 3)
hist(w, breaks = 100, main = 'Percent water')

w01 <- ceiling(w / 100) # convert anything above 0 to a 1
plot(w01, main = 'Yes/No water')
hist(w01, main = 'Yes/No water')

w1 <- classify(w01, cbind(0, NA))
plot(w1, main = 'Water only', col = 'blue')
unique(w1, main = 'Water only')

dist <- distance(w1)
plot(dist, main = 'Distance from nearest water')
hist(dist, main = 'Distance from nearest water')
writeRaster(dist, 'data/resource-rasters/distance-from-water.tif')

# create a raster of distance to water for all of BC ----
source('data/bc-shapefile.R')

# raster of tree cover
raster('data/resource-rasters/uncropped-rasters/consensus_full_class_1.tif') %>%
  crop(bc_unproj) %>%
  mask(bc_unproj) %>%
  writeRaster('data/resource-rasters/bc-forest.tif')
plot(raster('data/resource-rasters/bc-forest.tif'))

# raster of open water
layout(matrix(1:8, ncol = 2, byrow = TRUE))
w <- rast('data/resource-rasters/uncropped-rasters/consensus_full_class_12.tif') %>%
  crop(bc_unproj) %>%
  aggregate(8)
plot(w, main = 'Percent water')
hist(w, breaks = 100, main = 'Percent water')

w01 <- ceiling(w / 100) %>% # convert anything above 0 to a 1
  mask(st_as_sf(bc_unproj)) # remove sea
plot(w01, main = 'Yes/No water')
hist(w01, main = 'Yes/No water')

w1 <- classify(w01, cbind(0, NA))
plot(w1, main = 'Water only', col = 'blue')
hist(w1, main = 'Water only')

dist <- distance(w1)
dist <- mask(dist, st_as_sf(bc_unproj))
plot(dist, main = 'Distance from nearest water')
hist(dist, main = 'Distance from nearest water')
writeRaster(dist,
            'data/resource-rasters/bc-distance-from-water-coarse.tif')
