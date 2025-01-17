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


# import digital elevation raster (zoom > 6 gives excessively high peaks)
#' see `analysis/00-download-bc-dem.R` for downloading
dem <- rast('data/resource-rasters/bc-buffered-dem-z6.tif')
plot(dem)
plot(locs_ext, add = TRUE, lwd = 3)

# forest raster (original rasters from http://www.earthenv.org/landcover)
f1 <- rast('data/resource-rasters/uncropped-rasters/consensus_full_class_1.tif') %>%
  crop(dem)
f2 <- rast('~/../Downloads/consensus_full_class_2.tif') %>%
  crop(dem)
f3 <- rast('~/../Downloads/consensus_full_class_3.tif') %>%
  crop(dem)
f4 <- rast('~/../Downloads/consensus_full_class_4.tif') %>%
  crop(dem)

layout(matrix(1:4, ncol = 2))
plot(f1, main = 'layer 1: Evergreen/Deciduous Needleleaf Trees')
plot(f2, main = 'layer 2: Evergreen Broadleaf Trees')
plot(f3, main = 'layer 3: Deciduous Broadleaf Trees')
plot(f4, main = 'layer 4: Mixed/Other Trees')
layout(1)

f <- f1 + f2 + f3 + f4
layout(1:2)
plot(f1, main = 'Evergreen/Deciduous Needleleaf Trees')
plot(f, main = 'all trees')
layout(1)

max(values(f))
values(f) <- if_else(values(f) > 100, 100, values(f))
max(values(f))

writeRaster(f, 'data/resource-rasters/forest.tif')
plot(rast('data/resource-rasters/forest.tif'))
plot(locs_ext, add = TRUE, lwd = 3)

# raster of open water (original rasters from http://www.earthenv.org/landcover)
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
f %>%
  crop(bc_unproj) %>%
  mask(bc_unproj) %>%
  writeRaster('data/resource-rasters/bc-forest.tif')
plot(rast('data/resource-rasters/bc-forest.tif'))

# raster of open water
dist %>%
  crop(bc_unproj) %>%
  mask(bc_unproj) %>%
  writeRaster('data/resource-rasters/bc-distance-from-water.tif')
plot(rast('data/resource-rasters/bc-distance-from-water.tif'))
