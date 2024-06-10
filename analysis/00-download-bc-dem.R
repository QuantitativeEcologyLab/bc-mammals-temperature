library('sf')           # for spatial data
library('dplyr')        # for data wrangling
library('canadianmaps') # for shapefile of BC
library('elevatr')      # for downloading DEM (using version 0.99.0)
library('terra')        # for working with rasters
library('purrr')        # for functional programming
library('ctmm')         # for movement modeling
source('data/bc-shapefile.R')

ZOOM <- 3 #' resolution level for the DEM (see `?get_elev_raster`)

bc_dem <-
  get_elev_raster(locations = bc_unproj, z = ZOOM, clip = 'locations') %>%
  rast() #' `elevatr` v 0.99.0 still returs a `RasterLayer` object

plot(bc_dem)
plot(bc_unproj, add = TRUE)

if(FALSE) {
  # highest point in BC is mount fairweather at ~ 4e3 m
  #' `ZOOM > 6` results in excessively high elevation values
  max(values(bc_dem), na.rm = TRUE)
  get_elev_raster(locations = bc_unproj, z = 7, clip = 'locations') %>%
    rast() %>% #' `elevatr` v 0.99.0 still returs a `RasterLayer` object
    values() %>%
    max(na.rm = TRUE)
  
  # res is close enough
  res(bc_dem)
  res(rast('data/resource-rasters/bc-forest.tif'))
  res(rast('data/resource-rasters/bc-distance-from-water.tif'))
}

terra::writeRaster(bc_dem, paste0('data/bc-dem-z', ZOOM, '.tif'))

# download dem of buffered BC
bc_unproj_buff <- bc_unproj %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_as_sf() %>%
  st_buffer(300e3) # contains 99.9995% UDs

# check ranges
plot(bc_unproj_buff)
plot(bc_unproj, add = TRUE)
if(file.exists('models/movement-models-akdes-2024-06-06.rds')) {
  akdes <-
    readRDS('models/movement-models-akdes-2024-06-06.rds') %>%
    # remove animals that are not range resident (invalid weights)
    filter(! animal %in% c('SCEK014', 'SCEK014b', 'BW028')) %>%
    pull(akde) %>%
    map_dfr(\(.a) {
    SpatialPolygonsDataFrame.UD(.a, level.UD = 0.9995, level = 0) %>%
      st_as_sf() %>%
      st_geometry() %>%
      st_union() %>% # unite estimate with lower and upper 95% CIs
      st_set_crs(.a@info$projection) %>%
      st_transform('EPSG:4326') %>%
      st_as_sf()
  }) %>%
    st_make_valid() %>%
    st_union() %>%
    st_as_sf()
  
  ggplot() +
    geom_sf(data = bc_unproj_buff) +
    geom_sf(data = bc_unproj, fill = '#654321') +
    geom_sf(data = akdes, fill = 'forestgreen') +
    theme_map()
}

get_elev_raster(locations = bc_unproj_buff, z = ZOOM, clip = 'bbox') %>%
  rast() %>% #' `elevatr` v 0.99.0 still returs a `RasterLayer` object
  terra::writeRaster(paste0('data/resource-rasters/bc-buffered-dem-z',
                            ZOOM, '.tif'))

terra::plot(rast('data/resource-rasters/bc-buffered-dem-z3.tif'))
if(exists('akdes')) plot(akdes, add = TRUE, col = '#00000080')
