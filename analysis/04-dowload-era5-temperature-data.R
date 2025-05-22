library('dplyr')     # for data wrangling
library('ncdf4')     # for working with nc rasters
library('terra')     # for working with rasters
library('sf')        # for working with simple features
library('lubridate') # for working with dates
library('purrr')     # for functional programming
library('ctmm')      # for working with akdes
source('data/bc-shapefile.R')

#' using `{furrr}` causes: `Error: external pointer is not valid`

# find extent of the data for the download requests ----
# downloaded 1998-2023 data from https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels
# url above now redirects to the one below:
# downloaded 2024-2025 data from https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels
d <- readRDS('models/movement-models-akdes-2024-06-06.rds')

# find bounding box for the UDs for the download requests
shp <- d %>%
  pull(akde) %>%
  map(\(.a) {
    # extract a very wide UD with no CIs
    SpatialPolygonsDataFrame.UD(.a, level = 0, level.UD = 0.9999) %>%
      st_as_sf() %>%
      st_set_crs(.a@info$projection) %>% # using tpeqd projection
      st_transform('+proj=longlat') %>% # unproject
      bind_rows() %>%
      # take estimate only
      filter(grepl(pattern = 'est', x = name)) %>%
      # convert to a shapefile of the bounding box
      st_bbox() %>%
      st_as_sfc() %>%
      st_as_sf()
  }) %>%
  # join all bboxes together
  bind_rows() %>%
  st_combine() %>%
  # take a bbox from all of them
  st_bbox() %>%
  st_as_sfc() %>%
  st_as_sf() %>%
  # buffer by 10 kilometers
  st_buffer(1e4)

st_bbox(shp)

# rounded to 2 decimals:
#    xmin    ymin    xmax    ymax
# -141.61   46.49 -103.88   70.14

# check range of years
d %>%
  pull(tel) %>%
  map_dfr(data.frame) %>% #' `as.data.frame()` fails
  pull(timestamp) %>%
  year() %>%
  range()

# raster requests took <1 minute to fill and 5-15 minutes to download ----

# copy the downloads to the repo ----
downloaded_nc_files <- list.files(path = 'C:/Users/mezzinis/Downloads/',
                                  pattern = '\\.nc', full.names = TRUE)

# check which files are missing (old downloads included a time attribute)
c(1998:2023)[! 1998:2023 %in% map_int(downloaded_nc_files, \(.fn) {
  year(as.POSIXct(rast(.fn)[[1]]@pntr$time, tz = 'UTC'))
})]

# check which files are missing (2025 downloads don't include a time attr)
c(2024, 2025)[! c(2024, 2025) %in% map_int(downloaded_nc_files, \(.fn) {
  year(as.POSIXct(as.numeric(gsub('t2m_valid_time=', '', names(rast(.fn)[[1]]))), tz = 'UTC'))
})]

# add times to the 2024 and 2025 rasters
map(downloaded_nc_files, \(.fn) {
  r <- rast(.fn)
  
  y <- gsub('t2m_valid_time=', '', names(r)[1]) %>%
    as.numeric() %>%
    as.POSIXct(tz = 'UTC') %>%
    year()
  
  if(y %in% 2024:2025) {
    r <- rast(.fn)
    terra::time(r) <-
      gsub('t2m_valid_time=','',names(r)) %>%
      as.numeric() %>%
      as.POSIXct(tz = 'UTC')
    writeRaster(r, filename = gsub('.nc', '-w-time.nc', .fn))
  } # close if
  
  return(paste('Saved', .fn))
})

# copy the downloads to the repo
dir.create('data/ecmwf-era5-2m-temperature/original-downloads',
           recursive = TRUE)
file.copy(
  from = list.files(path = 'C:/Users/mezzinis/Downloads',
                    pattern = 'w-time\\.nc', full.names = TRUE),
  to = paste0('data/ecmwf-era5-2m-temperature/original-downloads/',
              list.files(path = 'C:/Users/mezzinis/Downloads',
                         pattern = 'w-time\\.nc', full.names = FALSE)))

# rename based on year and move global files to the "no-crs" folder ----
nc_files <- list.files(
  path = 'data/ecmwf-era5-2m-temperature/original-downloads/',
  pattern = '\\.nc$', full.names = TRUE)

dir.create('data/ecmwf-era5-2m-temperature/no-crs', recursive = TRUE)

#' **HERE**
map(nc_files, function(.file) {
  y <- year(as.POSIXct(rast(.file)[[1]]@pnt$time, tz = 'UTC'))
  
  file.copy(
    from = .file,
    to = paste0('data/ecmwf-era5-2m-temperature/no-crs/ecmwf-era5-2m-temp-',
                y, '.nc'))
})

# set the CRS for each of the moved rasters ----
# get the names of the renamed files
nc_files <- list.files(path = 'data/ecmwf-era5-2m-temperature/no-crs/',
                       pattern = '\\.nc', full.names = FALSE)

dir.create('data/ecmwf-era5-2m-temperature/epsg-4326', recursive = TRUE)

map(nc_files, function(.file) {
  # import the raster
  r <- rast(paste0('data/ecmwf-era5-2m-temperature/no-crs/', .file))
  
  # add the CRS (unprojected lat/long)
  set.crs(r, 'EPSG:4326')
  
  # create new filename
  new_filename <- paste0('data/ecmwf-era5-2m-temperature/epsg-4326/',
                         substr(.file, start = 1, stop = nchar(.file) - 3),
                         '-epsg-4326.nc')
  
  # save the raster (max compression does not reduce below 100 MB)
  writeCDF(x = r, filename = new_filename)
}, .progress = TRUE)

# check that all rasters have a layer for each hour
nc_files <- list.files(path = 'data/ecmwf-era5-2m-temperature/epsg-4326/',
                       pattern = '\\.nc', full.names = TRUE)

unique(map_dbl(nc_files, \(.fn) nlyr(rast(.fn)))) / 24
