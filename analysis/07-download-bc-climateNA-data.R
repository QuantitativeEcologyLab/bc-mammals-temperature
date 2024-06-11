library('climatenaR') # for downloading climate data and projections
library('purrr')      # for functional programming

#' if necessary, install the `climatenaR` package with
#' `remotes::install_github('burnett-m/climatenaR', build_vignettes = TRUE)`
#' *NOTE:* the `climatenaR` package requires the ClimateNA or climateBC
#'         software. See https://register.climatena.ca/ to download it.

#' metadata for `climateNA` data at:
#' https://s3-us-west-2.amazonaws.com/www.cacpd.org/documents/ClimateNAv7_manual.pdf

#' change the working directory as required by `climatenaR`
setwd('ClimateNA_v742')

if(! file.exists('bc-buffered-dem-z3.csv')) {
  #' convert the bc DEM to a csv as required by `climatenaR`
  demToCSV(file = '../data/resource-rasters/bc-buffered-dem-z3.tif',
           outdir = '.', # save in climateBC folder
           srs = NULL) # keep NULL if in lat/long
  
  # check the csv
  read.csv('bc-buffered-dem-z3.csv', nrows = 5)
}

if(! dir.exists('bc-buffered-dem-z3')) dir.create('bc-buffered-dem-z3')

# download climate data projections ----
if(! dir.exists('bc-buffered-dem-z3/projection-data')) {
  dir.create('bc-buffered-dem-z3/projection-data')
}

map(2011:2100,
    \(y) {
      cat(paste0('Downloading climate projections for ', y, '...\n'))
      projClimateNA(file = 'bc-buffered-dem-z3.csv',
                    tFrame = 'M', # monthly averages
                    exe = 'ClimateNA_v7.42.exe', # must be in wd
                    scen = '8GCM', # 8GCMs_ensemble General Circulation Model
                    ssp = c('S1', 'S2', 'S3', 'S5'), # SSP 4 not available
                    years = as.character(y)) # can only do one year at a time
    })

# check extent
if(FALSE) {
  library('ggplot2')
  library('dplyr')
  library('sf')
  theme_set(cowplot::theme_map())
  
  read.csv('bc-buffered-dem-z3/8GCMs_ensemble_ssp126@2011.csv') %>%
    ggplot() +
    geom_raster(aes(Longitude, Latitude, fill = Tave01)) +
    geom_sf(data = readRDS('../data/movement-models-speed-weights-temperature-2024-06-10.rds') %>%
              select(longitude, latitude) %>%
              st_as_sf(coords = c('longitude', 'latitude')) %>%
              st_bbox() %>%
              st_as_sfc() %>%
              st_as_sf(), fill = '#BBBBBB60') +
    scale_fill_distiller(type = 'div', palette = 5, limits = c(-40, 40))
}

# download historical climate data ----
if(! dir.exists('bc-buffered-dem-z3/historical-data')) {
  dir.create('bc-buffered-dem-z3/historical-data')
}

# tracking data ranges from 1998 to 2021
map(1998:2023,
    \(y) {
      cat(paste0('Downloading estimated historical data for ', y, '...\n'))
      histClimateNA(file = 'bc-buffered-dem-z3.csv',
                    dateR = as.character(y), # year
                    tFrame = 'M', # monthly averages
                    exe = 'ClimateNA_v7.42.exe', # must be in wd
                    outdir = 'bc-buffered-dem-z3/historical-data')
    })
