library('dplyr') # for data wrangling
library('sf')    # for spatial data

# import shapefile of bc
# download a shapefile of British Columbia
bc_shp <-
  sf::st_as_sf(canadianmaps::PROV) %>% # convert to spatial features (sf) object
  filter(PRENAME == 'British Columbia') %>% # filter to BC only
  sf::st_geometry() %>% # extract boundaries only
  # use projection used by the BC governemnt
  st_transform(crs = '+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs') %>%
  st_as_sf()
