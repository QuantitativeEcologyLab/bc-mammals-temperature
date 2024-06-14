library('sf')      # for spatial objects
library('dplyr')   # for data wrangling

prov_unproj <- canadianmaps::PROV %>%
  st_geometry() %>%
  st_as_sf()

prov <- st_transform(prov_unproj, 'EPSG:3005')

bc_unproj <- filter(canadianmaps::PROV, PRENAME == 'British Columbia') %>%
  st_geometry() %>%
  st_as_sf()

bc <- st_transform(bc_unproj, 'EPSG:3005')
