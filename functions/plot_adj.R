# function to plot telemetry points near problematic points
plot_adj <- function(id, n_adj = 10,
                     min_speed = NA, max_speed = NA,
                     min_distance = NA, max_distance = NA,
                     min_angle = NA, max_angle = NA,
                     min_dt = NA, max_dt = NA,
                     many = FALSE, reset_layout = TRUE, map = FALSE,
                     map_type = c('Esri.WorldImagery',
                                  'OpenStreetMap',
                                  'Esri.WorldStreetMap',
                                  'CartoDB.Positron',
                                  'CartoDB.DarkMatter')) {
  if(is.numeric(map_type)) {
    map_type <- c('Esri.WorldImagery', 'OpenTopoMap',
                  'OpenStreetMap.Mapnik')[map_type]
  }
  tel <- filter(d, animal == id)$tel[[1]] # find the telemetry
  tel <- filter(tel, ! outlier) # drop outliers
  
  # find problematic locations
  falses <- rep(FALSE, nrow(out))
  problematic <-
    which(
      # speeds (condition, yes, no, missing)
      dplyr::if_else(out$speed < min_speed, TRUE, FALSE, TRUE) &
        dplyr::if_else(out$speed > max_speed, TRUE, FALSE, TRUE) &
        # distances
        dplyr::if_else(out$distance < min_distance, TRUE, FALSE, TRUE) &
        dplyr::if_else(out$distance > max_distance, TRUE, FALSE, TRUE) &
        # angles
        dplyr::if_else(out$angle < min_angle, TRUE, FALSE, TRUE) &
        dplyr::if_else(out$angle > max_angle, TRUE, FALSE, TRUE) &
        # sampling time
        dplyr::if_else(out$dt < min_dt, TRUE, FALSE, TRUE) &
        dplyr::if_else(out$dt > max_dt, TRUE, FALSE, TRUE)) %>%
  sort()

if(length(problematic) > 50 & ! many) {
  stop(paste(
    paste(length(problematic), ' problematic points detected.'),
    '  Check your max and min values, and set `many = TRUE` if necessary.',
    '  (It may take a while to plot all the points...)',
    sep = '\n'))
}

# find adjacent indices
adj <- sapply(problematic, \(x) (x - n_adj):(x + n_adj))
adj[adj < 1] <- 1 # only allow positive indices 
adj[adj > nrow(tel)] <- nrow(tel) # no indices greater than n(rows)

if(map) {
  adj_long <- adj %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = 'group')# %>%
  # mutate(group = substr(group, 2, nchar(group)))
  
  pal <- khroma::color('romaO')(n_distinct(adj_long$group) + 1)
  pal <- pal[- length(pal)] # last and first are the same color
  
  if(many) {
    tel[adj_long$value, ] %>%
      data.frame() %>%
      bind_cols(adj_long) %>%
      sf::st_as_sf(coords = c('location.long', 'location.lat'),
                   crs = 4326) %>%
      # view data with a map background
      mapview::mapview(col = 'black', lwd = 3, col.regions ='darkorange',
                       alpha.regions = 1, map.types = map_type)
  } else {
    tel[adj_long$value, ] %>%
      data.frame() %>%
      bind_cols(adj_long) %>%
      mutate(problematic = factor(value == problematic)) %>%
      sf::st_as_sf(coords = c('location.long', 'location.lat'),
                   crs = 4326) %>%
      # view data with a map background
      mapview::mapview(zcol = 'group',  col = 'black', lwd = 3,
                       col.regions = pal, alpha.regions = 1,
                       map.types = map_type)
  }
} else {
  # change to telemetry format (note: outliers already removed above)
  tel <- as.telemetry(tel)
  layout(t(1:2))
  
  # plot of the telemetry with errors
  plot(tel[adj, ], units = FALSE, error = FALSE, col = 'black',
       type = 'l', main = paste0(id, ': locations with DOP'))
  plot(tel[adj, ], units = FALSE, add = TRUE, col = 'black',
       transparency.error = 0)
  plot(tel[problematic, ], units = FALSE, add = TRUE, col = 'red',
       transparency.error = 0, lwd = 2)
  
  # plot of the locations with lines colored by set of points
  plot(latitude ~ longitude, data.frame(tel)[adj, ], type = 'p',
       main = paste0(id, ': locations colored by set of points'))
  apply(X = adj, MARGIN = 2, FUN = \(p) {
    lines(latitude ~ longitude, data.frame(tel)[p, ], col = p,
          lwd = 2)
  })
  points(latitude ~ longitude, data.frame(tel)[adj, ], type = 'p')
  points(latitude ~ longitude, data.frame(tel)[problematic, ],
         col = 'red', pch = 19)
  if(reset_layout) layout(1)
}
}
