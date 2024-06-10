library('dplyr')     # for data wrangling (mutate(), %>%, etc.)
library('tidyr')     # for data wrangling (nest(), unnest(), pivot_*, etc.)
library('purrr')     # for functional programming (map_***(), etc.)
library('ctmm')      # for movement models
library('lubridate') # for working with dates
library('mapview')   # for interactive maps
library('ggplot2')   # for fancy plots
theme_set(theme_bw())

# source custom functions
source('functions/outlier_plots.R') # to plot outlier diagnostic plots
source('functions/check_animal.R') # to run diagnostic plots
source('functions/plot_adj.R') # to plot 20 adjacent locations
source('functions/flag_outlier.R') # to mark outliers
source('functions/remove_outlier_flags.R') # to start over with an animal

# dataset checks ----
# some NA timestamps (not enough to make a difference)
readRDS('data/tracking-data/all-tracking-data-not-cleaned.rds') %>%
  filter(is.na(timestamp)) %>%
  group_by(dataset_name) %>%
  summarize(n = n())

# ensure all only have one column of HDOP
readRDS('data/tracking-data/all-tracking-data-not-cleaned.rds') %>%
  filter(! is.na(timestamp)) %>%
  mutate(outlier = if_else(outlier, 1, 0)) %>%
  filter(! is.na(hdop) + ! is.na(HDOP) > 1)

# should ony have one column of DOP
readRDS('data/tracking-data/all-tracking-data-not-cleaned.rds') %>%
  filter(! is.na(timestamp)) %>%
  mutate(outlier = if_else(outlier, 1, 0)) %>%
  filter(! is.na(DOP) + ! is.na(gps.dop) > 1)

# import the full dataset ----
# import the dataset after dropping NAs
d <- readRDS('data/tracking-data/all-tracking-data-not-cleaned.rds') %>%
  filter(! is.na(timestamp)) %>%
  mutate(hdop = if_else(is.na(hdop), HDOP, hdop), # merge error columns
         dop = if_else(is.na(DOP), gps.dop, DOP),
         original_outliers = outlier,
         outlier = if_else(outlier, 1, 0)) %>% # make numeric
  select(! c(HDOP, DOP, gps.dop)) %>% # remove duplicate columns
  relocate(dop, .after = pdop) %>%
  nest(tel = ! c(species, dataset_name, animal))

# ensure animal IDs are unique
sum(duplicated(d$animal))

# check which species have DOP values
d %>%
  unnest(tel) %>%
  group_by(species) %>%
  summarise(across(c(dop, hdop, pdop), ~ sum(.x, na.rm = TRUE))) %>%
  mutate(has_dop = dop + hdop + pdop > 0) %>%
  arrange(desc(has_dop))

# some southern mountain caribou have both hdop and pdop hdop and dop, but
# different values. other datasets are ok
d %>%
  unnest(tel) %>%
  mutate(has_hdop = ! is.na(hdop),
         has_dop = ! is.na(dop),
         has_pdop = ! is.na(pdop)) %>%
  group_by(dataset_name, has_hdop, has_dop, has_pdop) %>%
  summarise(hdop = mean(hdop, na.rm = TRUE),
            pdop = mean(pdop, na.rm = TRUE),
            dop = mean(dop, na.rm = TRUE)) %>%
  mutate(finite = map_int(1:n(),
                          ~ sum(! is.nan(c(hdop[.x], pdop[.x], dop[.x])))))

d %>%
  filter(dataset_name == 'Rangifer_tarandus_southern_mountain') %>%
  unnest(tel) %>%
  filter(! is.na(hdop) & ! is.na(dop)) %>%
  plot(hdop ~ dop, ., main = paste('n =', nrow(.)))
abline(a = 0, b = 1, col = 'red', lwd = 2)

# 0 == ok, 1 == outlier, 2 == unsure
#' *NOTE:** `as.telemetry()` removes all rows with `outlier > 0`
unique(unnest(d, tel)$outlier)

# print outlier diagnostic plots
if(FALSE) { # initial diagnostic plots
  N <- nrow(d)
  
  for(i in 1:N) {
    cat('Running animal', i, 'of', N, paste0(d$species[i], '.\n'))
    png(filename = paste0('figures/outlier-diagnostics/',
                          stringr::str_replace_all(d$species[i], '_', '-'),
                          '-', d$animal[i], '.png'),
        width = 10, height = 10, units = 'in', res = 300)
    check_animal(id = d$animal[[i]], return_out = FALSE)
    dev.off()
  }
}

# Canis lupus (clean) ----
# BW008
out <- check_animal('BW008')
plot_adj('BW008', max_speed = 3.5)
plot_adj('BW008', max_speed = 3.5, map = TRUE)
plot_adj('BW008', max_angle = 178, max_speed = 1) # using a linear feature
plot_adj('BW008', max_angle = 90, max_speed = 2) # realistic movement
plot_adj('BW008', max_angle = 170, max_speed = 0.5) # fine sampling
out[out$angle > 170 & out$speed > 0.5, ] # 15-minute sampling

# BW009
out <- check_animal('BW009')
plot_adj('BW009', max_angle = 179, max_speed = 1) # potential outlier
plot_adj('BW009', max_angle = 179.7, max_speed = 1, n_adj=200) # uses path
plot_adj('BW009', max_angle = 175, max_speed = 1, n_adj = 12) # ok
plot_adj('BW009', max_angle = 170, min_angle = 175, max_speed = 1,
         n_adj = 20) # ok
plot_adj('BW009', max_angle = 170, max_speed = 0.5) # reasonable

# BW010
out <- check_animal('BW010')
plot_adj('BW010', max_speed = 3.5) # ok
plot_adj('BW010', max_speed = 1, max_angle = 170) # ok
plot_adj('BW010', max_speed = 2, max_angle = 135) # ok

# BW012
out <- check_animal('BW012')
plot_adj('BW012', max_speed = 3)
plot_adj('BW012', max_speed = 1, max_angle = 170)
# location is ok
tel <- as.telemetry(d$tel[[which(d$animal == 'BW012')]])
i <- which(out$speed > 2 & out$angle > 170)
plot(tel, error = FALSE, type = 'l', xlim = c(-2600, 0),
     ylim = c(300, 2200), col = 'black')
plot(tel[i + -10:10, ],
     error = FALSE, type = 'l', col = 'blue', lwd = 2, add = TRUE)
plot(tel[i, ],
     error = FALSE, col = 'blue', pch = 19, add = TRUE)
plot(tel, add = TRUE)
rm(tel, i)

# BW013
out <- check_animal('BW013')
plot_adj('BW013', max_speed = 1, max_angle = 170) # uses linear features

# BW014 has a clear GPS error
out <- check_animal('BW014')
plot_adj('BW014', max_speed = 4)
flag_outlier(id = 'BW014', max_speed = 4, value = 1)
out <- check_animal('BW014')
plot_adj('BW014', max_speed = 2) # ok
plot_adj('BW014', max_speed = 0.5, min_speed = 0.6, max_angle = 170,
         n_adj = 30) # seems near a linear feature
plot_adj('BW014', max_speed = 0.5, min_speed = 0.6, max_angle = 170,
         map = TRUE) # moving along a river

# B027
out <- check_animal('BW027')
plot_adj('BW027', max_speed = 3, max_angle = 90) # not problematic

# BW028
out <- check_animal('BW028')
plot_adj('BW028', max_speed = 3) # not problematic

# B029
out <- check_animal('BW029')
plot_adj('BW029', max_speed = 2, max_angle = 175) # need a finer scale
plot_adj('BW029', max_speed = 2, max_angle = 175, n_adj = 2) # ok

# B031
out <- check_animal('BW031')
plot_adj('BW031', max_speed = 1, max_angle = 170) # need a finer scale
plot_adj('BW031', max_speed = 1, max_angle = 170, n_adj = 5) # ok

# B033: sharp and fast turns are near linear features
out <- check_animal('BW033')
plot_adj('BW033', max_speed = 1.5, max_angle = 90, n_adj = 5)

# BW044nex
out <- check_animal('BW044nex')
plot_adj('BW044nex', max_speed = 1, max_angle = 170)

# BW046nex
out <- check_animal('BW046nex')
plot_adj('BW046nex', max_speed = 1, max_angle = 175, min_speed = 3) # ok
plot_adj('BW046nex', max_speed = 1, min_speed = 3, # problematic
         max_angle = 135, min_angle = 175, reset_layout = FALSE)
i <- which(out$speed > 1 & out$speed < 3 & # likely hunting
             out$angle > 135 & out$angle < 175) + 0:5
points(location.lat ~ location.long,
       d$tel[[which(d$animal == 'BW046nex')]][i, ], col = 'blue', cex = 2)
d$tel[[which(d$animal == 'BW046nex')]][i, 'outlier'] <- 1
plot_adj('BW046nex', max_speed = 1, max_angle = 170) # ok

# Cervus elaphus (clean) ----

#' *NOTE:* assuming movement at ~< 0.25 m/s is realistic, even if angle >
#' 170 degrees since it's too hard to distinguish between movement and gps
#' error at that scale, and the models account for GPS error anyway (which
#' I assume to be 10 m since I don't have calibration data. Mike said most
#' collars have an error ~10 m anyway).

# E002
out <- check_animal('E002')
plot_adj('E002', max_speed = 0.4)
plot_adj('E002', max_angle = 170, max_speed = 0.2)

# E003
out <- check_animal('E003')
plot_adj('E003', max_speed = 0.5)

# E006
out <- check_animal('E006')
plot_adj('E006', max_speed = 0.8) # ok
plot_adj('E006', max_speed = 0.4, max_angle = 170) # outlier
plot_adj('E006', max_speed = 0.4, max_angle = 170, map = TRUE) # see river
# dt is much larger: collar struggled to receive signal with tree cover
'hours' %#% out[which(out$speed > 0.4 & out$angle > 170) + (-4:4), 'dt']
'hours' %#% c(43282, 7254)
flag_outlier('E006', max_speed = 0.4, max_angle = 170, value = 1)
out <- check_animal('E006')
plot_adj('E006', max_speed = 0.2, max_angle = 170) # ok

# E007
out <- check_animal('E007')
plot_adj('E007', max_speed = 0.2, max_angle = 170) # ok

# E008
out <- check_animal('E008')
plot_adj('E008', max_speed = 0.3) # ok

# E011: last point is outlier (keeping previous points bc many 0 speeds)
out <- check_animal('E011')
plot_adj('E011', max_speed = 0.4, n_adj = 10)
plot_adj('E011', max_speed = 0.4, n_adj = 10, map = TRUE)
plot_adj('E011', max_speed = 0.4, n_adj = 9, reset_layout = FALSE)
which(out$speed > 0.4) == nrow(out)
i_tel <- which(d$animal == 'E011')
i <- (nrow(d$tel[[i_tel]]) - 7):nrow(d$tel[[i_tel]]) # drop last 8 points
points(location.lat ~ location.long, d$tel[[i_tel]][i, ], col = 'blue',
       cex = 2, pch = 19)
d$tel[[i_tel]][i, 'outlier'] <- 1
points(location.lat ~ location.long, d$tel[[i_tel]][i, ],
       col = as.numeric(d$tel[[i_tel]]$outlier[i]) + 1, cex = 2)
layout(1)

out <- check_animal('E011')
tel <- as.telemetry(d$tel[[i_tel]], mark.rm = TRUE)
plot(tel[nrow(tel) + -10:0, ], error = FALSE, type = 'b',
     pch = as.character(0:9))
plot_adj('E011', max_speed = 0.25)

# E013
out <- check_animal('E013')
plot_adj('E013', max_speed = 0.2, max_angle = 170) # escaped something?

# E015 has multiple GPS malfunctions
out <- check_animal('E015')
plot_adj('E015', max_speed = 1, reset_layout = FALSE)
# can't filter easily because the outliers don't have the max velocity
i <- d %>%
  pull(tel) %>%
  nth(which(d$animal == 'E015')) %>%
  mutate(i = 1:n()) %>%
  slice(which(timestamp ==
                as.POSIXct(out[which(out$speed > 0.6)[2], 't'])) +
          0:1) %>%
  pull(i) %>%
  first() + 0:10
points(location.lat ~ location.long, cex = 2, lwd = 3, col = 'blue',
       filter(d, animal == 'E015')$tel[[1]][i, ])
lines(location.lat ~ location.long,
      filter(d, animal == 'E015')$tel[[1]][i, ], lwd = 6, col = 'red3')
d$tel[[which(d$animal == 'E015')]][i, 'outlier'] <- 1
out <- check_animal('E015')
plot_adj('E015', max_speed = 0.60) # track looks ok now

# E016
out <- check_animal('E016')
plot_adj('E016', max_speed = 0.2, max_angle = 170) # one outlier
plot_adj('E016', max_speed = 0.2, max_angle = 179) # use higher upper bound
flag_outlier('E016', max_speed = 0.2, max_angle = 179, value = 1)
out <- check_animal('E016')
plot_adj('E016', max_speed = 0.2, max_angle = 170) # ok

# E017
out <- check_animal('E017')
plot_adj('E017', max_speed = 0.6) # ok
plot_adj('E017', max_angle = 170, max_speed = 0.2)

# E019
out <- check_animal('E019') # many high speeds at angle > 170
plot_adj('E019', max_speed = 0.65, max_angle = 135) # some outliers
plot_adj('E019', max_speed = 0.65, max_angle = 160, reset_layout = FALSE)
i <- which(out$speed > 0.65 & out$angle > 160)[c(1, 3)] # some outliers
points(location.lat ~ location.long,
       d$tel[[which(d$animal == 'E019')]][i, ], cex = 2, col = 'blue')
d$tel[[which(d$animal == 'E019')]][i, 'outlier'] <- 1

out <- check_animal('E019') # check other high speeds with angle > 160
plot_adj('E019', max_speed = 0.4, max_angle = 135, reset_layout = FALSE)
tel <- d$tel[[which(d$animal == 'E019')]]
i <- which(out$speed > 0.4 & out$angle > 135)
points(location.lat ~ location.long, cex = 1.5, col = 'black',
       tel[! tel$outlier, ][i, ],
       pch = as.character(0:(length(i) - 1)))
'hour' %#% out[i, 'dt'] # consistent sampling

# 0 is likely an error (erring on the side of caution)
# 1 is consistent with other movement
# 2 is likely an error
# 3 is consistent with other movement
# 4 is likely an error
# 5 is realistic
# 6 is likely an error
# 7 is realistic
# 8 is likely an error
# 9 is likely an error
i <- i[1 + c(0, 2, 4, 6, 8, 9)]
i_id <- which(d$animal == 'E019')
d$tel[[i_id]][! d$tel[[i_id]]$outlier, ][i, 'outlier'] <- 1
# ensure correct locations were removed
out <- check_animal('E019') # check other high speeds with angle > 160
plot_adj('E019', max_speed = 0.4, max_angle = 135, reset_layout = FALSE)

# E023
out <- check_animal('E023')
plot_adj('E023', max_speed = 1) # ok
plot_adj('E023', max_speed = 0.6, min_speed = 1) # ok
plot_adj('E023', max_speed = 0.3, max_angle = 170) # ok

# E024
out <- check_animal('E024')
plot_adj('E024', max_speed = 0.4) # ok
plot_adj('E024', max_speed = 0.2, max_angle = 170) # ok

# E025
out <- check_animal('E025')
plot_adj('E025', max_speed = 0.3, max_angle = 170) # ok

# E026
#' this elk does not have an unusual amount of data
hist(map_int(filter(d, species == 'Cervus elaphus')$tel, nrow))
abline(v = nrow(d$tel[[which(d$animal == 'E026')]]), lwd = 2)
d$tel[[which(d$animal == 'E026')]] %>%
  filter(! outlier) %>%
  sf::st_as_sf(coords = c('location.long', 'location.lat'),
               crs = 4326) %>%
  mapview(map.types = 'Esri.WorldImagery', col.regions = 'darkorange',
          col = 'black', lwd = 2)
out <- check_animal('E026') # many sharp turns at high speeds
plot_adj('E026', max_speed = 0.6, max_angle = 135)
plot_adj('E026', max_speed = 0.6, max_angle = 135, map = TRUE) # look north
flag_outlier('E026', max_speed = 0.6, max_angle = 135, value = 1)

out <- check_animal('E026')
plot_adj('E026', max_speed = 0.4, max_angle = 170) # NW one is ok
plot_adj('E026', max_speed = 0.4, max_angle = 170, map = TRUE)
plot_adj('E026', max_speed = 0.42, max_angle = 170) # W one may also be ok
plot_adj('E026', max_speed = 0.42, max_angle = 175) # W one may also be ok
flag_outlier('E026', max_speed = 0.42, max_angle = 175, value = 1)

out <- check_animal('E026')
plot_adj('E026', max_speed = 0.35, max_angle = 170) # NW, middle, & NE ok?
plot_adj('E026', max_speed = 0.35, min_speed = 0.4, max_angle = 170,
         reset_layout = FALSE)
i <- with(out, which(speed > 0.35 & speed < 0.4 & angle > 160))[c(3:4)]
points(location.lat ~ location.long,
       filter(d$tel[[which(d$animal == 'E026')]], ! outlier)[i, ],
       cex = 2, col = 'blue', lwd = 2)
#' values in `i` assume outliers have been removed
i_tel <- which(d$animal == 'E026')
d$tel[[i_tel]][d$tel[[i_tel]]$outlier == 0, ][i, 'outlier'] <- 1

out <- check_animal('E026')
plot_adj('E026', max_speed = 0.30, max_angle = 160) # likely ok

#' using an error-informed, continuous-time movement model substantially
#' reduces the issues caused by GPS error. Additionally, there are so many
#' elk data (for each individual and the number of elk, too) that a few
#' values from a single individual will have little to no effect.

# E027
out <- check_animal('E027')
plot_adj('E027', max_speed = 0.3, max_angle = 170, n_adj = 4) # not all
plot_adj('E027', max_speed = 0.3, max_angle = 175, n_adj = 4,
         reset_layout = FALSE)
i <- with(out, which(speed > 0.3 & angle > 175))[c(2, 3, 4, 5)]
points(location.lat ~ location.long,
       d$tel[[which(d$animal == 'E027')]][i, ],
       cex = 2, col = 'blue', lwd = 2)
'hours' %#% out$dt[i]
d$tel[[which(d$animal == 'E027')]][i, 'outlier'] <- 2 # not fully sure

# E030
out <- check_animal('E030')
plot_adj('E030', max_speed = 0.45) # ok
plot_adj('E030', max_speed = 0.3, max_angle = 135) # generally ok
plot_adj('E030', max_speed = 0.3, max_angle = 170) # possibly problematic
flag_outlier('E030', max_speed = 0.3, max_angle = 170, value = 1)

out <- check_animal('E030')
plot_adj('E030', max_speed = 0.2, max_angle = 170) # ok

# E032
out <- check_animal('E032')
plot_adj('E030', max_speed = 0.3, max_angle = 170) # ok

# E033
out <- check_animal(id = 'E033')
plot_adj('E033', max_speed = 0.6)

# E036
out <- check_animal(id = 'E036')
# the next point is along the same line 
plot_adj('E036', max_speed = 0.7, reset_layout = FALSE)
points(location.lat ~ location.long, pch = 19,
       filter(d, animal == 'E036')$tel[[1]][0:1 + which(out$speed > 0.7),])
plot_adj('E036', max_speed = 0.6, reset_layout = FALSE) # ok

# E037
out <- check_animal(id = 'E037')
plot_adj(id = 'E037', max_speed = 0.3, max_angle = 170) # ok

# E040
out <- check_animal(id = 'E040')
plot_adj('E040', max_speed = 0.35, max_angle = 160) # ok

# E042: range shift with highly stationary behavior
out <- check_animal(id = 'E042')
tel <- as.telemetry(d$tel[[which(d$animal == 'E042')]])
range(tel$timestamp)
plot(tel)
plot_adj(id = 'E042', many = TRUE, map = TRUE, max_speed = 0)
ctmm.guess(tel)

# E049
out <- check_animal(id = 'E049')
plot_adj('E049', max_speed = 0.4, max_angle = 170) # ok
plot_adj('E049', max_speed = 0.3, min_speed = 0.4, max_angle = 170) # ok

# E067
out <- check_animal(id = 'E067')
plot_adj('E067', max_speed = 0.6, max_angle = 170) # ok

# E068
out <- check_animal(id = 'E068')
plot_adj('E068', max_speed = 0.65, reset_layout = FALSE)
out[out$speed > 0.65, ] # 1st 2 points are right after capture (rownames)
points(location.lat ~ location.long,
       filter(d, animal == 'E068')$tel[[1]][1:2, ], col = 'blue', cex = 2)
d[d$animal == 'E068', ]$tel[[1]]$outlier[1:2] <- 1
out <- check_animal(id = 'E068')
plot_adj('E068', max_speed = 0.4, max_angle = 170)
plot_adj('E068', max_speed = 0.4, max_angle = 170, map = TRUE)
flag_outlier('E068', max_speed = 0.4, max_angle = 170, value = 1)
out <- check_animal(id = 'E068')
plot_adj('E068', max_speed = 0.25, max_angle = 170) # ok

# E070
out <- check_animal('E070')
plot_adj('E070', max_speed = 0.4, max_angle = 135) # ok

# E074
out <- check_animal(id = 'E074')
plot_adj('E074', max_speed = 0.4, max_angle = 170) # ok

# E075
out <- check_animal(id = 'E075')
plot_adj('E075', max_speed = 0.3, max_angle = 170) # ok

# E079
out <- check_animal(id = 'E079')
plot_adj('E079', max_speed = 0.6) # ok

# E080
out <- check_animal(id = 'E080')
plot_adj('E080', max_speed = 0.8) # ok
plot_adj('E080', max_speed = 0.5, max_angle = 170, reset_layout = FALSE)
i_tel <- which(d$animal == 'E080')
i <- which(out$speed > 0.5 & out$angle > 170)[1:2]
points(location.lat ~ location.long,
       filter(d$tel[[i_tel]], ! outlier)[i, ], col = 'blue',
       cex = 2)
d$tel[[i_tel]][! d$tel[[i_tel]]$outlier, ]$outlier[i] <- 1
out <- check_animal(id = 'E080')
plot_adj('E080', max_speed = 0.5, max_angle = 170) # ok

# E081
out <- check_animal(id = 'E081')
plot_adj('E081', max_speed = 0.2, max_angle = 170) # ok

# E082
out <- check_animal(id = 'E082')
plot_adj('E082', max_speed = 0.2, max_angle = 170) # ok

# E083
out <- check_animal(id = 'E083')
plot_adj('E083', max_speed = 0.3, max_angle = 135) # ok

# E084 is not worth keeping
out <- check_animal(id = 'E084')
plot(location.lat ~ location.long, filter(d, animal == 'E084')$tel[[1]],
     type = 'b', pch = 19, col = '#00000010')
nrow(filter(d, animal == 'E084')$tel[[1]])
range(filter(d, animal == 'E084')$tel[[1]]$timestamp)
d <- filter(d, animal != 'E084')

# E085
out <- check_animal(id = 'E085')
plot_adj('E085', max_speed = 1) # definitely outliers (> 20 km!)
flag_outlier(id = 'E085', max_speed = 1, value = 1)

out <- check_animal(id = 'E085')
plot_adj('E085', max_speed = 0.6) # maybe ok
with(out[which(out$speed > 0.6) + -10:10, ], # likely outlier
     plot('hours' %#% dt, dt * speed / 1e3,
          col = as.numeric(-10:10 == 0) + 1,
          pch = 19, ylab = 'SLD (dt * min speed, km)'))
flag_outlier('E085', max_speed = 0.6, value = 1)

out <- check_animal(id = 'E085')
plot_adj('E085', max_speed = 0.2, max_angle = 135) # ok

# E096
out <- check_animal(id = 'E096')
plot_adj('E096', max_speed = 0.4, max_angle = 170) # ok

# E100
out <- check_animal(id = 'E100')
plot_adj('E100', max_speed = 0.3, max_angle = 170) # ok

# E101
out <- check_animal(id = 'E101')
plot_adj('E101', max_speed = 0.4, max_angle = 170) # ok

# E106
d$tel[[which(d$animal == 'E106')]] %>%
  group_by(gps.fix.type.raw) %>%
  summarise(hdop = mean(dop, na.rm = TRUE))
d$tel[[which(d$animal == 'E106')]] %>%
  group_by(gps.fix.type.raw) %>%
  summarise(hdop = mean(! is.na(dop), na.rm = TRUE))
out <- check_animal(id = 'E106')
plot_adj('E106', max_speed = 0.4, max_angle = 170) # large DOP values; ok
d$tel[[which(d$animal == 'E106')]] %>%
  as.telemetry() %>%
  data.frame() %>%
  slice(which(out$speed > 0.4 & out$angle > 170)) %>%
  pull(HDOP)
plot_adj('E106', max_speed = 0.3, max_angle = 170) # ok
plot_adj('E106', max_speed = 0.2, max_angle = 170) # ok

# E115
# collar is very error-prone in SW forest
out <- check_animal(id = 'E115')
plot_adj('E115', many = TRUE, map = TRUE)
plot_adj('E115', max_speed = 0.4, max_angle = 170) # improbable
plot_adj('E115', max_speed = 0.4, max_angle = 170, n_adj = 50)#but possible
plot_adj('E115', max_speed = 0.4, max_angle = 175) # unlikely
slice(d$tel[[which(d$animal == 'E115')]],
      which(out$speed > 0.4 & out$angle > 175))$dop # low DOP
flag_outlier(id = 'E115', max_speed = 0.4, max_angle = 175, value = 2)

out <- check_animal('E115')
plot_adj('E115', max_speed = 0.3, max_angle = 170) # unlikely
flag_outlier('E115', max_speed = 0.3, max_angle = 170, value = 2)

out <- check_animal('E115')
plot_adj('E115', max_speed = 0.1, max_angle = 170, many = TRUE) # ok

# E116
out <- check_animal(id = 'E116')
plot_adj('E116', max_speed = 0.4, max_angle = 170)
flag_outlier(id = 'E116', max_speed = 0.4, max_angle = 170, value = 1)

out <- check_animal('E116')
plot_adj('E116', max_speed = 0.8)
flag_outlier(id = 'E116', max_speed = 0.8, value = 1)

out <- check_animal('E116')
plot_adj('E116', max_speed = 0.3, max_angle = 170) # improbable
plot_adj('E116', max_speed = 0.3, max_angle = 175) # but not impossible

# E117 has similar error patterns to E115
out <- check_animal(id = 'E117')
plot_adj('E117', max_speed = 0.7) # three outliers
plot_adj('E117', max_speed = 0.7, map = TRUE)
flag_outlier(id = 'E117', max_speed = 0.7, value = 1)

out <- check_animal('E117')
plot_adj('E117', many = TRUE, map = TRUE) # more errors when more forested
plot_adj('E117', max_speed = 0.4, max_angle = 170)

plot_adj('E117', max_speed = 0.4, max_angle = 170, map = TRUE)
flag_outlier(id = 'E117', max_speed = 0.4, max_angle = 170, value = 1)

out <- check_animal('E117')
plot_adj('E117', max_speed = 0.3, max_angle = 170) # maybe not all
plot_adj('E117', max_speed = 0.3, max_angle = 175) # likely all
plot_adj('E117', max_speed = 0.3, max_angle = 175, map = TRUE) # in forest
flag_outlier(id = 'E117', max_speed = 0.3, max_angle = 175, value = 1)

out <- check_animal('E117')
# still an excessive amount of sharp turns with high speeds (compare E118)
# dropping some of the data is better than dropping the whole animal
flag_outlier('E117', max_speed = 0.1, max_angle = 130, value = 2)
out <- check_animal('E117')

# E118 (compare with E117 above)
check_animal('E118', return_out = FALSE)

# E122
out <- check_animal(id = 'E122')
plot_adj('E122', max_speed = 0.45) # ok

# E125
out <- check_animal(id = 'E125')
plot_adj('E125', max_speed = 0.4, max_angle = 150)
flag_outlier('E125', max_speed = 0.4, max_angle = 150, value = 1)
out <- check_animal(id = 'E125') # ok

# E129
out <- check_animal(id = 'E129')
plot_adj('E129', max_speed = 0.6) # ok

# E131
out <- check_animal(id = 'E131')
plot_adj('E131', max_speed = 0.3, max_angle = 170) # N one is ok
plot_adj('E131', max_speed = 0.45, max_angle = 170) # outlier
flag_outlier('E131', max_speed = 0.45, max_angle = 170, value = 1)

out <- check_animal(id = 'E131') # ok

# E132
#' some high minimum speeds at long sampling intervals for this individual
out <- check_animal(id = 'E132', cap_dt = FALSE)
layout(matrix(c(0, 0, 0, 1), ncol = 2)) #
abline(v = 20, col = 'red')
abline(h = 0.08, col = 'red')
i <- which(out$dt > 20 %#% 'hours' & out$speed > 0.08)
out[i, ]
# small range shift at the daily scale within home range
plot_adj('E132', max_speed = 0.08, max_dt = 20 %#% 'hours')

# look at points near values with excessively high speeds
i_adj <- sapply(i, function(.i) .i + -10:10) %>%
  as.data.frame() %>% # convert matrix to a data frame
  pivot_longer(! c()) # use long format for a column of all values
tel <- as.telemetry(d$tel[[which(d$animal == 'E132')]])
adj <- left_join(slice(data.frame(out), i_adj$value),
                 slice(data.frame(tel), unique(i_adj$value)),
                 by = c('t')) %>%
  bind_cols(i_adj, .) %>%
  mutate(date = as.POSIXct(t),
         high_speed = speed > 0.08)
layout(matrix(1:4, ncol = 2))
for(g in seq(unique(i_adj$name))) {
  g_name <- unique(i_adj$name)[g]
  plot(y ~ x, filter(adj, name == g_name), type = 'l')
  points(y ~ x, data.frame(tel[i, ][g, ]), col = 'blue', pch = 19)
  rm(g, g_name)
} # ok
layout(1)

# large dt with minor "range shifts" ==> ok
plot('days' %#% dt ~ as.POSIXct(t), out[i_adj$value, ])
points('days' %#% dt ~ as.POSIXct(t), out[i, ], col = 'blue', pch = 19)

rm(i, i_adj, tel, adj)

# E134
out <- check_animal(id = 'E134')
plot_adj('E134', max_speed = 0.5, max_angle = 170) # ok

# E138
out <- check_animal(id = 'E138')
plot_adj('E138', max_speed = 0.8) # ok
plot_adj('E138', max_speed = 0.3, max_angle = 170) # ok

# E139
out <- check_animal(id = 'E139')
plot_adj('E139', max_speed = 0.3, max_angle = 170) # ok

# E140
out <- check_animal(id = 'E140')
plot_adj('E140', max_speed = 0.6)
flag_outlier('E140', max_speed = 0.6, value = 1)

out <- check_animal(id = 'E140') # ok

# E142
out <- check_animal(id = 'E142')
flag_outlier('E142', max_speed = 10, value = 1) # definitely an outlier

out <- check_animal(id = 'E142')
plot_adj('E142', max_speed = 0.8) # ok

# E143
out <- check_animal(id = 'E143')
plot_adj('E143', max_speed = 0.2, max_angle = 170) # ok

# E144
out <- check_animal(id = 'E144')
plot_adj('E144', max_speed = 0.3, max_angle = 170) # all ok

# E147
out <- check_animal(id = 'E147')
plot_adj('E147', max_speed = 0.3, max_angle = 170) # ok

# E161
out <- check_animal(id = 'E161')
plot_adj('E161', max_speed = 0.3, max_angle = 170) # ok

# E164: first point is definitely an outlier
out <- check_animal(id = 'E164')
tel <- d$tel[[which(d$animal == 'E164')]]
plot(location.lat ~ location.long, tel, type = 'b')
points(location.lat ~ location.long, tel[1, ], col = 'red', pch = 19)
rm(tel)
d$tel[[which(d$animal == 'E164')]]$outlier[1] <- 1

out <- check_animal(id = 'E164')
plot_adj('E164', max_speed = 0.3, max_angle = 160) # ok
rm(tel)

# E168
out <- check_animal(id = 'E168')
plot_adj('E168', max_speed = 0.3, max_angle = 170)
flag_outlier('E168', max_speed = 0.3, max_angle = 170, value = 1)

out <- check_animal(id = 'E168')
plot_adj('E168', max_speed = 0.2, max_angle = 170) # ok

# E170
out <- check_animal(id = 'E170')
plot_adj('E170', max_speed = 1)
plot_adj('E170', max_speed = 0.3, max_angle = 170) # ok

# E171
out <- check_animal(id = 'E171')
plot_adj('E171', max_speed = 0.6) # ok
plot_adj('E171', max_speed = 0.3, max_angle = 170) # ok

# E172
out <- check_animal(id = 'E172')
plot_adj('E172', max_speed = 0.3, max_angle = 150) # ok

# E178
out <- check_animal(id = 'E178')
plot_adj('E178', max_speed = 0.3, max_angle = 150) # ok

# E179
out <- check_animal(id = 'E179')
plot_adj('E179', max_speed = 0.8) # ok

# E181
out <- check_animal(id = 'E181')
plot_adj('E181', max_speed = 0.5, max_angle = 170) # ok

# E181
out <- check_animal(id = 'E182')
plot_adj('E182', max_speed = 0.3, max_angle = 170) # ok

# Oreamnos americanus (clean) ----
d %>%
  filter(species == 'Oreamnos americanus') %>%
  unnest(tel) %>%
  pull(hdop) %>%
  hist(main = 'Goat HDOP before as.telemetry()')

d %>%
  filter(species == 'Oreamnos americanus') %>%
  unnest(tel) %>%
  pull(pdop) %>%
  hist(main = 'Goat PDOP before as.telemetry()')

d %>%
  filter(species == 'Oreamnos americanus') %>%
  mutate(tel = map(tel, \(.t) data.frame(as.telemetry(.t)))) %>%
  unnest(tel) %>%
  pull(HDOP) %>%
  hist(main = 'Goat HDOP after as.telemetry()')

# pdop is not in the telemetries
d %>%
  filter(species == 'Oreamnos americanus') %>%
  mutate(tel = map(tel, \(.t) data.frame(as.telemetry(.t)))) %>%
  unnest(tel) %>%
  colnames()

# substitute missing HDOP values with PDOP values
d <- d %>%
  mutate(
    tel = if_else(species == 'Oreamnos americanus',
                  map(tel, \(.t) {
                    .t %>%
                      mutate(hdop = if_else(is.na(hdop), pdop, hdop)) %>%
                      return()
                  }),
                  tel))

# ensure values are ok
all((d %>%
       filter(species == 'Oreamnos americanus') %>%
       unnest(tel) %>%
       pull(hdop)) ==
      (d %>%
         filter(species == 'Oreamnos americanus') %>%
         mutate(tel = map(tel, \(.t) data.frame(as.telemetry(.t)))) %>%
         unnest(tel) %>%
         pull(HDOP)))

# all goats have some extreme outliers in Berlin
filter(d, species == 'Oreamnos americanus') %>%
  unnest(tel) %>%
  ggplot() +
  facet_wrap(~ animal) +
  geom_hex(aes(location.long, location.lat), bins = 10) +
  geom_vline(xintercept = 0, color = 'darkorange') +
  scale_fill_continuous(limits = c(1, NA), trans = 'log10') +
  theme(legend.position = c(0.8, 0.15))

filter(d, species == 'Oreamnos americanus') %>%
  unnest(tel) %>%
  filter(location.long < 0) %>%
  pull(location.long) %>%
  hist()
abline(v = -119.8, col = 'darkorange')

# set all locations with long > 0 to outliers
d <-
  mutate(
    d,
    tel = if_else(species == 'Oreamnos americanus',
                  true =
                    map(tel, \(.d) {
                      mutate(
                        .d,
                        outlier = if_else(location.long > -119.8, 1, 0),
                        original_outliers = outlier == 1)
                    }),
                  false = tel))

# 30548
out <- check_animal(id = '30548')
hist(out$angle) # a lot of sharp turns at low speeds
plot_adj('30548', max_speed = 0.4) # ok
plot_adj('30548', max_speed = 0.3, max_angle = 135) # ok
plot_adj('30548', max_speed = 0.3, max_angle = 135, map = TRUE) # ridges
plot_adj('30548', max_speed = 0.2, max_angle = 135) # ok
plot_adj('30548', max_speed = 0.1, max_angle = 170, n_adj = 3) # ok

# 30551
out <- check_animal(id = '30551') # another 124, many with speed = 0
# start and end are definitely not part of the goat's movement
d$tel[[which(d$animal == '30551')]] %>%
  filter(! outlier) %>%
  ggplot(aes(location.long, location.lat, col = as.POSIXct(timestamp))) +
  geom_path() +
  geom_point() +
  scale_color_viridis_c()
flag_outlier('30551', max_distance = 2e4, value = 1)

out <- check_animal(id = '30551') # last point is still a clear outlier
plot_adj('30551', max_speed = 0.6)
which(out$speed > 0.6) == nrow(out)
flag_outlier('30551', max_speed = 0.6, value = 1)

out <- check_animal(id = '30551') # all other speeds are reasonable
plot_adj('30551', max_speed = 0.4, n_adj = 20) # ok

# 30561
out <- check_animal(id = '30561')
plot_adj('30561', max_distance = 15e3, map = TRUE)
flag_outlier('30561', max_distance = 15e3, value = 1)

out <- check_animal(id = '30561') # some high speeds at sharp turns
plot_adj('30561', max_speed = 0.2, max_angle = 160)
plot_adj('30561', max_speed = 0.2, max_angle = 160, map = TRUE) # outliers
flag_outlier('30561', max_speed = 0.2, max_angle = 160, value = 1)

out <- check_animal(id = '30561') # speeds are reasonable
plot_adj('30561', max_speed = 0, map = TRUE, many = TRUE) # ok

# 30567
out <- check_animal(id = '30567')
plot_adj('30567', max_speed = 0.3, max_angle = 160)
flag_outlier('30567', max_speed = 0.3, max_angle = 160, value = 1)

out <- check_animal(id = '30567')
plot_adj('30567', max_speed = 0.3) # ok
plot_adj('30567', max_speed = 0.2, max_angle = 160)
flag_outlier('30567', max_speed = 0.2, max_angle = 160, value = 1)

out <- check_animal(id = '30567')
plot_adj('30567', max_speed = 0.15, max_angle = 160) # NE is ok
plot_adj('30567', max_speed = 0.15, max_angle = 170)
plot_adj('30567', max_speed = 0.15, max_angle = 170, n_adj = 3)
plot_adj('30567', max_speed = 0.15, max_angle = 170, n_adj = 3, map = TRUE)
flag_outlier('30567', max_speed = 0.15, max_angle = 170, value = 1)

out <- check_animal(id = '30567') # some possible outliers, but unclear
plot_adj('30567', max_speed = 0, many = TRUE, map = TRUE)
plot_adj('30567', max_speed = 0.1, max_angle = 170) # ok
out[out$speed > 0.1 & out$angle > 170, ]
plot_adj('30567', max_speed = 0.1, max_angle = 170, max_dt = 2e4) # ok

# 30575
out <- check_animal(id = '30575')
plot_adj('30575', max_speed = 0.8)
flag_outlier('30575', max_speed = 0.8, value = 1)

out <- check_animal(id = '30575')
plot_adj('30575', max_speed = 0.5)
flag_outlier('30575', max_speed = 0.5, value = 1)

out <- check_animal(id = '30575')
plot_adj('30575', max_speed = 0.3, max_angle = 170)
flag_outlier('30575', max_speed = 0.3, max_angle = 170, value = 1)

out <- check_animal(id = '30575')
plot_adj('30575', max_speed = 0.2, max_angle = 170)
flag_outlier('30575', max_speed = 0.2, max_angle = 170, value = 1)

out <- check_animal(id = '30575')
plot_adj('30575', max_speed = 0.175, max_angle = 170)
flag_outlier('30575', max_speed = 0.175, max_angle = 170, value = 1)

out <- check_animal(id = '30575')
plot_adj('30575', max_speed = 0.15, max_angle = 170)
plot_adj('30575', max_speed = 0.15, max_angle = 170, map = TRUE)
flag_outlier('30575', max_speed = 0.15, max_angle = 170, value = 1)

out <- check_animal(id = '30575')
plot_adj('30575', max_speed = 0.12, max_angle = 170)
plot_adj('30575', max_speed = 0.12, max_angle = 170, max_dt = 2e4)
flag_outlier('30575', max_speed = 0.12, max_angle = 170, max_dt = 2e4,
             value = 1)

out <- check_animal(id = '30575') # some points in forest, but speed is low
plot_adj('30575', many = TRUE, map = TRUE)
plot(as.telemetry(d$tel[[which(d$animal == '30575')]], mark.rm = TRUE))

# 30599
out <- check_animal(id = '30599')
plot_adj(id = '30599', max_speed = 0.3, max_angle = 170)
flag_outlier(id = '30599', max_speed = 0.3, max_angle = 170, value = 1)

out <- check_animal(id = '30599')
plot_adj(id = '30599', max_speed = 0.2, max_angle = 170)
flag_outlier(id = '30599', max_speed = 0.2, max_angle = 170, value = 1)

out <- check_animal(id = '30599')
plot_adj(id = '30599', max_speed = 0.3)
flag_outlier(id = '30599', max_speed = 0.3, value = 1)

out <- check_animal(id = '30599')
plot_adj(id = '30599', max_speed = 0.15, max_angle = 170)
flag_outlier(id = '30599', max_speed = 0.15, max_angle = 170, value = 1)

out <- check_animal(id = '30599')
plot_adj(id = '30599', max_speed = 0.1, max_angle = 175)
out[out$angle > 175 & out$speed > 0.1, ]
plot_adj(id = '30599', max_speed = 0.1, max_angle = 175, min_dt = 2e4)
i <- which(out$speed > 0.1 & out$angle > 175 & out$dt < 2e4)
tel <- d$tel[[which(d$animal == '30599')]]
d$tel[[which(d$animal == '30599')]][! tel$outlier, ][i, 'outlier'] <- 1

out <- check_animal(id = '30599') # ok

# 30613
out <- check_animal(id = '30613')
plot_adj('30613', max_speed = 0.4, max_angle = 170)
plot_adj('30613', max_speed = 0.4, max_angle = 170, map = TRUE)
flag_outlier('30613', max_speed = 0.4, max_angle = 170, value = 1)

out <- check_animal(id = '30613')
plot_adj('30613', max_speed = 0.2, max_angle = 170)
flag_outlier('30613', max_speed = 0.2, max_angle = 170, value = 1)

out <- check_animal(id = '30613')
plot_adj('30613', max_speed = 0.2, max_angle = 135)
plot_adj('30613', max_speed = 0.2, max_angle = 135, map = TRUE)
flag_outlier('30613', max_speed = 0.2, max_angle = 135, value = 1)

out <- check_animal(id = '30613')
plot_adj('30613', max_speed = 0.12, max_angle = 170)
flag_outlier('30613', max_speed = 0.12, max_angle = 170, value = 1)

out <- check_animal(id = '30613')
plot_adj('30613', max_speed = 0.1, max_angle = 170) # ok
plot_adj('30613', max_speed = 0.1, max_angle = 150) # ok
plot_adj('30613', max_speed = 0.05, max_angle = 170) # ok

# 30636
out <- check_animal(id = '30636')
flag_outlier('30636', max_distance = 15e3, value = 1)

out <- check_animal(id = '30636')
plot_adj('30636', max_speed = 0.4)
flag_outlier('30636', max_speed = 0.4, value = 1)

out <- check_animal(id = '30636')
plot_adj('30636', max_speed = 0.1, max_angle = 170) # ok

# 30642
out <- check_animal(id = '30642')
plot_adj('30642', max_speed = 1)
flag_outlier('30642', max_speed = 1, value = 1)

out <- check_animal(id = '30642')
plot_adj('30642', max_speed = 0.3, max_angle = 170)
out[out$speed > 0.3 & out$angle > 170, ]
plot_adj('30642', max_speed = 0.3, max_angle = 170, max_dt = 8100)
flag_outlier('30642', max_speed = 0.3, max_angle = 170, max_dt = 8100,
             value = 1)

out <- check_animal(id = '30642')
plot_adj('30642', max_speed = 0.2, max_angle = 170, min_speed = 0.3)
i_tel <- which(d$animal == '30642')
i <- which(out$speed > 0.2 & out$speed < 0.3 & out$angle > 170)
d$tel[[i_tel]][! d$tel[[i_tel]]$outlier, ][i, 'outlier'] <- 1
plot_adj('30642', max_speed = 0.2, max_angle = 170, min_speed = 0.3) # ok

out <- check_animal(id = '30642') # ok

# 30648
out <- check_animal(id = '30648')
flag_outlier('30648', max_speed = 1, value = 1)

out <- check_animal(id = '30648')
plot_adj('30648', max_speed = 0.3, max_angle = 170)
flag_outlier('30648', max_speed = 0.3, max_angle = 170, value = 1)

out <- check_animal(id = '30648')
plot_adj('30648', max_speed = 0.2, max_angle = 170) # ok
plot_adj('30648', max_speed = 0.4, max_angle = 130)
flag_outlier('30648', max_speed = 0.4, max_angle = 130, value = 1)

out <- check_animal(id = '30648')
plot_adj('30648', max_speed = 0.2, max_angle = 170)
plot_adj('30648', max_speed = 0.2, max_angle = 170, map = TRUE)
plot_adj('30648', max_speed = 0.2, max_angle = 175)
flag_outlier('30648', max_speed = 0.2, max_angle = 175, value = 1)

out <- check_animal(id = '30648')
plot_adj('30648', max_speed = 0.2, max_angle = 170) # ok
plot_adj('30648', max_speed = 0.1, max_angle = 175) # ok

# save goat data as a separate file
goats <- d %>%
  filter(species == 'Oreamnos americanus') %>%
  unnest(tel) %>%
  select(! c(pdop, dop, gps.fix.type.raw, NAV, gps.satellite.count))
write.csv(goats, 'data/tracking-data/goat-data.csv', row.names = FALSE)

# Puma concolor (clean) ----
# 246
out <- check_animal(id = '246')
plot_adj('246', max_speed = 3)
flag_outlier('246', max_speed = 3, value = 1)

out <- check_animal(id = '246')
plot_adj('246', max_speed = 0.5, n_adj = 3) # ok
plot_adj('246', max_speed = 0.2, max_angle = 170) # ok
plot_adj('246', max_speed = 0.2, max_angle = 170, n_adj = 4, map = TRUE)
# most southern point may be an outlier, but it looks ok
plot_adj('246', max_speed = 0.2, max_angle = 177, min_angle = 178) # ok
plot_adj('246', max_speed = 0.2, max_angle = 177, min_angle = 178,
         n_adj = 4, map = TRUE)
i <- which(out$speed > 0.2 & out$angle > 177 & out$angle < 178) + -5:5
tel <- d$tel[[which(d$animal == '246')]]
tel <- tel[! tel$outlier, ]
plot(location.lat ~ location.long, tel[i, ], type = 'b',
     pch = as.character(seq(length(i))))
'hours' %#% out[i, 'dt']
rm(i, tel, cols)

# 272 has a clear speed outlier, but locations are ok
out <- check_animal(id = '272')
plot_adj('272', max_speed = 1) # data seem ok
plot_adj('272', max_speed = 0.2, max_angle = 135) # ok

# C1
out <- check_animal(id = 'C1')
plot_adj('C1', max_speed = 0.2, max_angle = 170) # ok

# C2
out <- check_animal(id = 'C2')
plot_adj('C2', max_speed = 0.2, max_angle = 170, n_adj = 20) # gps error
plot_adj('C2', max_speed = 0.2, max_angle = 170, n_adj = 20, map = TRUE)
flag_outlier('C2', max_speed = 0.2, max_angle = 170, value = 1)

# C3
out <- check_animal(id = 'C3')
plot_adj('C3', max_speed = 1)
plot_adj('C3', max_speed = 1, map = TRUE)
flag_outlier('C3', max_speed = 1, value = 1)
out <- check_animal(id = 'C3') # remaining data are ok
plot_adj('C3', max_speed = 0.4, max_angle = 170, n_adj = 5) # ok
plot_adj('C3', max_speed = 0.3, max_angle = 170, n_adj = 5) # ok

# C4
out <- check_animal(id = 'C4')
plot_adj('C4', max_speed = 0.64)
plot_adj('C4', max_speed = 0.64, map = TRUE)
flag_outlier('C4', max_speed = 0.64, value = 1)
out <- check_animal(id = 'C4') # remaining data are ok
plot_adj('C4', max_speed = 0.2, max_angle = 170) # ok

# C5
out <- check_animal(id = 'C5')
plot_adj('C5', max_speed = 0.6) # ok
plot_adj('C5', max_speed = 0.4, max_angle = 170) # ok
plot_adj('C5', max_speed = 0.3, max_angle = 170) # ok

# C6
out <- check_animal(id = 'C6')
plot_adj('C6', max_speed = 0.4, max_angle = 170)
flag_outlier('C6', max_speed = 0.4, max_angle = 170, value = 1)
out <- check_animal(id = 'C6') # ok

# C7
out <- check_animal(id = 'C7')
plot_adj('C7', max_speed = 0.6, max_angle = 135)
flag_outlier('C7', max_speed = 0.6, max_angle = 135, value = 1)
out <- check_animal(id = 'C7') # ok
plot_adj('C7', max_speed = 0.4, max_angle = 135) # ok

# C8
out <- check_animal(id = 'C8')
plot_adj('C8', max_speed = 0.55) # ok
plot_adj('C8', max_speed = 0.55, map = TRUE)
'hours' %#% out[which(out$speed > 0.55) + -5:5, 'dt']
plot_adj('C8', max_speed = 0.2, max_angle = 170) # ok

# C9
out <- check_animal(id = 'C10')
plot_adj('C10', max_speed = 0.2, max_angle = 170) # ok

# C10
out <- check_animal(id = 'C9')
plot_adj('C9', max_speed = 0.2, max_angle = 170) # ok

# C11
out <- check_animal(id = 'C11')
plot_adj('C11', max_speed = 2)
flag_outlier('C11', max_speed = 2, value = 1)

out <- check_animal(id = 'C11')
plot_adj('C11', max_speed = 0.4, max_angle = 170) # NE ok
plot_adj('C11', max_speed = 0.4, max_angle = 175, map = TRUE) # outlier
'hours' %#% out[which(out$speed > 0.4 & out$angle > 175) + -3:3, 'dt']
flag_outlier('C11', max_speed = 0.4, max_angle = 175, value = 1)

out <- check_animal(id = 'C11')
plot_adj('C11', max_speed = 0.3, max_angle = 170) # NE ok

# C12
out <- check_animal(id = 'C12')
plot_adj('C12', max_speed = 0.8)
i <- which(out$speed > 0.8) + c(-1:0)
layout(t(0:1))
points(location.lat ~ location.long, d$tel[[which(d$animal == 'C12')]][i,],
       cex = 2, col = 'blue')
out$dt[i] # distance seems excessive for 2-hour sampling
d$tel[[which(d$animal == 'C12')]][i, 'outlier'] <- 2
out <- check_animal(id = 'C12')
plot_adj(id = 'C12', max_speed = 0.2, max_angle = 170) # ok
rm(i)

# c13
out <- check_animal(id = 'C13')
plot_adj('C13', max_speed = 0.3, max_angle = 170) # ok

# C14
out <- check_animal('C14')
plot_adj('C14', max_speed = 0.4, max_angle = 160) # ok

# C16
out <- check_animal('C16')
plot_adj('C16', max_speed = 0.3, max_angle = 170) # all outliers
flag_outlier('C16', max_speed = 0.3, max_angle = 170, value = 1)
out <- check_animal('C16') # ok

# C17
out <- check_animal('C17')
plot_adj('C17', max_speed = 0.2, max_angle = 170)
plot_adj('C17', max_speed = 0.35, max_angle = 170)
# one of the two is the last point of the telemetry
out[which((out$angle > 170 | is.na(out$angle)) & out$speed > 0.35), ]
i <- which(out$angle > 170 & out$speed > 0.35)
out[i, ]
# don't want to flag last point
d$tel[[which(d$animal == 'C17')]]$outlier[i] <- 1
out <- check_animal('C17')
plot_adj('C17', max_speed = 0.35, max_angle = 179) # ok

# C18
out <- check_animal('C18')
plot_adj('C18', max_speed = 0.6) # the two speeds are too similar
flag_outlier('C18', max_speed = 0.6, value = 2)
out <- check_animal('C18')
plot_adj('C18', max_speed = 0.2, max_angle = 170) # ok

# C19
out <- check_animal('C19')
plot_adj('C19', max_angle = 160, max_speed = 0.4) # NW is ok
plot_adj('C19', max_speed = 0.4, max_angle = 170)
flag_outlier('C19', max_speed = 0.4, max_angle = 170, value = 1)

out <- check_animal('C19') # ok

# C23
out <- check_animal('C23')
plot_adj('C23', max_speed = 0.3, max_angle = 160) # S is ok
plot_adj('C23', max_speed = 0.3, max_angle = 170) # unsure
plot_adj('C23', max_speed = 0.3, max_angle = 170, n_adj = 3)
out[which(out$angle > 170 & out$speed > 0.3) + -3:3, 'dt'] # too regular
plot_adj('C23', max_speed = 0.3, max_angle = 170, map = TRUE) # outlier
flag_outlier('C23', max_speed = 0.3, max_angle = 170, value = 1)
out <- check_animal('C23')

# C28
out <- check_animal('C28')
plot_adj('C28', max_speed = 0.3, max_angle = 170)
flag_outlier('C28', max_speed = 0.3, max_angle = 170, value = 1)
out <- check_animal('C28')
plot_adj('C28', max_speed = 0.25, max_angle = 170) # ok

# C29
out <- check_animal('C29')
plot_adj('C29', max_speed = 0.5, max_angle = 160)
flag_outlier('C29', max_speed = 0.5, max_angle = 160, value = 1)
out <- check_animal('C29')

# C30
out <- check_animal('C30')
plot_adj('C30', max_speed = 0.6)
flag_outlier('C30', max_speed = 0.6, value = 1)

out <- check_animal('C30')
plot_adj('C30', max_speed = 0.3, max_angle = 170)
flag_outlier('C30', max_speed = 0.3, max_angle = 170, value = 1)

out <- check_animal('C30')
plot_adj('C30', max_speed = 0.3, max_angle = 135) # ok

# C31
out <- check_animal('C31')
plot_adj('C31', max_speed = 0.37, max_angle = 160)
out[which(out$speed > 0.37 & out$angle > 160) + -3:3, 'dt'] # too regular
plot_adj('C31', max_speed = 0.37, max_angle = 160, map = TRUE) # outlier
flag_outlier('C31', max_speed = 0.37, max_angle = 160, value = 1)
out <- check_animal('C31') # ok

# Boreal caribou (clean) ----

# HSCEK063
out <- check_animal('HSCEK063')
plot_adj('HSCEK063', max_speed = 0.6) # ok
plot_adj('HSCEK063', max_speed = 0.2, max_angle = 135) # ok

# SCEK007b
out <- check_animal('SCEK007b')
plot_adj('SCEK007b', max_speed = 0.6, n_adj = 25) # ok

# SCEK009
out <- check_animal('SCEK009')
plot_adj('SCEK009', max_speed = 0.4, max_angle = 135)
plot_adj('SCEK009', max_speed = 0.4, max_angle = 135, n_adj = 5)
plot_adj('SCEK009', max_speed = 0.4, max_angle = 135, n_adj = 30) # ok

# SCEK020b migrates
out <- check_animal('SCEK020b')
plot(as.telemetry(d$tel[[which(d$animal == 'SCEK020b')]]))
plot_adj('SCEK020b', max_speed = 0, many = TRUE, map = TRUE)

# SCEK026b
out <- check_animal('SCEK026b')
plot_adj('SCEK026b', max_speed = 0.5, n_adj = 20) # ok

# SCEK033
out <- check_animal('SCEK033')
plot_adj('SCEK033', max_speed = 1) # ok

# SCEK036
out <- check_animal('SCEK036')
plot_adj('SCEK036', max_speed = 0.6) # ok

# SCEK044
out <- check_animal('SCEK044')
plot_adj('SCEK044', max_speed = 1) # ok

# SCEK047
out <- check_animal('SCEK047')
plot_adj('SCEK047', max_speed = 1) # ok

# SCEK050
out <- check_animal('SCEK050')
plot_adj('SCEK050', max_speed = 0.2, max_angle = 135) # ok

# SCEK079
out <- check_animal('SCEK079')
plot_adj('SCEK079', max_speed = 0.7) # ok

# SCEK079b
out <- check_animal('SCEK079b', cap_dt = FALSE)
plot_adj('SCEK079b', max_speed = 0.4) # ok
plot_adj('SCEK079b', max_speed = 0.15, max_angle = 170) # 2/4 outliers
plot_adj('SCEK079b', max_speed = 0.23, max_angle = 170) # possible outliers
flag_outlier('SCEK079b', max_speed = 0.23, max_angle = 170, value = 2)
out <- check_animal('SCEK079b', cap_dt = FALSE)

# SCEK114
out <- check_animal('SCEK114', cap_dt = FALSE)
plot_adj('SCEK114', max_speed = 0.4, max_angle = 160) # ok

# SCEK119
out <- check_animal('SCEK119', cap_dt = FALSE)
plot_adj('SCEK119', max_speed = 0.6, max_angle = 160) # outlier
flag_outlier('SCEK119', max_speed = 0.6, max_angle = 160, value = 2)
out <- check_animal('SCEK119', cap_dt = FALSE)

# SCEK133
out <- check_animal('SCEK133', cap_dt = FALSE)
plot_adj('SCEK133', max_speed = 0.9, max_angle = 170) # ok

# SCEK136
out <- check_animal('SCEK136', cap_dt = FALSE)
plot_adj('SCEK136', max_speed = 0.6, max_angle = 170) # unsure
plot_adj('SCEK136', max_speed = 0.6, max_angle = 170, n_adj = 40) # outlier
plot_adj('SCEK136', max_speed = 0.6, max_angle = 170, n_adj = 40, map = T)
flag_outlier('SCEK136', max_speed = 0.6, max_angle = 170, value = 1)
out <- check_animal('SCEK136', cap_dt = FALSE)

# SCEK161b
out <- check_animal('SCEK161b', cap_dt = FALSE)
plot_adj('SCEK161b', max_speed = 0.4, max_angle = 135) # ok

# SCEK163
out <- check_animal('SCEK163', cap_dt = FALSE)
plot_adj('SCEK163', max_speed = 0.8) # ok

# SCEK169
out <- check_animal('SCEK169', cap_dt = FALSE)
hist(out$speed)
plot_adj('SCEK169', max_speed = 10, n_adj = 1)

# SCEK170
out <- check_animal('SCEK170', cap_dt = FALSE)
hist(out$speed)
plot_adj('SCEK170', max_speed = 150, n_adj = 2)
plot_adj('SCEK170', max_speed = 2)
flag_outlier('SCEK170', max_speed = 2, value = 1)

out <- check_animal('SCEK170', cap_dt = FALSE)
plot_adj('SCEK170', max_speed = 2)
flag_outlier('SCEK170', max_speed = 2, value = 1)

out <- check_animal('SCEK170', cap_dt = FALSE)
plot_adj('SCEK170', max_speed = 2)
flag_outlier('SCEK170', max_speed = 2, value = 1)

out <- check_animal('SCEK170', cap_dt = FALSE)
plot_adj('SCEK170', max_speed = 1)
flag_outlier('SCEK170', max_speed = 1, value = 1)

out <- check_animal('SCEK170', cap_dt = FALSE) # ok

# SCEK176
out <- check_animal('SCEK176')
plot_adj('SCEK176', max_speed = 0.6) # ok

# SCEK177
out <- check_animal('SCEK177')
plot_adj('SCEK177', max_speed = 0.2, max_angle = 170)
plot_adj('SCEK177', max_speed = 0.2, max_angle = 170, n_adj = 5)
flag_outlier('SCEK177', max_speed = 0.2, max_angle = 170, value = 1)
out <- check_animal('SCEK177', cap_dt = FALSE)
which.max(out$dt) # ok
nrow(out)

# SCEK181
out <- check_animal('SCEK181', cap_dt = FALSE)
hist(out$speed)
plot_adj('SCEK181', max_speed = 2, many = TRUE)
flag_outlier('SCEK181', max_speed = 2, value = 1)

out <- check_animal('SCEK181', cap_dt = FALSE)
plot_adj('SCEK181', max_speed = 2)
flag_outlier('SCEK181', max_speed = 2, value = 1)

out <- check_animal('SCEK181', cap_dt = FALSE)
plot_adj('SCEK181', max_speed = 2)
flag_outlier('SCEK181', max_speed = 2, value = 1)

out <- check_animal('SCEK181', cap_dt = FALSE)
plot_adj('SCEK181', max_speed = 0.5) # ok

# SCEK185
out <- check_animal('SCEK185')
plot_adj('SCEK185', max_speed = 0.3) # ok

# SCEK194
out <- check_animal('SCEK194')
plot_adj('SCEK194', max_speed = 0.15, max_angle = 170, # one outlier
         reset_layout = FALSE)
i <- which(out$angle > 170 & out$speed > 0.15)[2]
points(location.lat ~ location.long,
       d$tel[[which(d$animal == 'SCEK194')]][i, ],
       col = 'blue', cex = 2)
d$tel[[which(d$animal == 'SCEK194')]]$outlier[i] <- 1
out <- check_animal('SCEK194')

# SCEK201
out <- check_animal('SCEK201')
range('hours' %#% out$dt, na.rm = TRUE)
plot_adj('SCEK201', max_speed = 0.3, max_angle = 170, n = 50)
flag_outlier('SCEK201', max_speed = 0.3, max_angle = 170, value = 1)

out <- check_animal('SCEK201')
plot_adj('SCEK201', max_speed = 0.15, max_angle = 170, n_adj = 30)
flag_outlier('SCEK201', max_speed = 0.15, max_angle = 170, value = 1)

out <- check_animal('SCEK201')
plot_adj('SCEK201', max_speed = 0.1, max_angle = 170, n_adj = 30) # ok

# SCEK204
out <- check_animal('SCEK204')
plot_adj('SCEK204', max_angle = 170, max_speed = 0.1) # ok

# SCEK206
out <- check_animal('SCEK206')
plot_adj('SCEK206', max_distance = 15000)
out[out$distance > 15000, ] # last two points are outliers
nrow(out)
flag_outlier('SCEK206', max_distance = 15000, value = 1)
out <- check_animal('SCEK206')

# SCEK212
out <- check_animal('SCEK212')
plot_adj('SCEK212', max_speed = 0.25, max_angle = 150) # ok

# SCEK213
out <- check_animal('SCEK213')
plot_adj('SCEK213', max_speed = 0.25, max_angle = 150) # ok

# SCEK214
out <- check_animal('SCEK214')
plot_adj('SCEK214', max_speed = 0.25, max_angle = 150) # ok

# SCEK216
out <- check_animal('SCEK216')
plot_adj('SCEK216', max_speed = 0.4, max_angle = 135) # ok

# SCEK238
out <- check_animal('SCEK238')
plot_adj('SCEK238', max_speed = 0.5) # ok

# Southern mountain caribou (clean) ----
filter(d, dataset_name == 'Rangifer_tarandus_southern_mountain') %>%
  unnest(tel) %>%
  pull(hdop) %>%
  hist(main = 'Mountain caribou HDOP')

#' some rowns have both`hdop` and `dop` values
filter(d, dataset_name == 'Rangifer_tarandus_southern_mountain') %>%
  unnest(tel) %>%
  mutate(missing_hdop = is.na(hdop),
         missing_dop = is.na(dop)) %>%
  group_by(dataset_name, missing_hdop, missing_dop) %>%
  summarise(hdop = mean(hdop, na.rm = TRUE),
            pdop = mean(pdop, na.rm = TRUE),
            dop = mean(dop, na.rm = TRUE))

#' some have substantially different non-NA values for `hdop` and `dop`
#' assuming `ctmm` takes `hdop` and that `hdop` is better than `dop`
filter(d, dataset_name == 'Rangifer_tarandus_southern_mountain') %>%
  unnest(tel) %>%
  filter(! is.na(hdop) & ! is.na(dop)) %>%
  transmute(diff = hdop - dop) %>%
  pull(diff) %>%
  hist()

# check date ranges
filter(d, dataset_name == 'Rangifer_tarandus_southern_mountain') %>%
  unnest(tel) %>%
  group_by(animal) %>%
  summarise(range = diff(range(as.POSIXct(timestamp)))) %>%
  pull(range) %>%
  as.numeric() %>%
  hist(xlab = 'Days')

# 92_prepenned
ID <- '92_prepenned'
out <- check_animal(ID)
plot_adj(ID, max_speed = 0.3)

# 101_prepenned
ID <- '101_prepenned'
out <- check_animal(ID)
plot_adj(ID, max_speed = 0.5)
which(out$speed > 0.5)
d$tel[[which(d$animal == ID)]][1, 'outlier'] <- 1 # ran after collaring
out <- check_animal(ID)
plot_adj(ID, max_speed = 0.2, max_angle = 170) # ok
plot_adj(ID, max_speed = 0.1, min_speed = 0.2, max_angle = 170) # ok

# 104_prepenned
ID <- '104_prepenned'
out <- check_animal(ID)
plot_adj(ID, max_speed = 0.4) # ok

# 105_prepenned
ID <- '105_prepenned'
out <- check_animal(ID)
plot_adj(ID, max_speed = 0.1, max_angle = 150) # ok

# 109_prepenned
ID <- '109_prepenned'
out <- check_animal(ID)
plot_adj(ID, max_speed = 0.07, max_angle = 170) # ok

# 117_prepenned
ID <- '117_prepenned'
out <- check_animal(ID)
plot_adj(ID, max_speed = 0.1) # ok

# 132_prepenned
ID <- '132_prepenned'
out <- check_animal(ID)
diff(range(as.POSIXct(out$t))) # only 15.5 days of data
plot_adj(ID, max_speed = 0.08) # ok
plot('hours' %#% out$dt)
plot_adj(ID, max_dt = 12 %#% 'hours') # ok

# 152_prepenned
ID <- '152_prepenned'
out <- check_animal(ID)
plot_adj(ID, max_speed = 0.1, max_angle = 135) # ok

# Ursus arctos horribilis (clean) ----
# 197
out <- check_animal('197')
plot_adj('197', max_speed = 0.3, max_angle = 170) # ok

# 242: first point is an outlier
out <- check_animal('242')
plot_adj('242', max_speed = 1, n_adj = 10)
which(out$speed > 1)
d$tel[[which(d$animal == '242')]][1, 'outlier'] <- 1
out <- check_animal('242')

# 259
out <- check_animal('259')
plot_adj('259', max_speed = 0.4, max_angle = 170, n_adj = 10)
flag_outlier('259', max_speed = 0.4, max_angle = 170, value = 1)
out <- check_animal('259')
plot_adj('259', max_speed = 0.2, max_angle = 135, n_adj = 5) # ok

# 260
out <- check_animal('260')
plot_adj('260', max_speed = 1, max_angle = 170, n_adj = 10)
flag_outlier('260', max_speed = 1, max_angle = 170, value = 1)
out <- check_animal('260')
plot_adj('260', max_speed = 0.4, max_angle = 170, n_adj = 10) # ok
plot_adj('260', max_speed = 0.3, max_angle = 170, n_adj = 10, # unsure
         min_speed = 0.4)
filter(data.frame(out), round(speed, 1) == 0.4 & angle > 170)
0.3509840 * 3645 / 1e3 # reasonable distance in an hour

# 262
out <- check_animal('262')
plot_adj('262', max_speed = 0.3, max_angle = 170) # ok
plot_adj('262', max_speed = 0.3, max_angle = 135) # ok

# 291
out <- check_animal('291')
plot_adj('291', max_speed = 0.3, max_angle = 170) # ok

# 292
out <- check_animal('292')
plot_adj('292', max_speed = 1.2, n_adj = 10)
flag_outlier(id = '292', max_speed = 1.2, value = 1)
out <- check_animal('292')
plot_adj('292', max_speed = 0.8) # ok
plot_adj('292', max_speed = 0.8, map = TRUE) # moving along a river

# 297
out <- check_animal('297')
plot_adj('297', max_speed = 0.5) # ok

# 317
out <- check_animal('317')
plot_adj('317', max_speed = 1.5, reset_layout = FALSE) # previous = outlier
i <- which(out$speed > 1.5) - 1
points(location.lat ~ location.long,d$tel[[which(d$animal == '317')]][i, ],
       cex = 2, col = 'blue')
d$tel[[which(d$animal == '317')]][i, 'outlier'] <- 1
out <- check_animal('317')
plot_adj('317', max_speed = 0.8) # ok
plot_adj('317', max_speed = 0.3, max_angle = 170) # possible outliers...
plot_adj('317', max_speed = 0.3, max_angle = 160) # but common at 0.3 m/s

# 318
out <- check_animal('318')
plot_adj('318', max_speed = 0.5) # ok
plot_adj('318', max_speed = 0.4, min_speed = 0.5) # ok

# to save progress ----
saveRDS(object = d,
        file = paste0('data/tracking-data/all-tracking-data-cleaned-',
                      format(Sys.time(), '%Y-%m-%d-%H-%M'),
                      '.rds'))
