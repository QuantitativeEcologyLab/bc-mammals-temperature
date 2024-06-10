library('dplyr') # for data wrangling (mutate(), %>%, etc.)
library('tidyr') # for data wrangling (unnest(), pivot_*, etc.)
library('purrr') # for functional programming (map_***(), etc.)
library('furrr') # for futures and parallel computing
library('ctmm')  # for movement models

#' `https://github.com/ctmm-initiative/ctmmlearn/blob/main/ctmm_error.R`
# goats are the only species with calibration data
goat_uere <-
  as.telemetry('data/tracking-data/goat-calibration-data.csv') %>%
  uere.fit()
summary(goat_uere) # only one error class

# some goats had mortality events
goat_info <-
  read.csv('data/tracking-data/oreamnos-americanus/goat_info.csv') %>%
  transmute(collar_id = as.character(collar_id),
            goat_id,
            capture_date = as.Date(capture_date),
            mortality_date = as.Date(mortality_date)) %>%
  as_tibble()

# running models separately to avoid having to re-run previous code
plan(multisession, workers = 10)
if(FALSE) {
  d <-
    readRDS('data/tracking-data/all-tracking-data-cleaned-2024-02-22-13-49.rds') %>%
    mutate(dataset_name =
             factor(dataset_name,
                    levels = c('Rangifer_tarandus_boreal',
                               'Canis_lupus_boreal',
                               'Rangifer_tarandus_southern_mountain',
                               'Elk in southwestern Alberta',
                               'Oreamnos_americanus',
                               'Puma_concolor_2',
                               'Puma_concolor_4',
                               'Ursus_arctos_horribilis'))) %>%
    arrange(dataset_name) %>% # to fit boreal caribou first
    left_join(goat_info, by = c('animal' = 'collar_id')) %>%
    # one collar used for two separate goats, but the correct animal ID was lost
    mutate(animal = case_when(
      species != 'Oreamnos americanus' ~ animal,
      animal == '30613' & mortality_date == Sys.Date() ~ 'CA12',
      animal == '30613' & ! is.na(mortality_date) ~ 'CA03',
      species == 'Oreamnos americanus' ~ goat_id),
      tel = imap(tel, \(.tel, .i) {
        if(species[.i] == 'Oreamnos americanus') {
          # remove data prior to capture
          .tel <- filter(.tel, as.Date(timestamp) > capture_date[.i])
          
          # remove data post mortality (if occurred)
          .date <- mortality_date[.i]
          .tel <- filter(.tel, as.Date(timestamp) < mortality_date[.i])
        }
        return(.tel)
      })) %>%
    mutate(tel = imap(tel, \(.tel, .i) {
      if(is.na(mortality_date[.i])) {
        .tel
      } else {
        filter(.tel, as.Date(timestamp) < mortality_date[.i])
      }
    })) %>%
    select(- goat_id, - mortality_date) %>%
    mutate(
      # drop unnecessary columns for goat data
      tel = imap(tel, \(x, .i) {
        if(species[.i] == 'Oreamnos americanus') {
          x <- select(x, - c(fix, pdop, dop, gps.fix.type.raw))
        }
        
        return(x)
      }),
      # convert data frames to telemetry objects
      tel = map(tel, \(x) as.telemetry(x, mark.rm = TRUE)),
      # add calibration-informed error or a reasonable guess
      tel = if_else(condition = species == 'Oreamnos americanus',
                    true = map(tel, \(x) {
                      uere(x) <- goat_uere
                      return(x)
                    }),
                    false = map(tel, \(x) {
                      # not changing default errors: assuming errors of 10 m
                      prior <- uere(x) # extract calibration object
                      prior$DOF[] <- 2 # low DOF for large uncertainty
                      uere(x) <- prior # assign prior to the data
                      return(x)
                    })),
      variogram = future_map(tel, \(x) ctmm.guess(data = x,
                                                  CTMM = ctmm(error = TRUE),
                                                  interactive = FALSE),
                             .progress = TRUE,
                             .options = furrr_options(seed = TRUE)),
      # movement model file name
      mm_file_name = paste0('models/movement-models/ctmm-', dataset_name,
                            '-', animal, '-', Sys.Date(), '.rds'))
} else {
  d <- readRDS('models/movement-models-2024-03-03.rds')
}

# fit movement models ----
if(! dir.exists('models/movement-models')) {
  dir.create('models/movement-models/')
}

# import models to see which ones are still missing
d <- mutate(d,
            movement_model = future_map(mm_file_name, \(.fn) {
              if(file.exists(.fn)) {
                return(readRDS(.fn))
              } else {
                return()
              }
            }))

# find missing models
missing_models <- which(map_chr(d$movement_model, class) == 'NULL')

# number of models left to fit
length(missing_models)

# proportion of models fit
round(1 - length(missing_models) / nrow(d), 2)

# fit missing models
future_map(.x = missing_models,
           \(i) ctmm.select(data = d$tel[[i]], CTMM = d$variogram[[i]]) %>%
             saveRDS(d$mm_file_name[[i]]),
           .progress = TRUE,
           .options = furrr_options(seed = TRUE))

# add new models
d <- mutate(d,
            movement_model = future_map(mm_file_name, \(.fn) {
              if(file.exists(.fn)) {
                return(readRDS(.fn))
              } else {
                return()
              }
            }))

# save progress
saveRDS(d, paste0('models/movement-models-', Sys.Date(), '.rds'))

# fit AKDEs ----
if(! dir.exists('models/akdes')) dir.create('models/akdes')

d <- mutate(d,
            akde_file_name = paste0('models/akdes/akde-', dataset_name,
                                    '-', d$animal, '-', Sys.Date(), '.rds'))

# find missing models
missing_akdes <- which(! file.exists(d$akde_file_name))

# number of akdes left to fit
length(missing_akdes)

# proportion of akdes fit
round(1 - length(missing_akdes) / nrow(d), 2)

future_map(.x = missing_akdes,
           \(i) akde(data = d$tel[[i]], CTMM = d$movement_model[[i]],
                     weights = TRUE) %>%
             saveRDS(d$akde_file_name[[i]]),
           .progress = TRUE)

d <- mutate(d,
            akde = future_map(akde_file_name, readRDS))

saveRDS(d, paste0('models/movement-models-akdes-', Sys.Date(), '.rds'))

# check computation time ----
T_END <- Sys.time()

T_END - T_START

# save goat data for Aimee Chhen ----
filter(d, species == 'Oreamnos americanus') %>%
  select(animal, species, tel, variogram, movement_model, akde) %>%
  saveRDS(paste0('models/movement-models-akdes-goats-', Sys.Date(), '.rds'))
