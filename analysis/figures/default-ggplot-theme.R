library('dplyr')   #' for `%>%`
library('ggplot2') #' for fancy figures

theme_set(theme_bw() +
            theme(panel.grid = element_blank(),
                  text = element_text(face = 'bold')))

PAL <- unname(c(khroma::color('bright')(6), brown = '#654321'))

SPECIES <- unique(readRDS('models/gamma-gam.rds')$model$species)
SPECIES_LABS <- gsub(' ', '~', SPECIES) %>%
  gsub('~mountain', ' mountain', .) %>%
  gsub('~\\(', '\\)~bold\\("(', x = .) %>%
  gsub(')$', ')"', .) %>%
  paste0('bolditalic(', ., ')')
N_SPECIES <- length(SPECIES)

COMMON_NAMES <- c('Wolves', 'Elk', 'Mountain goats', 'Cougars',
                  'Caribou (boreal)', 'Caribou (southern mountain)',
                  'Grizzly bears')

PAL <- PAL[order(COMMON_NAMES)] # rearrange based on order of common names
SPECIES <- SPECIES[order(SPECIES)]
COMMON_NAMES <- factor(COMMON_NAMES) # keep the vector in the species order
