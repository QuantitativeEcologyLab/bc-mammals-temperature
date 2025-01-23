library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('ggplot2')   # for fancy plots
source('analysis/figures/default-ggplot-theme.R') # bold text and no grids

ANIMALS <- readr::read_csv('data/tracking-data/telemetry-metadata.csv',
                           show_col_types = FALSE) %>%
  filter(keep_for_RSF == 'yes') %>%
  pull(animal)

d <- readRDS('data/tracking-data/rsf-data.rds') %>%
  filter(animal %in% ANIMALS)

a_u <- d %>%
  mutate(elevation_m = elevation_m / 1e3,
         dist_water_m = dist_water_m / 1e3) %>%
  pivot_longer(c(forest_perc, elevation_m, dist_water_m),
               names_to = 'variable') %>%
  mutate(
    lab = case_when(species == SPECIES[1] ~ SPECIES_LABS[1],
                        species == SPECIES[2] ~ SPECIES_LABS[2],
                        species == SPECIES[3] ~ SPECIES_LABS[3],
                        species == SPECIES[4] ~ SPECIES_LABS[4],
                        species == SPECIES[5] ~ SPECIES_LABS[5],
                        species == SPECIES[6] ~ SPECIES_LABS[6],
                        species == SPECIES[7] ~ SPECIES_LABS[7]),
    variable =
      case_when(variable == 'forest_perc' ~ "bold(Forest~cover~'(%)')",
                variable == 'elevation_m' ~ "bold(Elevation~'(km)')",
                variable == 'dist_water_m' ~ "bold(Distance~from~water~'(km)')") %>%
      factor(levels = c("bold(Forest~cover~'(%)')",
                        "bold(Elevation~'(km)')",
                        "bold(Distance~from~water~'(km)')"))) %>%
  ggplot() +
  facet_grid(lab ~ variable, scales = 'free', labeller = label_parsed,
             switch = 'x') +
  geom_histogram(aes(value, fill = factor(detected)),
                 position = 'identity', alpha = 0.5, bins = 10) +
  scale_fill_brewer(NULL, type = 'qual', palette = 6,
                    labels = c('Available', 'Used')) +
  scale_y_continuous(transform = 'sqrt') +
  labs(x = NULL, y = 'Count (square-root scale)') +
  theme(legend.position = 'top', strip.background.x = element_blank(),
        strip.placement = 'outside', strip.text.x = element_text(size = 11))

ggsave('figures/available-vs-used.png',
       a_u, width = 16, height = 17, dpi = 600)

a_u_w <- d %>%
  mutate(elevation_m = elevation_m / 1e3,
         dist_water_m = dist_water_m / 1e3) %>%
  pivot_longer(c(forest_perc, elevation_m, dist_water_m),
               names_to = 'variable') %>%
  mutate(
    species = case_when(species == SPECIES[1] ~ SPECIES_LABS[1],
                        species == SPECIES[2] ~ SPECIES_LABS[2],
                        species == SPECIES[3] ~ SPECIES_LABS[3],
                        species == SPECIES[4] ~ SPECIES_LABS[4],
                        species == SPECIES[5] ~ SPECIES_LABS[5],
                        species == SPECIES[6] ~ SPECIES_LABS[6],
                        species == SPECIES[7] ~ SPECIES_LABS[7]),
    variable =
      case_when(variable == 'forest_perc' ~ "bold(Forest~cover~'(%)')",
                variable == 'elevation_m' ~ "bold(Elevation~'(km)')",
                variable == 'dist_water_m' ~ "bold(Distance~from~water~'(km)')") %>%
      factor(levels = c("bold(Forest~cover~'(%)')",
                        "bold(Elevation~'(km)')",
                        "bold(Distance~from~water~'(km)')"))) %>%
  ggplot() +
  facet_grid(species ~ variable, scales = 'free', labeller = label_parsed,
             switch = 'x') +
  geom_histogram(aes(value, fill = factor(detected), weight = weight),
                 position = 'identity', alpha = 0.5, bins = 10) +
  scale_fill_brewer(NULL, type = 'qual', palette = 6,
                    labels = c('Available', 'Used')) +
  scale_y_continuous(transform = 'sqrt') +
  labs(x = NULL, y = 'Weighted count (square-root scale)') +
  theme(legend.position = 'top', strip.background.x = element_blank(),
        strip.placement = 'outside', strip.text.x = element_text(size = 11))

ggsave('figures/available-vs-used-weighted.png',
       a_u_w, width = 16, height = 17, dpi = 600)
