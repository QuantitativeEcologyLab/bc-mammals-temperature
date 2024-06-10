library(dplyr)
library(ggplot2)
theme_set(theme_bw())

wolf <- readxl::read_xlsx('data/tracking-data/webexportwolf.xlsx') %>%
  mutate(timestamp = as.POSIXct(paste(FixDate, FixTime),
                                format = "%Y-%m-%d %H:%M:%S"))
wolf

caribou <- readxl::read_xlsx('data/tracking-data/webexportcaribou.xlsx') %>%
  mutate(timestamp = as.POSIXct(paste(FixDate, FixTime),
                                format = "%Y-%m-%d %H:%M:%S"))
caribou

n_distinct(wolf$FieldID)
n_distinct(caribou$FieldID)

cowplot::plot_grid(
  ggplot() +
    geom_line(aes(FixDateTime, Collar, group = FieldID), wolf) +
    ggtitle('Wolves'),
  ggplot() +
    geom_line(aes(FixDateTime, Collar, group = FieldID), caribou) +
    ggtitle('Caribou'))
