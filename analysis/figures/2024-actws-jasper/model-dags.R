library('dagitty') # for DAGs
library('ggdag')   # for DAGs with ggplot2
library('ggplot2') # for fancy plots
theme_set(theme_dag() +
            theme(text = element_text(size = 10),
                  legend.position = 'none'))

# caribou only ----
#' also see `ggdag_status()` and `ggdag_adjustment_set(shadow = TRUE)`

# P(moving)
dag_car_mov <- dagify(moving ~ temp + DOY + TOD,
                      temp ~ DOY + TOD,
                      exposure = 'temp',
                      outcome = 'moving')
set.seed(2)
ggdag_status(dag_car_mov, node_size = 20, text_size = 3) +
  scale_color_manual(values = c('black', 'darkorange'))

ggsave('figures/actws-2024-jasper/dag-caribou-moving.png',
       width = 6, height = 4, dpi = 600, bg = 'transparent')

# speed | moving
dag_car_sp <- dagify(speed ~ temp + DOY + TOD + moving,
                     # moving ~ temp + DOY + TOD,
                     temp ~ DOY + TOD,
                     exposure = 'temp',
                     outcome = 'speed')
set.seed(7)
ggdag_status(dag_car_sp, node_size = 20, text_size = 3) +
  scale_color_manual(values = c('black', 'darkorange'))

ggsave('figures/actws-2024-jasper/dag-caribou-speed.png',
       width = 6, height = 4, dpi = 600, bg = 'transparent')

# habitat selection
dag_car_hs <- dagify(selected ~ forest + water + elevation,
                     forest ~ temp,
                     water ~ temp,
                     elevation ~ temp,
                     exposure = 'temp',
                     outcome = 'selected',
                     coords = data.frame(
                       x = c(10, 0, 1, 9, 5),
                       y = c(10, 0, 9, 1, 5),
                       name = c('selected', 'temp', 'forest', 'water',
                                'elevation')))
ggdag_status(dag_car_hs, node_size = 20, text_size = 3) +
  scale_color_manual(values = c('black', 'darkorange'))

ggsave('figures/actws-2024-jasper/dag-caribou-rsf.png',
       width = 6, height = 4, dpi = 600, bg = 'transparent')
