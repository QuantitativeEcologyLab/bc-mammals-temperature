library('ggplot2')
library('cowplot')

#' `get_legend()` for `{cowplot}` version 1.1.3 returns an empty plot 
get_legend <- function(.plot, position = 'top') {
  get_plot_component(.plot + theme(legend.position = position),
                     pattern = paste0('guide-box-', position),
                     return_all = TRUE)
}
