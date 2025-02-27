library('ggplot2')
library('cowplot')

#' `get_legend()` for `{cowplot}` version 1.1.3 returns an empty plot 
get_legend <- function(.plot) {
  get_plot_component(.plot + theme(legend.position = 'top'),
                     pattern = 'guide-box-top', return_all = TRUE)
}
