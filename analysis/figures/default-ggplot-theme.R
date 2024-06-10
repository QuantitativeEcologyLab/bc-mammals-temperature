library('ggplot2') # for fancy figures

theme_set(theme_bw() +
            theme(panel.grid = element_blank(),
                  text = element_text(face = 'bold')))

PAL <- unname(c(khroma::color('bright')(6), brown = '#654321'))
