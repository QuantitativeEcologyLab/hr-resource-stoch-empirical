library('ggplot2') # for fancy figures

theme_set(theme_bw() +
            theme(legend.position = 'none',
                  text = element_text(size = 15, face = 'bold'),
                  panel.grid = element_blank()))

pal <- c('#F09D00', '#3388E1', '#1C7734', '#992266')
