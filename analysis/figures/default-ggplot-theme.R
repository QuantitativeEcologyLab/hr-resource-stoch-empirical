library('ggplot2') # for fancy figures

theme_set(theme_bw() +
            theme(legend.position = 'none',
                  text = element_text(size = 15, face = 'bold'),
                  panel.grid = element_blank()))

pal <- c('#F09D00', '#3388E1', '#1C7734', '#992266')

ndvi_pal <- color('bukavu')(100) # wrong order
ndvi_pal <- c(ndvi_pal[c(seq(1, 40, length.out = (-0.1 + 1) * 10),
                         seq(100, 71, length.out = (0.2 + 0.1) * 10),
                         seq(70, 51, length.out = (0.9 - 0.2) * 10))],
              '#112A12')
create_ndvi_pal <- colorRampPalette(ndvi_pal)
ndvi_pal <- create_ndvi_pal(1e3)
