library('readr')
library('dplyr')
library('tidyr')
library('terra')

r_mean_ndvi <- rast('data/ndvi-rasters/marcus-estimated-mean.tif')
r_var_ndvi <- rast('data/ndvi-rasters/marcus-estimated-variance.tif')

layout(1:2)
plot(r_mean_ndvi)
plot(r_var_ndvi)
layout(1)

elton_traits <-
  read_delim('data/elton-traits-1.0/MamFuncDat.txt', delim = '\t',
             show_col_types = FALSE) %>%
  rename_with(stringr::str_to_snake) %>%
  select(scientific, body_mass_value, diet_inv:diet_plant_o) %>%
  mutate(body_mass_value = body_mass_value / 1e3) %>%
  rename(mass_kg = body_mass_value) %>%
  # add a row for elk, assuming it's identical to Cervus elaphus (red deer)
  bind_rows(.,
            filter(., scientific == 'Cervus elaphus') %>%
              mutate(scientific = 'Cervus canadensis'))

d <- readRDS('data/akde-tibble-2026-02-24.rds') %>%
  # product of UD with mean NDVI and variance
  # add raster of UD to calculate exposure to mean and variance of NDVI
  mutate(ud_raster = map(akde, function(.akde) {
    if(is.null(.akde)) {
      return(NULL)
    } else {
      rast(raster(.akde, DF = 'PMF')) %>%
        project('EPSG:4326') %>%
        return()
    }
  }),
  mean_ndvi = map_dbl(ud_raster, function(.r) {
    if(is.null(.r)) {
      NA_real_
    } else {
      # NOTE: projecting after normalizing to avoid bad interpolation
      values(.r) <- values(.r) / sum(values(.r), na.rm = TRUE)
      .r <- project(.r, r_mean_ndvi)
      .r <- .r * r_mean_ndvi
      return(sum(values(.r), na.rm = TRUE))
    }
  }),
  var_ndvi = map_dbl(ud_raster, function(.r) {
    if(is.null(.r)) {
      return(NA_real_)
    } else {
      # NOTE: projecting after normalizing to avoid bad interpolation
      values(.r) <- values(.r) / sum(values(.r), na.rm = TRUE)
      .r <- project(.r, r_var_ndvi)
      .r <- .r * r_var_ndvi
      return(sum(values(.r), na.rm = TRUE))
    }
  })) %>%
  # add mass and diet composition
  left_join(elton_traits, by = join_by('species' == 'scientific'))

if(sum(is.na(d)) == 0) {
  d %>%
    select(! c(tel, movement_model, akde, ud_raster)) %>%
    saveRDS(paste0('data/gam-dataset-', Sys.Date(), '.rds'))
} else {
  stop('There are NA values in the dataset!')
}
