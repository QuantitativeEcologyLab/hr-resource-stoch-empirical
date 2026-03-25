library('ctmm')
library('readr')
library('dplyr')
library('tidyr')
library('purrr')
library('terra')
library('sf')

r_mean_ndvi <- rast('../ndvi-stochasticity/output/long-term-preds.tif')[[1]]
r_var_ndvi <- rast('../ndvi-stochasticity/output/long-term-preds.tif')[[2]]
r_hfi <- rast('data/hfi-raster/hfp_2021_100m_v1-2_cog.tif')

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

# HFI raster captures roads well
regina <- st_as_sf(tibble(x = -104 - 36 / 60 - 24 / 60^2,
                y = 50 + 27 /60 + 17 / 60^2), coords = c('x', 'y')) %>%
  st_set_crs("EPSG:4326") %>%
  st_buffer(15e3)
mapview::mapview(regina, col.region = 'transparent', alpha.region = 0)
mapview::mapshot(regina)

st_transform(regina, crs(r_hfi)) %>%
  crop(r_hfi, .) %>%
  project(., "EPSG:4326") %>%
  crop(regina) %>%
  plot()

d <- readRDS('data/akde-tibble-2026-02-24.rds') %>%
  # product of UD with mean NDVI and variance
  # add raster of UD to calculate exposure to mean and variance of NDVI
  mutate(ud_raster = map(akde, function(.akde) {
    if (is.null(.akde)) {
      return(NULL)
    } else {
      .r <- rast(raster(.akde, DF = 'PMF'))
      values(.r) <- values(.r) / sum(values(.r), na.rm = TRUE)
      return(.r)
    }
  }),
  raster_covariates = map(ud_raster, function(.r) {
    if (is.null(.r)) {
      tibble(mean_ndvi = NA_real_,
             var_ndvi = NA_real_,
             hfi = NA_real_)
    } else {
      box <- st_bbox(.r) %>%
        st_as_sfc() %>%
        st_as_sf() %>%
        st_buffer(1e3) %>%
        st_transform(st_crs(r_hfi))
      
      r_hfi_cropped <- crop(r_hfi, box)
      
      #' projecting `.r` biases towards 0
      tibble(mean_ndvi = sum(values(.r * project(r_mean_ndvi, .r)),
                             na.rm = TRUE),
             var_ndvi = sum(values(.r * project(r_var_ndvi, .r)),
                            na.rm = TRUE),
             hfi = sum(values(.r * project(r_hfi_cropped, .r)),
                       na.rm = TRUE)) %>%
        return()         
    }
  }, .progress = 'Calculating raster covariates...')) %>%
  unnest(raster_covariates) %>%
  # add mass and diet composition
  left_join(elton_traits, by = join_by('species' == 'scientific'))

# some rivers have excessively high HFI
st_bbox(d$ud_raster[[1]]) %>%
  st_as_sfc() %>%
  st_as_sf() %>%
  st_transform(st_crs(r_hfi)) %>%
  crop(r_hfi, .) %>%
  project("EPSG:4326") %>%
  plot()
d$tel[[1]] %>% ctmm::as.sf() %>%
  st_transform(crs(r_hfi)) %>%
  st_geometry() %>%
  st_transform("EPSG:4326") %>%
  plot(add = TRUE, col = 'white')
d$tel[[1]] %>% ctmm::as.sf() %>% mapview::mapview()

if (sum(is.na(d)) == 0) {
  d %>%
    select(! c(tel, movement_model, akde, ud_raster)) %>%
    saveRDS(paste0('data/gam-dataset-', Sys.Date(), '.rds'))
} else {
  stop(paste('There are NA values in the dataset!',
             'Save the dataset manually to confirm.'))
}
