library('ctmm')
library('readr')
library('dplyr')
library('tidyr')
library('purrr')
library('terra')
library('sf')

r_mean_ndvi <- rast('../ndvi-stochasticity/output/long-term-estimates.tif')[[1]]
r_var_ndvi <- rast('../ndvi-stochasticity/output/long-term-estimates.tif')[[2]]
r_var_ndvi <- ifel(r_var_ndvi > 0.04, 0.04, r_var_ndvi)

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
              mutate(scientific = 'Cervus canadensis')) %>%
  mutate(scientific = if_else(scientific == "Martes pennanti", # update fisher
                              "Pekania pennanti", scientific))

# TODO: update to the newest dataset
d <- readRDS('data/akde-tibble-2026-09-13.rds') %>%
  # TODO: drop this filter after updating the dataset (but rm maned wolves!)
  filter(!grepl("Chrysocyon brachyurus", animal)) %>%
  # product of UD with mean NDVI and variance
  # add raster of UD to calculate exposure to mean and variance of NDVI
  mutate(ud_raster = map(akde, function(.akde) {
    if (is.null(.akde)) return(NULL)
    .r <- rast(raster(.akde, DF = 'PMF'))
    values(.r) <- values(.r) / sum(values(.r), na.rm = TRUE)
    .r
  }),
  raster_covariates = map(ud_raster, function(.r) {
    if (is.null(.r)) {
      tibble(mean_ndvi = NA_real_,
             var_ndvi = NA_real_)
    } else {
      tibble(mean_ndvi = sum(values(.r * project(r_mean_ndvi, .r)),
                             na.rm = TRUE),
             var_ndvi = sum(values(.r * project(r_var_ndvi, .r)),
                            na.rm = TRUE))
    }
  }, .progress = 'Calculating raster covariates...')) %>%
  unnest(raster_covariates) %>%
  mutate(species = case_when(
    species == "Dasyprocta" ~ "Dasyprocta punctata",
    species == "Bettongia penicillata ogilbyi" ~ "Bettongia penicillata",
    .default = species
  )) %>%
  # add mass and diet composition
  left_join(elton_traits, by = join_by('species' == 'scientific')) %>%
  mutate(species = factor(species))

# ensure all species have data
if (sum(is.na(d$mass_kg)) > 0) {
  filter(d, is.na(mass_kg)) %>%
    summarize(n = n(),
              .by = c(species, dataset_name)) %>%
    print()
  stop('At least one row has `NA` mass!')
}

if (sum(is.na(d)) > 0) {
  cli::cli_alert_warning(
    paste('There are NA values in columns',
          paste0(colnames(d)[colSums(is.na(d)) > 0], collapse = ', '), '. ',
          'Save the dataset manually to confirm.'))
}
d %>%
  select(! c(tel, movement_model, akde, ud_raster)) %>%
  saveRDS(paste0('data/gam-dataset-', Sys.Date(), '.rds'))
