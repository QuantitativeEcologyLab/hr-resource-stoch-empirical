library('readr') # for importing data
library('dplyr') # for data wrangling (mutate(), %>%, etc.)
library('tidyr') # for data wrangling (unnest(), pivot_*, etc.)
library('purrr') # for functional programming (map_***(), etc.)
library('furrr') # for futures and parallel computing
library('ctmm')  # for movement models

#' **UPDATE CTMM BEFORE RUNNING ANY MODELS**

# running models separately to avoid having to re-run previous code
plan(multisession, workers = future::availableCores(logical = FALSE) - 2)

d <-
  readRDS('data/telemetry-data/cleaned-telemetry-data-2025-07-22.rds') %>%
  mutate(
    animal_id = factor(animal_id),
    # # add calibration-informed error or a reasonable guess
    # tel = case_when(TRUE ~ map(tel, \(x) {
    #   # not changing default errors: assuming errors of 10 m
    #   prior <- uere(x) # extract calibration object
    #   prior$DOF[] <- 2 # low DOF for large uncertainty
    #   uere(x) <- prior # assign prior to the data
    #   return(x)
    # })),
    variogram = future_map(tel, \(x) ctmm.guess(data = x,
                                                CTMM = ctmm(error = TRUE),
                                                interactive = FALSE) %>%
                             list(),
                           .progress = TRUE,
                           .options = furrr_options(seed = TRUE)),
    # movement model file name
    mm_file_name = paste0('models/movement-models/ctmm-', dataset_name,
                          '-', animal_id, '.rds'))

# fit movement models ----
if(! dir.exists('models/movement-models')) {
  dir.create('models/movement-models')
}

# import models to see which ones are still missing
d <- mutate(d,
            movement_model = future_map(mm_file_name, \(.fn) {
              if(file.exists(.fn)) {
                return(readRDS(.fn))
              } else {
                return()
              }
            }))

# find missing models
missing_models <- which(map_chr(d$movement_model, class) == 'NULL')
length(missing_models) # number of missing models
round(length(missing_models) / nrow(d), 2) # proportion of missing models

# fit missing models
future_map(.x = missing_models,
           \(i) ctmm.select(data = d$tel[[i]], CTMM = d$variogram[[i]]) %>%
             saveRDS(d$mm_file_name[[i]]),
           .progress = TRUE,
           .options = furrr_options(seed = TRUE))

# add new models
d <- mutate(d,
            movement_model = future_map(mm_file_name, \(.fn) {
              if(file.exists(.fn)) {
                return(readRDS(.fn))
              } else {
                return()
              }
            }))

# save progress
saveRDS(d, paste0('models/movement-models-', Sys.Date(), '.rds'))

# fit AKDEs ----
if(! dir.exists('models/akdes')) dir.create('models/akdes')

d <- mutate(d,
            akde_file_name = paste0('models/akdes/akde-', dataset_name,
                                    '-', d$animal, '.rds'))

# find missing AKDEs
missing_akdes <- which(! file.exists(d$akde_file_name))
length(missing_akdes)
round(1 - length(missing_akdes) / nrow(d), 2)

future_map(.x = missing_akdes,
           \(i) akde(data = d$tel[[i]], CTMM = d$movement_model[[i]],
                     weights = TRUE) %>%
             saveRDS(d$akde_file_name[[i]]),
           .progress = TRUE)

d <- mutate(d, akde = future_map(akde_file_name, readRDS))

saveRDS(d, paste0('models/movement-models-akdes-', Sys.Date(), '.rds'))

# check computation time ----
T_END <- Sys.time()

T_END - T_START

#' **temporary file****************************************************
d <- readRDS('../bc-mammals-temperature/models/movement-models-akdes-2024-06-06.rds') %>%
  rename(animal_id = animal) %>%
  mutate(species = gsub(' horribilis', '', species))

d <- readRDS('models/movement-models-akdes-2025-07-XX.rds')

# make a spreadsheet of metadata ----
# EltonTraits 1.0 dataset: https://esajournals.onlinelibrary.wiley.com/doi/10.1890/13-1917.1
#' `diet_vend` = endotherm, `diet_vect` = ectotherm; see `metadata.htm`
traits <- read_delim('data/elton-traits-1.0/MamFuncDat.txt', delim = '\t',
                     show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  select(scientific, body_mass_value, diet_inv:diet_plant_o) %>%
  mutate(body_mass_value = body_mass_value / 1e3) %>%
  rename(mass_kg = body_mass_value) %>%
  # add a row for elk, assuming it's identical to red deer
  bind_rows(.,
            filter(., scientific == 'Cervus elaphus') %>%
              mutate(scientific = 'Cervus canadensis'))

d %>%
  arrange(paste0(dataset_name, '-', animal_id)) %>%
  transmute(
    animal_id, species, dataset_name, # keep some columns
    n_locations,
    median_sampling_freq_hours = map_dbl(tel, \(.t) {
      median(diff(.t$timestamp, units = 'hours'))
    }),
    sampling_duration_days = n_days,
    model_class = map_chr(movement_model, \(.m) {
      return(summary(.m)$name)
    }),
    has_speed = grepl('OUF', toupper(model_class)) |
      grepl('IOU', model_class),
    range_resident = '') %>%
  left_join(traits, by = join_by('species' == 'scientific')) %>%
  write.csv('data/telemetry-data/telemetry-metadata.csv', row.names = FALSE)

# ensure all species have a mass estimate
read_csv('data/telemetry-data/telemetry-metadata.csv') %>%
  filter(is.na(mass_kg))

# check range residency for each animal
imap(d$animal_id, \(.a, .i) {
  date_min <- as.Date(min(d$tel[[.i]]$timestamp))
  date_max <- as.Date(max(d$tel[[.i]]$timestamp))
  
  png(paste0('figures/movement-model-diagnostics/', d$dataset_name[.i],
             '-', .a, '-movement-model-diagnostics.png'),
      width = 12, height = 6, units = 'in', bg = 'white', res = 600)
  layout(t(1:2))
  plot(variogram(d$tel[[.i]]), d$movement_model[[.i]], fraction = 0.5)
  title(paste0(d$dataset_name[.i], '; ', d$animal_id[.i]))
  plot(d$tel[[.i]], d$akde[[.i]],
       col = color(d$tel[[.i]],
                   col.fn = \(i, alpha) {
                     dt <- min(diff(i), na.rm = TRUE)
                     t <- i / dt
                     t <- as.integer(t / (1 %#% 'day'))
                     t <- t - min(t) + 1
                     return(khroma::color('smoothrainbow')(length(t)))
                   }), error = FALSE, cex = 1, pch = 19)
  title(paste0('Start: ', date_min, '; ', 'End: ', date_max, '; ',
               'Duration: ', difftime(date_max, date_min, units = 'days'),
               ' days'))
  dev.off()
}, .progress = TRUE)
