library('readr')     # for importing data
library('sf')        # for shapefiles 
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('ctmm')      # for converting data to telemetry format
library('lutz')      # for getting time zones based on coordinates
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
source('analysis/figures/default-ggplot-theme.R')

d <-
  bind_rows(
    # clean data from other projects
    # Mezzini et al. 2026; https://github.com/QuantitativeEcologyLab/bc-mammals-temperature
    readRDS('data/hr-estimates/mezzini-bc-tempterature-2026/movement-models-akdes-2024-06-06.rds') %>%
      # drop non-range-resident individuals
      left_join(read_csv('data/hr-estimates/mezzini-bc-tempterature-2026/telemetry-metadata.csv',
                         show_col_types = FALSE) %>%
                  select(animal, keep_for_RSF),
                by = 'animal') %>%
      filter(keep_for_RSF == 'yes') %>%
      filter(animal != 'E035') %>% # has a HR of almost 3,000 km^2
      mutate(species = if_else(species == 'Ursus arctos horribilis',
                               'Ursus arctos', species)) %>%
      select(animal, species, dataset_name, tel, movement_model, akde),
    # Medici et al. 2022; https://doi.org/10.1186/s40462-022-00313-w
    readRDS('data/hr-estimates/medici-tapirs/medici-tapirs-final.rds') %>%
      rename(animal = name, tel = data, movement_model = model) %>%
      mutate(dataset_name = paste('Medici et al. (2022)', region),
             species = 'Tapirus terrestris') %>%
      select(animal, species, dataset_name, tel, movement_model, akde),
    # DeNicola white-tailed deer
    readRDS('data/hr-estimates/denicola-deer/denicola-ny-2025-full-telemetry-movement-models-2024-04-20.rds') %>%
      # drop data after mortality event
      filter(! grepl('Mortality', mortality.status)) %>%
      rename(movement_model = model) %>%
      mutate(dataset_name = case_when(
        study.area == 'staten_island' ~ 'DeNicola et al. (2025) (Staten Island)',
        study.area == 'rockefeller' ~ 'DeNicola et al. (2025) (Rockefeller)'),
        species = 'Odocoileus virginianus',
        animal = paste0(animal, ' (year ', study.year, ')')) %>%
      select(animal, species, dataset_name, tel, movement_model, akde),
    readRDS('data/hr-estimates/denicola-deer/denicola-se-2025-full-telemetry-movement-models.rds') %>%
      rename(movement_model = mm,
             akde = ud) %>%
      mutate(dataset_name = 'DeNicola et al. (2025) (South Euclid)',
             species = 'Odocoileus virginianus') %>%
      select(animal, species, dataset_name, tel, movement_model, akde),
    # Noonan
    tibble(species = 'Myrmecophaga tridactyla',
           dataset_name = 'Noonan giant anteaters',
           tel = list(NULL),
           movement_model = map(list.files('data/hr-estimates/noonan-anteater_movement/movement_models',
                                           full.names = TRUE),
                                function(.file) {
                                  load(.file)
                                  return(FIT)
                                }),
           animal = map_chr(movement_model, \(.mm) .mm@info$identity),
           akde = map(list.files('data/hr-estimates/noonan-anteater_movement/akdes/',
                                 full.names = TRUE),
                      function(.file) {
                        load(.file)
                        return(AKDE)
                      })),
    tibble(
      map('data/hr-estimates/noonan-anteater-2/Myrmecophaga_tridactyla.Rda',
          function(.fn) { # to convert Rda to an unnamed object
            load(.fn)
            data %>%
              select(! c(X.11:X)) %>%
              return()
          })[[1]]) %>%
      mutate(animal = ID) %>%
      nest(tel = ! animal) %>%
      filter(animal != 'Yoki') %>%
      mutate(species = 'Myrmecophaga tridactyla',
             dataset_name = 'Noonan giant anteaters 2',
             tel = map(tel, \(.tel) suppressMessages(as.telemetry(.tel))),
             movement_model = map(animal, function(.name) {
               load(paste0('data/hr-estimates/noonan-anteater-2/Fits_',
                           .name, '.Rda'))
               return(FIT)
             }),
             akde = list(NULL)) %>%
      select(animal, species, dataset_name, tel, movement_model, akde),
    tibble(species = 'Panthera onca',
           dataset_name = 'Noonan Maraca jaguars',
           tel = list(NULL),
           movement_model = map(list.files('data/hr-estimates/noonan-maraca-jaguars/movement_models',
                                           full.names = TRUE),
                                function(.file) {
                                  load(.file)
                                  return(FIT)
                                }),
           animal = map_chr(movement_model, \(.mm) .mm@info$identity),
           akde = map(animal, function(.animal) {
             load(paste0('data/hr-estimates/noonan-maraca-jaguars/akdes/',
                         'AKDE_', .animal, '.rda'))
             return(AKDE)
           })),
    tibble(species = 'Priodontes maximus',
           dataset_name = 'Noonan giant armadillo',
           tel = list(NULL),
           movement_model = map(list.files('data/hr-estimates/noonan-armadillo_movement/movement_models',
                                           full.names = TRUE),
                                function(.file) {
                                  load(.file)
                                  return(FIT)
                                }),
           animal = map_chr(movement_model, \(.mm) .mm@info$identity),
           akde = map(list.files('data/hr-estimates/noonan-armadillo_movement/akdes',
                                 full.names = TRUE),
                      function(.file) {
                        load(.file)
                        return(AKDE)
                      })),
    tibble(species = 'Tapirus terrestris',
           dataset_name = 'Noonan tapirs',
           tel = list(NULL),
           movement_model = map(list.files('data/hr-estimates/noonan-tapir_movement/movement_models',
                                           full.names = TRUE),
                                function(.file) {
                                  load(.file)
                                  return(FIT)
                                }),
           animal = map_chr(movement_model, \(.mm) .mm@info$identity),
           akde = map(list.files('data/hr-estimates/noonan-tapir_movement/akdes',
                                 full.names = TRUE),
                      function(.file) {
                        load(.file)
                        return(AKDE)
                      })),
    # giant anteaters (Issy Perry)
    readRDS('data/hr-estimates/perry-Anteater_data/Anteater_data/Anteater_tibble.rds') %>%
      # only keep range-resident anteaters
      filter(animal %in% c('Arya', 'Bella', 'Colete', 'Dom', 'Dumbo_1',
                           'Dumbo_2', 'George', 'Heather', 'Juju_2',
                           'Mulan', 'Nancy', 'Peter', 'Rita', 'Tim_2',
                           'Tim_3', 'Alexander', 'Annie', 'Beto',
                           'Bumpus', 'Cate', 'Christoffer', 'Elaine',
                           'Hannah', 'Jackson', 'Jane', 'Larry', 'Luigi',
                           'Makao', 'Margaret', 'Maria', 'Puji', 'Rodolfo',
                           'Sheron', 'Thomas', 'Anthony', 'Kyle',
                           'Little_Rick', 'Reid')) %>%
      # only keep non-orphaned anteaters
      filter(! animal %in% names(readRDS('data/hr-estimates/perry-Anteater_data/Anteater_data/FITS_UDs/FITS_Orphan.rds'))) %>%
      mutate(across(c(tel, movement_model, akde), unname),
             species = 'Myrmecophaga tridactyla'),
    # Northern Mountain Caribou (Lauren Mills)
    # TODO: missing tels and movement models
    readRDS('data/hr-estimates/mills-northern-mountain-caribou/akde_df.RDS') %>%
      as_tibble() %>%
      rename(animal = animal_id, akde = AKDE_object) %>%
      mutate(species = 'Rangifer tarandus',
             dataset_name = 'Mills Nothern Mountain Caribou',
             tel = list(NULL),
             movement_model = list(NULL),
             akde = unname(akde)) %>%
      select(animal, species, dataset_name, tel, movement_model, akde)
    # TODO: add pronghorn data
  ) %>%
  # remove duplicated individuals from Mike's datasets
  filter(! ((duplicated(animal) | duplicated(animal, fromLast = TRUE)) &
              grepl('Noonan giant anteaters', dataset_name))) %>%
  mutate(species = factor(species),
         animal = factor(paste0(animal, ' (', species, ')')),
         n_locations = map_int(tel, \(.tel) {
           if(is.null(.tel)) return(NA_integer_) else return(nrow(.tel))
         }),
         n_days = map_dbl(tel, \(.tel) {
           if(is.null(.tel)) {
             return(NA_real_)
           } else {
             return(diff(range(.tel$timestamp), units = 'days'))
           }
         }) %>%
           round()) %>%
  # # add missing movement models and AKDEs
  # mutate(movement_model = map2(movement_model, function(.tel, .mm) {
  #   if (class(.mm) == 'ctmm') {
  #     return(.mm)
  #   } else {
  #     .vg <- ctmm.guess(.tel, interactive = FALSE)
  #     .mm <- ctmm.select(data = .tel, CTMM = .vg)
  #     return(.mm)
  #   }
  # }, .progress = 'Fitting missing movement models...'),
  # akde = map2(.tel, .mm, function(.tel, .mm) {
  #   return(akde(data = .tel, CTMM = .mm, debias = TRUE, weights = TRUE))
  # }, .progress = 'Fitting missing AKDEs...')) %>%
  # estimate movement parameters
  mutate(hr_95_km2 = map_dbl(akde, function(.akde) {
    if (is.null(.akde)) {
      return(NA_real_)
    } else {
      # convert HR to km^2
      ci <- summary(.akde, units = FALSE, level.UD = 0.95)$CI
      return('km^2' %#% ci['area (square meters)', 'est'])
    }
  }, .progress = 'Estimating HRs...'),
  akde_ess_area = map_dbl(akde, function(.akde) {
    if(is.null(.akde)) {
      return(NA_real_)
    } else {
      return(summary(.akde)$DOF['area'])
    }
  }),
  keep_for_gam = akde_ess_area >= 4,
  #' add normalized weights based on `sqrt(ESS)`
  weight = if_else(keep_for_gam, sqrt(akde_ess_area), NA_real_),
  weight = weight / mean(weight, na.rm = TRUE),
  speed_km_day = map_dbl(movement_model, function(.mm) {
    if (is.null(.mm)) {
      return(NA_real_)
    } else {
      return('km/day' %#% speed(.mm, units = FALSE)$CI[, 'est'])
    }
  }, .progress = 'Estimating speeds...'),
  diffusion_km2_day =  'km^2/day' %#%
    map_dbl(movement_model, function(.mm) {
      if (is.null(.mm)) {
        return(NA_real_)
      } else {
        if(any(grepl('diffusion', rownames(summary(.m)$CI)))) {
          return(summary(.m, units = FALSE)$
                   CI['diffusion (square meters/second)', 'est'])
        } else {
          return(NA_real_)
        }
      }
    }))

# ensure all species have data
if(sum(is.na(d$mass_kg)) > 0) stop('At least one row has `NA` mass!')

# ensure there are no duplicated animal names across studies (should be 0)
if(sum(duplicated(d$animal)) > 0) {
  d %>%
    select(animal, species, dataset_name, n_days) %>%
    filter(duplicated(animal) | duplicated(animal, fromLast = TRUE)) %>%
    arrange(animal, dataset_name)
  
  stop('There are duplicate animal names!')
}

# check number of individuals per species and fixes per individual
d %>%
  group_by(species) %>%
  summarise(n_animals = n(),
            thousand_locations = sum(n_locations) / 1e3) %>%
  janitor::adorn_totals(where = 'row', na.rm = TRUE, fill = 'non-numeric')

#' TODO:
#' we might need to drop individuals with too little data. for example,
#' tracking a mongolian gazzelle for 10 days may give an unrealistically
#' small home range despite the low resources because the data doesn't
#' show the movement between patches
ggplot(d, aes(n_days)) +
  geom_histogram(bins = 50, fill = 'grey', color = 'black') +
  geom_vline(xintercept = 30, lty = 'dashed') +
  scale_x_log10('Number of days tracked', breaks = c(10, 30, 100, 1000)) +
  ylab('Count')

# currently dropping individuals with HR ESS < 10
ggplot(d, aes(akde_ess_area)) +
  geom_histogram(bins = 50, fill = 'grey', color = 'black') +
  geom_vline(xintercept = 4, lty = 'dashed') +
  scale_x_log10('Area ESS (log scale)',
                breaks = c(1, 4, 10, 100, 1000)) +
  ylab('Count')

# TODO: drop animals with excessively large HRs (might be fixed now)
# filter based on UD ESS?
hist(d$hr_95, breaks = 100)
hist(filter(d, keep_for_gam)$hr_95, breaks = 100)
filter(d, keep_for_gam, hr_95 > 8000) %>% pull(animal)

large_hrs <- filter(d, keep_for_gam, hr_95 > 8000)
DIM <- ceiling(sqrt(nrow(large_hrs)))

layout(matrix(1:(DIM)^2, ncol = DIM))
for(i in 1:nrow(large_hrs)) {
  plot(large_hrs$akde[[i]])
  plot(large_hrs$tel[[i]], add= TRUE, error = FALSE, type = 'l')
  title(paste0(large_hrs$animal[i], '\n95% HR:\t\t\t', round(d$hr_95[i]),
               ' km\U00B2'))
}
layout(1)

map_dbl(large_hrs$akde, function(.akde) {
  round(summary(.akde)$DOF['area'], 2)
})

# save the final dataset
saveRDS(d, paste0('data/akde-tibble-', Sys.Date(), '.rds'))

if(FALSE) {
  # Aikens at al. (collected by H. Sawyer); https://doi.org/10.1016/j.cub.2025.03.010
  # not using because there are no movement models fit and some are not range resident
  pronghorn <-
    st_read('data/telemetry-data/aikens-pronghorn/GPSdata.shp', quiet = TRUE) %>%
    st_transform(crs = 'EPSG:4326') %>%
    # assuming the data uses the local timezone
    mutate(time_zone = tz_lookup(., crs = st_crs(.), method = 'accurate')) %>%
    bind_cols(st_coordinates(.) %>%
                data.frame() %>%
                rename(long = X, lat = Y)) %>%
    st_drop_geometry() %>%
    transmute(animal = as.character(animalID),
              species = 'Antilocapra americana',
              dataset_name = 'Aikens et al. 2025 (Collected by Hall Sawyer)',
              individual.local.identifier = animal,
              long, lat, time_zone,
              timestamp = paste0(date, '-', hour) %>%
                ymd_h(tz = unique(time_zone)),
              outlier = FALSE, clean = FALSE) %>%
    arrange(animal, timestamp) %>% # original data are not sorted
    nest(tel = ! c(animal, species, dataset_name, clean)) %>%
    mutate(tel = map(tel, as.telemetry),
           vg = map(tel, \(.tel) ctmm.guess(.tel, interactive = FALSE)),
           movement_model = map(tel, vg, \(.tel, .vg) {
             list(ctmm.select(data = .tel, CTMM = .vg))
           }, .progress = TRUE))
  nrow(pronghorn)
  
  for(i in 1:nrow(pronghorn)) {
    png(paste0('figures/pronghorn-variograms/pronghorn-', i, '.png'))
    ctmm.guess(pronghorn$tel[[i]])
    dev.off()
  }
}