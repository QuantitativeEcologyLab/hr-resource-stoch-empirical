library('readr')     # for importing data
library('sf')        # for shapefiles 
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('furrr')     # for parallelized functional programming
library('ctmm')      # for converting data to telemetry format
library('lutz')      # for getting time zones based on coordinates
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
source('analysis/figures/default-ggplot-theme.R')
source('functions/import_rda.R')
source('functions/outlier_plots.R') # to plot outlier diagnostic plots
source('functions/check_animal.R') # to run diagnostic plots
source('functions/plot_adj.R') # to plot 20 adjacent locations
source('functions/flag_outlier.R') # to mark outliers
source('functions/remove_outlier_flags.R') # to start over with an animal

data('buffalo')
data('coati')
data('wolf')

# set up plan to run processes in parallel
plan(multisession, workers = availableCores(logical = FALSE) - 2)

# clean pronghorn data and fit models ----
# Aikens at al. (collected by H. Sawyer)
# https://doi.org/10.1016/j.cub.2025.03.010
pronghorn <-
  st_read('data/telemetry-data/aikens-pronghorn/GPSdata.shp',
          quiet = TRUE) %>%
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
              ymd_h(tz = unique(time_zone))) %>%
  arrange(animal, timestamp) %>% # original data are not sorted
  nest(tel = ! c(animal, species, dataset_name)) %>%
  mutate(tel = map(tel, as.telemetry))

if(FALSE) {
  map(1, function(i) {
    png(paste0('figures/pronghorn-variograms/pronghorn-', i, '.png'),
        width = 12, height = 10, units = "in", res = 200)
    ctmm.guess(pronghorn$tel[[i]])
    dev.off()
  }) %>%
    invisible()
}

# drop animals that are not range-resident within the observation period
pronghorn <- pronghorn %>% filter(! animal %in% c(22, 23, 26, 29, 30))

if(FALSE) {
  # all pronghorn telemetry data look clean
  out <- check_animal(id = 1, .data = pronghorn) # good
  out <- check_animal(id = 2, .data = pronghorn) # good
  out <- check_animal(id = 3, .data = pronghorn) # has a bit of an excursion
  out <- check_animal(id = 4, .data = pronghorn) # has a bit of an excursion
  out <- check_animal(id = 5, .data = pronghorn) # explores along a fence
  out <- check_animal(id = 8, .data = pronghorn) # explores along a fence
  out <- check_animal(id = 9, .data = pronghorn) # good
  out <- check_animal(id = 10, .data = pronghorn) # explores along a fence
  out <- check_animal(id = 11, .data = pronghorn) # good
  out <- check_animal(id = 12, .data = pronghorn) # fence makes a hole
  out <- check_animal(id = 13, .data = pronghorn) # good
  out <- check_animal(id = 14, .data = pronghorn) # good
  out <- check_animal(id = 17, .data = pronghorn) # explores along a fence
  out <- check_animal(id = 18, .data = pronghorn) # good
  out <- check_animal(id = 27, .data = pronghorn) # good
  out <- check_animal(id = 28, .data = pronghorn) # good
  out <- check_animal(id = 31, .data = pronghorn) # good
  out <- check_animal(id = 33, .data = pronghorn) # good
  out <- check_animal(id = 34, .data = pronghorn) # good
  out <- check_animal(id = 35, .data = pronghorn) # explores along a fence
  out <- check_animal(id = 36, .data = pronghorn) # good
  out <- check_animal(id = 37, .data = pronghorn) # good
  out <- check_animal(id = 38, .data = pronghorn) # explores along a fence
  out <- check_animal(id = 39, .data = pronghorn) # explores along a fence
  out <- check_animal(id = 40, .data = pronghorn) # good
  out <- check_animal(id = 41, .data = pronghorn) # good
  out <- check_animal(id = 42, .data = pronghorn) # good
  out <- check_animal(id = 43, .data = pronghorn) # good
  out <- check_animal(id = 44, .data = pronghorn) # good
}

nrow(pronghorn)

movebank_data <- readRDS('data/telemetry-data/movebank-data-2026-08-17.rds')

# compile full dataset ----
d <-
  bind_rows(
    movebank_data %>%
      rename(animal = animal_id) %>%
      mutate(movement_model = list(NULL),
             akde = list(NULL)) %>%
      select(animal, species, dataset_name, tel, movement_model, akde),
    # clean data from other projects
    # Mezzini et al. 2026; https://github.com/QuantitativeEcologyLab/bc-mammals-temperature
    readRDS('data/hr-estimates/mezzini-bc-tempterature-2026/movement-models-akdes-2024-06-06.rds') %>%
      # drop cougars to simplify data permission process
      filter(species != 'Puma concolor') %>%
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
    # DeNicola white-tailed deer (NY)
    readRDS('data/hr-estimates/denicola-deer/denicola-ny-2025-full-telemetry-movement-models-2024-04-20.rds') %>%
      filter(! grepl('Mortality', mortality.status)) %>% # post mortality
      rename(movement_model = model) %>%
      mutate(dataset_name = case_when(
        study.area == 'staten_island' ~ 'DeNicola et al. (2025) (Staten Island)',
        study.area == 'rockefeller' ~ 'DeNicola et al. (2025) (Rockefeller)'),
        species = 'Odocoileus virginianus',
        animal = paste0(animal, ' (year ', study.year, ')')) %>%
      select(animal, species, dataset_name, tel, movement_model, akde),
    # DeNicola white-tailed deer (South Euclid)
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
                      })) %>%
      filter(animal != 'Juliana'), # not range-resident 
    tibble(
      map('data/hr-estimates/noonan-anteater-2/Myrmecophaga_tridactyla.Rda',
          function(.fn) { # to convert Rda to an unnamed object
            load(.fn)
            data %>% select(! c(X.11:X)) %>% return()
          })[[1]]) %>%
      mutate(animal = ID) %>%
      nest(tel = ! animal) %>%
      filter(animal != 'Little Rick') %>% # duplicate from Issy's dataset
      filter(animal != 'Yoki') %>% # not a useful telemetry
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
    readRDS('data/hr-estimates/mills-northern-mountain-caribou/mills-nmc.rds') %>%
      mutate(species = 'Rangifer tarandus',
             dataset_name = 'Mills Nothern Mountain Caribou',
             movment_model = unname(movement_model),
             akde = unname(akde)) %>%
      select(animal, species, dataset_name, tel, movement_model, akde),
    # buffalo data from the {ctmm} package
    # AKDES and movement models provided by Mike via a Dropbox link
    #' from https://github.com/ctmm-initiative/ctmmlearn/blob/main/data
    tibble(species = 'Syncerus caffer',
           dataset_name = 'ctmm buffalo',
           movement_model =
             import_rda('data/telemetry-data/ctmm-data/buffalo.rda',
                        object_name = 'FITS'),
           animal = map_chr(movement_model, function(.mm) {
             return(.mm@info$identity)
           }),
           tel = buffalo[animal],
           akde = import_rda('data/telemetry-data/ctmm-data/buffalo_akdes.rda',
                             object_name = 'AKDES')[animal]),
    # coati data from the {ctmm} package
    tibble(species = 'Nasua narica',
           dataset_name = 'ctmm coati',
           movement_model = list(NULL),
           animal = names(coati),
           tel = coati[animal],
           akde = list(NULL)),
    # maned wolf data from the {ctmm} package
    tibble(species = 'Chrysocyon brachyurus',
           dataset_name = 'ctmm maned wolf',
           movement_model = list(NULL),
           animal = names(wolf),
           tel = wolf[animal],
           akde = list(NULL)),
    # banded armadillo from Mike Noonan
    tibble(tel = import_rda('data/telemetry-data/noonan-banded-armadillo/Tolypeutes_matacus.Rda',
                            object_name = 'Tolypeutes_matacus') %>%
             as.telemetry() %>%
             suppressMessages(),
           animal = names(tel),
           species = 'Tolypeutes matacus',
           dataset_name = 'Noonan southern tree banded armadillo') %>%
      mutate(tel = unname(tel)),
    # pronghorn data from Aikens at al. (collected by H. Sawyer)
    # https://doi.org/10.1016/j.cub.2025.03.010
    pronghorn %>%
      mutate(animal = paste0("Pronghorn", animal)),
    # alpine ibex from Movebank (movement models from Mike)
    tibble(tel = import_rda("data/telemetry-data/capra-ibex/Capra_ibex.Rda", "Capra_ibex") %>%
             as.telemetry() %>%
             suppressMessages(),
           animal = names(tel),
           species = "Capra ibex",
           dataset_name = "Monitoring of Capra ibex (Bovidae) populations in the western alps (project ALCOTRA LEMED-IBEX)",
           movement_model = map(animal, function(.a) {
             import_rda(paste0("data/hr-estimates/noonan-capra-ibex/",
                               "Fits_", .a, ".rda"), object_name = "FIT")
           }),
           akde = list(NULL))
  ) %>%
  # remove duplicated individuals from Mike's datasets
  filter(! (grepl('Noonan giant anteaters', dataset_name) &
              (duplicated(animal) | duplicated(animal, fromLast = TRUE)))) %>%
  mutate(species = factor(species),
         animal = factor(paste0(animal, ' (', species, ')')),
         n_locations = map_int(tel, \(.tel) {
           if (is.null(.tel)) return(NA_integer_) else return(nrow(.tel))
         }),
         n_days = map_dbl(tel, \(.tel) {
           if (is.null(.tel)) {
             return(NA_real_)
           } else {
             return(diff(range(.tel$timestamp), units = 'days'))
           }
         }) %>%
           round())

# telemetries are not available for ctmm maned wolf
filter(d, n_days < 0)

d$i <- 1:nrow(d)

# plot variograms ----
d <-
  d %>%
  mutate(variogram = map(i, function(.i) {
    .id <- substr(d$animal[.i], 1, 30)
    .tel <- tel[[.i]]
    if (is.null(.tel)) return(NULL)
    
    png(paste0("figures/variograms/", .i, "-", .id, ".png"),
        width = 6, height = 6, units = "in", res = 300)
    plot(variogram(.tel))
    dev.off()
    
    ctmm.guess(.tel, interactive = FALSE)
  }, .progress = TRUE))

# add missing movement models ----
# not using a single pipeline to have a more realistic ETA
d <-
  bind_rows(
    d %>% filter(!map_lgl(movement_model, is.null)),
    # fit missing movement models or UDs
    d %>% filter(map_lgl(movement_model, is.null)) %>%
      mutate(movement_model = future_map(1:n(), \(.i) {
        .mm <- movement_model[[.i]]
        if (class(.mm) == 'ctmm') {
          return(.mm)
        }
        .a <- as.character(animal[.i])
        model_filename <- paste0("models/movement-models/", .a, ".rds")
        if (file.exists(model_filename)) {
          return(readRDS(model_filename))
        }
        .tel <- tel[[.i]]
        if(is.null(.tel)) {
          return(NULL)
        }
        .vg <- variogram[[.i]]
        .mm <- ctmm.select(data = .tel, CTMM = .vg)
        saveRDS(.mm, model_filename)
        .mm
      }, .progress = TRUE, .options = furrr_options(seed = NULL))))

# ensure all movement models used the correct data
with(filter(d, !map_lgl(tel, is.null), !map_lgl(movement_model, is.null)),
     all(map_chr(tel, \(.tel) .tel@info$identity) ==
           map_chr(movement_model, \(.mm) .mm@info$identity)))

table(map_chr(d$movement_model, class))

# add missing AKDEs ----
# flag animals that are not range resident (can still use speed & diffusion)
# assessed via variogram figures; if migrating at least one full migration
d$range_resident <- ! d$i %in% c(
  18, 31, 44, 45, 48, 49, 50, 57, 59, 66, 90, 100, 105, 115, 116, 117, 133, 134,
  152, 153, 155, 156, 157, 160, 164, 166, 168, 169, 173, 176, 177, 181, 183,
  188, 190, 191, 194, 195, 196, 199, 201, 205, 208, 210, 211, 212, 214, 220,
  221, 222, 227, 228, 229, 231, 234, 237, 238, 239, 242, 243, 251, 253, 254,
  255, 258, 259, 261, 263, 268, 273, 274, 276, 277, 279, 281, 293, 295, 305,
  306, 317, 319, 328, 329, 341, 353, 358, 361, 362, 363, 366, 367, 368, 369,
  370, 371, 373, 374, 375, 378, 379, 380, 381, 383, 384, 405, 410, 419, 421,
  422, 430, 432, 444, 445, 447, 448, 449, 459, 500, 508, 511, 512, 514, 516,
  523, 526, 528, 530, 532, 533, 537, 553, 564, 567, 569, 570, 573, 574, 575,
  576, 577, 578, 580, 581, 582, 583, 584, 586, 587, 588, 589, 593, 596, 597,
  598, 599, 604, 606, 611, 614, 615, 620, 621, 622, 623, 630, 642, 643, 655,
  658, 663, 664, 666, 667, 686, 703, 718, 721, 731, 742, 756, 789, 790, 797,
  798, 804, 826, 841, 850, 851, 852, 868, 927, 932, 962, 990, 991, 992, 993,
  994, 996, 1013, 1021, 1041, 1139, 1155, 1181, 1192, 1217, 1224, 1257, 1263,
  1275, 1330, 1331, 1332, 1337, 1338, 1343, 1346, 1347, 1348, 1349, 1350, 1351,
  1352, 1353, 1354, 1355, 1356, 1357, 1358, 1360, 1361, 1363, 1367, 1399, 1408,
  1410, 1413, 1417, 1418, 1421, 1425, 1426, 1427, 1431, 1436, 1449, 1451, 1457,
  1461, 1464, 1469, 1476, 1479, 1489, 1497, 1501, 1503, 1507, 1509, 1510, 1512,
  1515, 1520, 1522, 1526, 1530, 1538, 1560, 1564, 1567, 1571, 1577, 1586, 1587,
  1596, 1597)

if (!dir.exists("models/akdes")) dir.create("models/akdes")

mean(map_lgl(d$akde, is.null))

d <- d %>%
  mutate(
    # import existing AKDEs
    akde = future_map2(animal, akde, \(.a, .ud) {
      if (!is.null(.ud)) {
        return(.ud)
      }
      
      akde_filename <- paste0("models/akdes/", .a, ".rds")
      if (file.exists(akde_filename)) {
        return(readRDS(akde_filename))
      } else {
        return(NULL)
      }
    }, .progress = TRUE))

mean(map_lgl(d$akde, is.null))

d <- d %>%
  bind_rows(
    # separate individuals that already have or can't have AKDEs 
    d %>% filter(!map_lgl(akde, is.null) | !range_resident | is.null(tel) |
                   is.null(movement_model)),
    # from akdes that can still be fit
    #' *NOTE:* can't use i column because filtering breaks the indexing
    d %>%
      filter(map_lgl(akde, is.null), range_resident, !is.null(tel),
             !is.null(movement_model)) %>%
      mutate(
        # fit remaining AKDEs
        akde = future_map(1:n(), \(.i) {
          .tel <- tel[[.i]]
          .mm <- movement_model[[.i]]
          .a <- animal[[.i]]
          .akde <- akde[[.i]]
          
          if (!is.null(.akde)) {
            return(.akde)
          }
          
          if (is.null(.tel) | is.null(.mm)) {
            return(NULL)
          }
          akde_filename <- paste0("models/akdes/", .a, ".rds")
          temp_filename <- paste0(.tel@info$identity, '-temp.png')
          png(temp_filename)
          .akde <- akde(data = .tel, CTMM = .mm, debias = TRUE,
                        weights = TRUE)
          dev.off()
          unlink(temp_filename)
          saveRDS(.akde, akde_filename)
          .akde
        }, .progress = TRUE)))

# estimate movement parameters ----
d <- d %>%
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
    if (is.null(.akde)) {
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
        if(any(grepl('diffusion', rownames(summary(.mm)$CI)))) {
          return(summary(.mm, units = FALSE)$
                   CI['diffusion (square meters/second)', 'est'])
        } else {
          return(NA_real_)
        }
      }
    }))

# warnings about fractal movement are ok: no speed was estimated
unique(dplyr::last_dplyr_warnings())

# save the final dataset
saveRDS(d, paste0('data/akde-tibble-', Sys.Date(), '.rds'))

# ensure there are no duplicated animal names across studies (should be 0)
if (sum(duplicated(d$animal)) > 0) {
  d %>%
    select(animal, species, dataset_name, n_days) %>%
    mutate(animal = gsub('[-_]', ' ', animal)) %>%
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
#' tracking a nomadic animal for 10 days may give an unrealistically
#' small home range despite the low resources because the data doesn't
#' show the movement between patches
ggplot(d, aes(n_days)) +
  geom_histogram(bins = 50, fill = 'grey', color = 'black') +
  geom_vline(xintercept = 30, lty = 'dashed') +
  scale_x_log10('Number of days tracked', breaks = c(10, 30, 100, 1000)) +
  ylab('Count')

# currently dropping individuals with HR ESS < 4
ggplot(d, aes(akde_ess_area)) +
  geom_histogram(bins = 50, fill = 'grey', color = 'black') +
  geom_vline(xintercept = 4, lty = 'dashed') +
  scale_x_log10('Area ESS (log scale)',
                breaks = c(1, 4, 10, 100, 1000)) +
  ylab('Count')

# TODO: drop animals with excessively large HRs (might be fixed now)
layout(1:2)
hist(d$hr_95_km2, breaks = 100)
hist(filter(d, keep_for_gam)$hr_95_km2, breaks = 100, xlim = range(d$hr_95_km2))
layout(1)
filter(d, keep_for_gam, hr_95_km2 > 8000) %>% pull(animal)

large_hrs <- filter(d, keep_for_gam, hr_95_km2 > 8000)
DIM <- ceiling(sqrt(nrow(large_hrs)))

layout(matrix(1:(DIM)^2, ncol = DIM))
for(i in 1:nrow(large_hrs)) {
  plot(large_hrs$akde[[i]])
  plot(large_hrs$tel[[i]], add= TRUE, error = FALSE, type = 'l')
  title(paste0(large_hrs$animal[i], '\n95% HR:\t\t\t', round(d$hr_95_km2[i]),
               ' km\U00B2'))
}
layout(1)

map_dbl(large_hrs$akde, function(.akde) {
  round(summary(.akde)$DOF['area'], 2)
})
