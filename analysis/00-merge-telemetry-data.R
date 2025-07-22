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

setwd('data/telemetry-data')
d <-
  bind_rows(
    # clean data from other projects
    # Mezzini et al. 2025; https://github.com/QuantitativeEcologyLab/bc-mammals-temperature
    readRDS('mezzini-bc-temperature/all-tracking-data-cleaned-2024-02-22-13-49.rds') %>%
      unnest(tel) %>%
      transmute(animal_id = animal,
                species, dataset_name, individual.local.identifier,
                long = location.long, lat = location.lat,
                timestamp = utc_timestamp, hdop, time_zone = 'UTC',
                outlier, clean = TRUE),
    # Medici et al. 2022; https://doi.org/10.1186/s40462-022-00313-w
    read_csv('medici-2022-tapirs/atlantica.csv', show_col_types = FALSE) %>%
      transmute(animal_id = individual.local.identifier,
                species = 'Tapirus terrestris',
                dataset_name = 'Medici et al. 2022 (AF)',
                individual.local.identifier,
                long = location.long, lat = location.lat,
                time_zone = tz_lookup_coords(lat = lat, lon = long,
                                             method = 'accurate'),
                timestamp = as.POSIXct(timestamp),
                outlier = Outlier, clean = TRUE),
    read_csv('medici-2022-tapirs/cerrado.csv', show_col_types = FALSE) %>%
      transmute(animal_id = individual.local.identifier,
                species = 'Tapirus terrestris',
                dataset_name = 'Medici et al. 2022 (CE)',
                individual.local.identifier,
                long = Longitude, lat = Latitude,
                time_zone = tz_lookup_coords(lat = lat, lon = long,
                                             method = 'accurate'),
                timestamp = as.POSIXct(timestamp),
                outlier = Outlier, clean = TRUE),
    read_csv('medici-2022-tapirs/pantanal.csv', show_col_types = FALSE) %>%
      transmute(animal_id = individual.local.identifier,
                species = 'Tapirus terrestris',
                dataset_name = 'Medici et al. 2022 (PA)',
                individual.local.identifier,
                long = Longitude, lat = Latitude,
                time_zone = tz_lookup_coords(lat = lat, lon = long,
                                             method = 'accurate'),
                timestamp = as.POSIXct(timestamp),
                outlier = Outlier, clean = TRUE),
    # data from other contributors that i still need to clean
    # Aikens at al. (collected by H. Sawyer); https://doi.org/10.1016/j.cub.2025.03.010
    st_read('aikens-pronghorn/GPSdata.shp', quiet = TRUE) %>%
      st_transform(crs = 'EPSG:4326') %>%
      # assuming the data uses the local timezone
      mutate(time_zone = tz_lookup(., crs = st_crs(.), method = 'accurate')) %>%
      bind_cols(st_coordinates(.) %>%
                  data.frame() %>%
                  rename(long = X, lat = Y)) %>%
      st_drop_geometry() %>%
      transmute(animal_id = as.character(animalID),
                species = 'Antilocapra americana',
                dataset_name = 'Aikens et al. 2025 (Collected by Hall Sawyer)',
                individual.local.identifier = animal_id, long, lat,
                time_zone,
                timestamp = paste0(date, '-', hour) %>%
                  ymd_h(tz = unique(time_zone)),
                outlier = FALSE, clean = FALSE) %>%
      arrange(animal_id, timestamp) # original data are not sorted
  ) %>%
  # add correct timezone (perform for each time zone separately)
  group_by(time_zone) %>%
  mutate(timestamp = timestamp %>%
           `tz<-`(unique(time_zone))) %>%
  ungroup() %>%
  # convert all times to UTC and make animal IDs unique
  mutate(timestamp = with_tz(timestamp, 'UTC'),
         animal_id = paste(animal_id, 'from', dataset_name)) %>%
  select(! time_zone) %>%
  # collect the datasets for each individual animal
  nest(data = ! c(animal_id, species, dataset_name, clean)) %>%
  mutate(tel = map(data, \(.d) {
    # convert to telemetries quietly (will be doing data cleaning later)
    suppressMessages(as.telemetry(.d, mark.rm = TRUE))
  }, .progress = TRUE),
  n_locations = map_int(tel, nrow),
  n_days = map_dbl(data, \(.d) {
    diff(range(.d$timestamp), units = 'days')
  }))
d

# ensure there are no duplicated animal names across studies (should be 0)
d %>%
  pull(animal_id) %>%
  duplicated() %>%
  sum()

# check number of individuals per species and locations per individual
d %>%
  group_by(species) %>%
  summarise(n_animals = n(),
            mean_locations = mean(map_int(data, nrow)))

#' we might need to drop individuals with too little data. for example,
#' tracking a mongolian gazzelle for 10 days may give an unrealistically
#' small home range despite the low resources because the data doesn't
#' show the movement between patches
ggplot(d, aes(n_days)) +
  geom_histogram(bins = 50, fill = 'grey', color = 'black') +
  geom_vline(xintercept = 30, lty = 'dashed') +
  scale_x_log10('Number of days tracked', breaks = c(10, 30, 100, 1000)) +
  ylab('Count')

# add collar calibrations when possible?

# save the final dataset
saveRDS(d, paste0('cleaned-telemetry-data-', Sys.Date(), '.rds'))
