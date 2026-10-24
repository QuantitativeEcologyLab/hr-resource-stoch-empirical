library('readr')     # for importing data
library('sf')        # for shapefiles 
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('ctmm')      # for converting data to telemetry format
library('lutz')      # for getting time zones based on coordinates
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
source('analysis/default-ggplot-theme.R')
source('functions/check_animal.R')
source('functions/flag_outlier.R')
source('functions/remove_outlier_flags.R')
source('functions/plot_start.R')
source('functions/plot_adj.R')
options(plot_outlier_diagnostics = FALSE) # plotting diagnostics is much slower
options(viewer = NULL) # ensures mapview::mapvioew() works

# functions for this script ----
import_movebank <- function(file_name) {
  d_name <- gsub('data/telemetry-data/public-movebank-data/', '', file_name) %>%
    gsub('\\.csv', '', .)
  
  .d <- readr::read_csv(file_name, name_repair = 'universal_quiet',
                        show_col_types = FALSE)
  
  if ('outlier' %in% colnames(.d)) {
    .d <- filter(.d, ! outlier)
  } else {
    .d$outlier <- FALSE
  }
  
  .d %>%
    transmute(
      dataset_name = d_name,
      animal_id = paste(`individual.local.identifier`, 'from', d_name),
      species = `individual.taxon.canonical.name`,
      dataset_name = `study.name`,
      individual.local.identifier = as.character(individual.local.identifier),
      long = `location.long`, lat = `location.lat`,
      timestamp = as.POSIXct(timestamp),
      outlier = FALSE
    ) %>%
    arrange(animal_id, timestamp) %>%
    #' drop duplicated locations -- not worth modelling since:
    #' 1. we don't have calibration data for the collars;
    #' 2. they will have little effect on HR size, diffusion, and mean speed;
    #' 3. they belong to few individuals.
    slice(1, .by = c(timestamp, animal_id))
}

plot_outlie <- function(.tel, .id) {
  .id <- substr(.id, 1, 30) # some names are too long
  png(paste0(.id, '.png'), width = 8, height = 8, units = 'in', res = 300)
  withr::defer(dev.off())
  layout(matrix(c(1, 1, 2, 3), ncol = 2, byrow = TRUE))
  plot(.tel, error = FALSE, type = 'l', col = 'black')
  suppressMessages(plot(outlie(d$tel[[1]])))
  .id
}

# import, clean, and tidy the movebank data ----
d0 <-
  list.files('data/telemetry-data/public-movebank-data',
             pattern = '\\.csv',
             full.names = TRUE) %>%
  stringr::str_subset('^(?!.*-reference-data)') %>% # drop reference data
  map(import_movebank) %>%
  list_rbind() %>%
  # collect the datasets for each individual animal
  nest(data = ! c(animal_id, species, dataset_name, clean)) %>%
  mutate(
    n_locations = map_int(data, nrow),
    n_days = map_dbl(data, \(.d) {
      diff(range(.d$timestamp), units = 'days')
    }))

# add missing species names
d <- d0 %>%
  mutate(
    species = if_else(is.na(species),
                      # to catch issues with datasets with more than 1 species
                      paste(unique(species[! is.na(species)]), collapse = '&'),
                      species),
    .by = dataset_name)
stopifnot(nrow(filter(d, grepl('&', species))) == 0)

# fix bad locations (all from the Sloth Conservation Foundation?)
if (FALSE) { # to only run in interactive runs
  cr_sloths <-
    d %>%
    mutate(bad_locs = map_int(data, \(.d) sum(abs(.d$long) > 180 | abs(.d$lat) > 90))) %>%
    select(dataset_name, species, data, bad_locs) %>%
    filter(bad_locs > 0) %>%
    unnest(data)
  range(cr_sloths$long)
  range(cr_sloths$lat)
  View(cr_sloths)
  
  hist(cr_sloths$lat)
  cr_sloths$lat <- cr_sloths$lat / 10^floor(log10(cr_sloths$lat))
  hist(cr_sloths$lat)
  
  hist(cr_sloths$long)
  cr_sloths$long <- cr_sloths$long / 10^(floor(log10(abs(cr_sloths$long))) - 1)
  cr_sloths$long <- - abs(cr_sloths$long) # some longitudes were positive...
  hist(cr_sloths$long)
  
  ggplot(cr_sloths, aes(long, lat)) +
    geom_point(shape = 1, alpha = 0.1)
}

message('Fixing telemetries from The Sloth Conservation Foundation:\n',
        'Some are missing decimal separators, and some longitudes are positive.')
d <- d %>%
  mutate(
    data = map2(data, dataset_name, function(.d, .n) {
      #' `if_else()` runs `log10(lat)` for all rows: returns `NaN` if `lat < 0`
      if (grepl('The Sloth Conservation Foundation', .n)) {
        .d <-
          mutate(.d,
                 lat = lat / 10^floor(log10(lat)),
                 long = -abs(long / 10^(floor(log10(abs(long))) - 1)))
      }
      .d
    }))

d %>%
  mutate(bad_locs = map_int(data, \(.d) {
    sum(abs(.d$long) > 180 | abs(.d$lat) > 90)
  })) %>%
  filter(bad_locs > 0) %>%
  nrow() %>%
  `==`(0) %>%
  stopifnot()

stopifnot(nrow(filter(d, duplicated(animal_id))) == 0)
stopifnot(nrow(filter(d, is.na(animal_id))) == 0)
stopifnot(nrow(filter(d, is.na(species))) == 0)

message('Dropping 3 Argos arctic fox telemetries: ',
        'Keeping the GPS telemetries.')
dup_foxes <- filter(d, species == 'Vulpes lagopus') %>%
  pull(animal_id) %>%
  gsub(' .*', '', .) %>%
  table() %>%
  `[`(., . > 1) %>%
  names()

d <-
  d %>%
  filter(! (species == 'Vulpes lagopus' &
              grepl(paste0(dup_foxes, collapse = '|'), animal_id) &
              grepl('Argos', animal_id)))

# drop telemetries not worth keeping
message('Dropping telemetries with less than 7 days.')
quantile(d$n_locations, c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 1))
quantile(d$n_days / 7, c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 1))
range(filter(d, n_locations > 5)$n_locations)

d <- filter(d, n_days > 7)

# check number of individuals per species and locations per individual
d %>%
  group_by(species) %>%
  summarise(n_animals = n(),
            mean_locations = mean(map_int(data, nrow))) %>%
  ggplot() +
  geom_point(aes(n_animals, species, color = log2(mean_locations))) +
  theme(legend.position = 'top')

#' we might need to drop individuals with too little data. for example,
#' tracking a mongolian gazzelle for 10 days may give an unrealistically
#' small home range despite the low resources because the data doesn't
#' show the movement between patches
ggplot(d, aes(n_days)) +
  geom_histogram(bins = 50, fill = 'grey', color = 'black') +
  geom_vline(xintercept = 30, lty = 'dashed') +
  scale_x_log10('Number of days tracked', breaks = c(10, 30, 100, 1000)) +
  ylab('Count')

# create telemetries
d$tel <- map(d$data, \(.d) {
  # convert to telemetries quietly (will be doing data cleaning later)
  suppressMessages(as.telemetry(.d, mark.rm = TRUE))
})

# check telemetry diagnostics
setwd('figures/movebank-data-telemetry-diagnostics/')
saved <- map2(d$tel, d$animal_id, plot_outlie, .progress = TRUE)
setwd('../..')

# clean the telemetries ----
if(! exists('d')) d <- readRDS('data/telemetry-data/movebank-data-2026-08-01.rds')

#' `cat(paste0('out <- check_animal(d$animal_id[', 1:nrow(d)'])'), sep = '\n)`
# skipping telemetries that look clean
# ignoring ouliers that seem reasonable: they won't affect the averages much
# 25: example of reasonable outlier
out <- check_animal(d$animal_id[25])

# 26
out <- check_animal(d$animal_id[26])
flag_outlier(d$animal_id[26], max_speed = 1000, value = 1)
out <- check_animal(d$animal_id[26])
flag_outlier(d$animal_id[26], max_speed = 200, value = 1)
out <- check_animal(d$animal_id[26])
# the speeds are impossibly high due to error in locations 1 second apart
min(diff(d$data[[26]]$timestamp))
hist(log(out$speed), breaks = 100)
# dropping all locations less than two minutes away from the previous one
# ideally, we would model the GPS error and include in the telemetry dataset,
# but we don't have the data and dropping th rows will only affect the averages
# minimally (other than possibly removing the possibility for some large biases)
hist(log(out$speed[diff(d$tel[[26]]$timestamp) < 60]), breaks = 100)
hist(log(as.numeric(diff(d$tel[[26]]$timestamp))), breaks = 100)
abline(v = log(120), col = 'red')
d$data[[26]]$outlier[diff(d$data[[26]]$timestamp) < 120] <- 1
out <- check_animal(d$animal_id[26])

# 28
out <- check_animal(d$animal_id[28])
flag_outlier(d$animal_id[28], max_speed = 6, value = 1)
out <- check_animal(d$animal_id[28])
hist(log(as.numeric(diff(d$tel[[28]]$timestamp))), breaks = 100)
abline(v = log(120), col = 'red')
d$data[[28]]$outlier[diff(d$data[[28]]$timestamp) < 120] <- 1
out <- check_animal(d$animal_id[28])

# 29
out <- check_animal(d$animal_id[29])
flag_outlier(d$animal_id[29], max_speed = 6, value = 1)
out <- check_animal(d$animal_id[29])
flag_outlier(d$animal_id[29], max_speed = 6, value = 1)
out <- check_animal(d$animal_id[29])

# 32
out <- check_animal(d$animal_id[32])
hist(log(as.numeric(diff(d$tel[[32]]$timestamp))), breaks = 100)
abline(v = log(120), col = 'red')
d$data[[32]]$outlier[diff(d$data[[32]]$timestamp) < 120] <- 1
out <- check_animal(d$animal_id[32])

# drop all elephant locations less than 2 minutes apart
d$data <- map_if(d$data, grepl("elephants", d$animal_id), function(.d) {
  .d$outlier[diff(.d$timestamp) < 120] <- 1
  .d
})

# 54
out <- check_animal(d$animal_id[54])
plot_adj(d$animal_id[54], max_speed = 5)
plot_adj(d$animal_id[54], max_speed = 15)
flag_outlier(d$animal_id[54], max_speed = 5, value = 2)
out <- check_animal(d$animal_id[54])
plot_adj(d$animal_id[54], max_speed = 3) # ok

# 56
out <- check_animal(d$animal_id[56])
flag_outlier(d$animal_id[56], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[56])
plot_adj(d$animal_id[56], max_speed = 3)
flag_outlier(d$animal_id[56], max_speed = 3, value = 1)
out <- check_animal(d$animal_id[56])

# 62
out <- check_animal(d$animal_id[62])
plot(d$tel[[62]], type = "l", error = FALSE, col = "black")
plot(d$tel[[62]][1:100, ], type = "l", error = FALSE, add = TRUE, col = "blue")
plot(d$tel[[62]][nrow(d$tel[[62]]) - 175:0, ], type = "l", error = FALSE,
     add = TRUE, lwd = 2)
# excursion is at the end. would have excluded it if it was at the start

# 87
out <- check_animal(d$animal_id[87])
flag_outlier(d$animal_id[87], max_distance = 2000, value = 1)
out <- check_animal(d$animal_id[87])

# 98
out <- check_animal(d$animal_id[98])
plot_adj(d$animal_id[98], max_distance = 2700)
flag_outlier(d$animal_id[98], max_distance = 2700, value = 1)
out <- check_animal(d$animal_id[98])

# 99
out <- check_animal(d$animal_id[99])
plot_adj(d$animal_id[99], max_distance = 1700)
flag_outlier(d$animal_id[99], max_distance = 1700, value = 1)
out <- check_animal(d$animal_id[99])

# 100
out <- check_animal(d$animal_id[100])
plot_adj(d$animal_id[100], max_speed = 0, map = TRUE, many = TRUE) # ok

# 133
out <- check_animal(d$animal_id[133])
flag_outlier(d$animal_id[133], max_distance = 2e5, value = 1)
out <- check_animal(d$animal_id[133])

# 135
out <- check_animal(d$animal_id[135])
plot_adj(d$animal_id[135], max_distance = 3e4)
flag_outlier(d$animal_id[135], max_distance = 1e4, value = 2)
out <- check_animal(d$animal_id[135])

# 141
out <- check_animal(d$animal_id[141])
plot_start(141, 10)
flag_outlier(d$animal_id[141], max_distance = 50000, value = 1)
out <- check_animal(d$animal_id[141])
flag_outlier(d$animal_id[141], max_distance = 25000, value = 1)
out <- check_animal(d$animal_id[141])

# 152
out <- check_animal(d$animal_id[152])
plot_start(152)
flag_outlier(d$animal_id[152], max_distance = 1e5, value = 1)
out <- check_animal(d$animal_id[152])

# 159
out <- check_animal(d$animal_id[159])
plot_adj(d$animal_id[159], max_speed = 0, many = TRUE, map = TRUE)

# 160 -- arctic foxes travel incredible distances
out <- check_animal(d$animal_id[160])
plot_adj(d$animal_id[160], max_speed = 0, many = TRUE, map = TRUE)

# dropping some extreme points because the AKDE may struggle to fit otherwise
# 175
out <- check_animal(d$animal_id[175])
plot_adj(d$animal_id[175], max_distance = 1e5, map = TRUE)
flag_outlier(d$animal_id[175], max_distance = 1e5, value = 1)
out <- check_animal(d$animal_id[175])

# 179
out <- check_animal(d$animal_id[179])
plot_adj(d$animal_id[179], max_distance = 4e4, map = TRUE)
flag_outlier(d$animal_id[179], max_distance = 4e4, value = 1)
out <- check_animal(d$animal_id[179])

# 195
out <- check_animal(d$animal_id[195])
flag_outlier(d$animal_id[195], max_speed = 1, value = 1)
out <- check_animal(d$animal_id[195])

# 196
out <- check_animal(d$animal_id[196])
flag_outlier(d$animal_id[196], max_speed = 2, value = 1)
out <- check_animal(d$animal_id[196])

# 265
out <- check_animal(d$animal_id[265])
flag_outlier(d$animal_id[265], max_speed = 0.5, value = 1)
out <- check_animal(d$animal_id[265])

# 267
out <- check_animal(d$animal_id[267])
flag_outlier(d$animal_id[267], max_speed = 0.5, value = 1)
out <- check_animal(d$animal_id[267])

# 270
out <- check_animal(d$animal_id[270])
flag_outlier(d$animal_id[270], max_speed = 0.5, value = 1)
out <- check_animal(d$animal_id[270])

# 297
d$data[[297]] <- filter(d$data[[297]], ! is.na(long)) # missing some coordinates
out <- check_animal(d$animal_id[297])
flag_outlier(d$animal_id[297], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[297])

# 299
d$data[[299]] <- filter(d$data[[299]], ! is.na(long)) # missing some coordinates
out <- check_animal(d$animal_id[299])
flag_outlier(d$animal_id[299], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[299])

# filter out all other rows with missing coords
d$data <- map(d$data, function(.d) filter(.d, ! is.na(long)))

# 300
out <- check_animal(d$animal_id[300])
flag_outlier(d$animal_id[300], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[300])
plot_adj(d$animal_id[300], max_speed = 6, map = TRUE)

# 301
out <- check_animal(d$animal_id[301])
flag_outlier(d$animal_id[301], max_speed = 20, value = 1)
out <- check_animal(d$animal_id[301])
flag_outlier(d$animal_id[301], max_speed = 12, value = 1)
out <- check_animal(d$animal_id[301])

# 303
out <- check_animal(d$animal_id[303])
flag_outlier(d$animal_id[303], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[303])

# 304
out <- check_animal(d$animal_id[304])
flag_outlier(d$animal_id[304], max_speed = 200, value = 1)
out <- check_animal(d$animal_id[304])
flag_outlier(d$animal_id[304], max_speed = 11, value = 1)
out <- check_animal(d$animal_id[304])

# 305
out <- check_animal(d$animal_id[305])
flag_outlier(d$animal_id[305], max_speed = 30, value = 1)
out <- check_animal(d$animal_id[305])
flag_outlier(d$animal_id[305], max_speed = 15, value = 1)
out <- check_animal(d$animal_id[305])
flag_outlier(d$animal_id[305], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[305])

# 311
out <- check_animal(d$animal_id[311])
flag_outlier(d$animal_id[311], max_speed = 40, value = 1)
out <- check_animal(d$animal_id[311])
plot_adj(d$animal_id[311], max_speed = 20) # ok

# 312
out <- check_animal(d$animal_id[312])
flag_outlier(d$animal_id[312], max_distance = 1e4, value = 1)
out <- check_animal(d$animal_id[312])

# 313
out <- check_animal(d$animal_id[313])
flag_outlier(d$animal_id[313], max_speed = 40, value = 1)
out <- check_animal(d$animal_id[313])

# 315
out <- check_animal(d$animal_id[315])
flag_outlier(d$animal_id[315], max_speed = 40, value = 1)
out <- check_animal(d$animal_id[315])

# 317
out <- check_animal(d$animal_id[317])
plot_adj(d$animal_id[317], max_speed = 17)
flag_outlier(d$animal_id[317], max_speed = 17, value = 1)
out <- check_animal(d$animal_id[317])

# 318
out <- check_animal(d$animal_id[318])
flag_outlier(d$animal_id[318], max_distance = 5e4, value = 1)
out <- check_animal(d$animal_id[318])
flag_outlier(d$animal_id[318], max_speed = 30, value = 1)
out <- check_animal(d$animal_id[318])
plot_adj(d$animal_id[318], max_speed = 10)
flag_outlier(d$animal_id[318], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[318])

# 319
out <- check_animal(d$animal_id[319])
flag_outlier(d$animal_id[319], max_speed = 15, value = 1)
out <- check_animal(d$animal_id[319])
flag_outlier(d$animal_id[319], max_speed = 8, value = 1)
out <- check_animal(d$animal_id[319])

# 320
# 10 m/s looks like a good cutoff
out <- check_animal(d$animal_id[320])
flag_outlier(d$animal_id[320], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[319])

# 321
out <- check_animal(d$animal_id[321])
flag_outlier(d$animal_id[321], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[319])

# 323
out <- check_animal(d$animal_id[323])
flag_outlier(d$animal_id[323], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[323])

# 324
out <- check_animal(d$animal_id[324])
flag_outlier(d$animal_id[324], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[324])

# 326
out <- check_animal(d$animal_id[326])
plot_adj(d$animal_id[326], max_speed = 10) # ok

# 327
out <- check_animal(d$animal_id[327])
flag_outlier(d$animal_id[327], max_speed = 100, value = 1)
out <- check_animal(d$animal_id[327])
flag_outlier(d$animal_id[327], max_distance = 5e5, value = 1)
out <- check_animal(d$animal_id[327])
flag_outlier(d$animal_id[327], max_speed = 20, value = 1)
out <- check_animal(d$animal_id[327])
flag_outlier(d$animal_id[327], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[327])

# 328
out <- check_animal(d$animal_id[328])
flag_outlier(d$animal_id[328], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[328])

# 330
out <- check_animal(d$animal_id[330])
flag_outlier(d$animal_id[330], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[330])

# 331
out <- check_animal(d$animal_id[331])
flag_outlier(d$animal_id[331], max_distance = 1e6, value = 1)
out <- check_animal(d$animal_id[331])
flag_outlier(d$animal_id[331], max_speed = 20, value = 1)
out <- check_animal(d$animal_id[331])
flag_outlier(d$animal_id[331], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[331])

# 333
out <- check_animal(d$animal_id[333])
flag_outlier(d$animal_id[333], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[333])

# 334
out <- check_animal(d$animal_id[334])

# 335
out <- check_animal(d$animal_id[335])
flag_outlier(d$animal_id[335], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[335])

# 336
out <- check_animal(d$animal_id[336])

# 337
out <- check_animal(d$animal_id[337])
flag_outlier(d$animal_id[337], max_speed = 25, value = 1)
out <- check_animal(d$animal_id[337])
flag_outlier(d$animal_id[337], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[337])

# 338
out <- check_animal(d$animal_id[338])
flag_outlier(d$animal_id[338], max_speed = 60, value = 1)
out <- check_animal(d$animal_id[338])
plot_adj(d$animal_id[338], max_speed = 30)
flag_outlier(d$animal_id[338], max_speed = 30, value = 1)
out <- check_animal(d$animal_id[338])
plot_adj(d$animal_id[338], max_speed = 15) # ok

# 340
out <- check_animal(d$animal_id[340])
flag_outlier(d$animal_id[340], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[340])

# 342
out <- check_animal(d$animal_id[342])
flag_outlier(d$animal_id[342], max_speed = 20, value = 1)
out <- check_animal(d$animal_id[342])
flag_outlier(d$animal_id[342], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[342])
flag_outlier(d$animal_id[342], max_speed = 8, value = 1)
out <- check_animal(d$animal_id[342])
flag_outlier(d$animal_id[342], max_speed = 6, value = 1)
out <- check_animal(d$animal_id[342])

# 344
out <- check_animal(d$animal_id[344])
flag_outlier(d$animal_id[344], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[344])

# 345
out <- check_animal(d$animal_id[345])
flag_outlier(d$animal_id[345], max_speed = 20, value = 1)
out <- check_animal(d$animal_id[345])

# 346
out <- check_animal(d$animal_id[346])
flag_outlier(d$animal_id[346], max_speed = 20, value = 1)
out <- check_animal(d$animal_id[346])
flag_outlier(d$animal_id[346], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[346])
flag_outlier(d$animal_id[346], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[346])
# the dataset is so large that we can be a bit liberal with the trimming
flag_outlier(d$animal_id[346], max_speed = 8, value = 1)
out <- check_animal(d$animal_id[346])

# 347
out <- check_animal(d$animal_id[347])
flag_outlier(d$animal_id[347], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[347])
flag_outlier(d$animal_id[347], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[347])

# 349
out <- check_animal(d$animal_id[349])
flag_outlier(d$animal_id[349], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[349])

# 352
out <- check_animal(d$animal_id[352])
flag_outlier(d$animal_id[352], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[352])

# 353
out <- check_animal(d$animal_id[353])
flag_outlier(d$animal_id[353], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[353])

# 354
out <- check_animal(d$animal_id[354])
flag_outlier(d$animal_id[354], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[354])

# 356
out <- check_animal(d$animal_id[356])
flag_outlier(d$animal_id[356], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[356])

# 357
out <- check_animal(d$animal_id[357])
flag_outlier(d$animal_id[357], max_speed = 10, value = 1)
out <- check_animal(d$animal_id[357])

# 411
out <- check_animal(d$animal_id[411])
flag_outlier(d$animal_id[411], max_speed = 1, value = 1)
out <- check_animal(d$animal_id[411])

# 421
out <- check_animal(d$animal_id[421]) # all high speeds are at small deviations

# 424
out <- check_animal(d$animal_id[424])
flag_outlier(d$animal_id[424], max_speed = 20, value = 1)
out <- check_animal(d$animal_id[424])

# 429
out <- check_animal(d$animal_id[429])
flag_outlier(d$animal_id[429], max_speed = 5, value = 1)
out <- check_animal(d$animal_id[429])

# 456
out <- check_animal(d$animal_id[456])
flag_outlier(d$animal_id[456], max_speed = 2, value = 1)
out <- check_animal(d$animal_id[456])

# 477
out <- check_animal(d$animal_id[477])
flag_outlier(d$animal_id[477], max_speed = 5, value = 1)
out <- check_animal(d$animal_id[477])

# 481
out <- check_animal(d$animal_id[481])
flag_outlier(d$animal_id[481], max_speed = 4, value = 1)
out <- check_animal(d$animal_id[481])

# 510
out <- check_animal(d$animal_id[510])
flag_outlier(d$animal_id[510], max_speed = 0.01, value = 1)
out <- check_animal(d$animal_id[510])
flag_outlier(d$animal_id[510], max_distance = 1000, value = 1)
out <- check_animal(d$animal_id[510])

# 517
out <- check_animal(d$animal_id[517])
flag_outlier(d$animal_id[517], max_speed = 0.05, value = 1)
out <- check_animal(d$animal_id[517])
flag_outlier(d$animal_id[517], max_speed = 0.002, value = 1)
out <- check_animal(d$animal_id[517])

# 518
out <- check_animal(d$animal_id[518])
flag_outlier(d$animal_id[518], max_speed = 0.001, value = 1)
out <- check_animal(d$animal_id[518])

# 520
out <- check_animal(d$animal_id[520])
flag_outlier(d$animal_id[520], max_speed = 0.01, value = 1)
out <- check_animal(d$animal_id[520])

# 523
out <- check_animal(d$animal_id[523])
flag_outlier(d$animal_id[523], max_distance = 40, value = 1)
out <- check_animal(d$animal_id[523])

# 549
out <- check_animal(d$animal_id[549])
flag_outlier(d$animal_id[549], max_distance = 1e4, value = 1)
out <- check_animal(d$animal_id[549])

# 562
out <- check_animal(d$animal_id[562])
flag_outlier(d$animal_id[562], max_speed = 2, value = 1)

# 568
out <- check_animal(d$animal_id[568])
flag_outlier(d$animal_id[568], max_speed = 2, value = 1)
out <- check_animal(d$animal_id[568])

# 595
out <- check_animal(d$animal_id[595])
flag_outlier(d$animal_id[595], max_speed = 5, value = 1)
out <- check_animal(d$animal_id[595])

# 597
out <- check_animal(d$animal_id[597])
flag_outlier(d$animal_id[597], max_distance = 1e3, value = 1)
out <- check_animal(d$animal_id[597])

# 598
out <- check_animal(d$animal_id[598])
flag_outlier(d$animal_id[598], max_distance = 1e3, value = 1)
out <- check_animal(d$animal_id[598])

# 612
out <- check_animal(d$animal_id[612])
flag_outlier(d$animal_id[612], max_speed = 3, value = 1)
out <- check_animal(d$animal_id[612])

# 616
out <- check_animal(d$animal_id[616])
flag_outlier(d$animal_id[616], max_speed = 5, value = 1)
out <- check_animal(d$animal_id[616])

# 625
out <- check_animal(d$animal_id[625])
# distances with high speeds are very short: likely just GPS error
plot_adj(d$animal_id[625], max_speed = 40)
plot_adj(d$animal_id[625], max_speed = 40, map = TRUE, n_adj = 3)

# 649
out <- check_animal(d$animal_id[649])
flag_outlier(d$animal_id[649], max_distance = 3000, value = 1)
out <- check_animal(d$animal_id[649])

# 651
out <- check_animal(d$animal_id[651])
flag_outlier(d$animal_id[651], max_distance = 3000, value = 1)
out <- check_animal(d$animal_id[651])

# 652
out <- check_animal(d$animal_id[652])
flag_outlier(d$animal_id[652], max_distance = 3000, value = 1)
out <- check_animal(d$animal_id[652])

# 653
out <- check_animal(d$animal_id[653])

# 654
out <- check_animal(d$animal_id[654])
flag_outlier(d$animal_id[654], max_distance = 20000, value = 1)
out <- check_animal(d$animal_id[654])

# 655
out <- check_animal(d$animal_id[655])
flag_outlier(d$animal_id[655], max_distance = 30000, value = 1)
out <- check_animal(d$animal_id[655])

# 661
out <- check_animal(d$animal_id[661])
flag_outlier(d$animal_id[661], max_distance = 4e4, value = 1)
out <- check_animal(d$animal_id[661])
flag_outlier(d$animal_id[661], max_distance = 15e3, value = 1)
out <- check_animal(d$animal_id[661])

# 670
out <- check_animal(d$animal_id[670])
flag_outlier(d$animal_id[670], max_distance = 2e4, value = 1)
out <- check_animal(d$animal_id[670])

# save the movebank dataset ----
saveRDS(d, paste0('data/telemetry-data/movebank-data-', Sys.Date(), '.rds'))
