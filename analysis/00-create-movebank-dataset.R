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

import_movebank <- function(file_name) {
  d_name <- gsub('data/telemetry-data/public-movebank-data/', '', file_name) %>%
    gsub('\\.csv', '', .)
  
  .d <- readr::read_csv(file_name, name_repair = 'universal_quiet',
                        show_col_types = FALSE)
  
  if ('outier' %in% colnames(.d)) {
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
      outlier = FALSE, clean = FALSE
    )
}

d <-
  list.files('data/telemetry-data/public-movebank-data',
             pattern = '\\.csv',
             full.names = TRUE) %>%
  stringr::str_subset('^(?!.*-reference-data)') %>% # drop reference data
  map(import_movebank) %>%
  list_rbind() %>%
  # collect the datasets for each individual animal
  nest(data = ! c(animal_id, species, dataset_name, clean)) %>%
  mutate(tel = map(data, \(.d) {
    # convert to telemetries quietly (will be doing data cleaning later)
    suppressMessages(as.telemetry(.d, mark.rm = TRUE))
  }),
  n_locations = map_int(tel, nrow),
  n_days = map_dbl(data, \(.d) {
    diff(range(.d$timestamp), units = 'days')
  }))
d

# ensure there are no duplicated animal names across studies (should be 0)
stopifnot(sum(duplicated(d$animal_id)) == 0)

stop('need to add missing species names:')
filter(d, is.na(species))$dataset_name

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

# save the movebank dataset
saveRDS(d, paste0('movebank-data-', Sys.Date(), '.rds'))
