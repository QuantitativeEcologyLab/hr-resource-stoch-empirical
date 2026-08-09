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
source('functions/check_animal.R')
source('functions/flag_outlier.R')
source('functions/remove_outlier_flags.R')

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
      outlier = FALSE, clean = FALSE
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
if(! exists('d')) d <- read_rds('movebank-data-2026-08-01.rds')

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

#' *continue here*

# 40
out <- check_animal(d$animal_id[40])

# 41
out <- check_animal(d$animal_id[41])

# 42
out <- check_animal(d$animal_id[42])

# 43
out <- check_animal(d$animal_id[43])

# 44
out <- check_animal(d$animal_id[44])

# 45
out <- check_animal(d$animal_id[45])

# 46
out <- check_animal(d$animal_id[46])

# 47
out <- check_animal(d$animal_id[47])

# 48
out <- check_animal(d$animal_id[48])

# 49
out <- check_animal(d$animal_id[49])

# 50
out <- check_animal(d$animal_id[50])

# 51
out <- check_animal(d$animal_id[51])

# 52
out <- check_animal(d$animal_id[52])

# 53
out <- check_animal(d$animal_id[53])

# 54
out <- check_animal(d$animal_id[54])

# 55
out <- check_animal(d$animal_id[55])

# 56
out <- check_animal(d$animal_id[56])

# 57
out <- check_animal(d$animal_id[57])

# 58
out <- check_animal(d$animal_id[58])

# 59
out <- check_animal(d$animal_id[59])

# 60
out <- check_animal(d$animal_id[60])

# 61
out <- check_animal(d$animal_id[61])

# 62
out <- check_animal(d$animal_id[62])

# 63
out <- check_animal(d$animal_id[63])

# 64
out <- check_animal(d$animal_id[64])

# 65
out <- check_animal(d$animal_id[65])

# 66
out <- check_animal(d$animal_id[66])

# 67
out <- check_animal(d$animal_id[67])

# 68
out <- check_animal(d$animal_id[68])

# 69
out <- check_animal(d$animal_id[69])

# 70
out <- check_animal(d$animal_id[70])

# 71
out <- check_animal(d$animal_id[71])

# 72
out <- check_animal(d$animal_id[72])

# 73
out <- check_animal(d$animal_id[73])

# 74
out <- check_animal(d$animal_id[74])

# 75
out <- check_animal(d$animal_id[75])

# 76
out <- check_animal(d$animal_id[76])

# 77
out <- check_animal(d$animal_id[77])

# 78
out <- check_animal(d$animal_id[78])

# 79
out <- check_animal(d$animal_id[79])

# 80
out <- check_animal(d$animal_id[80])

# 81
out <- check_animal(d$animal_id[81])

# 82
out <- check_animal(d$animal_id[82])

# 83
out <- check_animal(d$animal_id[83])

# 84
out <- check_animal(d$animal_id[84])

# 85
out <- check_animal(d$animal_id[85])

# 86
out <- check_animal(d$animal_id[86])

# 87
out <- check_animal(d$animal_id[87])

# 88
out <- check_animal(d$animal_id[88])

# 89
out <- check_animal(d$animal_id[89])

# 90
out <- check_animal(d$animal_id[90])

# 91
out <- check_animal(d$animal_id[91])

# 92
out <- check_animal(d$animal_id[92])

# 93
out <- check_animal(d$animal_id[93])

# 94
out <- check_animal(d$animal_id[94])

# 95
out <- check_animal(d$animal_id[95])

# 96
out <- check_animal(d$animal_id[96])

# 97
out <- check_animal(d$animal_id[97])

# 98
out <- check_animal(d$animal_id[98])

# 99
out <- check_animal(d$animal_id[99])

# 100
out <- check_animal(d$animal_id[100])

# 101
out <- check_animal(d$animal_id[101])

# 102
out <- check_animal(d$animal_id[102])

# 103
out <- check_animal(d$animal_id[103])

# 104
out <- check_animal(d$animal_id[104])

# 105
out <- check_animal(d$animal_id[105])

# 106
out <- check_animal(d$animal_id[106])

# 107
out <- check_animal(d$animal_id[107])

# 108
out <- check_animal(d$animal_id[108])

# 109
out <- check_animal(d$animal_id[109])

# 110
out <- check_animal(d$animal_id[110])

# 111
out <- check_animal(d$animal_id[111])

# 112
out <- check_animal(d$animal_id[112])

# 113
out <- check_animal(d$animal_id[113])

# 114
out <- check_animal(d$animal_id[114])

# 115
out <- check_animal(d$animal_id[115])

# 116
out <- check_animal(d$animal_id[116])

# 117
out <- check_animal(d$animal_id[117])

# 118
out <- check_animal(d$animal_id[118])

# 119
out <- check_animal(d$animal_id[119])

# 120
out <- check_animal(d$animal_id[120])

# 121
out <- check_animal(d$animal_id[121])

# 122
out <- check_animal(d$animal_id[122])

# 123
out <- check_animal(d$animal_id[123])

# 124
out <- check_animal(d$animal_id[124])

# 125
out <- check_animal(d$animal_id[125])

# 126
out <- check_animal(d$animal_id[126])

# 127
out <- check_animal(d$animal_id[127])

# 128
out <- check_animal(d$animal_id[128])

# 129
out <- check_animal(d$animal_id[129])

# 130
out <- check_animal(d$animal_id[130])

# 131
out <- check_animal(d$animal_id[131])

# 132
out <- check_animal(d$animal_id[132])

# 133
out <- check_animal(d$animal_id[133])

# 134
out <- check_animal(d$animal_id[134])

# 135
out <- check_animal(d$animal_id[135])

# 136
out <- check_animal(d$animal_id[136])

# 137
out <- check_animal(d$animal_id[137])

# 138
out <- check_animal(d$animal_id[138])

# 139
out <- check_animal(d$animal_id[139])

# 140
out <- check_animal(d$animal_id[140])

# 141
out <- check_animal(d$animal_id[141])

# 142
out <- check_animal(d$animal_id[142])

# 143
out <- check_animal(d$animal_id[143])

# 144
out <- check_animal(d$animal_id[144])

# 145
out <- check_animal(d$animal_id[145])

# 146
out <- check_animal(d$animal_id[146])

# 147
out <- check_animal(d$animal_id[147])

# 148
out <- check_animal(d$animal_id[148])

# 149
out <- check_animal(d$animal_id[149])

# 150
out <- check_animal(d$animal_id[150])

# 151
out <- check_animal(d$animal_id[151])

# 152
out <- check_animal(d$animal_id[152])

# 153
out <- check_animal(d$animal_id[153])

# 154
out <- check_animal(d$animal_id[154])

# 155
out <- check_animal(d$animal_id[155])

# 156
out <- check_animal(d$animal_id[156])

# 157
out <- check_animal(d$animal_id[157])

# 158
out <- check_animal(d$animal_id[158])

# 159
out <- check_animal(d$animal_id[159])

# 160
out <- check_animal(d$animal_id[160])

# 161
out <- check_animal(d$animal_id[161])

# 162
out <- check_animal(d$animal_id[162])

# 163
out <- check_animal(d$animal_id[163])

# 164
out <- check_animal(d$animal_id[164])

# 165
out <- check_animal(d$animal_id[165])

# 166
out <- check_animal(d$animal_id[166])

# 167
out <- check_animal(d$animal_id[167])

# 168
out <- check_animal(d$animal_id[168])

# 169
out <- check_animal(d$animal_id[169])

# 170
out <- check_animal(d$animal_id[170])

# 171
out <- check_animal(d$animal_id[171])

# 172
out <- check_animal(d$animal_id[172])

# 173
out <- check_animal(d$animal_id[173])

# 174
out <- check_animal(d$animal_id[174])

# 175
out <- check_animal(d$animal_id[175])

# 176
out <- check_animal(d$animal_id[176])

# 177
out <- check_animal(d$animal_id[177])

# 178
out <- check_animal(d$animal_id[178])

# 179
out <- check_animal(d$animal_id[179])

# 180
out <- check_animal(d$animal_id[180])

# 181
out <- check_animal(d$animal_id[181])

# 182
out <- check_animal(d$animal_id[182])

# 183
out <- check_animal(d$animal_id[183])

# 184
out <- check_animal(d$animal_id[184])

# 185
out <- check_animal(d$animal_id[185])

# 186
out <- check_animal(d$animal_id[186])

# 187
out <- check_animal(d$animal_id[187])

# 188
out <- check_animal(d$animal_id[188])

# 189
out <- check_animal(d$animal_id[189])

# 190
out <- check_animal(d$animal_id[190])

# 191
out <- check_animal(d$animal_id[191])

# 192
out <- check_animal(d$animal_id[192])

# 193
out <- check_animal(d$animal_id[193])

# 194
out <- check_animal(d$animal_id[194])

# 195
out <- check_animal(d$animal_id[195])

# 196
out <- check_animal(d$animal_id[196])

# 197
out <- check_animal(d$animal_id[197])

# 198
out <- check_animal(d$animal_id[198])

# 199
out <- check_animal(d$animal_id[199])

# 200
out <- check_animal(d$animal_id[200])

# 201
out <- check_animal(d$animal_id[201])

# 202
out <- check_animal(d$animal_id[202])

# 203
out <- check_animal(d$animal_id[203])

# 204
out <- check_animal(d$animal_id[204])

# 205
out <- check_animal(d$animal_id[205])

# 206
out <- check_animal(d$animal_id[206])

# 207
out <- check_animal(d$animal_id[207])

# 208
out <- check_animal(d$animal_id[208])

# 209
out <- check_animal(d$animal_id[209])

# 210
out <- check_animal(d$animal_id[210])

# 211
out <- check_animal(d$animal_id[211])

# 212
out <- check_animal(d$animal_id[212])

# 213
out <- check_animal(d$animal_id[213])

# 214
out <- check_animal(d$animal_id[214])

# 215
out <- check_animal(d$animal_id[215])

# 216
out <- check_animal(d$animal_id[216])

# 217
out <- check_animal(d$animal_id[217])

# 218
out <- check_animal(d$animal_id[218])

# 219
out <- check_animal(d$animal_id[219])

# 220
out <- check_animal(d$animal_id[220])

# 221
out <- check_animal(d$animal_id[221])

# 222
out <- check_animal(d$animal_id[222])

# 223
out <- check_animal(d$animal_id[223])

# 224
out <- check_animal(d$animal_id[224])

# 225
out <- check_animal(d$animal_id[225])

# 226
out <- check_animal(d$animal_id[226])

# 227
out <- check_animal(d$animal_id[227])

# 228
out <- check_animal(d$animal_id[228])

# 229
out <- check_animal(d$animal_id[229])

# 230
out <- check_animal(d$animal_id[230])

# 231
out <- check_animal(d$animal_id[231])

# 232
out <- check_animal(d$animal_id[232])

# 233
out <- check_animal(d$animal_id[233])

# 234
out <- check_animal(d$animal_id[234])

# 235
out <- check_animal(d$animal_id[235])

# 236
out <- check_animal(d$animal_id[236])

# 237
out <- check_animal(d$animal_id[237])

# 238
out <- check_animal(d$animal_id[238])

# 239
out <- check_animal(d$animal_id[239])

# 240
out <- check_animal(d$animal_id[240])

# 241
out <- check_animal(d$animal_id[241])

# 242
out <- check_animal(d$animal_id[242])

# 243
out <- check_animal(d$animal_id[243])

# 244
out <- check_animal(d$animal_id[244])

# 245
out <- check_animal(d$animal_id[245])

# 246
out <- check_animal(d$animal_id[246])

# 247
out <- check_animal(d$animal_id[247])

# 248
out <- check_animal(d$animal_id[248])

# 249
out <- check_animal(d$animal_id[249])

# 250
out <- check_animal(d$animal_id[250])

# 251
out <- check_animal(d$animal_id[251])

# 252
out <- check_animal(d$animal_id[252])

# 253
out <- check_animal(d$animal_id[253])

# 254
out <- check_animal(d$animal_id[254])

# 255
out <- check_animal(d$animal_id[255])

# 256
out <- check_animal(d$animal_id[256])

# 257
out <- check_animal(d$animal_id[257])

# 258
out <- check_animal(d$animal_id[258])

# 259
out <- check_animal(d$animal_id[259])

# 260
out <- check_animal(d$animal_id[260])

# 261
out <- check_animal(d$animal_id[261])

# 262
out <- check_animal(d$animal_id[262])

# 263
out <- check_animal(d$animal_id[263])

# 264
out <- check_animal(d$animal_id[264])

# 265
out <- check_animal(d$animal_id[265])

# 266
out <- check_animal(d$animal_id[266])

# 267
out <- check_animal(d$animal_id[267])

# 268
out <- check_animal(d$animal_id[268])

# 269
out <- check_animal(d$animal_id[269])

# 270
out <- check_animal(d$animal_id[270])

# 271
out <- check_animal(d$animal_id[271])

# 272
out <- check_animal(d$animal_id[272])

# 273
out <- check_animal(d$animal_id[273])

# 274
out <- check_animal(d$animal_id[274])

# 275
out <- check_animal(d$animal_id[275])

# 276
out <- check_animal(d$animal_id[276])

# 277
out <- check_animal(d$animal_id[277])

# 278
out <- check_animal(d$animal_id[278])

# 279
out <- check_animal(d$animal_id[279])

# 280
out <- check_animal(d$animal_id[280])

# 281
out <- check_animal(d$animal_id[281])

# 282
out <- check_animal(d$animal_id[282])

# 283
out <- check_animal(d$animal_id[283])

# 284
out <- check_animal(d$animal_id[284])

# 285
out <- check_animal(d$animal_id[285])

# 286
out <- check_animal(d$animal_id[286])

# 287
out <- check_animal(d$animal_id[287])

# 288
out <- check_animal(d$animal_id[288])

# 289
out <- check_animal(d$animal_id[289])

# 290
out <- check_animal(d$animal_id[290])

# 291
out <- check_animal(d$animal_id[291])

# 292
out <- check_animal(d$animal_id[292])

# 293
out <- check_animal(d$animal_id[293])

# 294
out <- check_animal(d$animal_id[294])

# 295
out <- check_animal(d$animal_id[295])

# 296
out <- check_animal(d$animal_id[296])

# 297
out <- check_animal(d$animal_id[297])

# 298
out <- check_animal(d$animal_id[298])

# 299
out <- check_animal(d$animal_id[299])

# 300
out <- check_animal(d$animal_id[300])

# 301
out <- check_animal(d$animal_id[301])

# 302
out <- check_animal(d$animal_id[302])

# 303
out <- check_animal(d$animal_id[303])

# 304
out <- check_animal(d$animal_id[304])

# 305
out <- check_animal(d$animal_id[305])

# 306
out <- check_animal(d$animal_id[306])

# 307
out <- check_animal(d$animal_id[307])

# 308
out <- check_animal(d$animal_id[308])

# 309
out <- check_animal(d$animal_id[309])

# 310
out <- check_animal(d$animal_id[310])

# 311
out <- check_animal(d$animal_id[311])

# 312
out <- check_animal(d$animal_id[312])

# 313
out <- check_animal(d$animal_id[313])

# 314
out <- check_animal(d$animal_id[314])

# 315
out <- check_animal(d$animal_id[315])

# 316
out <- check_animal(d$animal_id[316])

# 317
out <- check_animal(d$animal_id[317])

# 318
out <- check_animal(d$animal_id[318])

# 319
out <- check_animal(d$animal_id[319])

# 320
out <- check_animal(d$animal_id[320])

# 321
out <- check_animal(d$animal_id[321])

# 322
out <- check_animal(d$animal_id[322])

# 323
out <- check_animal(d$animal_id[323])

# 324
out <- check_animal(d$animal_id[324])

# 325
out <- check_animal(d$animal_id[325])

# 326
out <- check_animal(d$animal_id[326])

# 327
out <- check_animal(d$animal_id[327])

# 328
out <- check_animal(d$animal_id[328])

# 329
out <- check_animal(d$animal_id[329])

# 330
out <- check_animal(d$animal_id[330])

# 331
out <- check_animal(d$animal_id[331])

# 332
out <- check_animal(d$animal_id[332])

# 333
out <- check_animal(d$animal_id[333])

# 334
out <- check_animal(d$animal_id[334])

# 335
out <- check_animal(d$animal_id[335])

# 336
out <- check_animal(d$animal_id[336])

# 337
out <- check_animal(d$animal_id[337])

# 338
out <- check_animal(d$animal_id[338])

# 339
out <- check_animal(d$animal_id[339])

# 340
out <- check_animal(d$animal_id[340])

# 341
out <- check_animal(d$animal_id[341])

# 342
out <- check_animal(d$animal_id[342])

# 343
out <- check_animal(d$animal_id[343])

# 344
out <- check_animal(d$animal_id[344])

# 345
out <- check_animal(d$animal_id[345])

# 346
out <- check_animal(d$animal_id[346])

# 347
out <- check_animal(d$animal_id[347])

# 348
out <- check_animal(d$animal_id[348])

# 349
out <- check_animal(d$animal_id[349])

# 350
out <- check_animal(d$animal_id[350])

# 351
out <- check_animal(d$animal_id[351])

# 352
out <- check_animal(d$animal_id[352])

# 353
out <- check_animal(d$animal_id[353])

# 354
out <- check_animal(d$animal_id[354])

# 355
out <- check_animal(d$animal_id[355])

# 356
out <- check_animal(d$animal_id[356])

# 357
out <- check_animal(d$animal_id[357])

# 358
out <- check_animal(d$animal_id[358])

# 359
out <- check_animal(d$animal_id[359])

# 360
out <- check_animal(d$animal_id[360])

# 361
out <- check_animal(d$animal_id[361])

# 362
out <- check_animal(d$animal_id[362])

# 363
out <- check_animal(d$animal_id[363])

# 364
out <- check_animal(d$animal_id[364])

# 365
out <- check_animal(d$animal_id[365])

# 366
out <- check_animal(d$animal_id[366])

# 367
out <- check_animal(d$animal_id[367])

# 368
out <- check_animal(d$animal_id[368])

# 369
out <- check_animal(d$animal_id[369])

# 370
out <- check_animal(d$animal_id[370])

# 371
out <- check_animal(d$animal_id[371])

# 372
out <- check_animal(d$animal_id[372])

# 373
out <- check_animal(d$animal_id[373])

# 374
out <- check_animal(d$animal_id[374])

# 375
out <- check_animal(d$animal_id[375])

# 376
out <- check_animal(d$animal_id[376])

# 377
out <- check_animal(d$animal_id[377])

# 378
out <- check_animal(d$animal_id[378])

# 379
out <- check_animal(d$animal_id[379])

# 380
out <- check_animal(d$animal_id[380])

# 381
out <- check_animal(d$animal_id[381])

# 382
out <- check_animal(d$animal_id[382])

# 383
out <- check_animal(d$animal_id[383])

# 384
out <- check_animal(d$animal_id[384])

# 385
out <- check_animal(d$animal_id[385])

# 386
out <- check_animal(d$animal_id[386])

# 387
out <- check_animal(d$animal_id[387])

# 388
out <- check_animal(d$animal_id[388])

# 389
out <- check_animal(d$animal_id[389])

# 390
out <- check_animal(d$animal_id[390])

# 391
out <- check_animal(d$animal_id[391])

# 392
out <- check_animal(d$animal_id[392])

# 393
out <- check_animal(d$animal_id[393])

# 394
out <- check_animal(d$animal_id[394])

# 395
out <- check_animal(d$animal_id[395])

# 396
out <- check_animal(d$animal_id[396])

# 397
out <- check_animal(d$animal_id[397])

# 398
out <- check_animal(d$animal_id[398])

# 399
out <- check_animal(d$animal_id[399])

# 400
out <- check_animal(d$animal_id[400])

# 401
out <- check_animal(d$animal_id[401])

# 402
out <- check_animal(d$animal_id[402])

# 403
out <- check_animal(d$animal_id[403])

# 404
out <- check_animal(d$animal_id[404])

# 405
out <- check_animal(d$animal_id[405])

# 406
out <- check_animal(d$animal_id[406])

# 407
out <- check_animal(d$animal_id[407])

# 408
out <- check_animal(d$animal_id[408])

# 409
out <- check_animal(d$animal_id[409])

# 410
out <- check_animal(d$animal_id[410])

# 411
out <- check_animal(d$animal_id[411])

# 412
out <- check_animal(d$animal_id[412])

# 413
out <- check_animal(d$animal_id[413])

# 414
out <- check_animal(d$animal_id[414])

# 415
out <- check_animal(d$animal_id[415])

# 416
out <- check_animal(d$animal_id[416])

# 417
out <- check_animal(d$animal_id[417])

# 418
out <- check_animal(d$animal_id[418])

# 419
out <- check_animal(d$animal_id[419])

# 420
out <- check_animal(d$animal_id[420])

# 421
out <- check_animal(d$animal_id[421])

# 422
out <- check_animal(d$animal_id[422])

# 423
out <- check_animal(d$animal_id[423])

# 424
out <- check_animal(d$animal_id[424])

# 425
out <- check_animal(d$animal_id[425])

# 426
out <- check_animal(d$animal_id[426])

# 427
out <- check_animal(d$animal_id[427])

# 428
out <- check_animal(d$animal_id[428])

# 429
out <- check_animal(d$animal_id[429])

# 430
out <- check_animal(d$animal_id[430])

# 431
out <- check_animal(d$animal_id[431])

# 432
out <- check_animal(d$animal_id[432])

# 433
out <- check_animal(d$animal_id[433])

# 434
out <- check_animal(d$animal_id[434])

# 435
out <- check_animal(d$animal_id[435])

# 436
out <- check_animal(d$animal_id[436])

# 437
out <- check_animal(d$animal_id[437])

# 438
out <- check_animal(d$animal_id[438])

# 439
out <- check_animal(d$animal_id[439])

# 440
out <- check_animal(d$animal_id[440])

# 441
out <- check_animal(d$animal_id[441])

# 442
out <- check_animal(d$animal_id[442])

# 443
out <- check_animal(d$animal_id[443])

# 444
out <- check_animal(d$animal_id[444])

# 445
out <- check_animal(d$animal_id[445])

# 446
out <- check_animal(d$animal_id[446])

# 447
out <- check_animal(d$animal_id[447])

# 448
out <- check_animal(d$animal_id[448])

# 449
out <- check_animal(d$animal_id[449])

# 450
out <- check_animal(d$animal_id[450])

# 451
out <- check_animal(d$animal_id[451])

# 452
out <- check_animal(d$animal_id[452])

# 453
out <- check_animal(d$animal_id[453])

# 454
out <- check_animal(d$animal_id[454])

# 455
out <- check_animal(d$animal_id[455])

# 456
out <- check_animal(d$animal_id[456])

# 457
out <- check_animal(d$animal_id[457])

# 458
out <- check_animal(d$animal_id[458])

# 459
out <- check_animal(d$animal_id[459])

# 460
out <- check_animal(d$animal_id[460])

# 461
out <- check_animal(d$animal_id[461])

# 462
out <- check_animal(d$animal_id[462])

# 463
out <- check_animal(d$animal_id[463])

# 464
out <- check_animal(d$animal_id[464])

# 465
out <- check_animal(d$animal_id[465])

# 466
out <- check_animal(d$animal_id[466])

# 467
out <- check_animal(d$animal_id[467])

# 468
out <- check_animal(d$animal_id[468])

# 469
out <- check_animal(d$animal_id[469])

# 470
out <- check_animal(d$animal_id[470])

# 471
out <- check_animal(d$animal_id[471])

# 472
out <- check_animal(d$animal_id[472])

# 473
out <- check_animal(d$animal_id[473])

# 474
out <- check_animal(d$animal_id[474])

# 475
out <- check_animal(d$animal_id[475])

# 476
out <- check_animal(d$animal_id[476])

# 477
out <- check_animal(d$animal_id[477])

# 478
out <- check_animal(d$animal_id[478])

# 479
out <- check_animal(d$animal_id[479])

# 480
out <- check_animal(d$animal_id[480])

# 481
out <- check_animal(d$animal_id[481])

# 482
out <- check_animal(d$animal_id[482])

# 483
out <- check_animal(d$animal_id[483])

# 484
out <- check_animal(d$animal_id[484])

# 485
out <- check_animal(d$animal_id[485])

# 486
out <- check_animal(d$animal_id[486])

# 487
out <- check_animal(d$animal_id[487])

# 488
out <- check_animal(d$animal_id[488])

# 489
out <- check_animal(d$animal_id[489])

# 490
out <- check_animal(d$animal_id[490])

# 491
out <- check_animal(d$animal_id[491])

# 492
out <- check_animal(d$animal_id[492])

# 493
out <- check_animal(d$animal_id[493])

# 494
out <- check_animal(d$animal_id[494])

# 495
out <- check_animal(d$animal_id[495])

# 496
out <- check_animal(d$animal_id[496])

# 497
out <- check_animal(d$animal_id[497])

# 498
out <- check_animal(d$animal_id[498])

# 499
out <- check_animal(d$animal_id[499])

# 500
out <- check_animal(d$animal_id[500])

# 501
out <- check_animal(d$animal_id[501])

# 502
out <- check_animal(d$animal_id[502])

# 503
out <- check_animal(d$animal_id[503])

# 504
out <- check_animal(d$animal_id[504])

# 505
out <- check_animal(d$animal_id[505])

# 506
out <- check_animal(d$animal_id[506])

# 507
out <- check_animal(d$animal_id[507])

# 508
out <- check_animal(d$animal_id[508])

# 509
out <- check_animal(d$animal_id[509])

# 510
out <- check_animal(d$animal_id[510])

# 511
out <- check_animal(d$animal_id[511])

# 512
out <- check_animal(d$animal_id[512])

# 513
out <- check_animal(d$animal_id[513])

# 514
out <- check_animal(d$animal_id[514])

# 515
out <- check_animal(d$animal_id[515])

# 516
out <- check_animal(d$animal_id[516])

# 517
out <- check_animal(d$animal_id[517])

# 518
out <- check_animal(d$animal_id[518])

# 519
out <- check_animal(d$animal_id[519])

# 520
out <- check_animal(d$animal_id[520])

# 521
out <- check_animal(d$animal_id[521])

# 522
out <- check_animal(d$animal_id[522])

# 523
out <- check_animal(d$animal_id[523])

# 524
out <- check_animal(d$animal_id[524])

# 525
out <- check_animal(d$animal_id[525])

# 526
out <- check_animal(d$animal_id[526])

# 527
out <- check_animal(d$animal_id[527])

# 528
out <- check_animal(d$animal_id[528])

# 529
out <- check_animal(d$animal_id[529])

# 530
out <- check_animal(d$animal_id[530])

# 531
out <- check_animal(d$animal_id[531])

# 532
out <- check_animal(d$animal_id[532])

# 533
out <- check_animal(d$animal_id[533])

# 534
out <- check_animal(d$animal_id[534])

# 535
out <- check_animal(d$animal_id[535])

# 536
out <- check_animal(d$animal_id[536])

# 537
out <- check_animal(d$animal_id[537])

# 538
out <- check_animal(d$animal_id[538])

# 539
out <- check_animal(d$animal_id[539])

# 540
out <- check_animal(d$animal_id[540])

# 541
out <- check_animal(d$animal_id[541])

# 542
out <- check_animal(d$animal_id[542])

# 543
out <- check_animal(d$animal_id[543])

# 544
out <- check_animal(d$animal_id[544])

# 545
out <- check_animal(d$animal_id[545])

# 546
out <- check_animal(d$animal_id[546])

# 547
out <- check_animal(d$animal_id[547])

# 548
out <- check_animal(d$animal_id[548])

# 549
out <- check_animal(d$animal_id[549])

# 550
out <- check_animal(d$animal_id[550])

# 551
out <- check_animal(d$animal_id[551])

# 552
out <- check_animal(d$animal_id[552])

# 553
out <- check_animal(d$animal_id[553])

# 554
out <- check_animal(d$animal_id[554])

# 555
out <- check_animal(d$animal_id[555])

# 556
out <- check_animal(d$animal_id[556])

# 557
out <- check_animal(d$animal_id[557])

# 558
out <- check_animal(d$animal_id[558])

# 559
out <- check_animal(d$animal_id[559])

# 560
out <- check_animal(d$animal_id[560])

# 561
out <- check_animal(d$animal_id[561])

# 562
out <- check_animal(d$animal_id[562])

# 563
out <- check_animal(d$animal_id[563])

# 564
out <- check_animal(d$animal_id[564])

# 565
out <- check_animal(d$animal_id[565])

# 566
out <- check_animal(d$animal_id[566])

# 567
out <- check_animal(d$animal_id[567])

# 568
out <- check_animal(d$animal_id[568])

# 569
out <- check_animal(d$animal_id[569])

# 570
out <- check_animal(d$animal_id[570])

# 571
out <- check_animal(d$animal_id[571])

# 572
out <- check_animal(d$animal_id[572])

# 573
out <- check_animal(d$animal_id[573])

# 574
out <- check_animal(d$animal_id[574])

# 575
out <- check_animal(d$animal_id[575])

# 576
out <- check_animal(d$animal_id[576])

# 577
out <- check_animal(d$animal_id[577])

# 578
out <- check_animal(d$animal_id[578])

# 579
out <- check_animal(d$animal_id[579])

# 580
out <- check_animal(d$animal_id[580])

# 581
out <- check_animal(d$animal_id[581])

# 582
out <- check_animal(d$animal_id[582])

# 583
out <- check_animal(d$animal_id[583])

# 584
out <- check_animal(d$animal_id[584])

# 585
out <- check_animal(d$animal_id[585])

# 586
out <- check_animal(d$animal_id[586])

# 587
out <- check_animal(d$animal_id[587])

# 588
out <- check_animal(d$animal_id[588])

# 589
out <- check_animal(d$animal_id[589])

# 590
out <- check_animal(d$animal_id[590])

# 591
out <- check_animal(d$animal_id[591])

# 592
out <- check_animal(d$animal_id[592])

# 593
out <- check_animal(d$animal_id[593])

# 594
out <- check_animal(d$animal_id[594])

# 595
out <- check_animal(d$animal_id[595])

# 596
out <- check_animal(d$animal_id[596])

# 597
out <- check_animal(d$animal_id[597])

# 598
out <- check_animal(d$animal_id[598])

# 599
out <- check_animal(d$animal_id[599])

# 600
out <- check_animal(d$animal_id[600])

# 601
out <- check_animal(d$animal_id[601])

# 602
out <- check_animal(d$animal_id[602])

# 603
out <- check_animal(d$animal_id[603])

# 604
out <- check_animal(d$animal_id[604])

# 605
out <- check_animal(d$animal_id[605])

# 606
out <- check_animal(d$animal_id[606])

# 607
out <- check_animal(d$animal_id[607])

# 608
out <- check_animal(d$animal_id[608])

# 609
out <- check_animal(d$animal_id[609])

# 610
out <- check_animal(d$animal_id[610])

# 611
out <- check_animal(d$animal_id[611])

# 612
out <- check_animal(d$animal_id[612])

# 613
out <- check_animal(d$animal_id[613])

# 614
out <- check_animal(d$animal_id[614])

# 615
out <- check_animal(d$animal_id[615])

# 616
out <- check_animal(d$animal_id[616])

# 617
out <- check_animal(d$animal_id[617])

# 618
out <- check_animal(d$animal_id[618])

# 619
out <- check_animal(d$animal_id[619])

# 620
out <- check_animal(d$animal_id[620])

# 621
out <- check_animal(d$animal_id[621])

# 622
out <- check_animal(d$animal_id[622])

# 623
out <- check_animal(d$animal_id[623])

# 624
out <- check_animal(d$animal_id[624])

# 625
out <- check_animal(d$animal_id[625])

# 626
out <- check_animal(d$animal_id[626])

# 627
out <- check_animal(d$animal_id[627])

# 628
out <- check_animal(d$animal_id[628])

# 629
out <- check_animal(d$animal_id[629])

# 630
out <- check_animal(d$animal_id[630])

# 631
out <- check_animal(d$animal_id[631])

# 632
out <- check_animal(d$animal_id[632])

# 633
out <- check_animal(d$animal_id[633])

# 634
out <- check_animal(d$animal_id[634])

# 635
out <- check_animal(d$animal_id[635])

# 636
out <- check_animal(d$animal_id[636])

# 637
out <- check_animal(d$animal_id[637])

# 638
out <- check_animal(d$animal_id[638])

# 639
out <- check_animal(d$animal_id[639])

# 640
out <- check_animal(d$animal_id[640])

# 641
out <- check_animal(d$animal_id[641])

# 642
out <- check_animal(d$animal_id[642])

# 643
out <- check_animal(d$animal_id[643])

# 644
out <- check_animal(d$animal_id[644])

# 645
out <- check_animal(d$animal_id[645])

# 646
out <- check_animal(d$animal_id[646])

# 647
out <- check_animal(d$animal_id[647])

# 648
out <- check_animal(d$animal_id[648])

# 649
out <- check_animal(d$animal_id[649])

# 650
out <- check_animal(d$animal_id[650])

# 651
out <- check_animal(d$animal_id[651])

# 652
out <- check_animal(d$animal_id[652])

# 653
out <- check_animal(d$animal_id[653])

# 654
out <- check_animal(d$animal_id[654])

# 655
out <- check_animal(d$animal_id[655])

# 656
out <- check_animal(d$animal_id[656])

# 657
out <- check_animal(d$animal_id[657])

# 658
out <- check_animal(d$animal_id[658])

# 659
out <- check_animal(d$animal_id[659])

# 660
out <- check_animal(d$animal_id[660])

# 661
out <- check_animal(d$animal_id[661])

# 662
out <- check_animal(d$animal_id[662])

# 663
out <- check_animal(d$animal_id[663])

# 664
out <- check_animal(d$animal_id[664])

# 665
out <- check_animal(d$animal_id[665])

# 666
out <- check_animal(d$animal_id[666])

# 667
out <- check_animal(d$animal_id[667])

# 668
out <- check_animal(d$animal_id[668])

# 669
out <- check_animal(d$animal_id[669])

# 670
out <- check_animal(d$animal_id[670])

# 671
out <- check_animal(d$animal_id[671])

# 672
out <- check_animal(d$animal_id[672])

# 673
out <- check_animal(d$animal_id[673])

# 674
out <- check_animal(d$animal_id[674])


# save the movebank dataset ----
saveRDS(d, paste0('movebank-data-', Sys.Date(), '.rds'))
