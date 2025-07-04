# running a quick test with an early version of the work by Rekha Marcus:
# https://github.com/QuantitativeEcologyLab/stochasticity_parks
library('ctmm')
library('dplyr')
library('tidyr')
library('terra')
library('sf')
library('purrr')
library('ggplot2')
library('mgcv')
library('gratia')
source('analysis/figures/default-ggplot-theme.R')

# import rasters of mean, variance, and coefficient of variation in NDVI 
mu_r <- rast('data/marcus-rasters/mean_predNDVI_raster.nc')
s2_r <- rast('data/marcus-rasters/var_residuals_raster.nc')
cv_r <- rast('data/marcus-rasters/coeff_variation_raster.nc')

# extract values for each raster and average across GPS locations
d <-
  readRDS('../bc-mammals-temperature/models/movement-models-akdes-2024-06-06.rds') %>%
  left_join(read.csv('../bc-mammals-temperature/data/tracking-data/telemetry-metadata.csv') %>%
              select(animal, mass_g, keep_for_RSF),
            by = 'animal') %>%
  transmute(
    species = if_else(grepl('Rangifer', dataset_name),
                      dataset_name,
                      species),
    species = stringr::str_replace_all(species, '_', ' ') %>%
      stringr::str_replace_all('boreal', '(boreal)') %>%
      stringr::str_replace_all('southern mountain',
                               '(s. mountain)') %>%
      factor(),
    mass_g,
    keep_for_RSF,
    mu = map_dbl(tel, \(t) {
      data.frame(t) %>%
        select(longitude, latitude) %>%
        extract(x = mu_r, y = .) %>%
        pull(2) %>%
        mean(na.rm = TRUE)
    }),
    sigma2 = map_dbl(tel, \(t) {
      data.frame(t) %>%
        select(longitude, latitude) %>%
        extract(x = s2_r, y = .) %>%
        pull(2) %>%
        mean(na.rm = TRUE)
    }),
    cv = map_dbl(tel, \(t) {
      data.frame(t) %>%
        select(longitude, latitude) %>%
        extract(x = cv_r, y = .) %>%
        pull(2) %>%
        mean(na.rm = TRUE)
    }),
    hr_95 = map_dbl(akde, \(a) {
      summary(a, units = FALSE, level = 0.95)$CI['area (square meters)', 'est'] / 1e6
    }))

# plot unconditional relationships
d %>%
  pivot_longer(c(cv, mu, sigma2)) %>%
  ggplot(aes(value, hr_95)) +
  facet_grid(. ~ name, scales = 'free', switch = 'x') +
  geom_point(alpha = 0.2) +
  geom_smooth(aes(color = name, fill = name),
              alpha = 0.3, method = 'gam', formula = y ~ x) +
  xlab(NULL) +
  scale_y_log10() +
  scale_color_manual(values = c('#992266', '#F09D00', '#3388E1'),
                     aesthetics = c('color', 'fill')) +
  theme(strip.placement = 'outside', strip.background = element_blank(),
        legend.position = 'none')

# plot unconditional relationships by species
d %>%
  pivot_longer(c(cv, mu, sigma2)) %>%
  ggplot(aes(value, hr_95)) +
  facet_grid(species ~ name, scales = 'free', switch = 'x') +
  geom_point(alpha = 0.2) +
  geom_smooth(aes(color = name, fill = name),
              alpha = 0.3, method = 'gam', formula = y ~ x) +
  xlab(NULL) +
  scale_y_log10() +
  scale_color_manual(values = c('#992266', '#F09D00', '#3388E1'),
                     aesthetics = c('color', 'fill')) +
  theme(strip.placement = 'outside', strip.background = element_blank(),
        legend.position = 'none')

d %>%
  mutate(keep_for_RSF = factor(keep_for_RSF, levels = c('yes', 'no', 'check'))) %>%
  select(mass_g, mu, sigma2, cv, hr_95, keep_for_RSF) %>%
  pairs()

# fit a test model that does not account for body size
m <- bam(hr_95 ~
           # global effects
           s(mu, k = 5) +
           s(sigma2, k = 5) +
           # species-level effects
           s(species, bs = 're') +
           s(mu, species, bs = 'fs', k = 5) +
           s(sigma2, species, bs = 'fs', k = 5),
         family = Gamma(link = 'log'),
         data = d,
         subset = keep_for_RSF == 'yes',
         method = 'fREML',
         discrete = TRUE)
draw(m, residuals = TRUE, scales = 'fixed')

preds_mu <-
  tibble(species = 'new',
         mu = gratia:::seq_min_max(d$mu, n = 400),
         sigma2 = 0) %>%
  bind_cols(.,
            predict(m, newdata = ., type = 'link', se.fit = TRUE,
                    terms = c('(Intercept)', 's(mu)'),
                    discrete = FALSE) %>%
              as.data.frame()) %>%
  mutate(lwr_95 = exp(fit - se.fit * 1.96),
         est = exp(fit),
         upr_95 = exp(fit + se.fit * 1.96))

preds_sigma2 <-
  tibble(species = 'new',
         mu = 0,
         sigma2 = gratia:::seq_min_max(d$sigma2, n = 400)) %>%
  bind_cols(.,
            predict(m, newdata = ., type = 'link', se.fit = TRUE,
                    terms = c('(Intercept)', 's(sigma2)'),
                    discrete = FALSE) %>%
              as.data.frame()) %>%
  mutate(lwr_95 = exp(fit - se.fit * 1.96),
         est = exp(fit),
         upr_95 = exp(fit + se.fit * 1.96))

cowplot::plot_grid(
  ggplot() +
    geom_point(aes(mu, hr_95), d, color = pal[3], alpha = 0.3) +
    geom_ribbon(aes(mu, ymin = lwr_95, ymax = upr_95), preds_mu,
                fill = pal[1], alpha = 0.2) +
    geom_line(aes(mu, est), preds_mu, color = pal[1], linewidth = 1) +
    geom_line() +
    xlab('Mean NDVI') +
    scale_y_log10(expression(bold('95%'~HR~(log['10']~scale)))),
  ggplot() +
    geom_point(aes(sigma2, hr_95), d, color = pal[3], alpha = 0.3) +
    geom_ribbon(aes(sigma2, ymin = lwr_95, ymax = upr_95), preds_sigma2,
                fill = pal[2], alpha = 0.2) +
    geom_line(aes(sigma2, est), preds_sigma2, color = pal[2], linewidth = 1) +
    geom_line() +
    xlab('Variance in NDVI') +
    scale_y_log10(expression(bold('95%'~HR~(log['10']~scale)))),
  ncol = 1)

ggsave('figures/bc-test.png', height = 12, width = 12, dpi = 600, bg = 'white')
