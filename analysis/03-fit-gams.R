# include human pressure around the HR?
# could restrict to HFI values below some amount
library('dplyr')
library('tidyr')
library('mgcv')
library('ggplot2')
library('khroma')
library('gratia')
source('analysis/figures/default-ggplot-theme.R')

d <- readRDS('data/gam-dataset-2026-09-22.rds') %>%
  filter(mean_ndvi != 0) # drop individuals outside the current rasters

if (! any(colnames(d) == 'speed')) {
  d$speed <- rgamma(nrow(d), shape = 10, rate = 0.2)
  warning('Missing speed! Filled value with random Gamma values.')
}

if (! any(colnames(d) == 'diff')) {
  d$diffusion <- rgamma(nrow(d), shape = 10, rate = 0.2)
  warning('Missing diffusion! Filled value with random Gamma values.')
}

d %>%
  summarize(mean_ndvi = mean(mean_ndvi),
            var_ndvi = mean(var_ndvi),
            .by = species)

d %>%
  slice(1, .by = species) %>%
  select(species, mass_kg)

# HR size and mass
ggplot(d, aes(log(mass_kg), log(hr_95_km2))) +
  geom_point(alpha = 0.4) +
  labs(x = 'Mass (kg)',
       y = expression(paste(bold('95% HR size'~'(km'^'2'), bold(')')))) +
  theme(legend.position = 'top')

# plot unconditional relationships
d %>%
  pivot_longer(c(mean_ndvi, var_ndvi)) %>%
  ggplot(aes(value, hr_95_km2)) +
  facet_grid(. ~ name, scales = 'free', switch = 'x') +
  geom_point(alpha = 0.2) +
  geom_smooth(aes(color = name, fill = name),
              alpha = 0.3, method = 'gam', formula = y ~ x) +
  xlab(NULL) +
  scale_y_log10() +
  scale_color_manual(values = c('#F09D00', '#3388E1'),
                     aesthetics = c('color', 'fill')) +
  theme(strip.placement = 'outside', strip.background = element_blank(),
        legend.position = 'none')

# plot unconditional relationships by species
d %>%
  pivot_longer(c(mean_ndvi, var_ndvi)) %>%
  ggplot(aes(value, hr_95_km2)) +
  facet_grid(~ name, scales = 'free', switch = 'x') +
  geom_point(alpha = 0.2) +
  geom_smooth(aes(color = name, fill = name, group = species),
              alpha = 0.3, method = 'gam', formula = y ~ x) +
  xlab(NULL) +
  scale_y_log10() +
  scale_color_manual(values = c('#F09D00', '#3388E1'),
                     aesthetics = c('color', 'fill')) +
  theme(strip.placement = 'outside', strip.background = element_blank(),
        legend.position = 'none')

d %>%
  select(mass_kg, mean_ndvi, var_ndvi, hr_95_km2) %>%
  pairs()

# to get exponent for historical context
m_hr_mass_only <- gam(list(hr_95_km2 ~ log(mass_kg) + s(species, bs = 're'),
                           ~ 1),
                      family = gammals(),
                      data = d,
                      subset = keep_for_gam)
summary(m_hr_mass_only)
exp(coef(m_hr_mass_only)[1:2])

# add NDVI variables to account for deviations
m_hr_0 <- gam(list(hr_95_km2 ~
                     # global effects
                     s(mean_ndvi, k = 5) +
                     # s(mass_kg, k = 5) +
                     # species-level effects
                     s(species, bs = 're') +
                     s(mean_ndvi, species, bs = 'fs', k = 5),
                   ~ s(species, bs = 're') +
                     s(mean_ndvi, k = 5)),
              family = gammals(),
              data = d,
              subset = keep_for_gam &
                # in too urban of an area
                dataset_name != 'DeNicola et al. (2025) (Staten Island)',
              method = 'REML')

m_hr <- gam(list(hr_95_km2 ~
                   # global effects
                   s(mean_ndvi, k = 5) +
                   s(var_ndvi, k = 5) +
                   ti(mean_ndvi, var_ndvi, k = 3) +
                   # s(mass_kg, k = 5) +
                   # species-level effects
                   s(species, bs = 're') +
                   s(mean_ndvi, species, bs = 'fs', k = 5) +
                   s(var_ndvi, species, bs = 'fs', k = 5),
                 ~ s(species, bs = 're') +
                   te(mean_ndvi, var_ndvi, k = 5)),
            family = gammals(),
            data = d,
            subset = keep_for_gam &
              # in too urban of an area
              dataset_name != 'DeNicola et al. (2025) (Staten Island)',
            method = 'REML')
draw(m_hr, dist = 1)
ggsave('figures/hr-hgam-terms.png', width = 8, height = 8, units = 'in',
       dpi = 300, bg = 'white', scale = 2)

appraise(m_hr, point_alpha = 0.2)
ggsave('figures/hr-hgam-diagnostics.png', width = 8, height = 8, units = 'in',
       dpi = 300, bg = 'white')

# generalized likelihood test for model comparison
anova(m_hr_mass_only, m_hr, test = 'Chisq')
anova(m_hr_0, m_hr, test = 'Chisq')

# change in deviance explained
((summary(m_hr$null.deviance)) -
    (m_hr_0$deviance / m_hr_0$null.deviance)) * 100

if (any(! grepl('mass_kg', smooths(m_hr)))) {
  stop("The model does not have a term of mass and is thus not comparable to m_hr_0!")
} else {
  AIC(m_hr_0, m_hr) %>% mutate(AIC = AIC - min(AIC)) %>% arrange(AIC)
}

# is there evidence of anthropogenic pressure?
d %>%
  filter(dataset_name != 'DeNicola et al. (2025) (Staten Island)') %>%
  ggplot(aes(hfi, resid(m_hr))) +
  geom_point(alpha = 0.5) +
  geom_smooth()

# perform the same analysis but for speed ----
m_speed_0 <- gam(speed ~ log(mass_kg) + s(species, bs = 're'),
                 family = Gamma(link = 'log'),
                 data = d,
                 subset = keep_for_gam)

m_speed <- bam(speed ~
                 # global effects
                 s(mean_ndvi, k = 5) +
                 s(var_ndvi, k = 5) +
                 ti(mean_ndvi, var_ndvi, k = 5) +
                 # s(mass_kg, k = 5) +
                 # species-level effects
                 s(species, bs = 're') +
                 s(mean_ndvi, species, bs = 'fs', k = 5) +
                 s(var_ndvi, species, bs = 'fs', k = 5),
               family = Gamma(link = 'log'),
               data = d,
               method = 'fREML',
               discrete = TRUE)
draw(m_speed, residuals = TRUE)
appraise(m_speed, point_alpha = 0.2)
summary(m_speed)
AIC(m_speed_0, m_speed) %>% mutate(AIC = AIC - min(AIC)) %>% arrange(AIC)

# perform the same analysis but for diffusion
m_diff_0 <- gam(diffusion ~ log(mass_kg) + s(species, bs = 're'),
                family = Gamma(link = 'log'),
                data = d,
                subset = keep_for_gam)

m_diff <- bam(diffusion ~
                # global effects
                s(mean_ndvi, k = 5) +
                s(var_ndvi, k = 5) +
                ti(mean_ndvi, var_ndvi, k = 5) +
                # s(mass_kg, k = 5) +
                # species-level effects
                s(species, bs = 're') +
                s(mean_ndvi, species, bs = 'fs', k = 5) +
                s(var_ndvi, species, bs = 'fs', k = 5),
              family = Gamma(link = 'log'),
              data = d,
              method = 'fREML',
              discrete = TRUE)
draw(m_diff, residuals = TRUE)
appraise(m_diff, point_alpha = 0.2)
summary(m_diff)
AIC(m_diff_0, m_diff) %>% mutate(AIC = AIC - min(AIC)) %>% arrange(AIC)

# save models ----
saveRDS(m_hr, file = 'models/hr-hgam.rds')
saveRDS(m_speed, file = 'models/diffusion-hgam.rds')
saveRDS(m_diff, file = 'models/speed-hgam.rds')

saveRDS(m_hr_0, file = 'models/hr-null-hgam.rds')
saveRDS(m_speed_0, file = 'models/diffusion-null-hgam.rds')
saveRDS(m_diff_0, file = 'models/speed-null-hgam.rds')
