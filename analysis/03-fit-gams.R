# include human pressure around the HR?
# could restrict to HFI values below some amount
library('dplyr')
library('tidyr')
library('mgcv')
library('ggplot2')
library('khroma')
library('gratia')
source('analysis/figures/default-ggplot-theme.R')
source('functions/deviance.R')

d <- readRDS('data/gam-dataset-2026-09-26.rds') %>%
  filter(mean_ndvi != 0) # drop individuals outside the current rasters

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

# mean ndvi and DENVar are correlated
d %>%
  select(mass_kg, mean_ndvi, var_ndvi, hr_95_km2) %>%
  pairs()

# to get exponent for historical context
m_hr_0 <- gam(list(hr_95_km2 ~ log(mass_kg) + s(species, bs = 're'),
                           ~ 1),
                      family = gammals(),
                      data = d,
                      subset = keep_for_gam & range_resident)
summary(m_hr_0)
exp(coef(m_hr_0)[1:2])

# add NDVI variables to account for deviations
m_hr_ndvi <- gam(list(hr_95_km2 ~
                        # global effects
                        s(mean_ndvi, k = 5) +
                        log(mass_kg) +
                        # species-level effects
                        s(species, bs = 're') +
                        s(mean_ndvi, species, bs = 'fs', k = 5),
                      ~ s(species, bs = 're') +
                        s(mean_ndvi, k = 5)),
                 family = gammals(),
                 data = d,
                 subset = keep_for_gam & range_resident,
                 method = 'REML')
summary(m_hr_ndvi)

m_hr <- gam(list(hr_95_km2 ~
                   # global effects
                   s(mean_ndvi, k = 5) +
                   s(var_ndvi, k = 5) +
                   ti(mean_ndvi, var_ndvi, k = 3) +
                   log(mass_kg) +
                   # species-level effects
                   s(species, bs = 're') +
                   s(mean_ndvi, species, bs = 'fs', k = 5) +
                   s(var_ndvi, species, bs = 'fs', k = 5),
                 ~ s(species, bs = 're') +
                   te(mean_ndvi, var_ndvi, k = 5)),
            family = gammals(),
            data = d,
            subset = keep_for_gam & range_resident,
            method = 'REML')
draw(m_hr, dist = 1) +
  scale_x_continuous(transform = "log", breaks = 10^seq(-2, 3, by = 2),
                     labels = as.character)
ggsave('figures/hr-hgam-terms.png', width = 8, height = 8, units = 'in',
       dpi = 300, bg = 'white', scale = 2)

appraise(m_hr, point_alpha = 0.2)
ggsave('figures/hr-hgam-diagnostics.png', width = 8, height = 8, units = 'in',
       dpi = 300, bg = 'white')

summary(m_hr_ndvi)

# generalized likelihood test for model comparison
anova(m_hr_0, m_hr, test = 'Chisq')
anova(m_hr_ndvi, m_hr, test = 'Chisq')

# change in deviance relative to model with mean NDVI
perc_dev_change(m_hr, m_hr_ndvi)
perc_dev_prev_unexpl(m_hr, m_hr_ndvi)

# relative to mass-only model
perc_dev_change(m_hr, m_hr_0)
perc_dev_prev_unexpl(m_hr, m_hr_0)

compare_aic(m_hr, m_hr_ndvi, m_hr_0)

# perform the same analysis but for speed ----
m_speed_0 <- gam(speed_km_day ~ log(mass_kg) + s(species, bs = 're'),
                 family = Gamma(link = 'log'),
                 data = d,
                 subset = is.finite(speed_km_day))
summary(m_speed_0)

m_speed_ndvi <- gam(list(speed_km_day ~
                           # global effects
                           s(mean_ndvi, k = 5) +
                           log(mass_kg) +
                           # species-level effects
                           s(species, bs = 're') +
                           s(mean_ndvi, species, bs = 'fs', k = 5),
                         ~ s(species, bs = 're') +
                           s(mean_ndvi, k = 5)),
                    family = gammals(),
                    data = d,
                    subset = is.finite(speed_km_day),
                    method = 'REML')
summary(m_speed_ndvi)
perc_dev_change(m_speed_ndvi, m_speed_0)
perc_dev_prev_unexpl(m_speed_ndvi, m_speed_0)

m_speed <- gam(list(speed_km_day ~
                      # global effects
                      s(mean_ndvi, k = 5) +
                      s(var_ndvi, k = 5) +
                      ti(mean_ndvi, var_ndvi, k = 3) +
                      log(mass_kg) +
                      # species-level effects
                      s(species, bs = 're') +
                      s(mean_ndvi, species, bs = 'fs', k = 5) +
                      s(var_ndvi, species, bs = 'fs', k = 5),
                    ~ s(species, bs = 're') +
                      te(mean_ndvi, var_ndvi, k = 5)),
               family = gammals(),
               data = d,
               subset = is.finite(speed_km_day),
               method = 'REML')
summary(m_speed)
perc_dev_change(m_speed, m_speed_ndvi)
perc_dev_prev_unexpl(m_speed, m_speed_ndvi)

perc_dev_change(m_speed, m_speed_0)
perc_dev_prev_unexpl(m_speed, m_speed_0)

draw(m_speed, dist = 1)
appraise(m_speed, point_alpha = 0.2)
compare_aic(m_speed, m_speed_ndvi, m_speed_0)

# perform the same analysis but for diffusion
m_diff_0 <- gam(diffusion_km2_day ~ log(mass_kg) + s(species, bs = 're'),
                family = Gamma(link = 'log'),
                data = d)

m_diff_ndvi <- gam(list(diffusion_km2_day ~
                        # global effects
                        s(mean_ndvi, k = 5) +
                        log(mass_kg) +
                        # species-level effects
                        s(species, bs = 're') +
                        s(mean_ndvi, species, bs = 'fs', k = 5),
                      ~ s(species, bs = 're') +
                        s(mean_ndvi, k = 5)),
                 family = gammals(),
                 data = d,
                 subset = keep_for_gam & range_resident,
                 method = 'REML')
summary(m_diff_ndvi)

m_diff <- gam(list(diffusion_km2_day ~
                     # global effects
                     s(mean_ndvi, k = 5) +
                     s(var_ndvi, k = 5) +
                     ti(mean_ndvi, var_ndvi, k = 3) +
                     log(mass_kg) +
                     # species-level effects
                     s(species, bs = 're') +
                     s(mean_ndvi, species, bs = 'fs', k = 5) +
                     s(var_ndvi, species, bs = 'fs', k = 5),
                   ~ s(species, bs = 're') +
                     te(mean_ndvi, var_ndvi, k = 5)),
              family = gammals(),
              data = d,
              method = 'REML')
draw(m_diff, dist = 1)
appraise(m_diff, point_alpha = 0.2)
summary(m_diff)
compare_aic(m_diff, m_diff_ndvi, m_diff_0)

# save models ----
saveRDS(m_hr_0, file = 'models/hr-hgam-mass-only.rds')
saveRDS(m_speed_0, file = 'models/diffusion-hgam-mass-only.rds')
saveRDS(m_diff_0, file = 'models/speed-hgam-mass-only.rds')

saveRDS(m_hr_ndvi, file = 'models/hr-ndvi-only-hgam.rds')
saveRDS(m_speed_ndvi, file = 'models/diffusion-ndvi-only-hgam.rds')
saveRDS(m_diff_ndvi, file = 'models/speed-ndvi-only-hgam.rds')

saveRDS(m_hr, file = 'models/hr-full-hgam.rds')
saveRDS(m_speed, file = 'models/diffusion-full-hgam.rds')
saveRDS(m_diff, file = 'models/speed-full-hgam.rds')
