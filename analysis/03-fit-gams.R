# include human pressure around the HR?
# could restrict to HFI values below some amount
library('dplyr')
library('tidyr')
library('mgcv')
library('ggplot2')
library('khroma')
library('gratia')
source('analysis/figures/default-ggplot-theme.R')

d <- readRDS('data/gam-dataset-2026-03-24.rds') %>%
  filter(keep_for_gam) %>%
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
ggplot(d, aes(mass_kg, hr_95, color = species)) +
  geom_point(alpha = 0.4) +
  labs(x = 'Mass (kg)',
       y = expression(paste(bold('95% HR size'~'(km'^'2'), bold(')')))) +
  scale_color_discreterainbow(name = 'Species') +
  theme(legend.position = 'top')

# plot unconditional relationships
d %>%
  pivot_longer(c(mean_ndvi, var_ndvi)) %>%
  ggplot(aes(value, hr_95)) +
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
  ggplot(aes(value, hr_95)) +
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
  select(mass_kg, mean_ndvi, var_ndvi, hr_95) %>%
  pairs()

# to get exponent for historical context
m_hr_0 <- gam(hr_95 ~ log(mass_kg) + s(species, bs = 're'),
              family = Gamma(link = 'log'),
              data = d,
              subset = keep_for_gam)
summary(m_hr_0)
exp(coef(m_hr_0))

warning("Can't include mass yet because it's strongly correlated with E(NDVI)")
ggplot(d, aes(mass_kg, mean_ndvi)) +
  geom_point(aes(color = species)) +
  geom_smooth(method = 'lm') +
  scale_color_discreterainbow()

# add NDVI variables to account for deviations
# compare deviance explained
m_hr <- gam(list(hr_95 ~
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
            # in too urban of an area
            subset = dataset_name != 'DeNicola et al. (2025) (Staten Island)',
            method = 'REML')
draw(m_hr, dist = 1)
ggsave('figures/hr-model-terms.png', width = 8, height = 8, units = 'in',
       dpi = 300, bg = 'white')

appraise(m_hr, point_alpha = 0.2)
ggsave('figures/hr-model-diagnostics.png', width = 8, height = 8, units = 'in',
       dpi = 300, bg = 'white')

summary(m_hr)

