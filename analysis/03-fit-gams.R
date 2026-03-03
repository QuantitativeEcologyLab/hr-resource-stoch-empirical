# include human pressure around the HR?
# could restrict to HFI values below some amount
library('dplyr')
library('tidyr')
library('mgcv')
library('ggplot2')
library('khroma')
library('gratia')
source('analysis/figures/default-ggplot-theme.R')

d <- readRDS('data/gam-dataset-2026-03-03.rds') %>%
  filter(keep_for_gam) %>%
  filter(mean_ndvi != 0) %>% # drop individuals outside the current rasters
  mutate(species = factor(species))

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
  facet_wrap(. ~ species + name, scales = 'free', ncol = 4) +
  geom_point(alpha = 0.2) +
  geom_smooth(aes(color = name, fill = name),
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
m_0 <- gam(hr_95 ~ log(mass_kg),
           family = gaussian(link = 'log'),
           data = d,
           subset = keep_for_gam)
summary(m_0)
exp(coef(m_0))

# add NDVI variables to account for deviations
# compare deviance explained
m <- gam(hr_95 ~
           # global effects
           s(log1p(mean_ndvi), k = 3) +
           s(log1p(var_ndvi), k = 3) +
           #ti(mean_ndvi, var_ndvi, k = 5) +
           # s(mass_kg, k = 5) +
           # species-level effects
           s(species, bs = 're') +
           s(log1p(mean_ndvi), species, bs = 'fs', k = 5) +
           s(log1p(var_ndvi), species, bs = 'fs', k = 5),
         family = gaussian(link = 'log'),
         data = d,
         method = 'REML')
draw(m)

# regress residuals against HFI

# TODO: perform the same analysis but for diffusion and/or speed?