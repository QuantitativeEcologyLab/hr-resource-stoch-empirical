# include human pressure around the HR?
# could restrict to HFI values below some amount
library('dplyr')
library('mgcv')
library('ggplot2')
library('khroma')
library('gratia')
source('analysis/figures/default-ggplot-theme.R')

d <- readRDS('data/gam-dataset-2026-02-24.rds') %>%
  filter(keep_for_gam) %>%
  mutate(species = factor(species))

d %>%
  slice(1, .by = species) %>%
  select(species, mass_kg)

ggplot(d, aes(mass_kg, hr_95, color = species)) +
  geom_point(alpha = 0.4) +
  scale_color_discreterainbow() +
  theme(legend.position = 'top')

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