library('dplyr')
library('mgcv')
library('ggplot2')
library('khroma')
library('gratia')
library('terra')
library('sf')
library('tidyterra')
library('rgeoboundaries')
source('analysis/figures/default-ggplot-theme.R')

d <- readRDS('data/gam-dataset-2026-09-26.rds') %>%
  filter(mean_ndvi != 0) # drop individuals outside the rasters

world <- geoboundaries(adm_lvl = 0, type = 'simplified') %>%
  st_geometry()
plot(world)
r_mean_ndvi <- rast('../ndvi-stochasticity/output/long-term-estimates.tif')[[1]]
r_mean_ndvi <- ifel(r_mean_ndvi < -1, -1, r_mean_ndvi)
r_var_ndvi <- rast('../ndvi-stochasticity/output/long-term-estimates.tif')[[2]]
r_var_ndvi <- ifel(r_var_ndvi < 0, 0, r_var_ndvi)
r_var_ndvi <- ifel(r_var_ndvi > 0.04, 0.04, r_var_ndvi)

m_hr <- readRDS('models/hr-model.rds')

# Fig. 1: like figure 1 of Tucker et al. 2018:
# a mean NDVI
p_1a <- ggplot() +
  geom_sf(data = world) +
  geom_spatraster(data = r_mean_ndvi) +
  geom_sf(data = world, fill = "transparent", color = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn("Mean NDVI", colors = ndvi_pal, na.value = "white",
                       limits = c(-1, 1)) +
  theme(legend.position = "bottom", legend.key.width = rel(2))
p_1a

# b DENVar
p_1b <- ggplot() +
  geom_sf(data = world) +
  geom_spatraster(data = r_var_ndvi) +
  geom_sf(data = world, fill = "transparent", color = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_devon(name = "DENVar", na.value = "white") +
  theme(legend.position = "bottom", legend.key.width = rel(2))
p_1b

# c boxplots of HR size
p_1c <- d %>%
  mutate(median_hr_95_km2 = median(hr_95_km2), .by = species) %>%
  arrange(median_hr_95_km2) %>%
  mutate(species = factor(species, levels = unique(species))) %>%
  ggplot() +
  geom_boxplot(aes(log(hr_95_km2), species), fill = "grey",
               alpha = 0.5) +
  labs(x = "Home range size (km\U00B2, log scale)", y = NULL) +
  theme(axis.text.y = element_text(face = "bold.italic"))
p_1c

# d zoomed in versions of areas (could also be figs from below)

# or Fig 1 of https://www.pnas.org/doi/10.1073/pnas.2401042122
# - same species, but in different areas with different levels of NDVI and DENVar

# fig 2 regression plots for HR model

# fig 3: same as for fig 2, but for speed, diffusion, or ballistic lengthscale
