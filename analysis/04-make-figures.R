library('dplyr')
library('tidyr')
library('purrr')
library('ctmm')
library('mgcv')
library('ggplot2')
library('khroma')
library('gratia')
library('terra')
library('sf')
library('tidyterra')
library('rgeoboundaries')
library('rphylopic')
library('patchwork')
source('analysis/figures/default-ggplot-theme.R')

# TODO: update with latest dataset
d <- readRDS('data/gam-dataset-2026-10-04.rds') %>%
  filter(keep_for_gam, range_resident)

# find animals for panels 1C and 1D
plot(var_ndvi ~ mean_ndvi, d)

d_sub <- slice(d,
               c(which.min((mean_ndvi - 0.2)^2 + 20 * (var_ndvi - 0.005)^2),
                 which.min((mean_ndvi - 0.4)^2 + 20 * (var_ndvi - 0.02)^2),
                 which.min((mean_ndvi - 0.6)^2 + 20 * (var_ndvi - 0.04)^2)))

points(var_ndvi ~ mean_ndvi, d_sub, col = "red", pch = 19)

# TODO: update with latest dataset
d_akde <-  readRDS('data/akde-tibble-2026-09-13.rds') %>%
  filter(grepl(paste(c("LA12 from African elephants in Etosha National Park (data from Tsalyuk et al. 2018) (Loxodonta africana)",
                       "paf_pne (Capra ibex)", as.character(d_sub$animal)) %>%
                       gsub("\\(", "\\\\(", .) %>%
                       gsub("\\)", "\\\\)", .),
                     collapse = "|"
  ), animal))
layout(t(1:3))
plot(d_akde$tel[[1]], main = d_akde$animal[1])
plot(d_akde$tel[[2]][d_akde$tel[[1]]$y > -2e4, ], main = d_akde$animal[2])
plot(d_akde$tel[[3]], main = d_akde$animal[3])
plot(d_akde$tel[[4]], main = d_akde$animal[4])
plot(d_akde$tel[[5]], main = d_akde$animal[5])

extract_raster <- function(.tel, .raster) {
  .bbox <- as.sf(.tel) %>% st_bbox() %>% st_as_sfc() %>% st_as_sf()
  crop(.raster,
       st_transform(.bbox) %>% st_buffer(1e4)%>%
         st_transform(crs(.raster))) %>%
    project(., .tel@info$projection) %>%
    crop(.bbox %>% st_buffer(1000))
}

layout(matrix(1:10, ncol = 2, byrow = TRUE))
for(i in 1:5) {
  plot(extract_raster(d_akde$tel[[i]], .raster = r_mean_ndvi), main = "NDVI")
  plot(d_akde$tel[[i]], add = TRUE, error = FALSE, col = "black")
  
  plot(extract_raster(d_akde$tel[[i]], .raster = r_var_ndvi), main = "DENVar")
  plot(d_akde$tel[[i]], add = TRUE, error = FALSE, col = "black")
}

m_hr <- readRDS('models/hr-hgam-full.rds')

# fig 1 maps and HR ~ mass model ----
world <- geoboundaries(adm_lvl = 0, type = 'simplified') %>%
  st_geometry()
r_mean_ndvi <- rast('../ndvi-stochasticity/output/long-term-estimates.tif')[[1]]
r_mean_ndvi <- ifel(r_mean_ndvi < -1, -1, r_mean_ndvi)
r_var_ndvi <- rast('../ndvi-stochasticity/output/long-term-estimates.tif')[[2]]
r_var_ndvi <- ifel(r_var_ndvi < 0, 0, r_var_ndvi)
r_var_ndvi <- ifel(r_var_ndvi > 0.04, 0.04, r_var_ndvi)

# TODO: update with latest dataset
mean_locations <- readRDS("data/akde-tibble-2026-09-13.rds") %>%
  filter_out(dataset_name == "DeNicola et al. (2025) (Staten Island)") %>%
  # TODO: drop after updating AKDE dataset
  filter(species != "Chrysocyon brachyurus") %>%
  filter(keep_for_gam, !map_lgl(movement_model, is.null)) %>%
  filter(animal %in% d$animal) %>%
  # some CRSs removed due to data sensitivity
  mutate(mean = map(movement_model, function(.mm) {
    as.data.frame(.mm$mu) %>%
      st_as_sf(coords = c('x', 'y')) %>%
      st_set_crs(.mm@info$projection) %>%
      st_transform("EPSG:4326")
  })) %>%
  select(animal, species, dataset_name, mean) %>%
  unnest(mean) %>%
  st_as_sf() %>%
  aggregate(by = list(.$dataset_name), 
            function(x) x = x[1]) %>%
  st_centroid() %>%
  select(!Group.1)

# a mean NDVI
p_1a <- ggplot() +
  geom_sf(data = world, fill = "grey") +
  geom_spatraster(data = r_mean_ndvi) +
  geom_sf(data = world, fill = "transparent", color = "black") +
  geom_sf(data = mean_locations) +
  geom_sf(data = mean_locations, size = 0.25, color = "white") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn("NDVI", colors = ndvi_pal,
                       na.value = "white", limits = c(-1, 1)) +
  theme(legend.position = "right", legend.key.height = rel(2))
ggsave("figures/p-1a.png", p_1a, width = 8, height = 4, dpi = 300)

# b NDVI variance
p_1b <- ggplot() +
  geom_sf(data = world, fill = "grey") +
  geom_spatraster(data = r_var_ndvi) +
  geom_sf(data = world, fill = "transparent", color = "black") +
  geom_sf(data = mean_locations) +
  geom_sf(data = mean_locations, size = 0.25, color = "white") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_devon(name = "NDVI variance", na.value = "white", reverse = TRUE) +
  theme(legend.position = "right", legend.key.height = rel(2))
ggsave("figures/p-1b.png", p_1b, width = 8, height = 4, dpi = 300)

# d zoomed in versions of areas with NDVI
p_1c <-
  expand_grid(x = 1:10, y = 1:10, species = paste("Species", 1:4)) %>%
  mutate(NDVI = rbeta(n(), 2, 2)) %>%
  ggplot() +
  coord_equal() +
  facet_wrap(~ species, nrow = 1) +
  geom_raster(aes(x, y, fill = NDVI), show.legend = FALSE) +
  scale_x_continuous(NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, expand = c(0, 0)) +
  scale_fill_gradientn("NDVI", colors = ndvi_pal,
                       na.value = "white", limits = c(-1, 1))
ggsave("figures/p-1c.png", p_1d, width = 8, height = 2.5, dpi = 300)

# e zoomed in versions of areas with NDVI variance
p_1d <-
  expand_grid(x = 1:10, y = 1:10, species = paste("Species", 1:4)) %>%
  mutate(var_ndvi = rbeta(n(), 2, 2)) %>%
  ggplot() +
  coord_equal() +
  facet_wrap(~ species, nrow = 1) +
  geom_raster(aes(x, y, fill = var_ndvi), show.legend = FALSE) +
  scale_x_continuous(NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, expand = c(0, 0)) +
  scale_fill_devon(name = "NDVI variance", reverse = TRUE)
ggsave("figures/p-1d.png", p_1e, width = 8, height = 2.5, dpi = 300)

p_1 <- ((p_1a | p_1b) / (p_1d / p_1e)) +
  plot_annotation(tag_levels = 'A')

ggsave("figures/p-1.png", p_1, width = 8 * 2, height = 5 * 2, dpi = 300)
cli::cli_alert_info("Jittered masses to aid in visualization, but each species has a single mass.")

# fig 2 regression plots for HR model ----
# c HR ~ mass
null_preds <- tibble(mass_kg = exp(seq(log(0.01), log(1e4), length.out = 400)),
                     species = d$species[1]) %>%
  bind_cols(fitted_values(object = m_hr_0, data = ., scale = "link",
                          exclude = c("s(species)")) %>%
              filter(.parameter == "location"))

unique(d$species[which.min(d$mass_kg)])
unique(d$species[which.min(abs(d$mass_kg - 10))])
unique(d$species[which.max(d$mass_kg)])

p_1c <- ggplot(null_preds) +
  # Carollia castanea not available
  add_phylopic(name = "Desmodus rotundus", x = log(0.01), y = log(3),
               height = 1, verbose = TRUE) +
  # to pick specific phylopic
  add_phylopic(uuid = "ab6cfd4f-aef7-40fa-b5a5-1b79b7d112aa", x = log(10),
               y = log(1), height = 2, verbose = TRUE) +
  add_phylopic(name = "Loxodonta africana", x = log(5e3), y = log(50),
               height = 3, verbose = TRUE) +
  geom_ribbon(aes(log(mass_kg), ymin = .lower_ci, ymax = .upper_ci),
              alpha = 0.25) +
  geom_line(aes(log(mass_kg), .fitted)) +
  geom_jitter(aes(log(mass_kg), log(hr_95_km2)), d, alpha = 0.3,
              size = 0.5, height = 0, width = 0.01) +
  scale_x_continuous("Mass (kg)",
                     breaks = log(10^(-2:4)),
                     labels = \(x) as.character(exp(x))) +
  scale_y_continuous("Home-range size (km\U00B2)",
                     breaks = log(10^seq(-2, 6, by = 2)),
                     labels = \(x) as.character(exp(x)))
ggsave("figures/p-1c.png", p_1c, width = 8, height = 5, dpi = 300)


p_2a <-
  data_slice(m_hr, data = d, mean_ndvi = evenly(mean_ndvi)) %>%
  mutate(mass_kg = 1e4) %>%
  bind_cols(fitted_values(
    object = m_hr, data = ., terms = c("(Intercept)", "s(mean_ndvi)")) %>%
      filter(.parameter == "location")) %>%
  ggplot() +
  coord_cartesian(ylim = c(0, 30)) +
  geom_ribbon(aes(mean_ndvi, ymin = .lower_ci, ymax = .upper_ci),
              fill = pal[1], alpha = 0.3) +
  geom_line(aes(mean_ndvi, .fitted), color = pal[1], lwd = 1) +
  labs(x = "Mean NDVI", y = "Home-range size (km\U00B2)")

p_2b <-
  data_slice(m_hr, data = d, var_ndvi = evenly(var_ndvi)) %>%
  bind_cols(fitted_values(
    m_hr, data = ., terms = c("(Intercept)", "s(var_ndvi)")) %>%
      filter(.parameter == "location")) %>%
  ggplot() +
  coord_cartesian(ylim = c(0, 30)) +
  geom_ribbon(aes(var_ndvi, ymin = .lower_ci, ymax = .upper_ci),
              fill = pal[2], alpha = 0.3) +
  geom_line(aes(var_ndvi, .fitted), color = pal[2], lwd = 1) +
  labs(x = "NDVI variance", y = "Home-range size (km\U00B2)")

p_2c <-
  data_slice(m_hr, data = d,
             mean_ndvi = evenly(mean_ndvi),
             var_ndvi = evenly(var_ndvi)) %>%
  mutate(too_far = too_far(mean_ndvi, var_ndvi, d$mean_ndvi, d$var_ndvi,
                           dist = 0.3)) %>%
  filter(!too_far) %>%
  bind_cols(fitted_values(
    m_hr, data = ., terms = c("(Intercept)", "s(mean_ndvi)",
                              "s(var_ndvi)", "ti(mean_ndvi,var_ndvi)")) %>%
      filter(.parameter == "location")) %>%
  ggplot() +
  geom_raster(aes(mean_ndvi, var_ndvi, fill = log(.fitted))) +
  geom_contour(aes(mean_ndvi, var_ndvi, z = log(.fitted)),
               color = "black") +
  labs(x = "Mean NDVI", y = "NDVI variance") +
  scale_fill_gradient("Home-range size (km\U00B2)",
                      low = 'grey90', high = pal[3],
                      labels = exp, breaks = log(10^c(-1, 0, 1, 2, 3))) +
  theme(legend.position = "top", legend.key.width = rel(2))

p_2 <- ((p_2a | p_2b) / p_2c) +
  plot_annotation(tag_levels = 'A')
ggsave("figures/p-2.png", p_2, width = 8, height = 8, dpi = 300)
