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
library('maptiles')
source('analysis/figures/default-ggplot-theme.R')

# robinson objects
robinson_crs <- 'ESRI:54030'

bounds <-
  dplyr::tibble(x = c(-180:180, rep(180, 181), 180:-180, rep(-180, 181)),
                y = c(rep(-90, 361), -90:90, rep(90, 361), 90:-90)) %>%
  st_as_sf(coords = c('x', 'y'), crs = 'EPSG:4326') %>%
  dplyr::summarise(geometry = st_combine(geometry)) %>%
  st_cast('POLYGON') %>%
  st_transform(robinson_crs)

# TODO: update with latest dataset
d <- readRDS('data/gam-dataset-2026-10-04.rds') %>%
  filter(keep_for_gam, range_resident)
m_hr <- readRDS('models/hr-hgam-full.rds')

# fig 1 maps and HR ~ mass model ----
world <- geoboundaries(adm_lvl = 0, type = 'simplified') %>%
  # filter(shapeName != "Antarctica") %>%
  st_geometry() %>%
  st_transform(robinson_crs)
r_mean_ndvi <- rast('../ndvi-stochasticity/output/long-term-estimates.tif')[[1]]
r_mean_ndvi <- ifel(r_mean_ndvi < -1, -1, r_mean_ndvi)
r_mean_ndvi <- project(r_mean_ndvi, robinson_crs)
r_var_ndvi <- rast('../ndvi-stochasticity/output/long-term-estimates.tif')[[2]]
r_var_ndvi <- ifel(r_var_ndvi < 0, 0, r_var_ndvi)
r_var_ndvi <- ifel(r_var_ndvi > 0.04, 0.04, r_var_ndvi)
r_var_ndvi <- project(r_var_ndvi, robinson_crs)

# row of examples
plot(var_ndvi ~ mean_ndvi, d)
d_sub <- slice(d,
               c(which.min((mean_ndvi - 0.2)^2 + 20 * (var_ndvi - 0.005)^2),
                 which.min((mean_ndvi - 0.4)^2 + 20 * (var_ndvi - 0.02)^2),
                 which.min((mean_ndvi - 0.6)^2 + 20 * (var_ndvi - 0.04)^2)))
points(var_ndvi ~ mean_ndvi, d_sub, col = "red", pch = 19)

# TODO: update with latest dataset
d_akde <-  readRDS('data/akde-tibble-2026-09-13.rds')

d_akde_sub <- d_akde %>%
  filter(grepl(paste(as.character(d_sub$animal) %>%
                       gsub("\\(", "\\\\(", .) %>%
                       gsub("\\)", "\\\\)", .),
                     collapse = "|"
  ), animal))
plot(d_akde_sub$tel[[1]], main = d_akde_sub$animal[1])
plot(d_akde_sub$tel[[2]], main = d_akde_sub$animal[2])
plot(d_akde_sub$tel[[3]], main = d_akde_sub$animal[3])

transmute(d_sub,
          animal,
          species,
          mean_ndvi,
          var_ndvi,
          latitude = map_dbl(d_akde_sub$tel, \(.t) mean(.t$latitude)),
          longitude = map_dbl(d_akde_sub$tel, \(.t) mean(.t$longitude)))

mean_locations <- d_akde %>%
  # TODO: drop after updating AKDE dataset
  filter(species != "Chrysocyon brachyurus") %>%
  filter(keep_for_gam, !map_lgl(movement_model, is.null)) %>%
  filter(animal %in% d$animal) %>%
  # some CRSs removed due to data sensitivity
  mutate(mean = map(movement_model, function(.mm) {
    as.data.frame(.mm$mu) %>%
      st_as_sf(coords = c('x', 'y')) %>%
      st_set_crs(.mm@info$projection) %>%
      st_transform(robinson_crs)
  })) %>%
  select(animal, species, dataset_name, mean) %>%
  unnest(mean) %>%
  st_as_sf() %>%
  aggregate(by = list(.$dataset_name), 
            function(x) x = x[1]) %>%
  st_centroid() %>%
  select(!Group.1) %>%
  mutate(has_panel = dataset_name %in% d_sub$dataset_name)

# a mean NDVI
p_1a <- ggplot() +
  geom_sf(data = bounds, fill = "white") +
  geom_sf(data = world, fill = "grey") +
  geom_spatraster(data = r_mean_ndvi) +
  geom_sf(data = world, fill = "transparent", color = "black") +
  geom_sf(data = mean_locations) +
  geom_sf(data = mean_locations, size = 0.25, color = "white") +
  scale_fill_gradientn("Mean NDVI", colors = ndvi_pal,
                       na.value = "transparent", limits = c(-1, 1)) +
  theme_void() +
  theme(text = element_text(face = "bold"),
        legend.position = "top", legend.key.height = rel(0.5),
        legend.key.width = rel(3), legend.title = element_text(vjust = 1),
        legend.key.justification = "center",
        legend.box.margin = margin(t = 10))
ggsave("figures/p-1a.png", p_1a, width = 7, height = 4, dpi = 300,
       bg = "white")

# b NDVI variance
p_1b <- ggplot() +
  geom_sf(data = bounds, fill = "white") +
  geom_sf(data = world, fill = "grey") +
  geom_spatraster(data = r_var_ndvi) +
  geom_sf(data = world, fill = "transparent", color = "black") +
  geom_sf(data = mean_locations) +
  geom_sf(data = mean_locations, size = 0.25, color = "white") +
  scale_fill_devon(name = "NDVI variance", na.value = "transparent",
                   reverse = TRUE) +
  theme_void() +
  theme(text = element_text(face = "bold"),
        legend.position = "top", legend.key.height = rel(0.5),
        legend.key.width = rel(3), legend.title = element_text(vjust = 1),
        legend.key.justification = "center",
        legend.box.margin = margin(t = 10))
ggsave("figures/p-1b.png", p_1b, width = 7, height = 4, dpi = 300,
       bg = "white")

# bottom panels of examples telemetries ----
plot_tel <- function(tel, zoom, buff, decimals = 3) {
  tel_sf <- as.sf(tel) %>% st_transform("EPSG:4326")
  tel_coords <- st_coordinates(tel_sf)
  box <- st_as_sf(st_as_sfc(st_bbox(tel_sf)))
  # box_coords <- st_coordinates(box) %>% as.data.frame()
  # buff <- max(c(abs(diff(range(box_coords$X))),
  #               abs(diff(range(box_coords$Y)))))
  box <- st_buffer(box, dist = buff)
  
  basemap <- get_tiles(x = box, provider = "Esri.WorldImagery",
                       zoom = zoom)
  basemap <- crop(basemap, box)
  
  box_box <- st_bbox(box)
  
  get_breaks <- function(axis) {
    seq(box_box[paste0(axis, "min")] * 0.75 +
          box_box[paste0(axis, "max")] * 0.25,
        box_box[paste0(axis, "min")] * 0.25 +
          box_box[paste0(axis, "max")] * 0.75,
        length.out = 2) %>%
      round(decimals)
  }
  
  ggplot() +
    tidyterra::geom_spatraster_rgb(data = basemap, mask_projection = TRUE) +
    # geom_path(aes(X, Y), tel_coords, lwd = 1, color = "white") +
    geom_sf(data = tel_sf, size = 1, color = "white") +
    # geom_path(aes(X, Y), tel_coords) +
    geom_sf(data = tel_sf, size = 0.25) +
    scale_x_continuous(NULL, expand = c(0, 0), get_breaks("x")) +
    scale_y_continuous(NULL, expand = c(0, 0), get_breaks("y")) +
    theme(axis.text = element_text(size = rel(0.5)))
}

p_bottom <- (
  # Black-backed jackal near the Etosha Pan in Namibia
  plot_tel(d_akde_sub$tel[[2]], 13, 10000) |
    plot_spacer() |
    # Ibex in the Alps
    plot_tel(d_akde_sub$tel[[1]], 13, 1000) | # ibex 1
    plot_spacer() |
    # Sloth in Costa Rica
    plot_tel(d_akde_sub$tel[[3]], 18, 100, 4)) +
  plot_layout(widths = c(NA, 0.05, NA, 0.05, NA))

# to add legends: add a guide_area() patch and another letter in the layout
p_1 <-
  free(p_1a +
         geom_sf(data = filter(mean_locations, has_panel), size = 0.25,
                 show.legend = FALSE, color = "darkorange1"),
       side = "l") +
  free(p_1b +
         geom_sf(data = filter(mean_locations, has_panel), size = 0.25,
                 show.legend = FALSE, color = "darkorange1"),
       side = "l") +
  p_bottom +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = "
  AAAA
  AAAA
  BBBB
  BBBB
  CCCC") & 
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures/p-1.png", p_1, width = 7.5, height = 10, dpi = 300)

# fig 2 regression plots for HR model ----
mass_preds <- data_slice(m_hr, data = d, mass_kg = evenly(mass_kg)) %>%
  bind_cols(fitted_values(object = m_hr, data = ., scale = "link", 
                          terms = c("(Intercept)", "log(mass_kg)")) %>%
              filter(.parameter == "location"))

# find species for phylopics
unique(d$species[which.min(d$mass_kg)])
unique(d$species[which.min(abs(d$mass_kg - 10))])
unique(d$species[which.max(d$mass_kg)])

partial_residuals(m_hr, select = "s(var_ndvi)")

p_2a <- ggplot(mass_preds) +
  # Carollia castanea not available
  add_phylopic(name = "Desmodus rotundus", x = log(0.01), y = log(3),
               height = 1, verbose = TRUE) +
  # to pick specific phylopic
  add_phylopic(uuid = "ab6cfd4f-aef7-40fa-b5a5-1b79b7d112aa", x = log(10),
               y = log(1), height = 2, verbose = TRUE) +
  add_phylopic(name = "Loxodonta africana", x = log(3e3), y = log(50),
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
ggsave("figures/p-2a.png", p_2a, width = 8, height = 5, dpi = 300)
cli::cli_alert_info("Jittered masses to aid in visualization, but each species has a single mass.")

# TODO: switch data to partial residuals on response scale?
p_2b <-
  data_slice(m_hr, data = d, mean_ndvi = evenly(mean_ndvi),
             var_ndvi = evenly(var_ndvi)) %>%
  bind_cols(fitted_values(object = m_hr, data = ., scale = "link",
                          terms = c("(Intercept)",
                                    "s(mean_ndvi)", "s(var_ndvi)",
                                    "ti(mean_ndvi,var_ndvi)")) %>%
              filter(.parameter == "location")) %>%
  select(-c(.lower_ci, .upper_ci)) %>%
  summarise(.fitted = mean(.fitted),
            .se = sqrt(mean(.se^2)),
            .by = mean_ndvi) %>%
  rename(.estimate = .fitted) %>%
  add_confint() %>%
  ggplot() +
  geom_point(aes(mean_ndvi, log(hr_95_km2)), d, alpha = 0.3) +
  geom_ribbon(aes(mean_ndvi, ymin = .lower_ci, ymax = .upper_ci),
              fill = pal[1], alpha = 0.3) +
  geom_line(aes(mean_ndvi, .estimate), color = pal[1], lwd = 1) +
  xlab("Mean NDVI") +
  scale_y_continuous("Home-range size (km\U00B2)",
                     breaks = log(10^seq(-2, 6, by = 1)),
                     labels = \(x) as.character(exp(x)))
ggsave("figures/p-2b.png", p_2b, width = 8, height = 5, dpi = 300)

p_2c <-
  data_slice(m_hr, data = d, mean_ndvi = evenly(mean_ndvi),
             var_ndvi = evenly(var_ndvi)) %>%
  bind_cols(fitted_values(object = m_hr, data = ., scale = "link",
                          terms = c("(Intercept)",
                                    "s(mean_ndvi)", "s(var_ndvi)",
                                    "ti(mean_ndvi,var_ndvi)")) %>%
              filter(.parameter == "location")) %>%
  select(-c(.lower_ci, .upper_ci)) %>%
  summarise(.fitted = mean(.fitted),
            .se = sqrt(mean(.se^2)),
            .by = var_ndvi) %>%
  rename(.estimate = .fitted) %>%
  add_confint() %>%
  ggplot() +
  geom_point(aes(var_ndvi, log(hr_95_km2)), d, alpha = 0.3) +
  geom_ribbon(aes(var_ndvi, ymin = .lower_ci, ymax = .upper_ci),
              fill = pal[2], alpha = 0.3) +
  geom_line(aes(var_ndvi, .estimate), color = pal[2], lwd = 1) +
  xlab("NDVI variance") +
  scale_y_continuous("Home-range size (km\U00B2)",
                     breaks = log(10^seq(-2, 6, by = 1)),
                     labels = \(x) as.character(exp(x)))

# check animals in bottom right corner
table(as.character(filter(d, var_ndvi > 0.036)$species))
ecdf(d$mass_kg)(max(filter(d, var_ndvi > 0.036)$mass_kg))
cli::cli_alert_info(
  paste0("The animals in the bottom right corner have masses in the ",
         "bottom ",
         round(
           ecdf(d$mass_kg)(max(filter(d, var_ndvi > 0.036)$mass_kg)) * 100
         ),
         "th percentile.")
)
ggsave("figures/p-2c.png", p_2c, width = 8, height = 5, dpi = 300)

p_2 <- (p_2a / p_2b / p_2c) + plot_annotation(tag_levels = 'A')
ggsave("figures/p-2.png", p_2, width = 8, height = 15, dpi = 300)

# fig 3 interaction effects of HR model ----
p_3 <-
  data_slice(m_hr, data = d,
             mean_ndvi = evenly(mean_ndvi),
             var_ndvi = evenly(var_ndvi),
  ) %>%
  mutate(too_far = too_far(mean_ndvi, var_ndvi, d$mean_ndvi, d$var_ndvi,
                           dist = 0.5)) %>%
  filter(!too_far) %>%
  bind_cols(fitted_values(
    m_hr, data = ., terms = c("(Intercept)", "s(mean_ndvi)", "log(mass_kg)",
                              "s(var_ndvi)", "ti(mean_ndvi,var_ndvi)")) %>%
      filter(.parameter == "location")) %>%
  ggplot() +
  geom_raster(aes(mean_ndvi, var_ndvi, fill = log(.fitted))) +
  geom_contour(aes(mean_ndvi, var_ndvi, z = log(.fitted)),
               color = "black") +
  scale_x_continuous("Mean NDVI", expand = c(0, 0)) +
  scale_y_continuous("NDVI variance", expand = c(0, 0)) +
  scale_fill_gradient("Home-range size (km\U00B2)",
                      low = 'grey90', high = pal[3],
                      labels = exp, breaks = log(10^c(-1:4))) +
  theme(legend.position = "top", legend.key.width = rel(2))
cli::cli_alert_info(paste0(
  "Interaction plot is for an animal with a mass of ",
  round(unique(p_3@data$mass_kg)), " kg."))
ggsave("figures/p-3.png", p_3, width = 6, height = 6, dpi = 300)

library("ggcube")
p_3_persp <-
  data_slice(m_hr, data = d,
           mean_ndvi = evenly(mean_ndvi),
           var_ndvi = evenly(var_ndvi),
) %>%
  mutate(too_far = too_far(mean_ndvi, var_ndvi, d$mean_ndvi, d$var_ndvi,
                           dist = 0.5)) %>%
  filter(!too_far) %>%
  bind_cols(fitted_values(
    m_hr, data = ., terms = c("(Intercept)", "s(mean_ndvi)", "log(mass_kg)",
                              "s(var_ndvi)", "ti(mean_ndvi,var_ndvi)")) %>%
      filter(.parameter == "location")) %>%
  ggplot(aes(mean_ndvi, var_ndvi, z = log(.fitted))) +
  geom_surface_3d(color = "black") +
  coord_3d(ratio = c(1.5, 2, 1), expand = FALSE, panels = "zmin",
           light = light(direction = c(1, 1, 0)),
           yaw = 45, ) +
  scale_x_continuous("Mean NDVI", expand = c(0, 0)) +
  scale_y_continuous("NDVI variance", expand = c(0, 0)) +
  scale_fill_gradient("Home-range size (km\U00B2)",
                      low = 'grey90', high = pal[3],
                      labels = exp, breaks = log(10^c(-1:4))) +
  theme(legend.position = "top", legend.key.width = rel(2))
ggsave("figures/p-3-3D.png", p_3_persp, width = 6, height = 6, dpi = 300)
