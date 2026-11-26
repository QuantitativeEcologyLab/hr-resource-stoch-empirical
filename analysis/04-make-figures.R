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
library("ggcube")
source('analysis/default-ggplot-theme.R')

set.seed(3)
d_sim <- expand_grid(x = 1:20,
                     y = 1:20,
                     mu = c(0.1, 0.35, 0.6),
                     sigma2 = c(0.01, 0.02, 0.03)) %>%
  mutate(phi = sigma2 / (mu * (1 - mu)),
         a = mu * (1 / phi - 1),
         b = (1 - mu) * (1 / phi - 1),
         r = rbeta(n(), a, b),
         mu = paste("Mean =", mu),
         sigma2 = paste("Variance =", sigma2))
d_sim %>%
  summarize(min = min(r),
            max = max(r),
            range = max - min,
            mu = unique(mu),
            var = var(r),
            .by = c(mu, sigma2))
p_patches <-
  ggplot(d_sim) +
  facet_grid(sigma2 ~ mu) +
  geom_raster(aes(x, y, fill = r)) +
  # theme(legend.position = "top") +
  scale_x_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_fill_gradientn("Resources", colours = ndvi_pal, limits = c(-1, 1))

ggsave("figures/p-patches.png", width = 5, height = 5, dpi = 300,
       bg = "white")

# robinson objects
robinson_crs <- 'ESRI:54030'

bounds <-
  dplyr::tibble(x = c(-180:180, rep(180, 181), 180:-180, rep(-180, 181)),
                y = c(rep(-90, 361), -90:90, rep(90, 361), 90:-90)) %>%
  st_as_sf(coords = c('x', 'y'), crs = 'EPSG:4326') %>%
  dplyr::summarise(geometry = st_combine(geometry)) %>%
  st_cast('POLYGON') %>%
  st_transform(robinson_crs)

d <- readRDS("data/gam-dataset-2026-11-01.rds")
m_hr <- readRDS('models/hr-hgam-full.rds')

ggplot(d, aes(mean_ndvi, var_ndvi, color = log(mass_kg))) +
  facet_wrap(~ round(log10(mass_kg))) +
  geom_point() +
  geom_point(shape = 1, color = "black") +
  scale_color_bam() +
  theme(legend.position = "top")

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

mean_locations <- d %>%
  select(species, dataset_name, long, lat) %>%
  st_as_sf(coords = c("long", "lat"), crs = "EPSG:4326") %>%
  mutate(id = paste(species, dataset_name)) %>%
  summarize(species = unique(species),
            dataset_name = unique(dataset_name),
            geometry = st_centroid(st_union(geometry)), .by = id) %>%
  mutate(has_panel = dataset_name %in% d_sub$dataset_name)

d_akde_sub <-  readRDS('data/akde-tibble-2026-10-24.rds') %>%
  filter(grepl(paste(as.character(d_sub$animal) %>%
                       gsub("\\(", "\\\\(", .) %>%
                       gsub("\\)", "\\\\)", .),
                     collapse = "|"
  ), animal))
plot(d_akde_sub$tel[[1]], main = d_akde_sub$animal[1])
plot(d_akde_sub$tel[[2]], main = d_akde_sub$animal[2])
plot(d_akde_sub$tel[[3]], main = d_akde_sub$animal[3])

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
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave("figures/p-1.png", p_1, width = 7.5, height = 10, dpi = 300)

# fig 2 regression plots for HR model ----
mass_preds <- data_slice(m_hr, data = d, mass_kg = evenly(mass_kg)) %>%
  fitted_values(m_hr, data = ., terms = c("(Intercept)", "log(mass_kg)"),
                quantity = "expectation") %>%
  mutate(across(c(.fitted, .lower_ci, .upper_ci), log))

# find species for phylopics
unique(d$species[which.min(d$mass_kg)])
unique(d$species[which.min(abs(d$mass_kg - 10))])
unique(d$species[which.max(d$mass_kg)])

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

p_2b <-
  data_slice(m_hr, data = d, mean_ndvi = evenly(mean_ndvi),
             var_ndvi = c(0.005, 0.015, 0.03), mass_kg = median(mass_kg)) %>%
  fitted_values(object = m_hr, data = ., quantity = "expectation",
                exclude = c("s(species)",
                            "s(log(akde_ess_area))")) %>%
  mutate(across(c(.fitted, .lower_ci, .upper_ci), log)) %>%
  ggplot(aes(group = var_ndvi, color = var_ndvi, fill = var_ndvi)) +
  geom_point(aes(mean_ndvi, log(hr_95_km2)), d, alpha = 0.5) +
  geom_point(aes(mean_ndvi, log(hr_95_km2)), d, alpha = 0.3, shape = 1,
             color = "black") +
  # geom_ribbon(aes(mean_ndvi, ymin = .lower_ci, ymax = .upper_ci),
  #             alpha = 0.3, fill = "grey", color = "black", lwd = 1) +
  # geom_ribbon(aes(mean_ndvi, ymin = .lower_ci, ymax = .upper_ci),
  #             alpha = 0.3) +
  geom_line(aes(mean_ndvi, .fitted), lwd = 1.5, color = "black") +
  geom_line(aes(mean_ndvi, .fitted), lwd = 1) +
  xlab("Mean NDVI") +
  scale_y_continuous("Home-range size (km\U00B2)",
                     breaks = log(10^seq(-2, 6, by = 1)),
                     labels = \(x) as.character(exp(x))) +
  scale_color_devon(name = "NDVI variance", reverse = TRUE,
                    limits = c(0, NA)) +
  scale_fill_devon(name = "NDVI variance", reverse = TRUE,
                   limits = c(0, NA)) +
  theme(text = element_text(face = "bold"),
        legend.position = "top", legend.key.height = rel(0.5),
        legend.key.width = rel(2.75), legend.title = element_text(vjust = 1),
        legend.key.justification = "center",
        legend.box.margin = margin(t = 10))
p_2b
cli::cli_alert_info(
  c("Variance values for Fig. 2b are ",
    paste(round(unique(p_2b@data$var_ndvi), 3), collapse = ", "), "."))
ggsave("figures/p-2b.png", p_2b, width = 8, height = 5, dpi = 300)

p_2c <-
  data_slice(m_hr, data = d, mean_ndvi = c(0, 0.3, 0.6),
             var_ndvi = evenly(var_ndvi), mass_kg = median(mass_kg)) %>%
  fitted_values(object = m_hr, data = ., quantity = "expectation",
                exclude = c("s(log(akde_ess_area))",
                            "s(species)")) %>%
  mutate(across(c(.fitted, .lower_ci, .upper_ci), log)) %>%
  ggplot(aes(group = mean_ndvi, color = mean_ndvi, fill = mean_ndvi)) +
  geom_point(aes(var_ndvi, log(hr_95_km2)), d, alpha = 0.5) +
  geom_point(aes(var_ndvi, log(hr_95_km2)), d, alpha = 0.3,
             shape = 1, color = "black") +
  # geom_ribbon(aes(var_ndvi, ymin = .lower_ci, ymax = .upper_ci),
  #             alpha = 0.3) +
  geom_line(aes(var_ndvi, .fitted), lwd = 1.5, color = "black") +
  geom_line(aes(var_ndvi, .fitted), lwd = 1) +
  xlab("NDVI variance") +
  scale_y_continuous("Home-range size (km\U00B2)",
                     breaks = log(10^seq(-2, 6, by = 1)),
                     labels = \(x) as.character(exp(x))) +
  scale_fill_gradientn("Mean NDVI", colors = ndvi_pal,
                       aesthetics = c("color", "fill"),
                       na.value = "transparent", limits = c(-1, 1)) +
  theme(text = element_text(face = "bold"),
        legend.position = "top", legend.key.height = rel(0.5),
        legend.key.width = rel(2.75), legend.title = element_text(vjust = 1),
        legend.key.justification = "center",
        legend.box.margin = margin(t = 10))
p_2c
cli::cli_alert_info(
  c("Mean values for Fig. 2c are ",
    paste(round(unique(p_2c@data$mean_ndvi), 3), collapse = ", "), "."))
ggsave("figures/p-2c.png", p_2c, width = 8, height = 5, dpi = 300)

p_2 <- (p_2a / p_2b / p_2c) + plot_annotation(tag_levels = 'A')
ggsave("figures/p-2.png", p_2, width = 6, height = 15, dpi = 300)

# fig 3 interaction effects of HR model ----
r_values <- cbind(as.data.frame(round(r_mean_ndvi, 1)),
                  as.data.frame(round(r_var_ndvi, 2))) %>%
  filter(mu_hat > -0.1) %>%
  summarize(n = n(), .by = c(mu_hat, s2_hat)) %>%
  mutate(prop = n / sum(n))

r_values %>%
  filter(prop > 0.01) %>%
  ggplot(aes(mu_hat, s2_hat, fill = prop)) +
  geom_raster() +
  scale_fill_lapaz(name = "Proportion", reverse = TRUE) +
  theme(legend.position = "top")

DIST <- 0.2
pred_3 <- data_slice(m_hr, data = d,
                     mean_ndvi = seq(0, 0.7, length.out = 400),
                     var_ndvi = seq(0, 0.04, length.out = 400)) %>%
  fitted_values(m_hr, data = ., quantity = "expectation",
                exclude = c("s(log(akde_ess_area))",
                            "s(species)")) %>%
  mutate(too_far_a = too_far(mean_ndvi, var_ndvi, # in animal data,
                             d$mean_ndvi, d$var_ndvi,
                             dist = DIST),
         # dist from raster values didn't make a difference, but keeping
         # for stability
         too_far_w = too_far(mean_ndvi, var_ndvi, # in world rasters
                             r_values$mu_hat, r_values$s2_hat,
                             dist = DIST)) %>%
  mutate(.fitted = if_else(too_far_a, NA_real_, .fitted),
         .fitted = if_else(too_far_w, NA_real_, .fitted)) %>%
  filter(any(!is.na(.fitted)), .by = mean_ndvi)

p_int <-
  ggplot(pred_3) +
  geom_raster(aes(mean_ndvi, var_ndvi, fill = log(.fitted))) +
  geom_contour(aes(mean_ndvi, var_ndvi, z = log(.fitted)),
               color = "black", na.rm = TRUE) +
  geom_point(aes(mean_ndvi, var_ndvi), d, alpha = 0.3) +
  scale_x_continuous("Mean NDVI", expand = c(0, 0)) +
  scale_y_continuous("NDVI variance", expand = c(0, 0)) +
  scale_fill_gradientn(name = "Home-range size (km\U00B2)",
                       colors = hr_pal, labels = exp,
                       na.value = "transparent",
                       breaks = log(10^seq(-2, 8, by = 2))) +
  theme(legend.position = "top", legend.key.width = rel(2))
p_int
cli::cli_alert_info(paste0(
  "Interaction plot is for an animal with a mass of ",
  round(unique(pred_3$mass_kg)), " kg."))
ggsave("figures/p-3-interaction.png", p_int, width = 6, height = 6,
       dpi = 300)

# find panels to wrap around the interaction plot ----
p_int + annotate("point", r_values$mu_hat, r_values$s2_hat, size = 0.1)

# cutoffs for mean
low_mean <- r_mean_ndvi < 0.2; plot(low_mean)
med_mean <- abs(r_mean_ndvi - 0.4) < 0.05; plot(med_mean)
high_mean <- r_mean_ndvi > 0.5; plot(high_mean)

# cutoffs for variance (not even because distribution is skewed)
round(quantile(values(r_var_ndvi), na.rm = TRUE), 3)
low_var <- r_var_ndvi < 0.005; plot(low_var)
med_var <- abs(r_var_ndvi - 0.011) < 0.0025; plot(med_var)
high_var <- r_var_ndvi > 0.02; plot(high_var)

plot(low_mean & high_var) # Northern Finland
plot(med_mean & high_var) # Gabon
plot(high_mean & high_var) # Indonesia
plot(low_mean & med_var) # Canada
plot(med_mean & med_var) # Ukraine
plot(high_mean & med_var) # Brazil (Mato Grosso)
plot(low_mean & low_var) # Namibia
plot(med_mean & low_var) # Burkina Faso
plot(high_mean & low_var) # Brazil (Pará)

ex_coords <- tribble(
  ~mean, ~variance, ~country, ~long, ~lat,
  "low", "high", "Finland", 26.75, 67.68,
  "medium", "high", "Gabon", 10.74, -0.83,
  "high", "high", "Indonesia", 115.27, 0.08,
  "low", "medium", "Canada", -129.19, 64.31,
  "medium", "medium", "Ukraine", 25.71, 50.71,
  "high", "medium", "Brazil (Mato Grosso)", -54.79, -10.58,
  "low", "low", "Namibia", 19.14, -27.45,
  "medium", "low", "Burkina Faso", -5.36, 10.85,
  "high", "low", "Brazil (Pará)", -57.12, 0.10
) %>%
  mutate(mean = factor(mean, c("low", "medium", "high")),
         variance = factor(variance, c("low", "medium", "high"))) %>%
  st_as_sf(coords = c("long", "lat"), crs = "EPSG:4326")

p_ex_locations <- ggplot() +
  geom_sf(data = bounds, fill = "white") +
  geom_sf(data = world, fill = "grey", color = "black") +
  geom_sf(aes(shape = variance), size = 2.5,
          st_transform(ex_coords, robinson_crs)) +
  geom_sf(aes(color = mean, shape = variance),
          st_transform(ex_coords, robinson_crs)) +
  scale_color_manual("Mean NDVI",
                     values = c("#B48D5E", "#627224", "#175614")) +
  scale_shape("NDVI Variance") +
  theme_void() +
  theme(text = element_text(face = "bold"),
        legend.position = "top",
        legend.key.justification = "center",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.box.margin = margin(t = 10))
ggsave("figures/p-ex-locations.png", p_ex_locations,
       width = 7, height = 4.5, dpi = 300, bg = "white")

p_biomes <-
  map(1:nrow(ex_coords), function(.i) {
    .d <- slice(ex_coords, .i)
    r <- .d %>%
      st_buffer(5550) %>% # squares are ~11 km * 11 km
      crop(get_tiles(., provider = "Esri.WorldImagery", zoom = 14), .)
    
    pad_default <- unit(.25, "cm")
    width_hint <- 0.2 # npc (normalized parent coords)
    width_label <- strwidth("60 km", units = "inches", cex = .7)
    width_label <- unit(width_label, "in")
    
    ggplot() +
      geom_spatraster_rgb(data = r) +
      scale_x_continuous(expand = c(0, 0),
                         breaks = \(x) round(quantile(x, c(0.2, 0.8)), 2)) +
      scale_y_continuous(expand = c(0, 0),
                         breaks = \(x) round(quantile(x, c(0.2, 0.8)), 2)) +
      labs(subtitle = .d$country) +
      theme(axis.text = element_text(size = rel(0.6)),
            plot.subtitle = element_text(size = rel(0.5)))
  }) %>%
  wrap_plots()

mean_axis <-
  ggplot() +
  annotate("segment", x = 0, xend = 1, y = 0, yend = 0, lwd = 1,
           arrow = arrow(length = unit(0.3, "inches"), type = "closed",
                         angle = 22.5)) +
  annotate("text", x = 0.5, y = -0.1, label = "Mean NDVI",
           fontface = "bold", size = 5) +
  ylim(c(-0.2, 0.2)) +
  theme_void()

var_axis <-
  ggplot() +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 1, lwd = 1,
           arrow = arrow(length = unit(0.3, "inches"), type = "closed",
                         angle = 22.5)) +
  annotate("text", x = -0.1, y = 0.5, label = "NDVI variance", angle = 90,
           fontface = "bold", size = 5) +
  xlim(c(-0.2, 0.2)) +
  theme_void()

p_ex <-
  wrap_plots(mean_axis, var_axis, p_biomes, design = "
#########
BCCCCCCCC
BCCCCCCCC
BCCCCCCCC
BCCCCCCCC
BCCCCCCCC
BCCCCCCCC
BCCCCCCCC
BCCCCCCCC
#AAAAAAAA")

cli::cli_alert_info("Raster cells have a resolution of 0.05x0.05, which is about 5.55x5.55 km^2. Satellite images of biomes are about 11x11 km^2, so each picture would be 2x2 cells in the rasters.")

ggsave("figures/p-biomes.png", p_ex, width = 6, height = 6,
       units = "in", dpi = 300, bg = "white")

p_3 <- wrap_plots(p_int, plot_spacer(), wrap_elements(full = p_biomes),
                  widths = c(4, 0.5, 5)) +
  plot_annotation(tag_levels = c("A", "", "B"))

ggsave("figures/p-3.png", p_3, width = 12, height = 6,
       units = "in", dpi = 300, bg = "white")

# 3D version of p_int ----
d_est <- tibble(mean_ndvi = d$mean_ndvi,
                var_ndvi = d$var_ndvi,
                mass_kg = unique(pred_3$mass_kg),
                akde_ess_area = 1e6,
                species = d$species[1]) %>%
  fitted_values(object = m_hr, data = ., quantity = "expectation",
                exclude = c("s(log(akde_ess_area))", "s(species)"))

p_int_persp <-
  data_slice(m_hr, data = d,
             mean_ndvi = evenly(mean_ndvi),
             var_ndvi = evenly(var_ndvi),
  ) %>%
  fitted_values(m_hr, data = ., quantity = "expectation",
                exclude = c("s(log(akde_ess_area))", "s(species)")) %>%
  # dropping too far values makes the plot take much longer to save
  ggplot(aes(mean_ndvi, var_ndvi, z = log(.fitted))) +
  geom_surface_3d(
    # color = "black",
    alpha = 1,
    annotate = 
      annotate_3d(
        type = "point",
        x = d_est$mean_ndvi,
        y = d_est$var_ndvi,
        z = log(d_est$.fitted) + 0.1,
        size = 0.5)) +
  coord_3d(ratio = c(1.5, 2, 1), expand = FALSE, panels = "zmin",
           light = "none", yaw = 45) +
  scale_x_continuous("Mean NDVI", expand = c(0, 0)) +
  scale_y_continuous("NDVI variance", expand = c(0, 0)) +
  scale_fill_gradientn(name = "Home-range size (km\U00B2)", colors = hr_pal,
                       labels = exp, breaks = log(10^seq(0, 8, by = 1))) +
  theme(legend.position = "top", legend.key.width = rel(2))
ggsave("figures/p-3-3D.png", p_int_persp, width = 6, height = 6, dpi = 300)
