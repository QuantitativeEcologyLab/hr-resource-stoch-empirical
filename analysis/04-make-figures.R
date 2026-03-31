library('dplyr')
library('mgcv')
library('ggplot2')
library('khroma')
library('gratia')
source('analysis/figures/default-ggplot-theme.R')

m_hr <- readRDS('models/hr-model.rds')
m_speed <- readRDS('models/diffusion-model.rds')
m_diff <- readRDS('models/speed-model.rds')

# Fig. 1: like figure 1 of Tucker et al. 2018:
# a mean NDVI
# b DENVar
# c boxplots of HR size
# d zoomed in versions of areas (could also be figs from below)
# use shapefile in the background for countries with no data
# 
# or Fig 1 of https://www.pnas.org/doi/10.1073/pnas.2401042122
# - same species, but in different areas with different levels of NDVI and DENVar

# fig 2 regression plots for HR model

# fig 3: same as for fig 2, but for speed, diffusion, or ballistic lengthscale
