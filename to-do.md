- plan to submit to PNAS, science, nature, ecology letters
- follow structure in Tucker anthropocene paper (https://doi.org/10.1126/science.aam9712; could also see Jetz et al. 2004; 10.1126/science.1102138)
- stick to one, simple message: both E(R) and Var(R) are drivers of E(H), even after accounting for mass
- 1-3 (likely 1) paragraph of intro
  - theory: optimal foraging theory and marginal value theorem
  - previous papers averaged within species to remove effects of body size, E(R), Var(R), and others
- ~30 refs
- ~3-4 pages total
- figure 1:
  - HR vs body size
  - HR vs mean NDVI
  - HR vs DENVar
  - interaction term?
- figure 2: as 1 but for diffusion? (related to marginal value theorem and optimal foraging theory)
- figure 3: some case examples showing importance of nonlinear effects, or urban environments

methods:

- `ctmm::raster(x,DF="PMF",...)` creates a raster of probability masses, *but it does not sum to 1*!
- need to normalize the raster by dividing by the sum
- only use static rasters of mean and DENVar to avoid issues with time
