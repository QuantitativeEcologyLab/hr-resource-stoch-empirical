- stick to open data only for the thesis to make the thesis faster
- plan to submit to PNAS, science, nature, ecology letters
- stick to one, simple message: both E(R) and Var(R) are drivers of E(H), even after accounting for mass
- see Jetz et al. 2004 for reference (10.1126/science.1102138)
- 2-3 paragraphs for the intro:
  - weave in theory: optimal foraging theory and marginal value theorem
  - list limitations of the previous papers frameworks:
    - averaged within species to remove effects of body size, E(R), Var(R), and others
  - what is deviance explained of the model with mass, E(R), and Var(R)?
- create figures of Var(HR) vs:
  - mean NDVI (for a single species) to show that other drivers are important?
  - body mass (all species) to show large range of variation in responses
- add an analysis for diffusion the same as for HR to show that animals explore more at low E(R) and high Var(R)
  - related to marginal value theorem and optimal foraging theory

**NOTE**:
- `ctmm::raster(x,DF="PMF",...)` creates a raster of probability masses, but it does not sum to 1!
- need to normalize the raster by dividing by the sum
