# gain in percentage of deviance explained (m_0 the is simpler model)
perc_dev_change <- function(m_1, m_0) {
  x <- (1 - m_1$deviance / m_1$null.deviance) * 100 -
    (1 - m_0$deviance / m_0$null.deviance) * 100
  paste0(round(x, 2), "%")
}

# percentage of previously unexplained deviance now explained
perc_dev_prev_unexpl <- function(m_1, m_0) {
  x <- ((1 - m_1$deviance / m_1$null.deviance) - (1 - m_0$deviance / m_0$null.deviance)) /
    (m_0$deviance / m_0$null.deviance) * 100
  paste0(round(x, 2), "%")
}

# AIC comparisons
compare_aic <- function(...) {
  model_names <- as.character(rlang::exprs(...))
  models <- list(...)
  names(models) <- model_names
  dataset_size <- purrr::map(list(...), \(x) nrow(x$model))
  if (any(!any(parametric_terms(m_hr) == 'log(mass_kg)'))) {
    stop("The model does not have a term of mass and is thus not comparable to m_hr_ndvi!")
  } else if (!all.equal(nrow(m_hr$model), nrow(m_hr_ndvi$model),
                        nrow(m_hr_mass_only$model))) {
    stop("Datasets are not the same size!")
  } else {
    AIC(...) %>%
      mutate(model = model_names,
             AIC = AIC - min(AIC)) %>%
      relocate(model, .before = 1) %>%
      arrange(AIC) %>%
      `rownames<-`(NULL)
  }
}
