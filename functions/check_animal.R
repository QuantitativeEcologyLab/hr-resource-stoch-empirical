source('functions/outlier_plots.R')

#' function to run `outlie_plots()` for a given animal
check_animal <- function(id, UERE = NULL, return_out = TRUE, .data = d, ...) {
  .d <- filter(.data, animal_id == id)
  if (nrow(.d) != 1) {
    warning("pattern does not match one animal exactly.")
    print(.d)
  }
  tel <- suppressMessages(as.telemetry(.d$data[[1]], mark.rm = TRUE))
  outlier_plots(tel, return = return_out, ...)
}
