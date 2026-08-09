# function to flag outliers:
#' `value == 0` => point should be kept for analyses
#' `value == 1` => point to be removed: - distance is excessive,
#'                                          - speed is excessive,
#'                                          - outside of DOP ellipses,
#'                                          - turning angle is too small...
#' `value == 2` => point *may* need to be removed (but, e.g., speed is ok)
#' **NOTE:** `as.telemetry()` removes all rows with `outlier > 0`
flag_outlier <- function(id, ..., max_speed = 0, max_distance = 0,
                         max_angle = 0, max_dt = 0, value) {
  #' `remotes::install_github(repo = "poissonconsulting/chk")`
  chk::chk_unused(...)
  o <- out$speed >= max_speed &
    out$distance >= max_distance &
    (out$angle >= max_angle | is.na(out$angle)) &
    (out$dt >= max_dt | is.na(out$dt))  
  .data <- d$data[d$animal_id == id][[1]] # extract data
  if (is.logical(.data$outlier)) {
    .data$outlier <- as.numeric(.data$outlier)
    d$data[d$animal_id == id][[1]]$outlier <- .data$outlier
  }
  if (! 'original_outliers' %in% colnames(.data)) {
    d$data[d$animal_id == id][[1]]$original_outliers <<-
      d$data[d$animal_id == id][[1]]$outlier
  }
  .data[! .data$outlier, ][o, 'outlier'] <- value # skip outliers and add flags
  d$data[d$animal_id == id][[1]] <<- .data
}
