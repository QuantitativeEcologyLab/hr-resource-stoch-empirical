remove_outlier_flags <- function(id) {
  d$data[[which(d$animal_id == id)]]$outlier <<-
    d$data[[which(d$animal_id == id)]]$original_outliers
}
