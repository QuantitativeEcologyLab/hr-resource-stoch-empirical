remove_outlier_flags <- function(id) {
  i <- which(d$animal_id == id)
  
  if("original_outliers)" %in% colnames(d$data[[i]])) {
    d$data[[i]]$outlier <<- d$data[[i]]$original_outliers
  } else {
    d$data[[i]]$outlier <<- FALSE
  }
}
