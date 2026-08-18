plot_start <- function(id, n = 100) {
  if(class(id) == "characer") {
    i <- which(d$animal_id == id)
  } else {
    i <- id
  }
  
  if(length(i) == 0) {
    stop("ID does not exist.")
  }
  plot(d$tel[[i]], type = "l", error = FALSE, col = "black")
  plot(d$tel[[i]][1:n, ], type = "l", error = FALSE, lwd = 2, add = TRUE)
}
