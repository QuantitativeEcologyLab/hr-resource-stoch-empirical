# function to import Rda file but with a custom name
import_rda <- function(file, object_name = NA, rm_rda = TRUE) {
  load(file)
  return(get(object_name))
}

# test the function
if(FALSE) {
  URL <- 'https://github.com/ctmm-initiative/ctmmlearn/raw/refs/heads/main/data/buffalo.rda'
  class(import_rda(url(URL), object_name = 'FITS'))
  class(import_rda(url(URL), object_name = 'FITS')[[1]])
}
