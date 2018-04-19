
lazy_series <- function(bvec, i) {
  structure(list(
    bvec=bvec,
    i=i),
    class="lazy_series")
}

as.vector.lazy_series <- function(x) {
  series(x$bvec, x$i)
}



to_tibble.fmri_mem_dataset <- function(x) {
  
}
