.onLoad <- function(...) {
  .get <<- memoise::memoise(.get)
  .read <<- memoise::memoise(.read)
}
