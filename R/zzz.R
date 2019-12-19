.onLoad <- function(...) {
  .get <<- memoise::memoise(.get)
}
