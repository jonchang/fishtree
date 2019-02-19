#' fishtree: Interface to the Fish Tree of Life API
#'
#' fishtree makes it easy to access phylogenetic information about ray-finned
#' fishes with R. Download taxonomies, phylogenies, sequence matrices, and
#' diversification rate data with a simple set of functions.
#'
#' Implementation note: this package makes calls over the network and caches
#' its (sometimes large) results for faster loading. Because of this,
#' long-running R instances could use a lot of memory.
#'
#' @docType package
#' @name fishtree-package
NULL

release_questions <- function() {
  c(
    "Have you regenerated the README with `make`?",
    "Have you committed any changes to the knit README?",
    "Have you set the NOT_CRAN environment variable to 'true' to build the vignette properly?"
  )
}
