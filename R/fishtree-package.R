#' fishtree: Interface to the Fish Tree of Life API
#'
#' fishtree makes it easy to access phylogenetic information about ray-finned
#' fishes with R. Download taxonomies, phylogenies, sequence matrices, and
#' diversification rate data with a simple set of functions.
#'
#' Implementation note: this package makes calls over the network. It keeps
#' downloaded files in a cache on disk (see \code{\link{fishtree_cache_info}})
#' and their (sometimes large) parsed results in memory, so long-running R
#' instances could use a lot of memory.
#'
#' @docType package
#' @keywords internal
"_PACKAGE"


release_questions <- function() {
  c(
    "Have you regenerated the README with `make`?",
    "Have you committed any changes to the knit README?",
    "Have you set the NOT_CRAN environment variable to 'true' to build the vignette properly?"
  )
}
