# Non-exported utilities for the fishtree package

#' The base URL for the Fish Tree of Life API endpoint
#' @noRd
.baseurl <- "https://fishtreeoflife.org/"

#' Gets data from a URL, or a cache keyed by the same
#'
#' Uses the `.cache` internal object to check if an object keyed by `url` exists.
#' If it does, returns it. Otherwise, it will download the `url` to a temporary
#' file and load it with the `reader` function.
#'
#' @param url The URL to retrieve, and the key for the cached data for the same.
#' @param reader A function responsible for parsing the downloaded data. Its
#'   first parameter should take a file name.
#' @param quiet Should we \code{\link{download.file}} quietly? Defaults to `TRUE`.
#' @param ... Additional arguments passed to `reader`.
#' @return Whatever `reader` returns.
#' @import memoise
#' @noRd
.get <- memoise::memoise(function(url, reader, quiet = TRUE, ...) {
  if (rlang::is_missing(reader)) {
    rlang::abort("reader` must be specified when `url` is not in the cache.")
  }
  tmp <- tempfile()
  res <- tryCatch(suppressWarnings(utils::download.file(url, tmp, quiet = quiet)),
                  error = function(e) {
                    rlang::abort(paste("Download for URL", url, "failed with error:\n  ", e, "\nCheck your network status and consider retrying your request."))
                  })
  if (res != 0L) {
    rlang::abort(paste("Download for URL", url, "failed with error code", res))
  }
  reader(tmp, ...)
})

#' Reconcile names against a known good set
#'
#' Ensures that two lists of species names matche up. Automatically accounts
#' for underscores and spaces.
#'
#' @param wanted_names A vector of names to check for validity.
#' @param valid_names A vector of names known to be valid. Defaults to the
#'   species in the phylogeny from \code{\link{fishtree_phylogeny}}.
#' @param warn Warn the user if we find a mismatch? Defaults to `TRUE`.
#' @return A vector of valid names, possibly smaller than `wanted_names`.
#' @seealso \code\link[geiger]{name.check}
#' @noRd
.name_check <- function(wanted_names, valid_names = fishtree_phylogeny()$tip.label, warn = TRUE) {
  missing <- setdiff(gsub("_", " ", wanted_names), gsub("_", " ", valid_names))
  if (!rlang::is_empty(missing)) {
    tmp <- missing
    if (length(missing) > 5) tmp <- c(missing[1:5], paste("...and", length(missing) - 5, "others"))
    missing_str <- paste("*", tmp, collapse = "\n")
    if (warn) rlang::warn(paste0("Requested ", length(wanted_names), " but only found ", length(wanted_names) - length(missing), " species. Missing names:\n", missing_str))
  }
  intersect(gsub("_", " ", wanted_names), gsub("_", " ", valid_names))
}

# Auto detects the rank from the name and downloads the relevant taxonomy file
.fetch_rank <- function(name) {
  if (!rlang::is_scalar_character(name)) {
    rlang::warn(paste0("`name` should be length 1, not ", length(name), ". Only the first element was used (", name[1], ")."))
    name <- name[1]
  }

  tax <- fishtree_taxonomy()
  what <- tax[tax$name == name, "rank"]

  if (length(what) == 1) {
    context <- fishtree_taxonomy(name)
  } else {
    rlang::abort(paste0("Can't find data for `", name, "`."))
  }
  list(context, what)
}

# Splits an object of class DNAbin into partitions based on a RAxML-style partitions description.
.split_seqs <- function(sequence) {
  url <- paste0(.baseurl, "downloads/final_alignment.partitions")
  partitions <- .get(url, readLines)
  tt <- gsub("DNA, ", "", partitions, fixed = TRUE)
  splat <- strsplit(tt, "[= -]+")
  part_names <- sapply(splat, `[`, 1)
  part_start <- as.integer(sapply(splat, `[`, 2))
  part_end <- as.integer(sapply(splat, `[`, 3))
  result = list()
  for (ii in 1:length(part_names)) {
    result[[part_names[ii]]] <- sequence[, part_start[ii]:part_end[ii]]
  }
  result
}
