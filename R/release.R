# Tree releases and the data manifest

# The release used when none is given. Fixed per fishtree version, so a new
# release on the website never changes results.
.default_release <- "R2018"

.default_manifest <- "https://fishtreeoflife.org/data/manifest.json"

#' List releases on the Fish Tree of Life
#'
#' Lists every release available on the Fish Tree of Life. Choose one
#' with the `release` argument of the other functions, or for the whole session
#' with `options(fishtree.release = "...")`.
#'
#' @return A data frame with one row per release, with columns `id`, `title`,
#'   `status`, `native_taxonomy`, `doi` and `citation`.
#' @seealso \code{\link{fishtree_release}} for one release's citation and files.
#' @export
#' @examples
#' \dontrun{
#' fishtree_releases()
#' }
fishtree_releases <- function() {
  releases <- .manifest()$releases
  data.frame(
    id = vapply(releases, `[[`, "", "id"),
    title = vapply(releases, `[[`, "", "title"),
    status = vapply(releases, `[[`, "", "status"),
    native_taxonomy = vapply(releases, `[[`, "", "native_taxonomy"),
    doi = vapply(releases, function(r) r$doi %||% NA_character_, ""),
    citation = vapply(releases, `[[`, "", "citation"),
    stringsAsFactors = FALSE
  )
}

#' Get the details of a Fish Tree of Life release
#'
#' Describes one tree release. Includes citation information and checksums.
#'
#' @param release The release's id, from \code{\link{fishtree_releases}}. Defaults to
#'   `getOption("fishtree.release")`, or `"R2018"` if that is unset.
#' @return A list with the release's `id`, `title`, `status`, `citation`,
#'   `doi`, `license` and `native_taxonomy`, and `files`, a data frame of its
#'   files.
#' @seealso \code{\link{fishtree_releases}} to list releases.
#' @export
#' @examples
#' \dontrun{
#' rel <- fishtree_release("R2018")
#' rel$citation
#' rel$files[, c("role", "url", "bytes")]
#' }
fishtree_release <- function(release = getOption("fishtree.release", .default_release)) {
  rel <- .release(release)
  files <- lapply(rel$files, .resolve_file)
  out <- rel[intersect(c("id", "title", "status", "citation", "doi", "license", "native_taxonomy"), names(rel))]
  out$files <- data.frame(
    role = vapply(files, `[[`, "", "role"),
    title = vapply(files, function(f) f$title %||% NA_character_, ""),
    url = vapply(files, `[[`, "", "url"),
    format = vapply(files, `[[`, "", "format"),
    bytes = vapply(files, function(f) as.numeric(f$bytes), 0),
    sha256 = vapply(files, `[[`, "", "sha256"),
    stringsAsFactors = FALSE
  )
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.manifest_url <- function() {
  url <- getOption("fishtree.manifest", .default_manifest)
  if (!grepl("^[A-Za-z][A-Za-z0-9+.-]*:", url)) url <- paste0("file://", normalizePath(url, winslash = "/", mustWork = TRUE))
  url
}

.manifest <- function() {
  url <- .manifest_url()
  manifest <- .get(url, jsonlite::read_json)
  if (!identical(manifest$schema_version, 1L)) {
    rlang::abort(paste0("The data manifest at ", url, " has schema version ", manifest$schema_version,
                        ", which this version of fishtree can't read. Please update fishtree."))
  }
  manifest
}

.release <- function(release) {
  if (!rlang::is_scalar_character(release)) rlang::abort("`release` must be a single release id, like \"R2018\".")
  releases <- .manifest()$releases
  ids <- vapply(releases, `[[`, "", "id")
  if (!release %in% ids) {
    rlang::abort(paste0("Unknown release \"", release, "\". Available releases: ", paste(ids, collapse = ", "), "."))
  }
  releases[[match(release, ids)]]
}

# Resolves a URL relative to `base`, as a browser would.
.resolve_url <- function(url, base) {
  if (grepl("^[A-Za-z][A-Za-z0-9+.-]*:", url)) return(url)
  origin <- regmatches(base, regexpr("^[A-Za-z][A-Za-z0-9+.-]*://[^/]*", base))
  dir <- sub("[^/]*$", "", substring(base, nchar(origin) + 1))
  out <- character()
  for (seg in strsplit(paste0(dir, url), "/", fixed = TRUE)[[1]]) {
    if (seg == "..") out <- out[-length(out)]
    else if (seg != "." && seg != "") out <- c(out, seg)
  }
  paste0(origin, "/", paste(out, collapse = "/"), if (endsWith(url, "/") && length(out)) "/")
}

.resolve_file <- function(file) {
  file$url <- .resolve_url(file$url, .manifest_url())
  file
}

# The release's file with this role. complete_trees prefers the native taxonomy.
.release_file <- function(release, role) {
  rel <- .release(release)
  files <- Filter(function(f) identical(f$role, role), rel$files)
  if (role == "complete_trees") {
    native <- Filter(function(f) identical(f$taxonomy, rel$native_taxonomy), files)
    if (length(native)) files <- native
  }
  if (!length(files)) rlang::abort(paste0("Release ", release, " has no ", role, " file."))
  .resolve_file(files[[1]])
}

.read_release_file <- function(release, role, reader, ...) {
  .read(.fetch_file(.release_file(release, role)), reader, ...)
}

# The legacy per-taxon JSON API for a release, if it has one.
.legacy_api <- function(release) {
  rel <- .release(release)
  views <- Filter(function(v) identical(v$taxonomy, rel$native_taxonomy) && !is.null(v$legacy_api), rel$views)
  if (!length(views)) {
    rlang::abort(paste0("Looking up a `taxon` isn't available for release ", release,
                        " in this version of fishtree."))
  }
  base <- .resolve_url(views[[1]]$legacy_api, .manifest_url())
  list(
    index = paste0(sub("/$", "", base), ".json"),
    base = base,
    # Download paths inside the legacy JSON are relative to the site root.
    root = .resolve_url("../../", base)
  )
}
