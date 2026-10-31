# The download cache

#' Manage the Fish Tree of Life download cache
#'
#' Files from a release are downloaded once, checked against the release's
#' checksums, and kept in a cache. `fishtree_cache_info()` lists the cached
#' files and `fishtree_cache_clear()` deletes them. Cached files that the data
#' manifest no longer lists are deleted whenever fishtree downloads a file.
#'
#' The cache is in `tools::R_user_dir("fishtree", "cache")`. Set
#' `options(fishtree.cache_dir = ...)` to use another directory.
#'
#' @return `fishtree_cache_info()` returns a data frame with the `file` name
#'   (its SHA-256), `bytes` and `modified` time of each cached file.
#'   `fishtree_cache_clear()` returns the number of files deleted, invisibly.
#' @export
#' @examples
#' \dontrun{
#' fishtree_cache_info()
#' fishtree_cache_clear()
#' }
fishtree_cache_info <- function() {
  paths <- .cache_files()
  info <- file.info(paths)
  data.frame(file = basename(paths), bytes = info$size, modified = info$mtime, stringsAsFactors = FALSE)
}

#' @rdname fishtree_cache_info
#' @export
fishtree_cache_clear <- function() {
  paths <- .cache_files()
  unlink(paths)
  invisible(length(paths))
}

.cache_dir <- function() getOption("fishtree.cache_dir", tools::R_user_dir("fishtree", "cache"))

.cache_files <- function() list.files(.cache_dir(), pattern = "^[0-9a-f]{64}$", full.names = TRUE)

# Deletes cached files that no release or taxonomy in the manifest lists.
.prune_cache <- function() {
  manifest <- .manifest()
  current <- unlist(lapply(c(manifest$releases, manifest$taxonomies), function(x) vapply(x$files, `[[`, "", "sha256")))
  paths <- .cache_files()
  unlink(paths[!basename(paths) %in% current])
}

.sha256 <- function(path) {
  # tools::sha256sum() is new in R 4.5.0; older versions check only the size.
  tools <- asNamespace("tools")
  if (!exists("sha256sum", envir = tools)) return(NA_character_)
  unname(get("sha256sum", envir = tools)(path))
}

# Path to a cached copy of a release file, downloading and checking it if needed.
.fetch_file <- function(file) {
  dir <- .cache_dir()
  path <- file.path(dir, file$sha256)
  if (file.exists(path) && file.size(path) == file$bytes) return(path)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dir)
  on.exit(unlink(tmp))
  .download(file$url, tmp)
  if (!.matches(tmp, file) && grepl("^https?://", file$url)) {
    # Just after the website updates, a CDN can serve an older copy than the manifest lists.
    .download(.cache_bust(file$url), tmp)
  }
  if (!.matches(tmp, file)) {
    rlang::abort(paste("Download for URL", file$url, "doesn't match the size or checksum in the data manifest.",
                       "The file may have been corrupted in transit; consider retrying your request."))
  }
  file.rename(tmp, path)
  .prune_cache()
  path
}

.matches <- function(path, file) {
  hash <- .sha256(path)
  file.size(path) == file$bytes && (is.na(hash) || hash == file$sha256)
}

# Adds a unique query so caches pass the request through.
.cache_bust <- function(url) {
  paste0(url, if (grepl("?", url, fixed = TRUE)) "&" else "?", "nocache=", sprintf("%.0f", as.numeric(Sys.time()) * 1e6))
}

# Cached paths are named by sha256, so this is memoised in .onLoad.
.read <- function(path, reader, ...) reader(path, ...)
