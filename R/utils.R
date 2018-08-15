# Non-exported utilities for the fishtree package

.baseurl <- "https://fishtreeoflife.org/"

# Cache so we don't have to redownload the huge files
.cache <- rlang::new_environment()

# Gets data from a URL, or from the cache if it's there.
.get <- function(url, reader = NULL, quiet = TRUE, ...) {
  if (exists(url, envir = .cache)) {
    return(get(url, envir = .cache))
  }

  if (is.null(reader)) rlang::abort("reader` must be specified when `url` is not in the cache.")
  tmp <- tempfile()
  utils::download.file(url, tmp, quiet = quiet)
  assign(url, reader(tmp, ...), envir = .cache)
  .get(url)
}

# Auto detects the rank from the name and downloads the relevant taxonomy file
.fetch_rank <- function(name) {
  if (length(name) > 1) {
    rlang::warn(paste0("`", what, "` should be length 1, not ", length(name), ". Only the first element was used (", name[1], ")."))
    name <- name[1]
  }

  if (endsWith(name, "idae")) {
    context <- fishtree_taxonomy(family = name)
    what <- "family"
  } else if (endsWith(name, "iformes")) {
    context <- fishtree_taxonomy(order = name)
    what <- "order"
  } else {
    rlang::abort(paste0("Can't find data for ", name, " (only families and orders are currently supported)."))
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
