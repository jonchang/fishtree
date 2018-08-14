.cache <- rlang::new_environment()

.baseurl <- "https://fishtreeoflife.org/"

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

#' Download a phylogeny
#'
#' Retrieves a phylogeny from the fish tree of life website.
#'
#' @param name Download phylogenies for this rank (only families and orders are currently supported)
#' @param type Either \code{chronogram} or \code{phylogram}
#' @return An object of class \code{phylo}
#' @examples
#' surgeons <- fishtree_phylogeny("Acanthuridae")
#' # Chronograms may not be ultrametric due to numerical precision issues
#' ape::is.ultrametric(surgeons)
#' ape::is.ultrametric(surgeons, tol = 0.00001)
#' @seealso [ape::read.tree()]
#' @export
fishtree_phylogeny <- function(name = NULL, type = c("chronogram", "phylogram")) {
  type <- rlang::arg_match(type)

  if (is.null(name)) {
    if (type == "chronogram")
      return(.get("https://fishtreeoflife.org/downloads/actinopt_12k_treePL.tre.xz", ape::read.tree))
    if (type == "phylogram")
      return(.get("https://fishtreeoflife.org/downloads/actinopt_12k_raxml.tre.xz"), ape::read.tree)
  }

  res <- .fetch_rank(name)
  context <- res[[1]][[1]]
  what <- res[[2]]

  if (length(context) < 1) rlang::abort(paste("Can't find data for", what, name))

  url <- context[[type]]
  if (is.null(url)) {
      msg <- paste("Can't retrieve a", type, "for", what, name, "because")
      if (length(context$species) == 1) msg <- paste(msg, "it is monotypic.")
      else if (length(context$species) < 3) msg <- paste0(msg, " it has too few species (", length(context$species), ").")
      else if (length(context$sampled_species) < 3) msg <- paste0(msg, " it has too few species sampled (", length(context$sampled_species), ").")
      rlang::abort(msg)
  }
  return(.get(paste0(.baseurl, url), ape::read.tree))
}

#' Download a taxonomy
#'
#' Retrieves taxonomic information from the fish tree of life website. One of
#' \code{family} or \code{order} must be specified.
#'
#' @param family retrieve one or more families
#' @param order retrieve one or more orders
#' @return A list, with components containing data on the specified family or order.
#' @export
fishtree_taxonomy <- function(family = NULL, order = NULL) {
  if (!is.null(family) && !is.null(order))
    rlang::abort("Either `family` or `order` must be specified, not both.")

  if (!is.null(family)) {
    js <- .get("https://fishtreeoflife.org/api/family.json", jsonlite::fromJSON)
    ff <- js[family]
    if (length(ff) < 1) rlang::abort(paste("No results found for", family))
  }

  if (!is.null(order)) {
    js <- .get("https://fishtreeoflife.org/api/order.json", jsonlite::fromJSON)
    ff <- js[order]
    if (length(ff) < 1) rlang::abort(paste("No results found for", order))
  }
  return(ff)
}

.fastDNAbin <- function(dnafile, nlines = 0) {
  scan(dnafile, what = list(character(), character()), quiet = TRUE, nlines = nlines, strip.white = TRUE, skip = 1)
}

#' Download aligned sequences
#'
#' Retrieves an aligned sequence
#'
#' @param name Download an alignment for this rank (only families and orders are currently supported), or NULL for the alignment for all species.
#' @param split Split the alignment into a list of sequences by gene?
#' @return An object of class \code{DNAbin}
#' @seealso [ape::DNAbin()]
#' @export
fishtree_alignment <- function(name = NULL, split = FALSE) {
  if (is.null(name)) {
    url <- "https://fishtreeoflife.org/downloads/final_alignment.phylip.xz"
    if(!exists(url, envir = .cache)) rlang::inform("Loading alignment for all species, this will take a while...")
    nlines <- 11650
  } else {
    res <- .fetch_rank(name)
    context <- res[[1]][[1]]
    what <- res[[2]]

    if (is.null(context$matrix_phylip)) rlang::abort(paste("Can't find sequences for", what, name, "because there are too few sampled species."))
    url <- paste0(.baseurl, context$matrix_phylip)
    nlines <- length(context$sampled_species) + 2
  }

  dna <- .get(url, ape::read.dna, format = "sequential", nlines = nlines)
  if (split) dna <- .split_seqs(dna)
  return(dna)
}

.split_seqs <- function(sequence, raxml_partition_file = paste0(.baseurl, "downloads/final_alignment.partitions")) {
  partitions <- readLines(raxml_partition_file)
  tt <- stringr::str_replace_all(partitions, "DNA, ", "")
  splat <- stringr::str_split_fixed(tt, "[= -]+", 3)
  part_names <- splat[, 1]
  part_start <- as.integer(splat[, 2])
  part_end <- as.integer(splat[, 3])
  result = list()
  for (ii in 1:length(part_names)) {
    result[[part_names[ii]]] <- sequence[, part_start[ii]:part_end[ii]]
  }
  result
}

