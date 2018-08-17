
#' Phylogenies from the Fish Tree of Life
#'
#' Retrieves a phylogeny via the Fish Tree of Life API
#'
#' @param name Download phylogenies for this rank (only families and orders are currently supported)
#' @param type Either \code{chronogram} or \code{phylogram}
#' @return An object of class \code{phylo}
#' @references Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L., Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro, M. E. (2018). An inverse latitudinal gradient in speciation rate for marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1
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

#' Taxonomies and other data from the Fish Tree of Life
#'
#' Retrieves taxonomic and other information from the Fish Tree of Life API. One of
#' \code{family} or \code{order} must be specified.
#'
#' @param family retrieve one or more families
#' @param order retrieve one or more orders
#' @return A list, with components containing data on the specified family or order.
#' @export
#' @references Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L., Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro, M. E. (2018). An inverse latitudinal gradient in speciation rate for marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1
#' @examples
#' test <- fishtree_taxonomy(family = "Labridae")
#' paste("There are ", length(test$sampled_species), "sampled species out of",
#'       length(test$species), "in wrasses.")
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

#' Aligned sequences from the Fish Tree of Life
#'
#' Retrieves an aligned sequence via the Fish Tree of Life API.
#'
#' @param name Download an alignment for this rank (only families and orders are currently supported), or NULL for the alignment for all species.
#' @param split Split the alignment into a list of sequences by gene?
#' @return An object of class \code{\link[ape]{DNAbin}}, or a named list of the same if \code{split = TRUE}
#' @references Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L., Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro, M. E. (2018). An inverse latitudinal gradient in speciation rate for marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1
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

#' Tip rates for the Fish Tree of Life
#'
#' Downloads tip rates for the entire Fish Tree of Life, or for a specified subset. Tip rates can be thought of as an
#' instantaneous speciation or extinction rate; for example, a higher tip-specific speciation rate might imply that
#' a lineage is more likely to split a new lineage at the present time.
#'
#' @param name optionally subset by a taxonomy rank (currently restricted to family and order)
#' @param sampled_only only include taxa actually present in the phylogeny?
#' @return a data.frame
#' @export
#' @references
#' DR rates (supplement, section 1.2.2): Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012). The global diversity of birds in space and time. Nature, 491(7424), 444–448. doi:10.1038/nature11631
#'
#' BAMM rates: Rabosky, D. L. (2014). Automatic Detection of Key Innovations, Rate Shifts, and Diversity-Dependence on Phylogenetic Trees. PLoS ONE, 9(2), e89543. doi:10.1371/journal.pone.0089543
#'
#' Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L., Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro, M. E. (2018). An inverse latitudinal gradient in speciation rate for marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1
#'
#' Enhanced polytomy resolution strengthens evidence for global gradient in speciation rate for marine fishes. \url{https://fishtreeoflife.org/rabosky-et-al-2018-update/}
#' @examples
#' # Get cichlid rates and trees
#' rates <- fishtree_tip_rates("Cichlidae")
#' tree <- fishtree_phylogeny("Cichlidae")
#'
#' # Plot tree and extract plotting data
#' plot(tree, show.tip.label = FALSE)
#' obj <- get("last_plot.phylo", ape::.PlotPhyloEnv)
#'
#' # Generate a color ramp
#' ramp <- grDevices::colorRamp(c("black", "red"), bias = 10)
#' tiporder <- match(rates$species, gsub("_", " ", tree$tip.label))
#' scaled_rates <- rates$lambda.tv / max(rates$lambda.tv, na.rm = TRUE)
#' tipcols <- apply(ramp(scaled_rates), 1, function(x) do.call(rgb, as.list(x / 255)))
#'
#' # Place colored bars
#' for (ii in 1:length(tiporder)) {
#'     tip <- tiporder[ii]
#'     lines(x = c(obj$xx[tip] + 0.5, obj$xx[tip] + 0.5 + scaled_rates[ii]),
#'           y = rep(obj$yy[tip], 2),
#'           col = tipcols[ii])
#' }
fishtree_tip_rates <- function(name = NULL, sampled_only = TRUE) {
  rates <- .get("https://fishtreeoflife.org/downloads/tiprates.csv.xz", utils::read.csv, row.names = NULL)
  if (is.null(name)) return(rates)
  res <- .fetch_rank(name)
  ctx <- res[[1]][[1]]
  what <- res[[2]]
  if (sampled_only) wanted <- ctx$sampled_species
  else wanted <- ctx$species
  if (length(wanted) < 1) rlang::abort(paste("Can't get tip rates for", what, name, "because there are not enough sampled species"))
  return(rates[rates$species %in% wanted, ])
}
