
#' Get a phylogeny from the Fish Tree of Life
#'
#' Retrieves a phylogeny via the Fish Tree of Life API. If neither `species` nor `rank` are specified, returns the entire phylogeny.
#'
#' For maximum interoperability, `species` considers spaces and underscores equivalently. Internally, the phylogenies use underscores.
#'
#' @param species (Optionally) subset the results based on a vector of species names.
#' @param rank (Optionally) subset the results based on the supplied taxonomic rank.
#' @param type Either `"chronogram"` or `"phylogram"`. A chronogram has branch lengths proportional to units of time, while a phylogram has branch lengths proportional to the amount of character change. When retrieving a phylogeny by rank, and that rank is not recovered as monophyletic, acceptable types also include `"chronogram_mrca"` and `"phylogram_mrca"`, which returns a tree with *all* species descending from the common ancestor of species in the specified rank.
#' @return An object of class `"phylo"`.
#' @references Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L., Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro, M. E. (2018). An inverse latitudinal gradient in speciation rate for marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1
#' @examples
#' \dontrun{
#' # Get a phylogeny for a taxonomic rank
#' surgeons <- fishtree_phylogeny(rank = "Acanthuridae")
#'
#' # Get a phylogeny for only certain species
#' genomic_fish <- c("Oryzias latipes", "Tetraodon nigroviridis",
#'                   "Gasterosteus aculeatus", "Danio rerio")
#' fishtree_phylogeny(species = genomic_fish)
#'
#' # Chronograms may not be ultrametric due to numerical precision issues
#' # Consider using phytools::force.ultrametric
#' ape::is.ultrametric(surgeons)
#' ape::is.ultrametric(surgeons, tol = 0.00001)
#'
#' # Difference between MRCA trees and regular trees
#' gobies_mrca <- fishtree_phylogeny(rank = "Gobiidae", type = "chronogram_mrca")
#' gobies <- fishtree_phylogeny(rank = "Gobiidae", type = "chronogram")
#' # MRCA trees will have more tips for non-monophyletic groups
#' length(gobies_mrca$tip.label) > length(gobies$tip.label)
#' # Drop rogue tips in the MRCA tree
#' rogue_gobies <- fishtree_rogues("Gobiidae")
#' pruned_gobies <- ape::drop.tip(gobies_mrca, rogue_gobies)
#' # Now the trees are identical
#' setequal(gobies$tip.label, pruned_gobies$tip.label)
#' }

#' @seealso \code{\link[fishtree]{fishtree_rogues}}, \code{\link[ape]{read.tree}}, \code{\link[phytools]{force.ultrametric}}
#' @export
fishtree_phylogeny <- function(species, rank, type = c("chronogram", "phylogram", "chronogram_mrca", "phylogram_mrca")) {
  if (!rlang::is_missing(species) && !rlang::is_missing(rank)) rlang::abort("Must supply at most one of either `species` or `rank`, not both")

  type <- rlang::arg_match(type)
  fullurl <- switch(type,
                    chronogram = "https://fishtreeoflife.org/downloads/actinopt_12k_treePL.tre.xz",
                    chronogram_mrca = "https://fishtreeoflife.org/downloads/actinopt_12k_treePL.tre.xz",
                    phylogram = "https://fishtreeoflife.org/downloads/actinopt_12k_raxml.tre.xz",
                    phylogram_mrca = "https://fishtreeoflife.org/downloads/actinopt_12k_raxml.tre.xz")

  if (rlang::is_missing(rank)) {
    if (rlang::is_missing(species)) return(.get(fullurl, ape::read.tree))
    if (length(species) < 2) rlang::abort("Must include at least 2 tips in `species`")
    valid_names <- .name_check(species)
    if (length(valid_names) < 2) rlang::abort("Must include at least 2 sampled tips in `species`")
    tree <- .get(fullurl, ape::read.tree)
    return(ape::keep.tip(tree, tree$tip.label[tree$tip.label %in% gsub(" ", "_", valid_names)]))
  }

  res <- .fetch_rank(rank)
  context <- res[[1]][[1]]
  what <- res[[2]]

  if (length(context) < 1) rlang::abort(paste("Can't find data for", what, rank))

  url <- context[[type]]

  # `type` is *_mrca
  if (rlang::is_empty(url) && type %in% c("chronogram_mrca", "phylogram_mrca")) {
    rlang::inform(paste0("Can't retrieve a `", type, "` for ", what, " ", rank, " because the family is monophyletic."))
    type <- sub("_mrca", "", type)
    url <- context[[type]]
  }

  # `type` is NOT *_mrca and we've already possibly corrected for asking for a
  # MRCA tree for a monophyletic group
  if (rlang::is_empty(url)) {
    msg <- paste("Can't retrieve a", type, "for", what, rank, "because")
    if (length(context$species) == 1) msg <- paste(msg, "it is monotypic.")
    else if (length(context$species) < 3) msg <- paste0(msg, " it has too few species (", length(context$species), ").")
    else if (length(context$sampled_species) < 3) msg <- paste0(msg, " it has too few species sampled (", length(context$sampled_species), ").")
    rlang::abort(msg)
  }

  # Family not monophyletic?
  if (length(fishtree_rogues(rank)) > 0 && type %in% c("chronogram", "phylogram")) {
    rlang::inform(paste0(what, " ", rank, ' is not monophyletic. To instead retrieve the phylogeny descending from the common ancestor of all species in ', rank, ', use `type = "', paste0(type, "_mrca"), '"`'))
  }
  return(.get(paste0(.baseurl, url), ape::read.tree))
}

#' Get rogue taxa that break the monophyly of defined taxa
#'
#' For groups that were recovered as paraphyletic in the phylogenetic analysis,
#' uses the Fish Tree of Life API to identify which species prevented that clade
#' from being recovered as monophyletic.
#'
#' @param rank the (possibly paraphyletic) rank for which rogue or intruder species should be identified.
#' @return A vector of species names, potentially empty.
#' @export
#' @references Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L., Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro, M. E. (2018). An inverse latitudinal gradient in speciation rate for marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1
#' @examples
#' \dontrun{
#' fishtree_rogues("Gobiidae")   # several rogue taxa!
#' fishtree_rogues("Labridae")   # nice and monophlyetic
#' }
fishtree_rogues <- function(rank) {
  if (rlang::is_missing(rank))
    rlang::abort("`rank` must be specified.")

  res <- .fetch_rank(rank)
  context <- res[[1]][[1]]
  what <- res[[2]]

  if (length(context) < 1) rlang::abort(paste("Can't find data for", what, rank))

  return(as.vector(context[["rogues"]], mode = "character"))
}

#' Get taxonomies and other data from the Fish Tree of Life
#'
#' Retrieves taxonomic and other information from the Fish Tree of Life API.
#'
#' @param ranks One or more taxonomic ranks to retrieve.
#' @return A list, with components containing data on the specified taxa. If `ranks`
#' is unspecified, a data frame with all valid taxa is returned instead.
#' @export
#' @references Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L., Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro, M. E. (2018). An inverse latitudinal gradient in speciation rate for marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1
#' @examples
#' \dontrun{
#' tax <- fishtree_taxonomy(rank = "Labridae")
#' n_total <- length(tax$Labridae$species)
#' n_sampl <- length(tax$Labridae$sampled_species)
#' paste("There are", n_sampl, "sampled species out of", n_total, "in wrasses.")
#' }
fishtree_taxonomy <- function(ranks = NULL) {
  tax <- .get("https://fishtreeoflife.org/api/taxonomy.json", jsonlite::fromJSON)
  tax_df <- utils::stack(tax)
  colnames(tax_df) <- c("name", "rank")
  tax_df <- tax_df[c("rank", "name")]
  if (is.null(ranks)) return(tax_df)

  wanted <- tax_df[tax_df$name %in% ranks, ]
  if (nrow(wanted) < 1) rlang::abort("No matching taxa found.")

  output <- list()
  for (idx in 1:nrow(wanted)) {
    row <- wanted[idx, ]
    url <- paste0("https://fishtreeoflife.org/api/taxonomy/", row$rank, "/", row$name, ".json")
    output[[row$name]] <- .get(url, jsonlite::fromJSON)
  }

  output
}

#' Get aligned sequences from the Fish Tree of Life
#'
#' Retrieves an aligned sequence via the Fish Tree of Life API. If neither `species` nor `rank` are specified, returns the entire sequence matrix.
#'
#' @inheritParams fishtree_phylogeny
#' @param split Splits the output into a list by gene locus.
#' @return An object of class `"DNAbin"`, or a named list of the same if `split = TRUE``
#' @references Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L., Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro, M. E. (2018). An inverse latitudinal gradient in speciation rate for marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1
#' @seealso \link[ape]{DNAbin}
#' @export
#' @examples
#' \dontrun{
#' surgeon_dna <- fishtree_alignment(rank = "Acanthuridae", split = TRUE)
#' surgeon_dna[[1]]
#' par(mfrow = c(9, 3), mar = c(0.5, 0.5, 1, 0.5), xaxt = "n", yaxt = "n")
#' for (gene in names(surgeon_dna)) {
#'   image(surgeon_dna[[gene]], legend = FALSE, show.labels = FALSE)
#'   title(gene)
#' }
#' }
fishtree_alignment <- function(species, rank, split = FALSE) {
  if (!rlang::is_missing(species) && !rlang::is_missing(rank))
    rlang::inform("Supplying both `species` and `rank` arguments may limit the number of results you see.")

  if (rlang::is_missing(rank)) {
    url <- "https://fishtreeoflife.org/downloads/final_alignment.phylip.xz"
    nlines <- 11650
  } else {
    res <- .fetch_rank(rank)
    context <- res[[1]][[1]]
    what <- res[[2]]

    if (is.null(context$matrix_phylip)) rlang::abort(paste("Can't find sequences for", what, rank, "because there are too few sampled species."))
    url <- paste0(.baseurl, context$matrix_phylip)
    nlines <- length(context$sampled_species) + 2
  }

  dna <- .get(url, ape::read.dna, format = "sequential", nlines = nlines)
  if (!rlang::is_missing(species)) {
    valid_names <- .name_check(species)
    dna <- dna[gsub(" ", "_", valid_names), ]
  }
  if (split) dna <- .split_seqs(dna)
  dna
}

#' Get tip rates for the Fish Tree of Life
#'
#' Downloads tip rates for the entire Fish Tree of Life, or for a specified subset. Tip rates can be thought of as an
#' instantaneous speciation or extinction rate; for example, a higher tip-specific speciation rate might imply that
#' a lineage is more likely to split a new lineage at the present time. See Title (2019) in references for details.
#' If neither `species` nor `rank` are specified, returns the entire set of tip-specific diversification rates.
#'
#' @inheritParams fishtree_phylogeny
#' @param sampled_only Restricts the returned dataset to only those species that have genetic data available. Defaults to `TRUE`.
#' @return A data frame. Columns ending with `.tv` indicate time-variable BAMM runs; those ending in `.tc` are time-constant runs. The `dr` column refers to the DR statistic, while `lambda` and `mu` are speciation and extinction, respectively.
#' @export
#' @references
#' DR rates (supplement, section 1.2.2): Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012). The global diversity of birds in space and time. Nature, 491(7424), 444–448. doi:10.1038/nature11631
#'
#' Interpreting tip rates: Title, P. O., & Rabosky, D. L. (2019). Tip rates, phylogenies and diversification: What are we estimating, and how good are the estimates? Methods in Ecology and Evolution, 10(6), 821–834. doi:10.1111/2041-210x.13153
#'
#' BAMM rates: Rabosky, D. L. (2014). Automatic Detection of Key Innovations, Rate Shifts, and Diversity-Dependence on Phylogenetic Trees. PLoS ONE, 9(2), e89543. doi:10.1371/journal.pone.0089543
#'
#' Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L., Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro, M. E. (2018). An inverse latitudinal gradient in speciation rate for marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1
#'
#' Enhanced polytomy resolution strengthens evidence for global gradient in speciation rate for marine fishes. \url{https://fishtreeoflife.org/rabosky-et-al-2018-update/}
#' @examples
#' \dontrun{
#' # Get cichlid rates and trees
#' rates <- fishtree_tip_rates(rank = "Cichlidae")
#' tree <- fishtree_phylogeny(rank = "Cichlidae")
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
#' }
fishtree_tip_rates <- function(species, rank, sampled_only = TRUE) {
  if (!rlang::is_missing(species) && !rlang::is_missing(rank)) rlang::abort("Must supply at most one of either `species` or `rank`, not both")

  rates <- .get("https://fishtreeoflife.org/downloads/tiprates.csv.xz", utils::read.csv, row.names = NULL)
  if (rlang::is_missing(species) && rlang::is_missing(rank)) return(rates)

  if (!rlang::is_missing(species)) {
    # Figure out which species are sampled
    tree <- fishtree_phylogeny()
    tips <- gsub("_", " ", tree$tip.label)
    if (sampled_only) wanted <- rates$species[rates$species %in% tips]
    else wanted <- rates$species

    names_to_get <- .name_check(species, valid_names = wanted)
    return(rates[rates$species %in% names_to_get, ])
  }

  if (!rlang::is_missing(rank)) {
    res <- .fetch_rank(rank)
    sampled_species <- res[[1]][[1]]$sampled_species
    species <- res[[1]][[1]]$species
    what <- res[[2]]
    if (sampled_only) wanted <- sampled_species
    else wanted <- species
    if (rlang::is_empty(wanted)) rlang::abort(paste("Can't get tip rates for", what, rank, "because there are not enough sampled species"))
    return(rates[rates$species %in% wanted, ])
  }
}







#' Get complete (stochastically-resolved) phylogenies from the Fish Tree of Life
#'
#' Retrieves a complete phylogeny generated by stochastic polytomy resolution via the Fish Tree of Life API. If neither `species` nor `rank` are specified, returns the entire phylogeny. See Rabosky et al. (2018) and Chang et al. (2019) for details on how these phylogenies were built using stochastic polytomy resolution. WARNING: These phylogenies should generally not be used for downstream analyses of trait evolution. See Rabosky (2015) for details.
#'
#' @inheritParams fishtree_phylogeny
#' @param mc.cores Number of cores to use in \link[parallel]{mclapply} when subsetting the tree (default `1`)
#' @return An object of class `"multiPhylo"` that should probably not be used for analyses of trait evolution, including (but not limited to) \link[ape]{pic}, \link[ape]{ace}, \link[ape]{corBrownian}, \link[diversitree]{make.bisse}, or \link[hisse]{hisse}.
#' @export
#' @references
#' Rabosky, D. L. (2015). No substitute for real data: A cautionary note on the use of phylogenies from birth-death polytomy resolvers for downstream comparative analyses. Evolution, 69(12), 3207–3216. doi:10.1111/evo.12817
#'
#' Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L., Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro, M. E. (2018). An inverse latitudinal gradient in speciation rate for marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1
#'
#' Chang, J., Rabosky, D. L., & Alfaro, M. E. (2019). Estimating diversification rates on incompletely-sampled phylogenies: theoretical concerns and practical solutions. Systematic Biology. doi:10.1093/sysbio/syz081
#'
#' Enhanced polytomy resolution strengthens evidence for global gradient in speciation rate for marine fishes. \url{https://fishtreeoflife.org/rabosky-et-al-2018-update/}
#' @examples
#' \dontrun{
#' tree <- fishtree_complete_phylogeny(rank = "Acanthuridae")
#' sampled_tips <- fishtree_phylogeny(rank = "Acanthuridae")$tip.label
#' all_tips <- tree[[1]]$tip.label
#' new_tips <- setdiff(all_tips, sampled_tips)
#' par(mfrow = c(2,2))
#' for (ii in 1:4) {
#'   plot(tree[[ii]], show.tip.label = FALSE, no.margin = TRUE)
#'   ape::tiplabels(pch = 19, col = ifelse(tree[[ii]]$tip.label %in% new_tips, "red", NA))
#' }
#' }
fishtree_complete_phylogeny <- function(species, rank, mc.cores = getOption("mc.cores", 1L)) {
  if (!rlang::is_missing(species) && !rlang::is_missing(rank)) rlang::abort("Must supply at most one of either `species` or `rank`, not both")

  trees <- .get("https://fishtreeoflife.org/downloads/actinopt_full.trees.xz", ape::read.tree)
  if (rlang::is_missing(species) && rlang::is_missing(rank)) return(trees)

  if (!rlang::is_missing(rank)) {
    res <- .fetch_rank(rank)
    if (length(res[[1]][[1]]$rogues) > 0) rlang::warn(paste(res[[2]], rank, "is not monophyletic; only including species in this taxon"))
    valid_spp <- .name_check(res[[1]][[1]]$species, trees[[1]]$tip.label)
  } else if (!rlang::is_missing(species)) {
    valid_spp <- .name_check(species, trees[[1]]$tip.label)
  }
  valid_spp <- gsub(" ", "_", valid_spp) # fix up tip names
  res <- parallel::mclapply(trees, ape::keep.tip, tip = valid_spp, mc.cores = mc.cores)
  class(res) <- "multiPhylo"
  res
}



#get information on how many genes are sampled by rank or species
fishtree_gene_sampling <- function(species, rank, matrix = FALSE) {

  if (!rlang::is_missing(species) && !rlang::is_missing(rank)) rlang::abort("Must supply at most one of either `species` or `rank`, not both")

  gene_name_list <- list('12s','16s', '4c4','coi','cytb','enc1','ficd','glyt','hoxc6a','kiaa1239','myh6','nd2','nd4','panx2','plag12','ptr','rag1','rag2','rhodopsin','ripk4','sh3px3','sidkey','sreb2','svep1','tbr1','vcpip','zic1')

  #return all gene info
  #  if (rlang::is_missing(species) && rlang::is_missing(rank)) {
  #    return(sampling)
  #}
  #do a name.check?

  #if species is given
  if (!rlang::is_missing(species))  rlang::abort("Haven't been able to program for a species input")
  #species_list <- .name_check(species)

  #get list of species from rank
  if (!rlang::is_missing(rank)) {
    rank_info <- .fetch_rank(rank)[[1]][[1]]

    #get character of species, split into list
    species_char <- rank_info$species
    species <-strsplit(species_char, split = "        ")

    gene_list <- rank_info$gene_sampling

  }

  if (!matrix) {
    return(gene_list)
  }

  if (matrix) {
    numcol <- length(species)
    gene_matrix <- matrix(nrow = 0,ncol = numcol)
    colnames(gene_matrix) <- species
    for (ii in gene_list) {
      #create a vector of true false for the gene and species
      new_row <- integer(0)
      #if each species is in the sublist for a gene, true, otherwise false
      for (xx in species){
        if (xx %in% ii) {
          new_row <- append(new_row, 1)
        }
        else {
          new_row <- append(new_row, 0)
        }
      }
      gene_matrix <- rbind(gene_matrix, new_row)
    }
    rownames(gene_matrix) <- gene_name_list
  }

}
