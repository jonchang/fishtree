# Get a phylogeny from the Fish Tree of Life

Retrieves a phylogeny via the Fish Tree of Life API. If neither
\`species\` nor \`rank\` are specified, returns the entire phylogeny.

## Usage

``` r
fishtree_phylogeny(
  species,
  rank,
  type = c("chronogram", "phylogram", "chronogram_mrca", "phylogram_mrca")
)
```

## Arguments

- species:

  (Optionally) subset the results based on a vector of species names.

- rank:

  (Optionally) subset the results based on the supplied taxonomic rank.

- type:

  Either \`"chronogram"\` or \`"phylogram"\`. A chronogram has branch
  lengths proportional to units of time, while a phylogram has branch
  lengths proportional to the amount of character change. When
  retrieving a phylogeny by rank, and that rank is not recovered as
  monophyletic, acceptable types also include \`"chronogram_mrca"\` and
  \`"phylogram_mrca"\`, which returns a tree with \*all\* species
  descending from the common ancestor of species in the specified rank.

## Value

An object of class \`"phylo"\`.

## Details

For maximum interoperability, \`species\` considers spaces and
underscores equivalently. Internally, the phylogenies use underscores.

## References

Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L.,
Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro,
M. E. (2018). An inverse latitudinal gradient in speciation rate for
marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1

## See also

[`fishtree_rogues`](https://fishtree.fishtreeoflife.org/reference/fishtree_rogues.md),
[`read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html),
[`force.ultrametric`](https://rdrr.io/pkg/phytools/man/force.ultrametric.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Get a phylogeny for a taxonomic rank
surgeons <- fishtree_phylogeny(rank = "Acanthuridae")

# Get a phylogeny for only certain species
genomic_fish <- c("Oryzias latipes", "Tetraodon nigroviridis",
                  "Gasterosteus aculeatus", "Danio rerio")
fishtree_phylogeny(species = genomic_fish)

# Chronograms may not be ultrametric due to numerical precision issues
# Consider using phytools::force.ultrametric
ape::is.ultrametric(surgeons)
ape::is.ultrametric(surgeons, tol = 0.00001)

# Difference between MRCA trees and regular trees
gobies_mrca <- fishtree_phylogeny(rank = "Gobiidae", type = "chronogram_mrca")
gobies <- fishtree_phylogeny(rank = "Gobiidae", type = "chronogram")
# MRCA trees will have more tips for non-monophyletic groups
length(gobies_mrca$tip.label) > length(gobies$tip.label)
# Drop rogue tips in the MRCA tree
rogue_gobies <- fishtree_rogues("Gobiidae")
pruned_gobies <- ape::drop.tip(gobies_mrca, rogue_gobies)
# Now the trees are identical
setequal(gobies$tip.label, pruned_gobies$tip.label)
} # }
```
