# Get taxonomies and other data from the Fish Tree of Life

Retrieves taxonomic and other information from the Fish Tree of Life
API.

## Usage

``` r
fishtree_taxonomy(
  taxon = NULL,
  release = getOption("fishtree.release", .default_release),
  ranks = deprecated()
)
```

## Arguments

- taxon:

  One or more taxon names to retrieve, such as \`"Labridae"\`.

- release:

  The tree release to use, from
  [`fishtree_releases`](https://fishtree.fishtreeoflife.org/reference/fishtree_releases.md).
  Defaults to \`getOption("fishtree.release")\`, or \`"R2018"\` if that
  is unset. Looking up a \`taxon\` is only available for \`"R2018"\`.

- ranks:

  Deprecated; use \`taxon\` instead.

## Value

If \`taxon\` is unspecified, a data frame of all valid taxa, with
columns \`rank\` and \`name\`. Otherwise, a list with one element per
taxon, named by taxon. Each element is a list with:

- \`species\`:

  All species in the taxon.

- \`sampled_species\`:

  The species in the taxon that are in the phylogeny.

- \`unsampled_species\`:

  The species in the taxon that are not in the phylogeny.

- \`taxonomy\`:

  A list, named by rank, of the names at each rank among the taxon's
  species. Look up ranks by name, as in \`taxonomy\$order\`; which ranks
  are present, and their order, depend on the taxonomy.

Future versions of fishtree will keep these fields and the columns of
the data frame. The other fields are deprecated and will be removed in
fishtree 1.0: use
[`fishtree_rogues`](https://fishtree.fishtreeoflife.org/reference/fishtree_rogues.md)
instead of \`rogues\`;
[`fishtree_tip_rates`](https://fishtree.fishtreeoflife.org/reference/fishtree_tip_rates.md)
instead of \`tiprates_dr\`, \`tiprates_bamm_lambda\`,
\`tiprates_bamm_mu\`, \`mean_dr\`, \`mean_bamm_lambda\` and
\`mean_bamm_mu\`; and
[`fishtree_phylogeny`](https://fishtree.fishtreeoflife.org/reference/fishtree_phylogeny.md)
and
[`fishtree_alignment`](https://fishtree.fishtreeoflife.org/reference/fishtree_alignment.md)
instead of the download paths \`chronogram\`, \`phylogram\`,
\`chronogram_mrca\`, \`phylogram_mrca\`, \`matrix_phylip\` and
\`matrix_nexus\`. \`gene_sampling\` has no replacement.

## References

Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L.,
Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro,
M. E. (2018). An inverse latitudinal gradient in speciation rate for
marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1

## Examples

``` r
if (FALSE) { # \dontrun{
tax <- fishtree_taxonomy("Labridae")
n_total <- length(tax$Labridae$species)
n_sampl <- length(tax$Labridae$sampled_species)
paste("There are", n_sampl, "sampled species out of", n_total, "in wrasses.")
} # }
```
