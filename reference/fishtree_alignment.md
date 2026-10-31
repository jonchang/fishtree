# Get aligned sequences from the Fish Tree of Life

Retrieves an aligned sequence via the Fish Tree of Life API. If neither
\`species\` nor \`taxon\` are specified, returns the entire sequence
matrix.

## Usage

``` r
fishtree_alignment(
  species,
  taxon,
  split = FALSE,
  release = getOption("fishtree.release", .default_release),
  rank = deprecated()
)
```

## Arguments

- species:

  (Optionally) subset the results based on a vector of species names.

- taxon:

  (Optionally) subset the results to the taxon with this name, such as a
  family or order (e.g., \`"Acanthuridae"\`).

- split:

  Splits the output into a list by gene locus.

- release:

  The tree release to use, from
  [`fishtree_releases`](https://fishtree.fishtreeoflife.org/reference/fishtree_releases.md).
  Defaults to \`getOption("fishtree.release")\`, or \`"R2018"\` if that
  is unset. Looking up a \`taxon\` is only available for \`"R2018"\`.

- rank:

  Deprecated; use \`taxon\` instead.

## Value

An object of class \`"DNAbin"\`, or a named list of the same if \`split
= TRUE“

## References

Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L.,
Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro,
M. E. (2018). An inverse latitudinal gradient in speciation rate for
marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1

## See also

[DNAbin](https://rdrr.io/pkg/ape/man/DNAbin.html)

## Examples

``` r
if (FALSE) { # \dontrun{
surgeon_dna <- fishtree_alignment(taxon = "Acanthuridae", split = TRUE)
surgeon_dna[[1]]
par(mfrow = c(9, 3), mar = c(0.5, 0.5, 1, 0.5), xaxt = "n", yaxt = "n")
for (gene in names(surgeon_dna)) {
  image(surgeon_dna[[gene]], legend = FALSE, show.labels = FALSE)
  title(gene)
}
} # }
```
