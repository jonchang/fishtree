# Get taxonomies and other data from the Fish Tree of Life

Retrieves taxonomic and other information from the Fish Tree of Life
API.

## Usage

``` r
fishtree_taxonomy(ranks = NULL)
```

## Arguments

- ranks:

  One or more taxonomic ranks to retrieve.

## Value

A list, with components containing data on the specified taxa. If
\`ranks\` is unspecified, a data frame with all valid taxa is returned
instead.

## References

Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L.,
Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro,
M. E. (2018). An inverse latitudinal gradient in speciation rate for
marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1

## Examples

``` r
if (FALSE) { # \dontrun{
tax <- fishtree_taxonomy(rank = "Labridae")
n_total <- length(tax$Labridae$species)
n_sampl <- length(tax$Labridae$sampled_species)
paste("There are", n_sampl, "sampled species out of", n_total, "in wrasses.")
} # }
```
