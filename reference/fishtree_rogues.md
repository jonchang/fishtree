# Get rogue taxa that break the monophyly of defined taxa

For groups that were recovered as paraphyletic in the phylogenetic
analysis, uses the Fish Tree of Life API to identify which species
prevented that clade from being recovered as monophyletic.

## Usage

``` r
fishtree_rogues(rank)
```

## Arguments

- rank:

  the (possibly paraphyletic) rank for which rogue or intruder species
  should be identified.

## Value

A vector of species names, potentially empty.

## References

Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L.,
Friedman, M., Kashner, K., Garilao, C., Near, T. J., Coll, M., Alfaro,
M. E. (2018). An inverse latitudinal gradient in speciation rate for
marine fishes. Nature, 559(7714), 392–395. doi:10.1038/s41586-018-0273-1

## Examples

``` r
if (FALSE) { # \dontrun{
fishtree_rogues("Gobiidae")   # several rogue taxa!
fishtree_rogues("Labridae")   # nice and monophlyetic
} # }
```
