# fishtree

The goal of `fishtree` is to provide an easy interface to the Fish Tree of Life API.

## Installation

You can install the released version of fishtree from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("fishtree")
```

## Example

``` r
library(fishtree)
library(ape)
phy <- fishtree_phylogeny("Acanthuridae")
branchin.times(phy)
```

