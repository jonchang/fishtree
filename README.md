[![Build Status](https://travis-ci.com/jonchang/fishtree.svg?token=CAAYReeKsDcnZM7jk2wY&branch=master)](https://travis-ci.com/jonchang/fishtree)
[![CRAN status](https://www.r-pkg.org/badges/version/fishtree)](https://cran.r-project.org/package=fishtree)

# fishtree

The goal of `fishtree` is to provide an easy interface to the Fish Tree of Life API, to download taxonomies, phylogenies, diversification rate information, and other data for ray-finned fishes.

## Example

Retrieve a phylogeny for the surgeonfishes and plot the phylogeny and lineage through time plot.

``` r
library(fishtree)
library(ape)
phy <- fishtree_phylogeny(rank = "Acanthuridae")

par(mfrow=c(2, 1))
plot(phy, show.tip.label = FALSE)
ltt.plot(phy)
```


## Installation

You can install the released version of fishtree from [CRAN](https://CRAN.R-project.org) with:

```r
install.packages("fishtree")
```

Alternatively, download the development version with devtools:

```r
devtools::install_github("jonchang/fishtree")
```

## References

The manuscript for this package is currently in review.
