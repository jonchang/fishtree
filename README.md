
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.com/jonchang/fishtree.svg?token=CAAYReeKsDcnZM7jk2wY&branch=master)](https://travis-ci.com/jonchang/fishtree)
[![CRAN
status](https://www.r-pkg.org/badges/version/fishtree)](https://cran.r-project.org/package=fishtree)

# fishtree

The goal of `fishtree` is to provide an easy interface in R to the Fish
Tree of Life API, to download taxonomies, phylogenies, diversification
rate information, and other data for ray-finned fishes. It should
seamlessly integrate with the rest of the R ecosystem, especially the
package `ape`.

## Example

See a list of taxa available to download.

``` r
library(fishtree)

tax <- fishtree_taxonomy()
head(tax)
#>         rank        name
#> 1      class Actinopteri
#> 2      class   Cladistia
#> 3   subclass Chondrostei
#> 4   subclass Neopterygii
#> 5 infraclass    Holostei
#> 6 infraclass   Teleostei
```

Retrieve a phylogeny for the surgeonfishes and plot the phylogeny and
lineage through time plot.

``` r
library(ape)
phy <- fishtree_phylogeny(rank = "Acanthuridae")
phy
#> 
#> Phylogenetic tree with 67 tips and 66 internal nodes.
#> 
#> Tip labels:
#>  Acanthurus_mata, Acanthurus_blochii, Acanthurus_xanthopterus, Acanthurus_bariene, Acanthurus_dussumieri, Acanthurus_leucocheilus, ...
#> 
#> Rooted; includes branch lengths.
```

``` r
par(mfrow=c(2, 1))
plot(phy, show.tip.label = FALSE)
ltt.plot(phy)
```

You can also check out the vignettes for more detailed examples.

## Installation

You can install the released version of fishtree from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("fishtree")
```

Alternatively, download the development version with devtools:

``` r
devtools::install_github("jonchang/fishtree")
```

## References

To cite `fishtree` in publications use:

Chang, J., Rabosky, D. L., Smith, S. A., & Alfaro, M. E. (2019). An R
package and online resource for macroevolutionary studies using the
ray‐finned fish tree of life. Methods in Ecology and Evolution. doi:
[10.1111/2041-210x.13182](https://doi.org/10.1111/2041-210x.13182)

The primary data source for `fishtree` was published as a part of:

Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L.,
Friedman, M., … Alfaro, M. E. (2018). An inverse latitudinal gradient in
speciation rate for marine fishes. Nature, 559(7714), 392–395. doi:
[10.1038/s41586-018-0273-1](https://doi.org/10.1038/s41586-018-0273-1)

## License

The `fishtree` package is licensed under a [2-clause BSD
license](https://opensource.org/licenses/BSD-2-Clause).

## Sponsorship

Please consider sponsoring the maintenance of `fishtree` via [GitHub
Sponsors](https://github.com/sponsors/jonchang).

## Releasing

    withr::with_envvar(c("NOT_CRAN" = "true"), devtools::release())
