## Test environments

* local macOS 27.2, R 4.6.2
* GitHub Actions, Ubuntu 24.04, R release and R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## This release

fishtree now finds its files through a data manifest published on
fishtreeoflife.org, checks downloads against their checksums, and keeps them in
a cache in `tools::R_user_dir("fishtree", "cache")`, which
`fishtree_cache_clear()` empties. The `rank` argument is deprecated in favor of
`taxon`, using lifecycle.

## Downstream dependencies

I checked the 4 reverse dependencies (aptg, FishPhyloMaker, phylosem,
prepR4pcm) with `tools::check_packages_in_dir()`. None has new problems:

* FishPhyloMaker has the NOTE about unused Imports that it has on CRAN. It
  passes `rank` to `fishtree_taxonomy()`, which still works with a
  deprecation warning.
* phylosem's vignette that uses fishtree builds. Its other vignette needs the
  Bioconductor package graph, which wasn't installed.
