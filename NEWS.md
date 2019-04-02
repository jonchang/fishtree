# fishtree 0.3.0

# fishtree (development version)

* New function: `fishtree_complete_phylogeny`, gets complete phylogenies
  where unsampled taxa were placed using stochastic polytomy resolution.

* Use package `memoise` instead of rolling our own caching functionality.

* New global option `fishtree.quiet`: when set to `FALSE`, be more verbose
  during downloads (can help with debugging issues).

* Some internal reorganization

# fishtree 0.2.0

* New function: `fishtree_rogues`, identifies rogue/intruder taxa that break
  the monophyly of a selected group.

* `fishtree_phylogeny` now permits downloads of MRCA trees for paraphyletic
  groups. Suggested by an anonymous reviewer.

* `fishtree_taxonomy` revamped. It now takes a single argument `ranks = ...`
  where you can retrieve taxonomic information for any valid taxonomic rank
  in our taxonomy. Calling it without arguments returns a data frame of all
  valid taxa names. Suggested by an anonymous reviewer.

# fishtree 0.1.0

Initial release.
