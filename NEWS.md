# fishtree (development version)

* The internal `.get()` function now avoids a build-time dependency on
  `memoise::memoise()`, per https://github.com/r-lib/memoise/issues/76

# fishtree 0.3.2

* `fishtree_alignment()` now works correctly with `species = ...` arguments.

* The internal `.name_check` function now emits structured warning output.
  This means that users supplying species names to e.g., `fishtree_phylogeny`
  that do not have matches in the taxonomy may use functions that capture
  warning messages to programatically determine which species are missing,
  without parsing error messages.

* Converted some warnings to messages.

* Minor documentation changes.

# fishtree 0.3.1

* `fishtree_tip_rates` now warns on missing / invalid species.

# fishtree 0.3.0

* New function: `fishtree_complete_phylogeny`, gets complete phylogenies
  where unsampled taxa were placed using stochastic polytomy resolution.

* Use package `memoise` instead of rolling our own caching functionality.

* New global option `fishtree.quiet`: when set to `FALSE`, be more verbose
  during downloads (can help with debugging issues).

* Some internal reorganization.

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

* Initial release.
