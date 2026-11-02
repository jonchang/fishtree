# fishtree 0.4.0

* `taxon` replaces `rank` as the argument that names a taxon, e.g.,
  `fishtree_phylogeny(taxon = "Acanthuridae")`, and replaces `ranks` in
  `fishtree_taxonomy()`. `rank` and `ranks` still work, with a warning, and
  will be removed in fishtree 1.0.

* `fishtree_taxonomy()` documents the fields that future versions will keep:
  `species`, `sampled_species`, `unsampled_species` and `taxonomy` for each
  taxon, and the `rank` and `name` columns of the table of taxa. Its other
  fields are deprecated and will be removed in fishtree 1.0. Look up ranks in
  `taxonomy` by name, not position.

* New `fishtree_releases()` and `fishtree_release()` list the tree releases on
  the Fish Tree of Life and describe one: how to cite it, its license, and its
  files with their URLs, sizes and SHA-256 checksums.

* Every function has a new `release` argument, which defaults to
  `getOption("fishtree.release", "R2018")`. The default is fixed for each
  version of fishtree, so a new release on the website never changes your
  results. Looking up a `taxon` is only available for `"R2018"`.

* fishtree now finds its files through the data manifest at
  https://fishtreeoflife.org/data/manifest.json. Set
  `options(fishtree.manifest = ...)` to use another copy.

* Downloads are checked against their size and, on R 4.5.0 and later, their
  SHA-256 checksum in the manifest. They're kept in a cache that lasts between
  sessions; see `fishtree_cache_info()` and `fishtree_cache_clear()`.
  Cached files that the manifest no longer lists are deleted whenever fishtree
  downloads a file.

* Downloads now wait at least 10 minutes before timing out.

* Taxa with spaces in their names (e.g., "Incertae sedis in Carangaria")
  now work with `fishtree_taxonomy()` and all functions that take a
  `taxon = ...` argument.

# fishtree 0.3.4

* This release fixes some tests, and has no new user-facing features.

# fishtree 0.3.3

* The internal `.get()` function now avoids a build-time dependency on
  `memoise::memoise()`, per https://github.com/r-lib/memoise/issues/76

* Vignettes and citations have been updated to use the latest versions,
  where relevant.

# fishtree 0.3.2

* `fishtree_alignment()` now works correctly with `species = ...` arguments.

* The internal `.name_check` function now emits structured warning output.
  This means that users supplying species names to e.g., `fishtree_phylogeny`
  that do not have matches in the taxonomy may use functions that capture
  warning messages to programmatically determine which species are missing,
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
