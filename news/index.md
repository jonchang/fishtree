# Changelog

## fishtree (development version)

- New
  [`fishtree_releases()`](https://fishtree.fishtreeoflife.org/reference/fishtree_releases.md)
  and
  [`fishtree_release()`](https://fishtree.fishtreeoflife.org/reference/fishtree_release.md)
  list the tree releases on the Fish Tree of Life and describe one: how
  to cite it, its license, and its files with their URLs, sizes and
  SHA-256 checksums.

- Every function has a new `release` argument, which defaults to
  `getOption("fishtree.release", "R2018")`. The default is fixed for
  each version of fishtree, so a new release on the website never
  changes your results. Looking up taxa by `rank` is only available for
  `"R2018"`.

- fishtree now finds its files through the data manifest at
  <https://fishtreeoflife.org/data/manifest.json>. Set
  `options(fishtree.manifest = ...)` to use another copy.

- Downloads are checked against their size and, on R 4.5.0 and later,
  their SHA-256 checksum in the manifest. They’re kept in a cache that
  lasts between sessions; see
  [`fishtree_cache_info()`](https://fishtree.fishtreeoflife.org/reference/fishtree_cache_info.md)
  and
  [`fishtree_cache_clear()`](https://fishtree.fishtreeoflife.org/reference/fishtree_cache_info.md).

- Downloads now wait at least 10 minutes before timing out.

- Taxa with spaces in their names (e.g., “Incertae sedis in Carangaria”)
  now work with
  [`fishtree_taxonomy()`](https://fishtree.fishtreeoflife.org/reference/fishtree_taxonomy.md)
  and all functions that take a `rank = ...` argument. Their download
  URLs are now percent-encoded.

## fishtree 0.3.4

CRAN release: 2021-01-30

- This release fixes some tests, and has no new user-facing features.

## fishtree 0.3.3

CRAN release: 2020-10-31

- The internal `.get()` function now avoids a build-time dependency on
  [`memoise::memoise()`](https://memoise.r-lib.org/reference/memoise.html),
  per <https://github.com/r-lib/memoise/issues/76>

- Vignettes and citations have been updated to use the latest versions,
  where relevant.

## fishtree 0.3.2

CRAN release: 2019-12-17

- [`fishtree_alignment()`](https://fishtree.fishtreeoflife.org/reference/fishtree_alignment.md)
  now works correctly with `species = ...` arguments.

- The internal `.name_check` function now emits structured warning
  output. This means that users supplying species names to e.g.,
  `fishtree_phylogeny` that do not have matches in the taxonomy may use
  functions that capture warning messages to programmatically determine
  which species are missing, without parsing error messages.

- Converted some warnings to messages.

- Minor documentation changes.

## fishtree 0.3.1

CRAN release: 2019-05-31

- `fishtree_tip_rates` now warns on missing / invalid species.

## fishtree 0.3.0

CRAN release: 2019-04-03

- New function: `fishtree_complete_phylogeny`, gets complete phylogenies
  where unsampled taxa were placed using stochastic polytomy resolution.

- Use package `memoise` instead of rolling our own caching
  functionality.

- New global option `fishtree.quiet`: when set to `FALSE`, be more
  verbose during downloads (can help with debugging issues).

- Some internal reorganization.

## fishtree 0.2.0

CRAN release: 2019-02-26

- New function: `fishtree_rogues`, identifies rogue/intruder taxa that
  break the monophyly of a selected group.

- `fishtree_phylogeny` now permits downloads of MRCA trees for
  paraphyletic groups. Suggested by an anonymous reviewer.

- `fishtree_taxonomy` revamped. It now takes a single argument
  `ranks = ...` where you can retrieve taxonomic information for any
  valid taxonomic rank in our taxonomy. Calling it without arguments
  returns a data frame of all valid taxa names. Suggested by an
  anonymous reviewer.

## fishtree 0.1.0

CRAN release: 2018-11-30

- Initial release.
