# Manage the Fish Tree of Life download cache

Files from a release are downloaded once, checked against the release's
checksums, and kept in a cache. \`fishtree_cache_info()\` lists the
cached files and \`fishtree_cache_clear()\` deletes them.

## Usage

``` r
fishtree_cache_info()

fishtree_cache_clear()
```

## Value

\`fishtree_cache_info()\` returns a data frame with the \`file\` name
(its SHA-256), \`bytes\` and \`modified\` time of each cached file.
\`fishtree_cache_clear()\` returns the number of files deleted,
invisibly.

## Details

The cache is in \`tools::R_user_dir("fishtree", "cache")\`. Set
\`options(fishtree.cache_dir = ...)\` to use another directory.

## Examples

``` r
if (FALSE) { # \dontrun{
fishtree_cache_info()
fishtree_cache_clear()
} # }
```
