# List releases on the Fish Tree of Life

Lists every release available on the Fish Tree of Life. Choose one with
the \`release\` argument of the other functions, or for the whole
session with \`options(fishtree.release = "...")\`.

## Usage

``` r
fishtree_releases()
```

## Value

A data frame with one row per release, with columns \`id\`, \`title\`,
\`status\`, \`native_taxonomy\`, \`doi\` and \`citation\`.

## See also

[`fishtree_release`](https://fishtree.fishtreeoflife.org/reference/fishtree_release.md)
for one release's citation and files.

## Examples

``` r
if (FALSE) { # \dontrun{
fishtree_releases()
} # }
```
