# Get the details of a Fish Tree of Life release

Describes one tree release. Includes citation information and checksums.

## Usage

``` r
fishtree_release(release = getOption("fishtree.release", .default_release))
```

## Arguments

- release:

  The release's id, from
  [`fishtree_releases`](https://fishtree.fishtreeoflife.org/reference/fishtree_releases.md).
  Defaults to \`getOption("fishtree.release")\`, or \`"R2018"\` if that
  is unset.

## Value

A list with the release's \`id\`, \`title\`, \`status\`, \`citation\`,
\`doi\`, \`license\` and \`native_taxonomy\`, and \`files\`, a data
frame of its files.

## See also

[`fishtree_releases`](https://fishtree.fishtreeoflife.org/reference/fishtree_releases.md)
to list releases.

## Examples

``` r
if (FALSE) { # \dontrun{
rel <- fishtree_release("R2018")
rel$citation
rel$files[, c("role", "url", "bytes")]
} # }
```
