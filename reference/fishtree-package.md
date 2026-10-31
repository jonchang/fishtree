# fishtree: Interface to the Fish Tree of Life API

fishtree makes it easy to access phylogenetic information about
ray-finned fishes with R. Download taxonomies, phylogenies, sequence
matrices, and diversification rate data with a simple set of functions.

## Details

Implementation note: this package makes calls over the network. It keeps
downloaded files in a cache on disk (see
[`fishtree_cache_info`](https://fishtree.fishtreeoflife.org/reference/fishtree_cache_info.md))
and their (sometimes large) parsed results in memory, so long-running R
instances could use a lot of memory.

## See also

Useful links:

- <https://fishtreeoflife.org/>

- <https://github.com/jonchang/fishtree>

- Report bugs at <https://github.com/jonchang/fishtree/issues>

## Author

**Maintainer**: Jonathan Chang <me@jonathanchang.org>
([ORCID](https://orcid.org/0000-0002-3811-1254))

Other contributors:

- Michael E Alfaro \[contributor\]

- Daniel L Rabosky \[contributor\]

- Stephen A Smith \[contributor\]
