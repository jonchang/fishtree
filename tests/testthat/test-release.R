test_that("relative URLs resolve like a browser's", {
  base <- "https://fishtreeoflife.org/data/manifest.json"
  expect_equal(fishtree:::.resolve_url("../downloads/a.tre.xz", base), "https://fishtreeoflife.org/downloads/a.tre.xz")
  expect_equal(fishtree:::.resolve_url("../api/taxonomy/", base), "https://fishtreeoflife.org/api/taxonomy/")
  expect_equal(fishtree:::.resolve_url("../../", "https://fishtreeoflife.org/api/taxonomy/"), "https://fishtreeoflife.org/")
  expect_equal(fishtree:::.resolve_url("https://zenodo.org/x", base), "https://zenodo.org/x")
  expect_equal(fishtree:::.resolve_url("../downloads/a", "file:///tmp/site/data/manifest.json"), "file:///tmp/site/downloads/a")
})

test_that("releases are listed from the manifest", {
  local_fixture()
  releases <- fishtree_releases()
  expect_equal(releases$id, c("R2018", "X2030"))
  expect_equal(releases$status, c("stable", "preview"))
  expect_equal(releases$doi, c("10.1234/test.2018", NA))
})

test_that("a release describes its files", {
  local_fixture()
  rel <- fishtree_release()
  expect_equal(rel$id, "R2018")
  expect_equal(rel$license, "CC0-1.0")
  expect_equal(rel$files$role, c("chronogram", "phylogram", "complete_trees", "tip_rates", "alignment", "partitions"))
  expect_true(all(startsWith(rel$files$url, "file://")))
  expect_true(all(file.exists(sub("^file://", "", rel$files$url))))
  expect_error(fishtree_release("R1999"), "Available releases: R2018, X2030")
})

test_that("whole-release files are read through the manifest", {
  local_fixture()
  expect_equal(ape::Ntip(fishtree_phylogeny()), 4)
  expect_equal(fishtree_phylogeny(type = "phylogram")$edge.length[1], 0.05)
  expect_equal(ape::Ntip(fishtree_phylogeny(species = c("Aus bus", "Bus_dus"))), 2)
  expect_equal(nrow(fishtree_tip_rates()), 5)
  expect_warning(rates <- fishtree_tip_rates(species = c("Aus bus", "Aus zus")), "only found 1 species")
  expect_equal(rates$species, "Aus bus")
  expect_equal(dim(fishtree_alignment()), c(4, 12))
  expect_equal(names(fishtree_alignment(split = TRUE)), c("gene1", "gene2"))
  expect_length(fishtree_complete_phylogeny(), 2)
  expect_equal(ape::Ntip(fishtree_complete_phylogeny(species = c("Aus bus", "Aus zus"))[[1]]), 2)
})

test_that("release can be chosen per call or per session", {
  local_fixture()
  expect_equal(ape::Ntip(fishtree_phylogeny(release = "X2030")), 5)
  withr::local_options(fishtree.release = "X2030")
  expect_equal(ape::Ntip(fishtree_phylogeny()), 5)
  expect_error(fishtree_tip_rates(), "Release X2030 has no tip_rates file")
})

test_that("taxon lookups use the release's legacy API", {
  local_fixture()
  expect_equal(fishtree_taxonomy()$name, c("Ausidae", "Busidae"))
  expect_equal(fishtree_taxonomy("Ausidae")$Ausidae$sampled_species, c("Aus bus", "Aus cus"))
  expect_equal(fishtree_phylogeny(taxon = "Ausidae")$tip.label, c("Aus_bus", "Aus_cus"))
  expect_equal(fishtree_rogues("Ausidae"), character())
  expect_equal(dim(fishtree_alignment(taxon = "Ausidae")), c(2, 12))
  expect_equal(fishtree_tip_rates(taxon = "Ausidae")$species, c("Aus bus", "Aus cus"))
  expect_error(fishtree_phylogeny(taxon = "Ausidae", release = "X2030"), "isn't available for release X2030")
})

test_that("an unknown manifest schema version is refused", {
  site <- local_fixture_copy()
  local_fixture(site)
  path <- file.path(site, "data", "manifest.json")
  writeLines(sub('"schema_version": 1', '"schema_version": 2', readLines(path)), path)
  expect_error(fishtree_releases(), "Please update fishtree")
})
