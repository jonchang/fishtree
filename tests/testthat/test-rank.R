context("test-rank")

test_that("path encoding keeps slashes", {
  expect_equal(fishtree:::.encode_path("downloads/taxonomy/order/Incertae sedis in Carangaria.tre"),
               "downloads/taxonomy/order/Incertae%20sedis%20in%20Carangaria.tre")
})

test_that("taxa with spaces work", {
  skip_if_offline()
  expect_s3_class(suppressMessages(fishtree_phylogeny(taxon = "Incertae sedis in Carangaria")), "phylo")
})

test_that("rank and ranks still work, with a deprecation warning", {
  local_fixture()
  lifecycle::expect_deprecated(tree <- fishtree_phylogeny(rank = "Ausidae"))
  expect_equal(tree, fishtree_phylogeny(taxon = "Ausidae"))
  lifecycle::expect_deprecated(rogues <- fishtree_rogues(rank = "Ausidae"))
  expect_equal(rogues, fishtree_rogues("Ausidae"))
  lifecycle::expect_deprecated(dna <- fishtree_alignment(rank = "Ausidae"))
  expect_equal(dna, fishtree_alignment(taxon = "Ausidae"))
  lifecycle::expect_deprecated(rates <- fishtree_tip_rates(rank = "Ausidae"))
  expect_equal(rates, fishtree_tip_rates(taxon = "Ausidae"))
  lifecycle::expect_deprecated(trees <- fishtree_complete_phylogeny(rank = "Ausidae"))
  expect_equal(trees, fishtree_complete_phylogeny(taxon = "Ausidae"))
  lifecycle::expect_deprecated(tax <- fishtree_taxonomy(ranks = "Ausidae"))
  expect_equal(tax, fishtree_taxonomy("Ausidae"))
  # `rank =` partially matches `ranks`
  lifecycle::expect_deprecated(fishtree_taxonomy(rank = "Ausidae"))
})

test_that("fishtree_taxonomy returns the documented fields", {
  local_fixture()
  expect_named(fishtree_taxonomy(), c("rank", "name"))
  tax <- fishtree_taxonomy(c("Ausidae", "Busidae"))
  expect_named(tax, c("Ausidae", "Busidae"))
  for (taxon in tax) {
    expect_true(all(c("species", "sampled_species", "unsampled_species", "taxonomy") %in% names(taxon)))
  }
  expect_equal(tax$Ausidae$taxonomy$family, "Ausidae")
})
