context("test-rank")

test_that("path encoding keeps slashes", {
  expect_equal(fishtree:::.encode_path("downloads/taxonomy/order/Incertae sedis in Carangaria.tre"),
               "downloads/taxonomy/order/Incertae%20sedis%20in%20Carangaria.tre")
})

test_that("rank with spaces works", {
  skip_if_offline()
  expect_s3_class(suppressMessages(fishtree_phylogeny(rank = "Incertae sedis in Carangaria")), "phylo")
})
