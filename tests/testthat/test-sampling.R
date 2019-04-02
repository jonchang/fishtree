context("test-sampling")

test_that("phylogeny works", {
  spp <- c("Acanthurus achilles","Acanthurus auranticavus","Acanthurus bahianus", "Nonexistant")
  expect_equal(3, ape::Ntip(fishtree_phylogeny(species = spp)))
  expect_warning(fishtree_phylogeny(species = spp), "Missing names")
  expect_is(fishtree_phylogeny(species = spp), "phylo")
})


