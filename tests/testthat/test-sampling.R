context("test-sampling")

test_that("phylogeny works", {
  spp <- c("Acanthurus achilles","Acanthurus auranticavus","Acanthurus bahianus")
  expect_equal(3, ape::Ntip(fishtree_phylogeny(species = spp)))
  expect_is(fishtree_phylogeny(species = spp), "phylo")
})

test_that("warn on extra species", {
  spp <- c("Acanthurus achilles","Acanthurus auranticavus","Acanthurus bahianus", "Nonexistant")
  expect_warning(fishtree_phylogeny(species = spp), "only found 3 species")
})


test_that("alignment single species works", {
    skip_on_cran() # Takes a long time
    spp <- "Abalistes_stellaris"
    expect_s3_class(fishtree_alignment(species = spp), "DNAbin")
})
