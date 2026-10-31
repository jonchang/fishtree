test_that("downloads are cached by sha256", {
  local_fixture()
  fishtree_phylogeny()
  info <- fishtree_cache_info()
  expect_equal(info$file, fishtree_release()$files$sha256[1])
  expect_equal(fishtree_cache_clear(), 1)
  expect_equal(nrow(fishtree_cache_info()), 0)
})

test_that("a file that doesn't match the manifest is rejected", {
  site <- local_fixture_copy()
  local_fixture(site)
  cat("(Aus_bus:1,Aus_cus:1);\n", file = file.path(site, "downloads", "phylogram.tre"))
  expect_error(fishtree_phylogeny(type = "phylogram"), "doesn't match the size or checksum")
  expect_equal(nrow(fishtree_cache_info()), 0)

  skip_if(getRversion() < "4.5.0", "sha256 checks need R 4.5.0")
  tree <- readLines(file.path(site, "downloads", "chronogram.tre"))
  writeLines(sub("Cus_eus", "Cus_fus", tree), file.path(site, "downloads", "chronogram.tre"))
  expect_error(fishtree_phylogeny(), "doesn't match the size or checksum")
})
