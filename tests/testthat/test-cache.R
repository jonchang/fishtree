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

test_that("a mismatched download over HTTP is retried once past caches", {
  local_fixture()
  good <- sub("^file://", "", fishtree:::.release_file("R2018", "chronogram")$url)
  file <- fishtree:::.release_file("R2018", "chronogram")
  file$url <- "https://example.org/downloads/chronogram.tre"
  urls <- character()
  local_mocked_bindings(.download = function(url, dest) {
    urls <<- c(urls, url)
    if (length(urls) == 1) writeLines("stale", dest) else file.copy(good, dest, overwrite = TRUE)
    invisible(dest)
  })
  path <- fishtree:::.fetch_file(file)
  expect_length(urls, 2)
  expect_match(urls[2], "^https://example.org/downloads/chronogram.tre\\?nocache=[0-9]+$")
  expect_equal(readLines(path), readLines(good))
})

test_that("cache busting keeps an existing query", {
  expect_match(fishtree:::.cache_bust("https://example.org/a"), "^https://example.org/a\\?nocache=[0-9]+$")
  expect_match(fishtree:::.cache_bust("https://example.org/a?b=1"), "^https://example.org/a\\?b=1&nocache=[0-9]+$")
})

test_that("files the manifest no longer lists are pruned", {
  local_fixture()
  stale <- file.path(getOption("fishtree.cache_dir"), strrep("0", 64))
  dir.create(dirname(stale), recursive = TRUE, showWarnings = FALSE)
  writeLines("old", stale)
  fishtree_phylogeny()
  expect_false(file.exists(stale))
  expect_equal(fishtree_cache_info()$file, fishtree_release()$files$sha256[1])
})
