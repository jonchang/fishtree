# Points at the offline fixture site in fixtures/site, with an empty cache.
local_fixture <- function(site = test_path("fixtures", "site"), env = parent.frame()) {
  withr::local_options(
    fishtree.manifest = file.path(site, "data", "manifest.json"),
    fishtree.cache_dir = withr::local_tempdir(.local_envir = env),
    fishtree.release = NULL,
    .local_envir = env
  )
}

# A copy of the fixture site that a test can change.
local_fixture_copy <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  file.copy(test_path("fixtures", "site"), dir, recursive = TRUE)
  file.path(dir, "site")
}
