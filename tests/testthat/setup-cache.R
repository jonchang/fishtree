# Keeps downloads out of the user's cache during tests.
withr::local_options(
  fishtree.cache_dir = withr::local_tempdir(.local_envir = testthat::teardown_env()),
  .local_envir = testthat::teardown_env()
)
