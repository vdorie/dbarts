if (requireNamespace("tinytest", quietly = TRUE)) {
  tinytest::test_package("dbarts")
} else if (nzchar(Sys.getenv("CI", ""))) {
  stop("'tinytest' is not available: install it to run the test suite under CI")
} else {
  cat("package 'tinytest' not available; cannot run unit tests\n")
}
