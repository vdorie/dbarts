if (requireNamespace("tinytest", quietly = TRUE)) {
  tinytest::test_package("dbarts")
} else {
  cat("package 'tinytest' not available; cannot run unit tests\n")
}
