Sys.unsetenv("R_TESTS")
if (require(testthat, quietly = TRUE)) {
  require(dbarts)
  test_check("dbarts")
} else {
  cat("package 'testthat' not available; cannot run unit tests\n")
}

