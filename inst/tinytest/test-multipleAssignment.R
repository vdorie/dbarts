massign <- dbarts:::massign
"[<-.lval" <- dbarts:::"[<-.lval"

source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

# test that multiple assignment works with missing arguments
rh <- c(2, 5)
massign[a, b] <- rh
expect_equal(a, 2)
expect_equal(b, 5)
rm(a, b)

massign[a, ] <- rh
expect_equal(a, 2)
expect_error(b, "object 'b' not found")
rm(a)
massign[, b] <- rh
expect_equal(b, 5)
expect_error(a, "object 'a' not found")
rm(b)

warnings.unnamedRhs <- captureWarnings(massign[a = b, ] <- rh)
expect_equal(length(warnings.unnamedRhs), 1L)
expect_match(
  conditionMessage(warnings.unnamedRhs[[1L]]),
  "right-hand-side of assignment is unnamed; using position only"
)
expect_inherits(warnings.unnamedRhs[[1L]], "dbartsPositionalArgsWarning")
rm(a)

rm(rh)

# test that multiple assignment works with named arguments
rh <- c(a = 2, b = 5)
massign[a, c = b] <- rh
expect_equal(a, 2)
expect_equal(c, 5)
expect_error(b, "object 'b' not found")
rm(a, c)

massign[c = a, ] <- rh
expect_equal(c, 2)

expect_error(
  massign[c = d, ] <- rh,
  "'d' not present in right-hand-side of assignment"
)

massign[c = a, d = a] <- rh
expect_equal(c, 2)
expect_equal(d, 2)

warnings.dupLhs <- captureWarnings(massign[c = a, c = b] <- rh)
expect_equal(length(warnings.dupLhs), 1L)
expect_match(
  conditionMessage(warnings.dupLhs[[1L]]),
  "names on left-hand-side of assignment appear more than once: c"
)
expect_inherits(warnings.dupLhs[[1L]], "dbartsDuplicateNameWarning")
rm(c)

rh <- c(a = 2, a = 5)
warnings.dupRhs <- captureWarnings(massign[b = a, c] <- rh)
expect_equal(length(warnings.dupRhs), 1L)
expect_match(
  conditionMessage(warnings.dupRhs[[1L]]),
  "'a' present multiple times in right-hand-side of assignment"
)
expect_inherits(warnings.dupRhs[[1L]], "dbartsDuplicateNameWarning")
expect_equal(b, 2)
expect_equal(c, 5)
rm(b, c)

rm(rh)

rm(list = c("[<-.lval", "massign"))
