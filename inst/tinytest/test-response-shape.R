source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# A multi-column response reaching codeResponse() (a Surv object, an n x 2
# matrix, or a wider n x K one) used to flatten column-major with as.double()
# before its row count was compared against 'x', producing a false "'x' must
# have the same number of observations as 'y'". Each shape now names itself
# and the route that does accept it, before codeResponse ever runs.
x <- testData$x
n <- nrow(x)

nx2 <- cbind(a = testData$y, b = testData$y + 1)
expect_error(
  dbarts::dbartsData(x, nx2),
  "'y' is an n x 2 matrix"
)
expect_error(
  dbarts::dbartsData(x, nx2),
  "family = \"aft\"/\"hazard\""
)

nxK <- cbind(a = testData$y, b = testData$y + 1, c = testData$y + 2)
expect_error(
  dbarts::dbartsData(x, nxK),
  "'y' is an n x 3 matrix"
)
expect_error(
  dbarts::dbartsData(x, nxK),
  "family = \"multinomial\""
)

if (requireNamespace("survival", quietly = TRUE)) {
  survY <- survival::Surv(abs(testData$y), rep(c(0L, 1L), length.out = n))

  expect_error(
    dbarts::dbartsData(x, survY),
    "'y' is a survival response \\(Surv\\)"
  )

  # of the inheriting surfaces, dbarts() itself is not one for a Surv
  # response specifically - family = "auto" auto-dispatches it to "aft"
  # before dbartsData's ingestion is ever reached, so it legitimately fits
  # rather than refusing; xbart has no such special case
  expect_error(dbarts::xbart(x, survY), "'y' is a survival response \\(Surv\\)")
  groupBy <- factor(rep(seq_len(4L), length.out = n))
  expect_error(
    dbarts::rbart_vi(x, survY, group.by = groupBy),
    "'y' is a survival response \\(Surv\\)"
  )
}

expect_error(dbarts::dbarts(x, nx2), "'y' is an n x 2 matrix")
expect_error(dbarts::xbart(x, nx2), "'y' is an n x 2 matrix")
groupBy <- factor(rep(seq_len(4L), length.out = n))
expect_error(
  dbarts::rbart_vi(x, nx2, group.by = groupBy),
  "'y' is an n x 2 matrix"
)

# --- the formula interface reaches the same guards ---------------------------
# A cbind(...) left-hand side is a multi-column response like any other. The
# aft routes rewrite a (time, status) left-hand side to a plain numeric BEFORE
# the data object is built, and bart2's multinomial cbind() route maps it to
# 'counts', so only a caller reaching dbartsData directly lands here.
frm <- data.frame(
  a = testData$y,
  b = testData$y + 1,
  c = testData$y + 2,
  p = x[, 1L],
  q = x[, 2L]
)
expect_error(dbarts::dbartsData(cbind(a, b) ~ p, frm), "'y' is an n x 2 matrix")
expect_error(
  dbarts::dbartsData(cbind(a, b, c) ~ p, frm),
  "'y' is an n x 3 matrix"
)

# a matrix offset is a per-category shift wherever it arrives, so the formula
# interface refuses it by shape too, not by a length mismatch
expect_error(
  dbarts::dbartsData(a ~ p, frm, offset = matrix(0.0, n, 3L)),
  "which only family = \"multinomial\" accepts"
)
# its test twin already shared the matrix branches' set-aside
expect_error(
  dbarts::dbartsData(a ~ p, frm, test = frm, offset.test = matrix(0.0, n, 3L)),
  "'offset.test' must be a numeric vector"
)

# the aft route's own bare (time, status) left-hand side is rewritten ahead of
# the guard, so it still fits (bart2's multinomial cbind() count response is
# pinned in test-multinomial-surface.R)
aftFrame <- data.frame(
  time = exp(abs(testData$y) + 1),
  status = rep_len(c(1L, 0L), n),
  p = x[, 1L],
  q = x[, 2L]
)
aftFit <- dbarts::rbart_vi(
  cbind(time, status) ~ p + q,
  data = aftFrame,
  group.by = groupBy,
  family = "aft",
  n.samples = 5L,
  n.burn = 2L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 5L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(aftFit$family, "aft")

rm(x, n, nx2, nxK, groupBy, frm, aftFrame, aftFit)
