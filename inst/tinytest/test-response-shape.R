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

rm(x, n, nx2, nxK, groupBy)
