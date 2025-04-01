# test that logical arguments are identified and 'bounded'
expect_error(
  dbarts::dbartsControl(verbose = NA),
  "'verbose' must be TRUE/FALSE"
)
expect_error(
  dbarts::dbartsControl(verbose = "not-a-logical"),
  "'verbose' must be TRUE/FALSE"
)

expect_error(
  dbarts::dbartsControl(keepTrainingFits = NA),
  "'keepTrainingFits' must be TRUE/FALSE"
)
expect_error(
  dbarts::dbartsControl(keepTrainingFits = "not-a-logical"),
  "'keepTrainingFits' must be TRUE/FALSE"
)

expect_error(
  dbarts::dbartsControl(useQuantiles = NA),
  "'useQuantiles' must be TRUE/FALSE"
)
expect_error(
  dbarts::dbartsControl(useQuantiles = "not-a-logical"),
  "'useQuantiles' must be TRUE/FALSE"
)

expect_error(
  dbarts::dbartsControl(updateState = NA),
  "'updateState' must be TRUE/FALSE"
)
expect_error(
  dbarts::dbartsControl(updateState = "not-a-logical"),
  "'updateState' must be TRUE/FALSE"
)


# test that integer arguments are identified and bounded
expect_error(
  dbarts::dbartsControl(n.samples = "not-an-integer"),
  "'n.samples' must be coercible to type: integer"
)
expect_error(
  dbarts::dbartsControl(n.samples = -1L),
  "'n.samples' must be a non-negative integer"
)

expect_error(
  dbarts::dbartsControl(n.burn = "not-an-integer"),
  "'n.burn' must be coercible to type: integer"
)
expect_error(
  dbarts::dbartsControl(n.burn = NA_integer_),
  "'n.burn' must be a non-negative integer"
)
expect_error(
  dbarts::dbartsControl(n.burn = -1L),
  "'n.burn' must be a non-negative integer"
)

expect_error(
  dbarts::dbartsControl(n.trees = "not-an-integer"),
  "'n.trees' must be coercible to type: integer"
)
expect_error(
  dbarts::dbartsControl(n.trees = NA_integer_),
  "'n.trees' must be a positive integer"
)
expect_error(
  dbarts::dbartsControl(n.trees = 0L),
  "'n.trees' must be a positive integer"
)

expect_error(
  dbarts::dbartsControl(n.chains = "not-an-integer"),
  "'n.chains' must be coercible to type: integer"
)
expect_error(
  dbarts::dbartsControl(n.chains = NA_integer_),
  "'n.chains' must be a positive integer"
)
expect_error(
  dbarts::dbartsControl(n.chains = 0L),
  "'n.chains' must be a positive integer"
)

expect_error(
  dbarts::dbartsControl(n.threads = "not-an-integer"),
  "'n.threads' must be coercible to type: integer"
)
expect_error(
  dbarts::dbartsControl(n.threads = NA_integer_),
  "'n.threads' must be a positive integer"
)
expect_error(
  dbarts::dbartsControl(n.threads = 0L),
  "'n.threads' must be a positive integer"
)

expect_error(
  dbarts::dbartsControl(n.thin = "not-an-integer"),
  "'n.thin' must be coercible to type: integer"
)
expect_error(
  dbarts::dbartsControl(n.thin = NA_integer_),
  "'n.thin' must be a non-negative integer"
)
expect_error(
  dbarts::dbartsControl(n.thin = -1L),
  "'n.thin' must be a non-negative integer"
)

expect_error(
  dbarts::dbartsControl(printEvery = "not-an-integer"),
  "'printEvery' must be coercible to type: integer"
)
expect_error(
  dbarts::dbartsControl(printEvery = NA_integer_),
  "'printEvery' must be a non-negative integer"
)
expect_error(
  dbarts::dbartsControl(printEvery = -1L),
  "'printEvery' must be a non-negative integer"
)

expect_error(
  dbarts::dbartsControl(printCutoffs = "not-an-integer"),
  "'printCutoffs' must be coercible to type: integer"
)
expect_error(
  dbarts::dbartsControl(printCutoffs = NA_integer_),
  "'printCutoffs' must be a non-negative integer"
)
expect_error(
  dbarts::dbartsControl(printCutoffs = -1L),
  "'printCutoffs' must be a non-negative integer"
)

expect_error(
  dbarts::dbartsControl(n.cuts = "not-an-integer"),
  "'n.cuts' must be coercible to type: integer"
)
expect_error(
  dbarts::dbartsControl(n.cuts = NA_integer_),
  "'n.cuts' must contain only positive integers"
)
expect_error(
  dbarts::dbartsControl(n.cuts = -1L),
  "'n.cuts' must contain only positive integers"
)

expect_error(
  dbarts::dbartsControl(rngKind = "not-an-rng"),
  "unrecognized rng kind 'not-an-rng"
)
expect_error(
  dbarts::dbartsControl(rngNormalKind = "not-an-rng"),
  "unrecognized rng normal kind 'not-an-rng'"
)

