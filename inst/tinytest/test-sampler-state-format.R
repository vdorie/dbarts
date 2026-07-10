# The additive-by-name state format (docs/design/public-surface.md 2,
# c-api-growth): setState reads per-chain blocks by name behind an encoding
# floor, so future additive versions still load, a genuinely older encoding is
# refused, and a missing REQUIRED block is named. States are opaque R lists
# with attributes, so these are pure attribute surgery (no C needed).

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

set.seed(99L)
bartFit <- dbarts::bart(
  testData$x,
  testData$y,
  ntree = 3L,
  ndpost = 7L,
  nskip = 0L,
  keeptrees = TRUE,
  verbose = FALSE
)
preds <- predict(bartFit, testData$x)

state <- bartFit$fit$state
expect_inherits(state, "bartcoreState")
expect_equal(attr(state, "formatVersion"), 3L)

# anti-orphan: a FUTURE additive version still loads. The floor is >=, and the
# reader looks blocks up by name, so an unknown future block would just be
# ignored - an additive release never orphans an older reader's states.
future <- state
attr(future, "formatVersion") <- 4L
bartFit$fit$setState(future)
expect_equal(predict(bartFit, testData$x), preds)

# floor: an encoding BELOW the floor is refused, naming both versions.
old <- state
attr(old, "formatVersion") <- 2L
expect_error(
  bartFit$fit$setState(old),
  pattern = "encoding version 2.*oldest this dbarts \\(3\\)"
)

# naming: a missing REQUIRED per-chain block is refused, naming the block.
missingSigma <- state
missingSigma[[1L]][["sigma"]] <- NULL
expect_error(
  bartFit$fit$setState(missingSigma),
  pattern = "missing required block 'sigma'"
)

# a REQUIRED block present but of the wrong type is named as malformed, not
# missing - the two-message convention.
badSigma <- state
badSigma[[1L]][["sigma"]] <- "not a number"
expect_error(
  bartFit$fit$setState(badSigma),
  pattern = "block 'sigma' is malformed"
)

# a missing REQUIRED forests block is likewise named.
missingForests <- state
missingForests[[1L]][["forests"]] <- NULL
expect_error(
  bartFit$fit$setState(missingForests),
  pattern = "missing required block 'forests'"
)

# default: removing an OPTIONAL block still loads. rng.state absence only
# forfeits bitwise continuation; prediction from the saved trees is unaffected.
noRng <- state
noRng[[1L]][["rng.state"]] <- NULL
bartFit$fit$setState(noRng)
expect_equal(predict(bartFit, testData$x), preds)
