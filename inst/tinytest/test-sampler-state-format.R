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
expect_equal(attr(state, "formatVersion"), 1L)

# anti-orphan: a FUTURE additive version still loads. The floor is >=, and the
# reader looks blocks up by name, so an unknown future block would just be
# ignored - an additive release never orphans an older reader's states.
future <- state
attr(future, "formatVersion") <- 2L
bartFit$fit$setState(future)
expect_equal(predict(bartFit, testData$x), preds)

# floor: an encoding BELOW the floor is refused, naming both versions. 0 is
# also what a state with no version attribute reads as.
old <- state
attr(old, "formatVersion") <- 0L
expect_error(
  bartFit$fit$setState(old),
  pattern = "encoding version 0.*oldest this dbarts \\(1\\)"
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

# leaf.scale: the newest OPTIONAL per-forest block (multiforest-mutation-gaps
# item 3). Pin its name, that an absent block still loads - the append-only
# contract, and what it now installs.
forest1 <- state[[1L]][["forests"]][[1L]]
expect_true("leaf.scale" %in% names(forest1))
expect_true(is.numeric(forest1[["leaf.scale"]]))
expect_equal(length(forest1[["leaf.scale"]]), 1L)

noLeafScale <- state
noLeafScale[[1L]][["forests"]][[1L]][["leaf.scale"]] <- NULL
bartFit$fit$setState(noLeafScale)
expect_equal(predict(bartFit, testData$x), preds)

# a present block is type/length checked and named when malformed
badLeafScale <- state
badLeafScale[[1L]][["forests"]][[1L]][["leaf.scale"]] <- c(1, 2)
expect_error(
  bartFit$fit$setState(badLeafScale),
  pattern = "block 'leaf.scale' is malformed"
)

# the leaf prior is mu ~ N(0, (scale / k)^2) and the state already restored k;
# it now restores the scale too, so a donor's node.scale survives a restore the
# way its k always has. CONSEQUENCE: a setModel(node.scale) issued after the
# last storeState() no longer survives a save/load re-creation.
control.ls <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  updateState = FALSE
)
makeSF <- function() {
  set.seed(5L)
  dbarts::dbarts(
    testData$x,
    testData$y,
    control = control.ls,
    resid.prior = dbarts::dbartsPriors$fixed(0.3)
  )
}
grabState <- function(s) {
  s$storeState()
  s$state
}

donor.sf <- makeSF()
donor.sf$model@node.scale <- 1.5
donor.sf$setModel(donor.sf$model)
invisible(donor.sf$run(25L, 5L))
state.sf <- grabState(donor.sf)
# nodeScale / sqrt(numTrees)
expect_equal(state.sf[[1L]][["forests"]][[1L]][["leaf.scale"]], 1.5 / 5)

dest.sf <- makeSF()
expect_equal(
  grabState(dest.sf)[[1L]][["forests"]][[1L]][["leaf.scale"]],
  0.5 / 5
)
dest.sf$setState(state.sf)
expect_equal(
  grabState(dest.sf)[[1L]][["forests"]][[1L]][["leaf.scale"]],
  1.5 / 5
)

# stripped, the destination keeps what it constructed
state.stripped <- state.sf
state.stripped[[1L]][["forests"]][[1L]][["leaf.scale"]] <- NULL
old.sf <- makeSF()
old.sf$setState(state.stripped)
expect_equal(
  grabState(old.sf)[[1L]][["forests"]][[1L]][["leaf.scale"]],
  0.5 / 5
)

# hostile VALUES match k's posture exactly - no new refusal class. A leaf scale
# is strictly positive, so non-finite and non-positive fall closed to the
# construction value; anything else flows through as k's does.
for (badValue in list(NaN, 0, -1)) {
  hostile <- state.sf
  hostile[[1L]][["forests"]][[1L]][["leaf.scale"]] <- badValue
  hostile.sf <- makeSF()
  hostile.sf$setState(hostile)
  expect_equal(
    grabState(hostile.sf)[[1L]][["forests"]][[1L]][["leaf.scale"]],
    0.5 / 5
  )
}

rm(
  forest1,
  noLeafScale,
  badLeafScale,
  control.ls,
  makeSF,
  grabState,
  donor.sf,
  state.sf,
  dest.sf,
  state.stripped,
  old.sf,
  badValue,
  hostile,
  hostile.sf
)
