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

# floor: an encoding BELOW the floor is refused, naming both versions. 0 is
# also what a state with no version attribute reads as.
old <- state
attr(old, "formatVersion") <- 0L
expect_error(
  bartFit$fit$setState(old),
  pattern = "encoding version 0.*oldest this dbarts \\(3\\)"
)

# the floor is what makes a BLOCK RENAME safe. Version 1 spelled the amplitude
# glue "bcf"; were such a state read here it would find "glue" absent, default
# it as an optional block, and leave the amplitudes at their construction
# values - a wrong answer, not an error. The version check runs BEFORE any
# block is read, so a state still carrying the old name is refused by version.
priorEncoding <- state
priorEncoding[[1L]][["bcf"]] <- c(1, 1, 1, 1)
attr(priorEncoding, "formatVersion") <- 1L
expect_error(
  bartFit$fit$setState(priorEncoding),
  pattern = "encoding version 1.*oldest this dbarts \\(3\\)"
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

# leaf.scale: the newest OPTIONAL per-forest block. Pin its name, that an
# absent block still loads - the append-only contract, and what it now
# installs.
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
makeSF <- function(x, y) {
  set.seed(5L)
  dbarts::dbarts(
    x,
    y,
    control = control.ls,
    resid.prior = dbarts::dbartsPriors$fixed(0.3)
  )
}
grabState <- function(s) {
  s$storeState()
  s$state
}

donor.sf <- makeSF(testData$x, testData$y)
donor.sf$model@node.scale <- 1.5
donor.sf$setModel(donor.sf$model)
invisible(donor.sf$run(25L, 5L))
state.sf <- grabState(donor.sf)
# nodeScale / sqrt(numTrees)
expect_equal(state.sf[[1L]][["forests"]][[1L]][["leaf.scale"]], 1.5 / 5)

dest.sf <- makeSF(testData$x, testData$y)
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
old.sf <- makeSF(testData$x, testData$y)
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
  hostile.sf <- makeSF(testData$x, testData$y)
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

# --- a refused install is transactional ------------------------------------
# getPointer re-creates the engine behind a dead pointer and installs the
# stored state into it. When that install is REFUSED - the version floor
# above, a malformed block, anything - the object must be left exactly as it
# was: no engine bound, the stored state untouched. Binding the freshly
# created engine before the install succeeds leaves a live but UNFITTED
# sampler, which the next run silently samples from (stumps, not the saved
# forests) and then stores over the fitted state with.

control.tx <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 3L,
  updateState = TRUE,
  verbose = FALSE
)
set.seed(7L)
sampler.tx <- dbarts::dbarts(testData$x, testData$y, control = control.tx)
invisible(sampler.tx$run(10L, 3L))
sampler.tx$storeState()

tempFile <- tempfile()
saveRDS(sampler.tx, file = tempFile)
rm(sampler.tx)

# a version stamp below the floor stands in for any refusal the install can
# raise; it is the one refusal reachable through pure attribute surgery
revived <- readRDS(tempFile)
stale <- revived$state
attr(stale, "formatVersion") <- 2L
revived$state <- stale

expect_error(revived$run(0L, 1L), pattern = "encoding version 2")
# the second run refuses IDENTICALLY rather than sampling a stump
expect_error(revived$run(0L, 1L), pattern = "encoding version 2")
# nothing was bound and nothing was overwritten: the stored state still holds
# the fitted forests, so a storeState()-less save still carries them
expect_false(.Call(dbarts:::C_dbarts_bartcore_isValidPointer, revived$pointer))
expect_identical(revived$state, stale)

# $setState's own dead-pointer branch is transactional the same way: after a
# refused install the next revival resumes from the state that is still
# stored, so the run after a refusal is bitwise the run without one
refused.tx <- readRDS(tempFile)
expect_error(refused.tx$setState(stale), pattern = "encoding version 2")
expect_identical(refused.tx$state, readRDS(tempFile)$state)
set.seed(11L)
draws.refused <- refused.tx$run(0L, 2L)

clean.tx <- readRDS(tempFile)
set.seed(11L)
draws.clean <- clean.tx$run(0L, 2L)
expect_identical(draws.refused$train, draws.clean$train)

unlink(tempFile)
rm(
  control.tx,
  tempFile,
  revived,
  stale,
  refused.tx,
  draws.refused,
  clean.tx,
  draws.clean
)
