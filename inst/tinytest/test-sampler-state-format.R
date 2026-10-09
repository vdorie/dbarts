# The additive-by-name state format: setState reads per-chain blocks by
# name behind an encoding
# floor, so future additive versions still load, a genuinely older encoding is
# refused, and a missing REQUIRED block is named. States are opaque R lists
# with attributes, so these are pure attribute surgery (no C needed).

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

set.seed(99L)
bartFit <- dbarts::bartBT(
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

# the floor is what makes a BLOCK RENAME safe. An older encoding spelled the
# amplitude glue "bcf"; were such a state read here it would find "glue"
# absent, default it as an optional block, and leave the amplitudes at their
# construction values - a wrong answer, not an error. The version check runs
# BEFORE any block is read, so a state still carrying the old name is refused
# by version.
priorEncoding <- state
priorEncoding[[1L]][["bcf"]] <- c(1, 1, 1, 1)
attr(priorEncoding, "formatVersion") <- 0L
expect_error(
  bartFit$fit$setState(priorEncoding),
  pattern = "encoding version 0.*oldest this dbarts \\(1\\)"
)

# sigma is present only where the saver drew it, so a state without it loads
# and the sampler keeps its own
expect_true(is.numeric(state[[1L]][["sigma"]]))
missingSigma <- state
missingSigma[[1L]][["sigma"]] <- NULL
bartFit$fit$setState(missingSigma)
expect_equal(predict(bartFit, testData$x), preds)

# a block present but of the wrong type is named as malformed, not missing -
# the two-message convention.
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

# leaf.scale: the leaf prior's scale is model, so a state holds no block for
# it, and k is written only where it is drawn. A state written before, which
# carries the block, still installs: the block is not read, whatever it holds.
forest1 <- state[[1L]][["forests"]][[1L]]
expect_false("leaf.scale" %in% names(forest1))
expect_null(forest1[["k"]])
for (oldValue in list(0.2, c(1, 2), "not a number")) {
  withLeafScale <- state
  withLeafScale[[1L]][["forests"]][[1L]][["leaf.scale"]] <- oldValue
  bartFit$fit$setState(withLeafScale)
  expect_equal(predict(bartFit, testData$x), preds)
}

# the donor's model does not ride its state: a destination keeps its own leaf
# prior and its own fixed sigma, and a state that carries the donor's leaf
# scale, k and sigma, as one written before did, installs the same way
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
    family = gaussian(sigma = dbarts::dbartsPriors$fixed(0.3))
  )
}
grabState <- function(s) {
  s$storeState()
  s$state
}

donor.sf <- makeSF(testData$x, testData$y)
donor.sf$model@leaf.scale <- 1.5
donor.sf$setModel(donor.sf$model)
invisible(donor.sf$run(25L, 5L))
state.sf <- grabState(donor.sf)
expect_null(state.sf[[1L]][["sigma"]])

dest.sf <- makeSF(testData$x, testData$y)
priorBefore <- dest.sf$getLeafPrior()
sigmaBefore <- dest.sf$getSigmas()
expect_true(priorBefore$k.scale != donor.sf$getLeafPrior()$k.scale)
dest.sf$setState(state.sf)
expect_identical(dest.sf$getLeafPrior(), priorBefore)
expect_equal(dest.sf$getSigmas(), sigmaBefore)

state.old <- state.sf
state.old[[1L]][["forests"]][[1L]][["leaf.scale"]] <- 1.5 / 5
state.old[[1L]][["forests"]][[1L]][["k"]] <- 4
state.old[[1L]][["sigma"]] <- 2
old.sf <- makeSF(testData$x, testData$y)
old.sf$setState(state.old)
expect_identical(old.sf$getLeafPrior(), priorBefore)
expect_equal(old.sf$getSigmas(), sigmaBefore)
expect_identical(old.sf$run(0L, 3L), dest.sf$run(0L, 3L))

rm(
  forest1,
  oldValue,
  withLeafScale,
  control.ls,
  makeSF,
  grabState,
  donor.sf,
  state.sf,
  dest.sf,
  priorBefore,
  sigmaBefore,
  state.old,
  old.sf
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
attr(stale, "formatVersion") <- 0L
revived$state <- stale

expect_error(revived$run(0L, 1L), pattern = "encoding version 0")
# the second run refuses IDENTICALLY rather than sampling a stump
expect_error(revived$run(0L, 1L), pattern = "encoding version 0")
# nothing was bound and nothing was overwritten: the stored state still holds
# the fitted forests, so a storeState()-less save still carries them
expect_false(.Call(dbarts:::C_dbarts_bartcore_isValidPointer, revived$pointer))
expect_identical(revived$state, stale)

# $setState's own dead-pointer branch is transactional the same way: after a
# refused install the next revival resumes from the state that is still
# stored, so the run after a refusal is bitwise the run without one
refused.tx <- readRDS(tempFile)
expect_error(refused.tx$setState(stale), pattern = "encoding version 0")
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
