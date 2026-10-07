# A Student-t state and the rows at weight zero. A scale is drawn every sweep
# at every row, and at a row of weight zero it is drawn without that row's
# residual. So a Student-t state names those rows, in "weights.zero", and an
# install under other weights redraws the scale of exactly the rows that
# enter the likelihood - at zero then, positive and active now - which is
# what the same setWeights call does on a sampler holding the state under the
# weights it was stored under. Every other row keeps its stored scale.

set.seed(9021L, sample.kind = "Rejection")
n <- 81L
x <- matrix(runif(n * 2L), n)
f <- 2.5 * x[, 1L] - 1.2
yBinary <- as.double(rbinom(n, 1L, plogis(f)))
yContinuous <- as.double(f + rnorm(n))
wA <- rep(c(1, 4, 4), length.out = n)
wB <- rep(c(2, 2, 5), length.out = n)
nu <- 5

controlFor <- function(chains = 1L, seed = 902L) {
  dbarts::dbartsControl(
    n.chains = chains,
    n.threads = 1L,
    n.trees = 20L,
    updateState = FALSE,
    seed = seed
  )
}
# zeros go in through $setWeights, creation warning about them; a NULL df is
# estimated
studentAt <- function(weights, control = controlFor(), df = nu) {
  sampler <- dbarts::dbarts(
    x,
    yContinuous,
    weights = wA,
    control = control,
    family = if (is.null(df)) dbarts:::student() else dbarts:::student(df)
  )
  sampler$setWeights(weights)
  sampler
}
storedFrom <- function(sampler) {
  invisible(sampler$run(10L, 5L))
  sampler$storeState()
  sampler$state
}
scalesOf <- function(sampler) as.matrix(sampler$getLatents())
storedScales <- function(state) {
  vapply(state, function(chain) chain[["latents"]], numeric(n))
}
without <- function(state, what) {
  for (name in what) {
    attr(state, name) <- NULL
  }
  state
}

# stored under: rows 1:20 at weight zero
wStored <- replace(wA, 1:20, 0)
destinations <- list(
  own = wStored,
  sameZeroRows = replace(wB, 1:20, 0),
  zeroRowsMoved = replace(wB, 62:81, 0), # 1:20 enter, 62:81 leave
  allPositive = wB, # 1:20 enter
  moreZeroRows = replace(wStored, 21:30, 0) # 21:30 leave
)
wMoved <- destinations$zeroRowsMoved
# rows 5:10 enter under wMoved, rows 30:35 are positive under both
mask <- replace(rep(1, n), c(5:10, 30:35), 0)

# --- what a state carries ---------------------------------------------------
# a byte per row beside the digest, 1 at a zero weight, on a Student-t state
# only; the mask is not in it and the encoding version does not move for it
state <- storedFrom(studentAt(wStored))
stored <- storedScales(state)
record <- attr(state, "weights.zero")
expect_true(is.raw(record))
expect_identical(record, as.raw(wStored == 0))
expect_identical(attr(state, "formatVersion"), 1L)
expect_identical(
  attr(storedFrom(studentAt(wB)), "weights.zero"),
  raw(n)
)
unweightedState <- storedFrom(dbarts::dbarts(
  x,
  yContinuous,
  control = controlFor(),
  family = dbarts:::student(nu)
))
expect_identical(attr(unweightedState, "weights.zero"), raw(n))
masked <- studentAt(wStored)
masked$setActiveRows(mask)
maskedState <- storedFrom(masked)
expect_identical(attr(maskedState, "weights.zero"), as.raw(wStored == 0))
gaussianAt <- function(weights) {
  dbarts::dbarts(x, yContinuous, weights = weights, control = controlFor())
}
gaussianState <- storedFrom(gaussianAt(wA))
for (family in c("logistic", "probit")) {
  other <- dbarts::dbarts(x, yBinary, family = family, control = controlFor())
  expect_null(attr(storedFrom(other), "weights.zero"))
}
expect_null(attr(gaussianState, "weights.zero"))

# --- the identity with setWeights -------------------------------------------
# The restored sampler takes the state under the destination's weights; the
# twin takes it under the stored ones, which draws nothing, and then the same
# setWeights call. The two are one sampler: the rows that enter are redrawn,
# the rest hold the stored scales, and the next draws agree.
for (chains in 1:2) {
  for (df in list(nu, NULL)) {
    donor <- storedFrom(studentAt(wStored, controlFor(chains), df))
    kept <- storedScales(donor)
    for (destination in destinations) {
      entering <- wStored == 0 & destination > 0
      restored <- studentAt(destination, controlFor(chains), df)
      expect_true(restored$setState(donor))
      twin <- studentAt(wStored, controlFor(chains), df)
      twin$setState(donor)
      expect_identical(scalesOf(twin), kept)
      twin$setWeights(destination)
      redrawn <- scalesOf(restored) != kept
      expect_identical(scalesOf(restored), scalesOf(twin))
      expect_identical(which(rowSums(redrawn) > 0), which(entering))
      expect_true(all(redrawn[entering, ]))
      expect_identical(restored$run(0L, 3L), twin$run(0L, 3L))
    }
  }
}

# an unweighted sampler's state is one stored under no zero rows, so under
# weights with zeros no row enters and nothing is redrawn
underZeros <- studentAt(wStored)
expect_true(underZeros$setState(unweightedState))
expect_identical(scalesOf(underZeros), storedScales(unweightedState))

# two restores of one state under the same other weights are one sampler
again <- studentAt(wMoved)
again$setState(state)
restored <- studentAt(wMoved)
restored$setState(state)
expect_identical(scalesOf(again), scalesOf(restored))
expect_identical(again$run(0L, 3L), restored$run(0L, 3L))

# --- under a mask -----------------------------------------------------------
# On a live sampler the mask is in force at the install: an entering row that
# is masked keeps its stored scale and is redrawn when the mask lifts, with
# every other masked row at positive weight, as on the twin.
restored <- studentAt(wMoved)
restored$setActiveRows(mask)
invisible(restored$run(3L, 1L))
expect_true(restored$setState(state))
twin <- studentAt(wStored)
twin$setActiveRows(mask)
twin$setState(state)
twin$setWeights(wMoved)
expect_identical(scalesOf(restored), scalesOf(twin))
expect_identical(
  which(scalesOf(restored)[, 1L] != stored[, 1L]),
  which(wStored == 0 & wMoved > 0 & mask == 1)
)
held <- scalesOf(restored)
restored$setActiveRows(NULL)
twin$setActiveRows(NULL)
expect_identical(scalesOf(restored), scalesOf(twin))
expect_identical(
  which(scalesOf(restored)[, 1L] != held[, 1L]),
  which(wMoved > 0 & mask == 0)
)

# A copy is made again from the state, the mask going back after it. Stored
# under the mask and copied after a swap between weights with the same zero
# rows, no row enters: the masked rows are not rows at weight zero.
masked$setWeights(destinations$sameZeroRows, updateState = FALSE)
expect_identical(scalesOf(masked$copy()), storedScales(maskedState))

# --- a copy and a reload after the weights changed without a store ----------
# Both are made from the state stored under the earlier weights and moved to
# the current ones as setWeights moved the live sampler.
for (chains in 1:2) {
  live <- studentAt(wStored, controlFor(chains))
  kept <- storedScales(storedFrom(live))
  live$setWeights(wMoved, updateState = FALSE)
  expect_identical(
    which(rowSums(scalesOf(live) != kept) > 0),
    which(wStored == 0)
  )
  copied <- live$copy()
  file <- tempfile()
  saveRDS(live, file = file)
  reloaded <- readRDS(file)
  unlink(file)
  expect_true(max(abs(scalesOf(copied) - scalesOf(live))) <= 1e-13)
  expect_true(max(abs(scalesOf(reloaded) - scalesOf(live))) <= 1e-13)
  ahead <- live$run(0L, 3L)$train
  expect_true(max(abs(copied$run(0L, 3L)$train - ahead)) <= 1e-13)
  expect_true(max(abs(reloaded$run(0L, 3L)$train - ahead)) <= 1e-13)
}

# --- undoing a weight change ------------------------------------------------
# Store, propose other weights, sweep, go back: by the old weights and the
# stored state in either order. A sampler that never proposed is the
# reference.
proposing <- function(proposed) {
  sampler <- studentAt(wStored)
  kept <- storedFrom(sampler)
  sampler$setWeights(proposed)
  invisible(sampler$run(0L, 3L))
  list(sampler = sampler, state = kept)
}
reference <- studentAt(wStored)
reference$setState(storedFrom(reference))
reference <- reference$run(0L, 5L)
for (proposed in destinations[c("sameZeroRows", "zeroRowsMoved")]) {
  # the old weights first: the stored chain
  undone <- proposing(proposed)
  expect_identical(storedScales(undone$state), stored)
  undone$sampler$setWeights(wStored)
  expect_true(undone$sampler$setState(undone$state))
  expect_identical(scalesOf(undone$sampler), stored)
  expect_identical(undone$sampler$run(0L, 5L), reference)
}
# the state first, the zero rows the same under the proposal: the stored chain
undone <- proposing(destinations$sameZeroRows)
expect_true(undone$sampler$setState(undone$state))
expect_identical(scalesOf(undone$sampler), stored)
undone$sampler$setWeights(wStored)
expect_identical(scalesOf(undone$sampler), stored)
expect_identical(undone$sampler$run(0L, 5L), reference)
# the state first, the zero rows moved: the rows positive under both hold the
# stored scales, and the rows the proposal took to zero are redrawn as the old
# weights bring them back
undone <- proposing(wMoved)
expect_true(undone$sampler$setState(undone$state))
undone$sampler$setWeights(wStored)
both <- wStored > 0 & wMoved > 0
back <- wStored > 0 & wMoved == 0
expect_identical(scalesOf(undone$sampler)[both, ], stored[both, ])
expect_true(all(scalesOf(undone$sampler)[back, ] != stored[back, ]))

# --- a state from before the record -----------------------------------------
# With no record nothing says which rows were out, so every row at positive
# weight and active is redrawn, each chain off its own generator; with no
# digest either, nothing is.
for (chains in 1:2) {
  donor <- storedFrom(studentAt(wStored, controlFor(chains)))
  kept <- storedScales(donor)
  restored <- studentAt(wMoved, controlFor(chains))
  restored$setActiveRows(mask)
  expect_true(restored$setState(without(donor, "weights.zero")))
  inLikelihood <- wMoved > 0 & mask == 1
  redrawn <- scalesOf(restored) != kept
  expect_true(all(redrawn[inLikelihood, ]))
  expect_false(any(redrawn[!inLikelihood, ]))
  if (chains == 2L) {
    redrawn <- scalesOf(restored)[inLikelihood, ]
    expect_true(all(redrawn[, 1L] != redrawn[, 2L]))
  }
  expect_true(all(is.finite(restored$run(0L, 3L)$train)))
  restored <- studentAt(wMoved, controlFor(chains))
  restored$setState(without(donor, c("weights.zero", "weights.digest")))
  expect_identical(scalesOf(restored), kept)
}

# --- refusals ---------------------------------------------------------------
# a record is a raw 0 or 1 for each row of this sampler, on any family, and a
# refused state touches nothing
restored <- studentAt(wMoved)
invisible(restored$run(3L, 1L))
held <- scalesOf(restored)
malformed <- state
attr(malformed, "weights.zero") <- wStored == 0
expect_error(
  restored$setState(malformed),
  pattern = "malformed zero-weight rows in bartcore state"
)
attr(malformed, "weights.zero") <- replace(record, 3L, as.raw(2L))
expect_error(
  restored$setState(malformed),
  pattern = "malformed zero-weight rows in bartcore state"
)
for (wrongLength in list(record[-1L], c(record, as.raw(0L)))) {
  attr(malformed, "weights.zero") <- wrongLength
  expect_error(
    restored$setState(malformed),
    pattern = "state is not consistent with this sampler"
  )
}
expect_identical(scalesOf(restored), held)
# the record is checked whether or not it will be used: under the state's own
# weights, where nothing is redrawn, and on a state with no digest
matched <- studentAt(wStored)
invisible(matched$run(3L, 1L))
twin <- studentAt(wStored)
invisible(twin$run(3L, 1L))
attr(malformed, "weights.zero") <- replace(record, 3L, as.raw(2L))
for (offered in list(malformed, without(malformed, "weights.digest"))) {
  expect_error(
    matched$setState(offered),
    pattern = "malformed zero-weight rows in bartcore state"
  )
}
expect_identical(scalesOf(matched), scalesOf(twin))
expect_identical(matched$run(0L, 3L), twin$run(0L, 3L))
# on a gaussian state a well-formed one is read and not used
withRecord <- gaussianState
attr(withRecord, "weights.zero") <- record
plain <- gaussianAt(wB)
expect_true(plain$setState(gaussianState))
marked <- gaussianAt(wB)
expect_true(marked$setState(withRecord))
expect_identical(marked$run(0L, 3L), plain$run(0L, 3L))
attr(withRecord, "weights.zero") <- replace(record, 3L, as.raw(2L))
expect_error(
  gaussianAt(wB)$setState(withRecord),
  pattern = "malformed zero-weight rows in bartcore state"
)

# --- an entering row does not keep its stored scale -------------------------
# A scale stored for a row at weight zero was drawn without the row's
# residual. Transformed by the conditional at the stored fit and sigma - gamma
# with shape (nu + 1) / 2 and rate (nu + w_i (y_i - f_i)^2 / sigma^2) / 2 - a
# redrawn scale is uniform, and a kept one has a mean near 0.62. Forty rounds
# of twenty entering rows. The bound on the mean, 0.05, is 4.9 standard errors
# at 800 draws, and the p-value is one a sample from the conditional falls
# under once in a million; both are wide because the draws are not the same on
# every platform, so each is a fresh sample there. This guards a kept scale
# and nothing finer: a wrong weight, fit or sigma in the conditional is caught
# by the bitwise identity with the setWeights twin above, not here.
source <- studentAt(wStored)
restored <- studentAt(wB)
entering <- which(wStored == 0)
transformed <- numeric(0)
for (round in seq_len(40L)) {
  draw <- source$run(3L, 1L)
  source$storeState()
  restored$setState(source$state)
  rate <- 0.5 *
    (nu +
      wB[entering] *
        (yContinuous[entering] - draw$train[entering, 1L])^2 /
        draw$sigma[1L]^2)
  transformed <- c(
    transformed,
    pgamma(restored$getLatents()[entering], shape = 0.5 * (nu + 1), rate = rate)
  )
}
expect_true(abs(mean(transformed) - 0.5) < 0.05)
expect_true(ks.test(transformed, "punif")$p.value > 1e-6)

rm(
  n,
  x,
  f,
  yBinary,
  yContinuous,
  wA,
  wB,
  nu,
  controlFor,
  studentAt,
  storedFrom,
  scalesOf,
  storedScales,
  without,
  wStored,
  destinations,
  wMoved,
  mask,
  state,
  stored,
  record,
  unweightedState,
  masked,
  maskedState,
  gaussianAt,
  gaussianState,
  family,
  other,
  chains,
  df,
  donor,
  kept,
  destination,
  entering,
  restored,
  twin,
  redrawn,
  underZeros,
  again,
  held,
  live,
  copied,
  file,
  reloaded,
  ahead,
  proposing,
  reference,
  proposed,
  undone,
  both,
  back,
  inLikelihood,
  malformed,
  wrongLength,
  matched,
  offered,
  withRecord,
  plain,
  marked,
  source,
  transformed,
  round,
  draw,
  rate
)
