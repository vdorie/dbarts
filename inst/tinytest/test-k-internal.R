# k is each chain's parameter on the internal scale, measured against the
# data's scale: k.scale is fixed by the family and the response mapping and
# is never redefined by how a prior is written. A prior written with an sd is
# sugar - an sd of s is k = k.scale / s, and invchi(df, s) is
# k ~ chi(df, k.scale / s) - restated whenever the mapping is re-derived, so
# the named sd keeps its meaning in response units while a drawn k, and so
# its current sd, lags with the leaves until the next draw. Every install -
# a restored state, a copy, a reload, a warm start - puts the chain in as
# stored on the internal scale, read against the receiving sampler's mapping,
# and converts nothing. One block per operation; the reads are the sampler's
# own readers and, beneath them, the engine's calibration.

for (name in c("normal", "chi", "invchi", "gp")) {
  assign(name, dbartsPriors[[name]])
}

set.seed(11)
n <- 200L
x <- matrix(runif(n * 3L), n, 3L)
colnames(x) <- paste0("x", 1:3)
y <- 2 * sin(pi * x[, 1L]) + x[, 2L] + rnorm(n, sd = 0.4)
z <- as.double(y > median(y))
wide <- 3 * y
# the data's scale: half the response range
kScale <- (max(y) - min(y)) / 2

kControl <- function(n.chains = 1L, keepTrees = FALSE, n.samples = 4L) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = 1L,
    n.trees = 20L,
    n.samples = n.samples,
    n.burn = 0L,
    keepTrees = keepTrees,
    updateState = FALSE,
    seed = 5L
  )
}
make <- function(leaf.prior, response = y, control = kControl(), ...) {
  dbarts(x, response, control = control, leaf.prior = leaf.prior, ...)
}
burned <- function(...) {
  sampler <- make(...)
  invisible(sampler$run(40L, 1L))
  sampler
}
stored <- function(sampler) {
  sampler$storeState()
  sampler$state
}
# a single chain's k, which the reader names
kOf <- function(sampler) unname(sampler$getK())
# the engine's own reading, one row per chain, and one column of it
engine <- function(sampler) {
  .Call(dbarts:::C_dbarts_bartcore_getLeafPrior, sampler$getPointer(), 0L)
}
engineOf <- function(sampler, column) unname(engine(sampler)[, column])
# k, k.scale and the spread in force
reading <- function(sampler) {
  kScale <- sampler$getLeafPrior()$k.scale
  k <- kOf(sampler)
  list(k = k, k.scale = kScale, sd = kScale / k)
}
expectReading <- function(sampler, k, kScale, info) {
  read <- reading(sampler)
  expect_identical(read$k, k, info = info)
  expect_equal(read$k.scale, kScale, tolerance = 1e-12, info = info)
  expect_equal(read$sd, kScale / k, tolerance = 1e-12, info = info)
}
# a state's numbers, less its record of the mapping it was read under
unlabelled <- function(state) {
  for (chain in seq_along(state)) {
    state[[chain]]$fit.scale <- NULL
  }
  state
}
internal <- function(sampler, values) {
  prior <- sampler$getLeafPrior()
  (values - prior$response.shift) / prior$response.scale
}

# --- creation: the k spellings are untouched, and an sd is k = k.scale / s
# against the data's scale, which every spelling reports ---
expectReading(make(normal(k = chi(1.5, 2))), 2, kScale, "k = chi(1.5, 2)")
expectReading(make(normal(k = chi(1.5, 4))), 2, kScale, "k = chi(1.5, 4)")
expectReading(make(normal(k = 3)), 3, kScale, "k = 3")
fixedSd <- make(normal(sd = 0.5))
expect_identical(kOf(fixedSd), fixedSd$getLeafPrior()$k.scale / 0.5)
expectReading(fixedSd, kOf(fixedSd), kScale, "sd = 0.5")
expect_equal(reading(fixedSd)$sd, 0.5, tolerance = 1e-12)
expect_identical(engineOf(fixedSd, "named.sd"), 0.5)
expect_true(is.na(engineOf(make(normal(k = 3)), "named.sd")))
# a prior on the sd: its chi scale is the translation, and the chain starts
# there, so at the named sd
drawnSd <- make(normal(sd = invchi(3, 0.5)))
expect_identical(kOf(drawnSd), drawnSd$getLeafPrior()$k.scale / 0.5)
expect_identical(
  engineOf(drawnSd, "k.prior.scale"),
  engineOf(drawnSd, "prior.scale") / 0.5
)
expectReading(drawnSd, kOf(drawnSd), kScale, "sd = invchi(3, 0.5)")
# the improper limit names no scale: chi(df, Inf), started at 2
improper <- make(normal(sd = invchi(3, 0)))
expectReading(improper, 2, kScale, "sd = invchi(3, 0)")
expect_identical(engineOf(improper, "k.prior.scale"), Inf)
expect_true(is.na(engineOf(improper, "named.sd")))
# the binary families' scales are fixed: 3 and pi sqrt(3)
expectReading(
  make(normal(k = chi(1.5, 2)), z, family = "probit"),
  2,
  3,
  "probit"
)
expectReading(make(normal(sd = 0.5), z, family = "probit"), 6, 3, "probit sd")
expectReading(
  make(normal(k = chi(1.5, 2)), z, family = "logistic"),
  2,
  pi * sqrt(3),
  "logistic"
)
# the model's encoding of the prior does not move
expect_identical(fixedSd$model@prior.scale, 1)
expect_identical(drawnSd$model@leaf.hyperprior, chi(3, 2))

# --- the identities: an sd draws what its k draws ---
draws <- function(sampler, n.samples = 60L) sampler$run(0L, n.samples)
# probit: sd = invchi(1.5, 1.5) is the k default's chain, bit for bit
expect_identical(
  draws(make(normal(sd = invchi(1.5, 1.5)), z, family = "probit")),
  draws(make(normal(k = chi(1.5, 2)), z, family = "probit"))
)
for (s in c(0.5, 0.7, 1.3)) {
  expect_identical(
    draws(make(normal(sd = s))),
    draws(make(normal(k = kScale / s))),
    info = s
  )
  # a fixed k turned into chi() keeps k, so it is the drawn sd's start; a
  # sampler created at a fixed k records no k draws, so k is read at the end
  started <- make(normal(k = kScale / s))
  started$setLeafPrior(normal(k = chi(3, kScale / s)))
  sdDrawn <- make(normal(sd = invchi(3, s)))
  expect_identical(
    draws(sdDrawn)[c("train", "sigma")],
    draws(started)[c("train", "sigma")],
    info = s
  )
  expect_identical(kOf(sdDrawn), kOf(started), info = s)
  expect_false(identical(kOf(sdDrawn), kScale / s), info = s)
}

# --- setLeafPrior: k.scale never moves, so k is untouched into a drawn
# prior, and a fixed value stated is taken as written ---
chain <- burned(normal(k = chi(1.5, 2)))
k0 <- kOf(chain)
expect_false(identical(k0, 2))
chain$setLeafPrior(normal(sd = invchi(3, 1)))
expectReading(chain, k0, kScale, "into sd = invchi(3, 1)")
expect_identical(
  engineOf(chain, "k.prior.scale"),
  engineOf(chain, "prior.scale")
)
chain$setLeafPrior(normal(sd = invchi(3, 2)))
expectReading(chain, k0, kScale, "into sd = invchi(3, 2)")
chain$setLeafPrior(normal(sd = 0.5))
expectReading(
  chain,
  chain$getLeafPrior()$k.scale / 0.5,
  kScale,
  "into sd = 0.5"
)
held <- kOf(chain)
chain$setLeafPrior(normal(k = chi(1.5, 2)))
expectReading(chain, held, kScale, "into k = chi(1.5, 2)")
chain$setLeafPrior(normal(k = 3))
expectReading(chain, 3, kScale, "into k = 3")

# --- setModel with another leaf prior keeps k, and so the spread ---
chain <- burned(normal(k = chi(1.5, 2)))
k0 <- kOf(chain)
chain$setModel(make(normal(sd = invchi(3, 1)))$model)
expectReading(chain, k0, kScale, "setModel into sd = invchi(3, 1)")
expect_identical(engineOf(chain, "named.sd"), 1)
chain$setModel(make(normal(k = chi(2, 3)))$model)
expectReading(chain, k0, kScale, "setModel into k = chi(2, 3)")
expect_identical(engineOf(chain, "k.prior.scale"), 3)

# --- xbart's warm sweep goes through the engine's setModel per cell: the
# previous cell's spread is kept into a drawn cell, whichever spelling ---
remodel <- function(sampler, leaf.prior) {
  dbarts:::bartcoreSetModel(
    list(ptr = sampler$getPointer()),
    make(leaf.prior)$model,
    sampler$data,
    sampler$control
  )
}
cell <- burned(normal(sd = 0.25))
remodel(cell, normal(sd = invchi(3, 1)))
expect_equal(reading(cell)$sd, 0.25, tolerance = 1e-12)
expect_identical(engineOf(cell, "k.prior.scale"), engineOf(cell, "prior.scale"))
cell <- burned(normal(sd = invchi(3, 0.5)))
spread <- reading(cell)$sd
expect_false(isTRUE(all.equal(spread, 0.5)))
remodel(cell, normal(sd = invchi(3, 2)))
expect_identical(reading(cell)$sd, spread)
expect_identical(
  engineOf(cell, "k.prior.scale"),
  engineOf(cell, "prior.scale") / 2
)
# a grid mixing a fixed and a drawn cell of different scales runs, and
# repeats under a seed
gridLoss <- function() {
  xbart(
    x,
    y,
    sd = list(0.25, invchi(3, 1)),
    n.samples = 20L,
    n.reps = 2L,
    n.burn = c(20L, 10L),
    n.trees = 10L,
    n.threads = 1L,
    seed = 5L,
    verbose = FALSE
  )
}
firstGrid <- gridLoss()
expect_true(all(is.finite(firstGrid)))
expect_identical(gridLoss(), firstGrid)

# --- an install across a change of prior, same mapping: k as stored, the
# spread kept; a sampler holding k fixed keeps its own ---
chain <- burned(normal(k = chi(1.5, 2)))
state <- stored(chain)
k0 <- kOf(chain)
chain$setLeafPrior(normal(sd = invchi(3, 1)))
chain$setState(state)
expectReading(chain, k0, kScale, "restore after a change of prior")
recipient <- make(normal(sd = invchi(3, 1)))
recipient$setState(state)
expectReading(recipient, k0, kScale, "state into an sd-named sampler")
recipient <- make(normal(k = 3))
recipient$setState(state)
expectReading(recipient, 3, kScale, "state into a fixed k")

# --- installs within one mapping: a sampler's own restore, a copy and a
# reload are one chain, bit for bit, with a held and a drawn named sd, before
# a re-anchor and after one whose state was stored after it ---
reloadOf <- function(sampler) {
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  saveRDS(sampler, file)
  readRDS(file)
}
for (leafPrior in list(normal(sd = 0.5), normal(sd = invchi(3, 0.5)))) {
  for (reanchor in c(FALSE, TRUE)) {
    info <- paste(format(leafPrior@prior.sd), reanchor)
    source <- burned(leafPrior)
    if (reanchor) {
      source$setResponse(wide + 1, updateScale = TRUE)
      invisible(source$run(5L, 1L))
    }
    state <- stored(source)
    copied <- source$copy()
    reloaded <- reloadOf(source)
    for (other in list(copied, reloaded)) {
      expect_identical(other$getLeafPrior(), source$getLeafPrior(), info = info)
      expect_identical(kOf(other), kOf(source), info = info)
      expect_identical(engine(other), engine(source), info = info)
      expect_identical(stored(other), state, info = info)
    }
    expect_true(source$setState(state), info = info)
    restored <- draws(source, 10L)
    expect_identical(draws(copied, 10L), restored, info = info)
    expect_identical(draws(reloaded, 10L), restored, info = info)
  }
}

# --- sigma is stored and installed on the internal scale, with no pass
# through response units: a value the response multiplier does not carry
# there and back, in either order, goes in and reads back bit for bit, in the
# sampler itself, a copy and a reload ---
sigmaSource <- burned(normal(k = 3))
sigmaState <- stored(sigmaSource)
range <- sigmaSource$getLeafPrior()$response.scale
expect_identical(
  unname(sigmaSource$getSigmas()),
  sigmaState[[1L]]$sigma * range
)
candidates <- sigmaState[[1L]]$sigma * (1 + seq_len(8191L) / 4096)
missed <- candidates[
  (candidates * range) / range != candidates &
    (candidates / range) * range != candidates
]
expect_true(length(missed) > 0L)
sigmaState[[1L]]$sigma <- missed[[1L]]
expect_true(sigmaSource$setState(sigmaState))
expect_identical(stored(sigmaSource)[[1L]]$sigma, missed[[1L]])
expect_identical(unname(sigmaSource$getSigmas()), missed[[1L]] * range)
sigmaSource$setState(sigmaState)
sigmaCopy <- sigmaSource$copy()
sigmaReload <- reloadOf(sigmaSource)
expect_identical(stored(sigmaCopy)[[1L]]$sigma, missed[[1L]])
expect_identical(stored(sigmaReload)[[1L]]$sigma, missed[[1L]])
expect_identical(sigmaCopy$getSigmas(), sigmaSource$getSigmas())
sigmaSource$setState(sigmaState)
sigmaDraws <- draws(sigmaSource, 10L)
expect_identical(draws(sigmaCopy, 10L), sigmaDraws)
expect_identical(draws(sigmaReload, 10L), sigmaDraws)

# --- installs across mappings: setState and a warm start onto a response 3
# times as wide put in the stored internal numbers. Nothing is converted, so
# the fit, the spread and sigma are all 3 times wider, a pure rescale of the
# response being no change on the internal scale ---
for (leafPrior in list(normal(k = chi(1.5, 2)), normal(sd = invchi(3, 1)))) {
  info <- if (is.null(leafPrior@prior.sd)) "k spelling" else "sd spelling"
  donor <- burned(leafPrior)
  state <- stored(donor)
  recipient <- make(leafPrior, wide)
  unitsBefore <- stored(recipient)[[1L]]$fit.scale
  expect_true(recipient$setState(state), info = info)
  installed <- stored(recipient)
  expect_identical(installed[[1L]]$fit.scale, unitsBefore, info = info)
  expect_identical(unlabelled(installed), unlabelled(state), info = info)
  expect_identical(kOf(recipient), kOf(donor), info = info)
  expect_equal(reading(recipient)$k.scale, 3 * kScale, tolerance = 1e-12)
  expect_equal(reading(recipient)$sd, 3 * reading(donor)$sd, tolerance = 1e-12)
  expect_equal(
    recipient$getSigmas(),
    3 * donor$getSigmas(),
    tolerance = 1e-12,
    info = info
  )
  expect_equal(
    recipient$getFitsWithoutOffset(),
    3 * donor$getFitsWithoutOffset(),
    tolerance = 1e-12,
    info = info
  )
  warmed <- make(leafPrior, wide)
  warmed$installTrees(donor)
  expect_identical(stored(warmed)[[1L]]$fit.scale, unitsBefore, info = info)
  expect_identical(kOf(warmed), kOf(donor), info = info)
  expect_identical(
    stored(warmed)[[1L]]$sigma,
    state[[1L]]$sigma,
    info = info
  )
  expect_identical(
    stored(warmed)[[1L]]$forests[[1L]]$tree.values,
    state[[1L]]$forests[[1L]]$tree.values,
    info = info
  )
  expect_equal(
    warmed$getFitsWithoutOffset(),
    3 * donor$getFitsWithoutOffset(),
    tolerance = 1e-12,
    info = info
  )
}
# no leaf model refuses a state from another shift
gpDonor <- burned(gp("x2"))
gpRecipient <- make(gp("x2"), wide + 10)
expect_true(gpRecipient$setState(stored(gpDonor)))
expect_equal(
  internal(gpRecipient, gpRecipient$getFitsWithoutOffset()),
  internal(gpDonor, gpDonor$getFitsWithoutOffset()),
  tolerance = 1e-11
)

# --- a re-anchor and an install differ on the kept draws: a re-anchor
# converts them, so predict on them does not move, and an install brings
# them in as stored, so predict reads them on the new scale ---
keepControl <- kControl(keepTrees = TRUE)
kept <- make(normal(k = 2), control = keepControl)
invisible(kept$run(20L, 4L))
keptState <- stored(kept)
keptPredictions <- kept$predict(x)
keptInternal <- internal(kept, keptPredictions)
kept$setResponse(wide + 1, updateScale = TRUE)
expect_equal(kept$predict(x), keptPredictions, tolerance = 1e-12)
installedKept <- make(normal(k = 2), wide + 1, control = keepControl)
expect_true(installedKept$setState(keptState))
expect_false(isTRUE(all.equal(installedKept$predict(x), keptPredictions)))
expect_equal(
  internal(installedKept, installedKept$predict(x)),
  keptInternal,
  tolerance = 1e-12
)

# --- a re-anchor: a k-named prior and a drawn k stay put, a held sd is
# restated so its spread is the sd, and a prior on the sd is restated while
# the k drawn so far lags with the leaves until the next draw ---
anchored <- burned(normal(k = chi(1.5, 2)))
k0 <- kOf(anchored)
anchored$setResponse(wide, updateScale = TRUE)
expectReading(anchored, k0, 3 * kScale, "re-anchor, k spelling")
anchored <- burned(normal(sd = 0.5))
anchored$setResponse(wide, updateScale = TRUE)
expect_identical(kOf(anchored), anchored$getLeafPrior()$k.scale / 0.5)
expect_equal(reading(anchored)$k.scale, 3 * kScale, tolerance = 1e-12)
expect_equal(reading(anchored)$sd, 0.5, tolerance = 1e-12)
anchored <- burned(normal(sd = invchi(3, 1)))
k0 <- kOf(anchored)
spread <- reading(anchored)$sd
anchored$setResponse(wide, updateScale = TRUE)
expectReading(anchored, k0, 3 * kScale, "re-anchor, drawn sd")
expect_equal(reading(anchored)$sd, 3 * spread, tolerance = 1e-12)
expect_identical(
  engineOf(anchored, "k.prior.scale"),
  engineOf(anchored, "prior.scale") / 1
)
# setOffset and setData re-anchor too
for (mutate in list(
  function(s) s$setOffset(2 * x[, 3L] + 0.5, updateScale = TRUE),
  function(s) s$setData(dbartsData(x, wide + 1))
)) {
  heldAnchor <- burned(normal(sd = 0.5))
  scaleBefore <- heldAnchor$getLeafPrior()$k.scale
  mutate(heldAnchor)
  expect_false(isTRUE(all.equal(
    heldAnchor$getLeafPrior()$k.scale,
    scaleBefore
  )))
  expect_identical(kOf(heldAnchor), heldAnchor$getLeafPrior()$k.scale / 0.5)
  drawnAnchor <- burned(normal(sd = invchi(3, 1)))
  k0 <- kOf(drawnAnchor)
  mutate(drawnAnchor)
  expect_identical(kOf(drawnAnchor), k0)
  expect_identical(
    engineOf(drawnAnchor, "k.prior.scale"),
    engineOf(drawnAnchor, "prior.scale")
  )
}
# without the scale update nothing moves
unmoved <- burned(normal(sd = invchi(3, 1)))
before <- engine(unmoved)
unmoved$setResponse(wide)
unmoved$setOffset(rep_len(2, n))
expect_identical(engine(unmoved), before)

# --- a re-creation after a re-anchor puts the chains at the recorded
# mapping before any state goes in: a copy, a reload and setState on a dead
# pointer, with a held sd, whose spread is still the sd, and a drawn one,
# whose chi scale is still k.scale / sd. With the response then swapped back
# under the pinned mapping the re-created sampler's own data name another
# mapping than the record, so the chains are there only because the
# re-creation moved them, and the sd was restated with the move ---
for (leafPrior in list(normal(sd = 0.5), normal(sd = invchi(3, 0.5)))) {
  for (pinned in c(FALSE, TRUE)) {
    drawn <- !is.numeric(leafPrior@prior.sd)
    info <- paste(if (drawn) "drawn sd" else "held sd", pinned)
    source <- burned(leafPrior)
    source$setResponse(wide + 1, updateScale = TRUE)
    invisible(source$run(5L, 1L))
    if (pinned) {
      source$setResponse(y, updateScale = FALSE)
      invisible(source$run(2L, 1L))
    }
    state <- stored(source)
    dead <- unserialize(serialize(source, NULL))
    expect_true(dead$setState(state), info = info)
    for (other in list(source$copy(), reloadOf(source), dead)) {
      expect_identical(engine(other), engine(source), info = info)
      read <- function(column) engineOf(other, column)
      expect_identical(
        stored(other)[[1L]]$fit.scale,
        state[[1L]]$fit.scale,
        info = info
      )
      expect_equal(read("prior.scale"), 3 * kScale, tolerance = 1e-12)
      if (drawn) {
        expect_identical(
          read("k.prior.scale"),
          read("prior.scale") / 0.5,
          info = info
        )
      } else {
        expect_equal(read("prior.sd"), 0.5, tolerance = 1e-12, info = info)
      }
    }
  }
}

# on the pinned fixture the move carries the rest of the model too: a held
# sigma keeps its response-units value, and a variance forest's calibration
# and surface are the recorded mapping's, so a copy, a reload and a setState
# on a dead pointer agree with the source
reanchoredAndPinned <- function(...) {
  source <- burned(normal(sd = 0.5), ...)
  source$setResponse(3 * y + 1, updateScale = TRUE)
  invisible(source$run(5L, 1L))
  source$setResponse(y, updateScale = FALSE)
  invisible(source$run(2L, 1L))
  source
}
heldSource <- reanchoredAndPinned(family = gaussian(sigma = fixed(4)))
heldSigma <- unname(heldSource$getSigmas())
expect_equal(heldSigma, 2, tolerance = 1e-12)
heldState <- stored(heldSource)
heldDead <- unserialize(serialize(heldSource, NULL))
expect_true(heldDead$setState(heldState))
for (other in list(heldSource$copy(), reloadOf(heldSource), heldDead)) {
  expect_identical(unname(other$getSigmas()), heldSigma)
}
varSource <- reanchoredAndPinned(variance = TRUE)
varState <- stored(varSource)
varCopy <- varSource$copy()
varReload <- reloadOf(varSource)
varDead <- unserialize(serialize(varSource, NULL))
expect_true(varDead$setState(varState))
expect_true(varSource$setState(varState))
varRestored <- draws(varSource, 10L)
expect_identical(draws(varCopy, 10L), varRestored)
expect_identical(draws(varReload, 10L), varRestored)
expect_identical(draws(varDead, 10L), varRestored)

# --- what a fit and the readers report under an sd spelling: the real k
# against the data's scale, the data's scale, the prior as named, and the sd
# itself ---
fitOf <- function(leaf.prior, response = y) {
  bart(
    x,
    response,
    leaf.prior = leaf.prior,
    n.trees = 20L,
    n.burn = 20L,
    n.samples = 20L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 3L,
    verbose = FALSE
  )
}
drawnFit <- suppressMessages(fitOf(normal(sd = invchi(3, 0.5))))
expect_equal(drawnFit$leaf.prior$k.scale, kScale, tolerance = 1e-12)
expect_identical(drawnFit$leaf.prior$leaf.prior, normal(sd = invchi(3, 0.5)))
expect_identical(
  extract(drawnFit, "leaf.prior.sd"),
  drawnFit$leaf.prior$k.scale / extract(drawnFit, "k")
)
expect_null(drawnFit$fixed$k)
# a held sd is read exactly. The held k is k.scale / s, and k.scale over it
# misses s by an ulp on some responses, so each sd is read on a response
# picked for that miss, and on the plain one
for (s in c(0.7, 1.3, 0.1)) {
  missing <- 0
  for (scale in 1 + seq_len(60L) / 7) {
    probed <- make(normal(sd = s), scale * y)$getLeafPrior()$k.scale
    if (probed / (probed / s) != s) {
      missing <- scale
      break
    }
  }
  expect_true(missing > 0, info = s)
  for (scale in c(1, missing)) {
    heldFit <- suppressMessages(fitOf(normal(sd = s), scale * y))
    info <- paste(s, scale)
    heldScale <- heldFit$leaf.prior$k.scale
    expect_identical(heldFit$fixed$k, heldScale / s, info = info)
    expect_identical(extract(heldFit, "leaf.prior.sd"), s, info = info)
    expect_identical(summary(heldFit)$fixed$leaf.prior.sd, s, info = info)
    expect_identical(heldFit$leaf.prior$leaf.prior, normal(sd = s), info = info)
    expect_null(heldFit$k, info = info)
  }
  expect_false(identical(heldScale / heldFit$fixed$k, s), info = s)
}
# a write-back of the prior a sampler reports moves no bit, under both forms
for (leafPrior in list(normal(sd = 0.7), normal(sd = invchi(3, 0.7)))) {
  written <- burned(leafPrior)
  twin <- burned(leafPrior)
  expect_identical(written$getLeafPrior()$leaf.prior, leafPrior)
  before <- engine(written)
  written$setLeafPrior(written$getLeafPrior()$leaf.prior)
  expect_identical(engine(written), before)
  expect_identical(draws(written, 10L), draws(twin, 10L))
}
