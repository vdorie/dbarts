# The $setForestWeights(forest, weights) surface over the engine channel
# test-forest-weights.R already pins end to end. Not re-pinned here: the
# capability probe, the sigma degrees of freedom, and the substituted-response
# invariance of excluded rows. Covered here instead: the 1-based boundary
# conversion at resolveForestIndex, the mirror-and-reapply field at all THREE
# re-creation sites (getPointer, setState, and copy), that the channel never
# moves the reported mean, and the two refusals reachable through the public
# route (a multinomial engine is never reachable this way - see
# test-active-rows-pins.R for its internal-route pin, which this does not
# duplicate).

set.seed(97)
n <- 220L
p <- 4L
x <- matrix(runif(n * p), n, p)
z <- rbinom(n, 1L, 0.5)
mu <- 2 * sin(pi * x[, 1L]) + x[, 2L]
tau <- 1 + 2 * x[, 3L]
y <- mu + z * tau + rnorm(n, sd = 0.2)
wt <- runif(n, 0.2, 2)

seededControl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 30L,
    n.samples = 6L,
    updateState = FALSE,
    seed = 71L,
    ...
  )
}
build <- function() {
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z))),
    control = seededControl()
  )
}

# --- 1-based boundary conversion: setForestWeights(1L, ) reaches the
# engine's forest 0 (prognostic) and setForestWeights(2L, ) reaches forest 1
# (treatment) - the same calls a raw bridge install at those 0-based indices
# makes on an identically seeded sampler - and 0L is refused by
# resolveForestIndex before any bridge call ---
viaR5.0 <- build()
viaR5.0$setForestWeights(1L, wt)
result.r5.0 <- viaR5.0$run(0L, 3L)

viaRaw.0 <- build()
.Call(
  dbarts:::C_dbarts_bartcore_setForestWeights,
  viaRaw.0$getPointer(),
  0L,
  wt
)
result.raw.0 <- viaRaw.0$run(0L, 3L)
expect_identical(result.r5.0$train, result.raw.0$train)

viaR5.1 <- build()
viaR5.1$setForestWeights(2L, wt)
result.r5.1 <- viaR5.1$run(0L, 3L)

viaRaw.1 <- build()
.Call(
  dbarts:::C_dbarts_bartcore_setForestWeights,
  viaRaw.1$getPointer(),
  1L,
  wt
)
result.raw.1 <- viaRaw.1$run(0L, 3L)
expect_identical(result.r5.1$train, result.raw.1$train)

# the two 1-based targets really do route to different engine forests
expect_false(identical(result.r5.0$train, result.r5.1$train))

expect_error(build()$setForestWeights(0L, wt), "positive integer")

# NA, NaN, Inf, and negative forest weights are refused at the R surface,
# not only downstream by the bridge's !R_FINITE(w) || w < 0.0 backstop
wt.bad <- wt
wt.bad[1L] <- NA_real_
expect_error(
  build()$setForestWeights(1L, wt.bad),
  "'weights' must be finite and non-negative"
)
wt.bad[1L] <- NaN
expect_error(
  build()$setForestWeights(1L, wt.bad),
  "'weights' must be finite and non-negative"
)
wt.bad[1L] <- Inf
expect_error(
  build()$setForestWeights(1L, wt.bad),
  "'weights' must be finite and non-negative"
)
wt.bad[1L] <- -1
expect_error(
  build()$setForestWeights(1L, wt.bad),
  "'weights' must be finite and non-negative"
)

# --- the mirror, both halves, at all THREE re-creation sites. The
# correct oracle is the SAME (reloaded or restated) object compared against
# itself with the mirror switched off, not a weighted build against a
# different plain build - that comparison is satisfied by the stored-state
# difference alone and stays green even with reapplyForestWeights deleted.
# Forcing forestWeights to list() right before the site under test fires its
# own reapply isolates exactly that call: identical seeds and identical
# starting state mean the only way the continuations can still differ is
# whether the weight was actually reinstalled ---
buildWeighted <- function() {
  s <- build()
  s$setForestWeights(2L, wt)
  invisible(s$run(0L, 3L))
  s$storeState()
  s
}
reloadRds <- function(s) {
  file <- tempfile(fileext = ".rds")
  saveRDS(s, file)
  reloaded <- readRDS(file)
  unlink(file)
  reloaded
}

# site 1: getPointer()'s own create + setState re-creation branch, forced by
# the dead pointer a save/load round trip leaves behind. Both copies come
# from the SAME serialized file, so they start from bit-identical state
weightedFile <- tempfile(fileext = ".rds")
saveRDS(buildWeighted(), weightedFile)
mirrorOn <- readRDS(weightedFile)
mirrorOff <- readRDS(weightedFile)
unlink(weightedFile)
mirrorOff$forestWeights <- list() # before either copy's first engine access
expect_false(identical(
  mirrorOn$run(0L, 1L)$train,
  mirrorOff$run(0L, 1L)$train
))

# site 2: setState()'s own reapply, reached with a pointer that never went
# invalid - no serialization here at all - which isolates that call from the
# getPointer branch exercised above. The field is written directly (as
# forestWeights stores it, 1-based via a hole at index 1) rather than through
# setForestWeights: that method's own C-level install would stick on an
# already-valid pointer regardless of the field, since setState touches
# trees/sigma/k/leaf scale/dispersion/latents/residual df and nothing about
# a forest weight buffer - so a real pre-install would confound this
# comparison rather than isolate setState's reapply call
donorState <- buildWeighted()$state
stateMirrorOn <- build()
stateMirrorOn$forestWeights <- list(NULL, wt)
stateMirrorOff <- build()
stateMirrorOn$setState(donorState)
stateMirrorOff$setState(donorState)
expect_false(identical(
  stateMirrorOn$run(0L, 1L)$train,
  stateMirrorOff$run(0L, 1L)$train
))

# site 3: $copy(), which builds the dupe through a plain $new + $setState
# that on its own would reapply the dupe's still-empty forestWeights field.
# A copy's dupe is always a FRESH pointer (state round trips restore tree
# structure, sigma, and the like, never a chain's rng stream position -
# see inst/common/stateContinuation.R - so a copy's continuation is not
# comparable to the live original's own continuation regardless of forest
# weights). The comparable pair is instead two sources sharing the SAME
# donor state and differing only in forestWeights, each then copied: both
# dupes are equally fresh pointers, so the only remaining difference is
# whichever weight $copy() reapplied
donor <- buildWeighted()
origOn <- build()
origOn$setState(donor$state)
origOn$forestWeights <- list(NULL, wt)
origOff <- build()
origOff$setState(donor$state)
dupeOn <- origOn$copy()
dupeOff <- origOff$copy()
expect_identical(dupeOn$forestWeights, origOn$forestWeights)
expect_false(identical(
  dupeOn$run(0L, 1L)$train,
  dupeOff$run(0L, 1L)$train
))

# --- the channel is a leaf-conditional precision, never a multiplier on
# the reported location - with s = 0 on every row of the treatment forest,
# train still reconstructs exactly from the forest fits and the glue, so
# forest 1's full (data-uninformed but nonzero) contribution is still summed
# in ---
zeroed <- build()
zeroed$setForestWeights(2L, rep(0, n))
result.zeroed <- zeroed$run(0L, 5L)
mu0 <- zeroed$getForestFits(1L)[, 1L]
tau0 <- zeroed$getForestFits(2L)[, 1L]
glue0 <- zeroed$getForestAmplitudes()
zeroed$storeState()
fitScale0 <- zeroed$state[[1L]]$fit.scale
scale0 <- fitScale0[2L] - fitScale0[1L]
shift0 <- scale0 * 0.5 + fitScale0[1L]
bz0 <- ifelse(z != 0, glue0[3L, 1L], glue0[2L, 1L])
expect_equal(
  scale0 * (glue0[1L, 1L] * mu0 + bz0 * tau0) + shift0,
  result.zeroed$train[, 5L],
  tolerance = 1e-10
)
expect_true(abs(mean(tau0)) > 0.1)

# --- the two refusals reachable through the public method: an
# ordinary (non-BCF) sampler, and a multinomial fit's $fit, whose per-forest
# leaf scales come from a log-sum-exp margin rather than an amplitude ---
plain <- dbarts(x, y, control = seededControl())
expect_error(
  plain$setForestWeights(1L, wt),
  "requires a sampler that carries forest amplitudes"
)

set.seed(83)
labels <- factor(rbinom(n, 1L, 0.5))
fitMN <- bart2(
  x,
  labels,
  family = "multinomial",
  keepTrees = TRUE,
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = 2L,
  n.samples = 3L,
  verbose = FALSE
)
expect_error(
  fitMN$fit$setForestWeights(1L, rep(1, n)),
  "log-sum-exp"
)
