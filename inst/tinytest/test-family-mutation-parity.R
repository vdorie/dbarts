# Setter parity for the family/conduit cells that carried no mutation test:
# the offset setter under ordinal and nbinom, and the discrete-time hazard's
# response and offset setters. The oracle is create-vs-mutate - a sampler
# built WITH the quantity and one built without it then given the setter are
# the same sampler, bitwise, on every channel the family reports. Each arm
# pins its engine RNG through dbartsControl(seed = ), so no arm rests on R's
# stream, and each family carries a negative control - the un-mutated sampler
# - so a parity cannot pass by the setter having done nothing.

n.burn <- 10L
n.samples <- 5L
n.trees <- 10L

# --- ordinal: setOffset ---
# The cumulative-probit latent is f(x) + o, so the offset enters every
# threshold Metropolis step and every truncated-normal z draw; the compared
# channels are the latent fits, the threshold draws, and the live z.

set.seed(99)
n <- 120L
x <- matrix(runif(n * 3L), n, 3L)
eta <- 2 * (x[, 1L] - 0.5)
z <- eta + rnorm(n)
codes <- 1L + (z > 0) + (z > 0.8)
offsetOrdinal <- 0.5 * (x[, 2L] - 0.5)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = n.trees,
  updateState = FALSE,
  seed = 13L
)

recordOrdinal <- function(sampler) {
  result <- sampler$run(n.burn, n.samples)
  list(
    train = result$train,
    thresholds = result$thresholds,
    latents = sampler$getLatents()
  )
}

ordinalSampler <- function(offset) {
  if (is.null(offset)) {
    dbarts(x, codes, family = "ordinal", control = control, verbose = FALSE)
  } else {
    dbarts(
      x,
      codes,
      family = "ordinal",
      offset = offset,
      control = control,
      verbose = FALSE
    )
  }
}

ordinal.created <- recordOrdinal(ordinalSampler(offsetOrdinal))

ordinalSet <- ordinalSampler(NULL)
expect_null(ordinalSet$data@offset)
ordinalSet$setOffset(offsetOrdinal)
expect_equal(ordinalSet$data@offset, offsetOrdinal)
ordinal.set <- recordOrdinal(ordinalSet)

expect_identical(ordinal.set$train, ordinal.created$train)
expect_identical(ordinal.set$thresholds, ordinal.created$thresholds)
expect_identical(ordinal.set$latents, ordinal.created$latents)

# negative control: without the setter the same sampler draws elsewhere, so
# the parity above is not the offset being ignored
ordinal.none <- recordOrdinal(ordinalSampler(NULL))
expect_false(isTRUE(all.equal(ordinal.none$train, ordinal.created$train)))

# --- nbinom: setOffset ---
# The count mean is r exp(psi + o), so a log-exposure offset moves the
# Polya-Gamma omega draws and the dispersion step as well as the latent fits.

set.seed(101)
yCount <- rnbinom(n, size = 5L, mu = exp(0.9 * (x[, 1L] - 0.5)))
offsetNbinom <- 0.4 * (x[, 2L] - 0.5)

recordNbinom <- function(sampler) {
  result <- sampler$run(n.burn, n.samples)
  list(
    train = result$train,
    dispersion = result$dispersion,
    latents = sampler$getLatents(),
    r = sampler$getDispersion()
  )
}

nbinomSampler <- function(offset) {
  if (is.null(offset)) {
    dbarts(x, yCount, family = "nbinom", control = control, verbose = FALSE)
  } else {
    dbarts(
      x,
      yCount,
      family = "nbinom",
      offset = offset,
      control = control,
      verbose = FALSE
    )
  }
}

nbinom.created <- recordNbinom(nbinomSampler(offsetNbinom))

nbinomSet <- nbinomSampler(NULL)
expect_null(nbinomSet$data@offset)
nbinomSet$setOffset(offsetNbinom)
expect_equal(nbinomSet$data@offset, offsetNbinom)
nbinom.set <- recordNbinom(nbinomSet)

expect_identical(nbinom.set$train, nbinom.created$train)
expect_identical(nbinom.set$dispersion, nbinom.created$dispersion)
expect_identical(nbinom.set$latents, nbinom.created$latents)
expect_identical(nbinom.set$r, nbinom.created$r)

nbinom.none <- recordNbinom(nbinomSampler(NULL))
expect_false(isTRUE(all.equal(nbinom.none$train, nbinom.created$train)))

# --- the hazard fixture ---
# A hazard sampler is an ordinary binary sampler over the person-period
# design, so its setters take EXPANDED vectors of length N' = sum(terminal
# period), not the n subject-level ones creation ingests. The fixture only has
# to be a valid (time, status) pair - nothing here recovers a hazard - and the
# two statuses share one time vector so the two expansions share one design.

set.seed(517L)
xHaz <- matrix(runif(n * 3L), n, 3L)
fHaz <- 0.9 * sin(pi * xHaz[, 1L]) + 0.6 * xHaz[, 2L]
K <- 5L
timeHaz <- as.double(sample.int(K, n, replace = TRUE))
statusA <- as.double(rbinom(n, 1L, pnorm(fHaz - 1)))
statusB <- statusA
flip <- sample.int(n, 30L)
statusB[flip] <- 1 - statusB[flip]
offsetHazard <- 0.5 * (xHaz[, 2L] - 0.5)

hazardSampler <- function(status, offset = NULL) {
  response <- cbind(timeHaz, status)
  if (is.null(offset)) {
    dbarts(
      xHaz,
      response,
      family = "hazard",
      control = control,
      verbose = FALSE
    )
  } else {
    dbarts(
      xHaz,
      response,
      family = "hazard",
      offset = offset,
      control = control,
      verbose = FALSE
    )
  }
}

recordHazard <- function(sampler) {
  result <- sampler$run(n.burn, n.samples)
  list(train = result$train, latents = sampler$getLatents())
}

# the expansion's own oracle, coded here rather than read off a sampler: the
# rows are subject-major with the period counter restarting at 1, the period
# column is appended LAST, and y'_ik = status_i * 1{k = terminal_i}
hazardLayout <- hazardSampler(statusA)
periodColumn <- hazardLayout$data@x[, ncol(hazardLayout$data@x)]
subjectOf <- cumsum(periodColumn == 1)
terminal <- match(timeHaz, sort(unique(timeHaz)))
expandStatus <- function(status) {
  as.double(periodColumn == terminal[subjectOf] & status[subjectOf] == 1)
}
expect_identical(hazardLayout$data@y, expandStatus(statusA))
yHazardB <- expandStatus(statusB)
expect_identical(hazardSampler(statusB)$data@y, yHazardB)

# --- hazard: setResponse ---
# Equalized arms, both reaching the compared state through the setter. A
# probit constructor cold-starts its latents at z = 2 y - 1 and draws nothing,
# where setResponse redraws z from the truncated normal against the current
# fits; a created sampler and a mutated one are therefore different states by
# construction, and the parity the setter owes is that its post-state depends
# only on the response installed, never on the one it replaced.

hazardSwapped <- hazardSampler(statusA)
hazardSwapped$setResponse(yHazardB)
expect_identical(hazardSwapped$data@y, yHazardB)
hazard.swapped <- recordHazard(hazardSwapped)

hazardSelf <- hazardSampler(statusB)
hazardSelf$setResponse(yHazardB)
hazard.self <- recordHazard(hazardSelf)

expect_identical(hazard.swapped$train, hazard.self$train)
expect_identical(hazard.swapped$latents, hazard.self$latents)

# negative control: the same slot given the OTHER response draws elsewhere
hazardOther <- hazardSampler(statusB)
hazardOther$setResponse(hazardLayout$data@y)
hazard.otherResponse <- recordHazard(hazardOther)
expect_false(isTRUE(all.equal(hazard.otherResponse$train, hazard.self$train)))

# and the cold start really is a different state, which is why the arms above
# are equalized rather than compared against a plain creation
hazard.createdB <- recordHazard(hazardSampler(statusB))
expect_false(isTRUE(all.equal(hazard.createdB$train, hazard.self$train)))

# --- hazard: setOffset ---
# Creation replicates a subject-level offset down the person-period rows; the
# setter takes that expanded vector, and installing it lands where creation's
# replication did. Nothing here is drawn, so the arms compare directly.

hazardOffsetCreated <- hazardSampler(statusA, offsetHazard)
offsetExpanded <- hazardOffsetCreated$data@offset
expect_identical(offsetExpanded, offsetHazard[subjectOf])
hazardOffset.created <- recordHazard(hazardOffsetCreated)

hazardOffsetSet <- hazardSampler(statusA)
expect_null(hazardOffsetSet$data@offset)
hazardOffsetSet$setOffset(offsetExpanded)
expect_equal(hazardOffsetSet$data@offset, offsetExpanded)
hazardOffset.set <- recordHazard(hazardOffsetSet)

expect_identical(hazardOffset.set$train, hazardOffset.created$train)
expect_identical(hazardOffset.set$latents, hazardOffset.created$latents)

hazardOffset.none <- recordHazard(hazardSampler(statusA))
expect_false(isTRUE(all.equal(
  hazardOffset.none$train,
  hazardOffset.created$train
)))
