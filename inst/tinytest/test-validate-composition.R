# dbartsValidateComposition: simulation-based calibration over a host-supplied
# one-sweep step.
#
# The oracle is ANALYTIC - a conjugate normal-normal model whose posterior is
# closed form - so no sampler runs here and a verdict is the validator's own.
# theta ~ N(0, 1) and y_i | theta ~ N(theta, 1), i = 1..n, give a posterior
# N(v * sum(y), v) with v = 1 / (1 + n). Three steps share it: the EXACT draw
# (known good), the posterior MEAN substituted for the draw (the measured
# defect class), and a draw at TWICE the posterior sd (over-dispersion, so the
# diagnostic is shown to catch both directions).
#
# The cells pass n.thin = 1L because the exact conjugate step draws the
# posterior outright, leaving consecutive states independent, and n.burn = 0L
# because their init is exact at theta0, leaving the chain stationary from the
# first state. The API's own defaults, 30 and 200, serve hosts whose step mixes
# and whose init is only approximate.

nObs <- 5L
postVar <- 1 / (1 + nObs)

drawPrior <- function() c(theta = rnorm(1L))
simulateData <- function(theta0) rnorm(nObs, theta0[["theta"]], 1)
initState <- function(theta0, data) {
  list(
    theta = theta0[["theta"]],
    mean = postVar * sum(data),
    sd = sqrt(postVar)
  )
}
stepExact <- function(state) {
  state$theta <- rnorm(1L, state$mean, state$sd)
  state
}
stepMean <- function(state) {
  state$theta <- state$mean
  state
}
stepWide <- function(state) {
  state$theta <- rnorm(1L, state$mean, 2 * state$sd)
  state
}
theta <- function(state) c(theta = state$theta)

validate <- function(
  step = stepExact,
  functionals = theta,
  n.replications = 200L,
  n.draws = 200L,
  n.thin = 1L,
  n.burn = 0L,
  ...
) {
  dbartsValidateComposition(
    drawPrior,
    simulateData,
    initState,
    step,
    functionals,
    n.replications = n.replications,
    n.draws = n.draws,
    n.thin = n.thin,
    n.burn = n.burn,
    ...
  )
}

# the se of the mean of R ranks uniform on {0, ..., L} is
# sqrt(((L + 1)^2 - 1)/12/R), 4.1 at R = L = 200, so this slack is ~5 se; a
# rank rule that counts ties on the wrong side of theta0 shifts the mean by
# tens of ranks
meanRankGap <- function(fit) {
  max(abs(colMeans(fit$ranks) - fit$mean.rank.target))
}

# --- the discriminating trio

set.seed(20L)
good <- validate(stepExact)
expect_true(good$pass)
expect_equal(good$verdicts$verdict, "PASS")
expect_equal(dim(good$ranks), c(200L, 1L))
expect_equal(colnames(good$ranks), "theta")
expect_true(all(good$ranks >= 0L & good$ranks <= 200L))
expect_equal(good$L, 200L)
expect_equal(good$alpha, 0.05)
expect_equal(good$alpha.adjusted, 0.05)
expect_equal(good$mean.rank.target, 100)
expect_true(good$verdicts$ecdf.diff <= good$verdicts$band)
expect_true(meanRankGap(good) < 20)

goodOutput <- capture.output(print(good))
expect_true(any(grepl("PASS", goodOutput, fixed = TRUE)))
expect_false(any(grepl("FLAG", goodOutput, fixed = TRUE)))
expect_true(any(grepl("200 replications", goodOutput, fixed = TRUE)))

set.seed(20L)
meanForDraw <- validate(stepMean)
expect_false(meanForDraw$pass)
expect_equal(meanForDraw$verdicts$verdict, "FLAG")
expect_true(meanForDraw$verdicts$ecdf.diff > meanForDraw$verdicts$band)
# the posterior mean is the same number every draw, so theta0 lands at one end
# of the L draws or the other and nowhere between
expect_true(all(meanForDraw$ranks %in% c(0L, 200L)))
expect_true(any(grepl(
  "FLAG",
  capture.output(print(meanForDraw)),
  fixed = TRUE
)))

set.seed(20L)
overDispersed <- validate(stepWide)
expect_false(overDispersed$pass)
expect_equal(overDispersed$verdicts$verdict, "FLAG")
expect_true(overDispersed$verdicts$ecdf.diff > overDispersed$verdicts$band)
# over-dispersed draws straddle theta0, so its rank piles up in the middle
expect_true(mean(abs(overDispersed$ranks - 100)) < 45)

# --- DERIVED functionals, neither of them a component of the prior draw
#
# 'absTheta' and 'positive' appear nowhere in what drawPrior returns, which is
# why the ranked quantity has to be the functional evaluated at the init state
# and never a name match against theta0. 'positive' is also ATOM-valued: every
# draw ties with 0 or with 1, the case the discrete rank rule exists for.

derived <- function(state) {
  c(absTheta = abs(state$theta), positive = as.numeric(state$theta > 0))
}

set.seed(31L)
derivedFit <- validate(stepExact, derived)
expect_equal(colnames(derivedFit$ranks), c("absTheta", "positive"))
expect_true(derivedFit$pass)
expect_equal(derivedFit$verdicts$verdict, c("PASS", "PASS"))
expect_true(meanRankGap(derivedFit) < 20)
# the tie-break is what keeps an atom-valued functional off the rank ends
expect_true(sum(derivedFit$ranks[, "positive"] %in% c(0L, 200L)) < 20L)

# --- the Bonferroni, over enough functionals to see it
#
# Eight INDEPENDENT conjugate scalars in one state: correlated functionals flag
# together, so only independent ones distinguish a per-call correction from
# none. M = 8 also pushes the band simulation past its default of 2000, the
# 1 - alpha quantile needing at least 20 null draws in its tail.

nParam <- 8L
paramNames <- paste0("t", seq_len(nParam))
drawPriorMany <- function() setNames(rnorm(nParam), paramNames)
simulateMany <- function(theta0) {
  matrix(rnorm(nParam * nObs, rep(theta0, each = nObs)), nObs, nParam)
}
initMany <- function(theta0, data) {
  list(
    theta = unname(theta0),
    mean = postVar * colSums(data),
    sd = sqrt(postVar)
  )
}
stepMany <- function(state) {
  state$theta <- rnorm(nParam, state$mean, state$sd)
  state
}
functionalsMany <- function(state) setNames(state$theta, paramNames)

set.seed(23L)
many <- dbartsValidateComposition(
  drawPriorMany,
  simulateMany,
  initMany,
  stepMany,
  functionalsMany,
  n.replications = 200L,
  n.draws = 200L,
  n.thin = 1L,
  n.burn = 0L
)
expect_equal(colnames(many$ranks), paramNames)
expect_true(many$pass)
expect_equal(many$alpha.adjusted, 0.05 / nParam)
expect_true(meanRankGap(many) < 20)
# the band each functional is judged against is the CORRECTED one, not the
# nominal band it would earn on its own; at this seed one functional's
# statistic sits between the two, so an uncorrected call would flag a
# composition that is calibrated
corrected <- dbarts:::rankUniformity(many$ranks[, 1L], 200L, alpha = 0.05 / 8)
nominal <- dbarts:::rankUniformity(many$ranks[, 1L], 200L, alpha = 0.05)
expect_equal(many$verdicts$band[1L], corrected$ecdfBand)
expect_true(corrected$ecdfBand > nominal$ecdfBand)
expect_true(any(many$verdicts$ecdf.diff > nominal$ecdfBand))
expect_true(all(many$verdicts$ecdf.diff <= many$verdicts$band))

# --- the random number stream
#
# rankUniformity fixes its own stream so that a band is reproducible from the
# ranks alone. Without the restore every call would leave the caller's stream
# where that fixed stream ended - the same state whatever the caller had set,
# the number of draws taken being fixed by R, L and nSim.

set.seed(101L)
firstRun <- validate(n.replications = 30L, n.draws = 40L)
firstSeed <- .Random.seed
set.seed(202L)
secondRun <- validate(n.replications = 30L, n.draws = 40L)
secondSeed <- .Random.seed
expect_false(identical(firstSeed, secondSeed))
expect_false(identical(firstRun$ranks, secondRun$ranks))

# 'seed' seeds the whole run, so two calls agree rank for rank
set.seed(101L)
seeded <- validate(n.replications = 30L, n.draws = 40L, seed = 7L)
set.seed(202L)
reseeded <- validate(n.replications = 30L, n.draws = 40L, seed = 7L)
expect_identical(seeded$ranks, reseeded$ranks)

# the guard restores the absence of a stream too, not only its state
rm(".Random.seed", envir = globalenv())
expect_true(dbarts:::rankUniformity(rep(0:9, 20L), 9L, nSim = 50L)$pass)
expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
set.seed(3L)

# --- input validation

expect_error(
  dbartsValidateComposition(1, simulateData, initState, stepExact, theta),
  "'drawPrior' must be a function"
)
expect_error(
  dbartsValidateComposition(drawPrior, simulateData, initState, stepExact, 1),
  "'functionals' must be a function"
)
expect_error(validate(n.replications = 1L), "'n.replications' must be a single")
expect_error(validate(n.draws = 2.5), "'n.draws' must be a single")
expect_error(validate(n.thin = 0L), "'n.thin' must be a single")
expect_error(validate(n.burn = -1L), "'n.burn' must be a single")
expect_error(validate(alpha = 0), "'alpha' must be a single number")
expect_error(validate(alpha = c(0.05, 0.1)), "'alpha' must be a single number")
expect_error(validate(seed = "one"), "'seed' must be a single number")

# a functional that changes length or name between the initial state and a
# stepped one ranks something other than what it drew
stepFlagged <- function(state) {
  state <- stepExact(state)
  state$stepped <- TRUE
  state
}
grown <- function(state) {
  if (is.null(state$stepped)) c(theta = state$theta) else c(1, 2)
}
renamed <- function(state) {
  if (is.null(state$stepped)) c(theta = state$theta) else c(other = state$theta)
}
sameOnes <- "the ranked functionals must be the same ones"
expect_error(validate(stepFlagged, grown, 2L, 2L), sameOnes)
expect_error(validate(stepFlagged, renamed, 2L, 2L), sameOnes)

numericVector <- "'functionals' must return a finite, non-empty numeric vector"
expect_error(
  validate(stepExact, function(state) "theta", 2L, 2L),
  numericVector
)
expect_error(
  validate(stepExact, function(state) c(theta = NA_real_), 2L, 2L),
  numericVector
)

# an unnamed functional is labelled positionally rather than refused
unnamedFit <- validate(stepExact, function(state) state$theta, 2L, 2L)
expect_equal(colnames(unnamedFit$ranks), "f1")
