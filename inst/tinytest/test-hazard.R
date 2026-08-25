# The discrete-time hazard R surface: person-period ingestion sugar over
# the binary families. The family
# adds no engine code - a hazard fit remaps to probit/logistic before any
# family-keyed switch, so $family reads the binary token and every link-keyed
# generic stays correct; the hazard provenance lives on the $periods marker,
# which survivalProbabilities dispatches on. The bitwise reduction gate lives
# in benchmarks/R/hazard-reduction.R; this file covers the surface.

set.seed(517L)
n <- 250L
p <- 3L
x <- matrix(runif(n * p), n, p)
f <- 0.9 * sin(pi * x[, 1L]) + 0.6 * x[, 2L] - 0.4 * x[, 3L]

# a genuinely discrete grid: a small integer number of periods, baseline
# hazard rising over time, covariate effect through f
K <- 6L
baseline <- seq(-2.0, -0.5, length.out = K)
simulate <- function(link) {
  invlink <- if (link == "logistic") plogis else pnorm
  event <- rep(K + 1L, n)
  for (i in seq_len(n)) {
    for (k in seq_len(K)) {
      if (runif(1L) < invlink(baseline[k] + f[i])) {
        event[i] <- k
        break
      }
    }
  }
  cens <- sample.int(K, n, replace = TRUE)
  status <- as.double(event <= cens & event <= K)
  time <- as.double(pmin(event, cens, K))
  list(time = time, status = status)
}

set.seed(11L)
d <- simulate("probit")

fitArgs <- list(
  n.trees = 40L,
  n.burn = 80L,
  n.samples = 150L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 13L
)

# ---- B2: named training predictors (the ordinary matrix-interface spelling,
# bart.Rd:75) are the unnamed fixture's counterpart. The re-expanded design's
# appended period column used to go unnamed on the matrix-training and
# matrix-newdata branches, so a named fit died in predict's column-name match;
# a data.frame training input hit the same branch once bart2 converted it to
# a named matrix. Small and cheap: only the naming path is under test.
set.seed(929L)
n.named <- 60L
x.named <- matrix(runif(n.named * p), n.named, p)
colnames(x.named) <- c("x1", "x2", "x3")
f.named <- 0.9 *
  sin(pi * x.named[, 1L]) +
  0.6 * x.named[, 2L] -
  0.4 * x.named[, 3L]
event.named <- rep(K + 1L, n.named)
for (i in seq_len(n.named)) {
  for (k in seq_len(K)) {
    if (runif(1L) < pnorm(baseline[k] + f.named[i])) {
      event.named[i] <- k
      break
    }
  }
}
cens.named <- sample.int(K, n.named, replace = TRUE)
d.named <- list(
  time = as.double(pmin(event.named, cens.named, K)),
  status = as.double(event.named <= cens.named & event.named <= K)
)
namedFitArgs <- modifyList(
  fitArgs,
  list(n.trees = 10L, n.burn = 10L, n.samples = 10L)
)

fit.namedMat <- do.call(
  bart2,
  c(
    list(x.named, cbind(d.named$time, d.named$status), family = "hazard"),
    namedFitArgs
  )
)
expect_equal(
  colnames(fit.namedMat$fit$data@x),
  c("x1", "x2", "x3", "period")
)
sp.named <- survivalProbabilities(fit.namedMat)
expect_equal(
  dim(sp.named),
  c(nrow(fit.namedMat$yhat.train), length(fit.namedMat$periods), n.named)
)

xn.named <- matrix(runif(4L * p), 4L, p)
colnames(xn.named) <- c("x1", "x2", "x3")
sp.named.mat <- survivalProbabilities(fit.namedMat, newdata = xn.named)
expect_equal(
  dim(sp.named.mat),
  c(nrow(fit.namedMat$yhat.train), length(fit.namedMat$periods), 4L)
)
sp.named.df <- survivalProbabilities(
  fit.namedMat,
  newdata = as.data.frame(xn.named)
)
expect_equal(sp.named.df, sp.named.mat)

# a data.frame training input is coerced to a named matrix by bart2, so it
# takes the same branch and must work identically
fit.namedDf <- do.call(
  bart2,
  c(
    list(
      as.data.frame(x.named),
      cbind(d.named$time, d.named$status),
      family = "hazard"
    ),
    namedFitArgs
  )
)
sp.namedDf <- survivalProbabilities(fit.namedDf)
expect_equal(dim(sp.namedDf), dim(sp.named))

rm(
  n.named,
  x.named,
  f.named,
  event.named,
  cens.named,
  d.named,
  namedFitArgs,
  fit.namedMat,
  sp.named,
  xn.named,
  sp.named.mat,
  sp.named.df,
  fit.namedDf,
  sp.namedDf,
  i,
  k
)

# ---- both tokens, marker, and the remap ($family reads the binary token) ----

fit.probit <- do.call(
  bart2,
  c(list(x, cbind(d$time, d$status), family = "hazard"), fitArgs)
)
expect_identical(fit.probit[["family"]], "probit")
expect_false(is.null(fit.probit[["periods"]]))
expect_equal(fit.probit$periods, sort(unique(d$time)))
# a hazard fit carries no sigma (it is a binary fit under the hood)
expect_true(is.null(fit.probit[["sigma"]]))

fit.logit <- do.call(
  bart2,
  c(list(x, cbind(d$time, d$status), family = "hazard.logistic"), fitArgs)
)
expect_identical(fit.logit[["family"]], "logistic")
expect_false(is.null(fit.logit[["periods"]]))

# hazard.probit is an accepted alias for the probit link
fit.alias <- do.call(
  bart2,
  c(list(x, cbind(d$time, d$status), family = "hazard.probit"), fitArgs)
)
expect_identical(fit.alias[["family"]], "probit")
# same seed/design, so byte-identical to the "hazard" token
expect_equal(fit.alias$yhat.train, fit.probit$yhat.train)

# ---- the plogis lock: the marker design keeps link dispatch correct ----
# type = "ev" transforms the latent by the fit's $family link. A logit-hazard
# fit must use plogis (not the default pnorm) - the wrong number with no error
# had the hazard token been recorded in $family instead of the binary one.
xn <- matrix(runif(4L * p), 4L, p)
xn.exp <- cbind(xn[rep(seq_len(4L), K), ], rep(seq_len(K), each = 4L))
lat.logit <- predict(fit.logit, xn.exp, type = "bart")
ev.logit <- predict(fit.logit, xn.exp, type = "ev")
expect_equal(ev.logit, plogis(lat.logit))
lat.probit <- predict(fit.probit, xn.exp, type = "bart")
ev.probit <- predict(fit.probit, xn.exp, type = "ev")
expect_equal(ev.probit, pnorm(lat.probit))

# ---- Surv ingestion gives an identical fit to the two-column form ----
surv <- structure(
  cbind(time = d$time, status = d$status),
  class = "Surv",
  type = "right"
)
fit.surv <- do.call(bart2, c(list(x, surv, family = "hazard"), fitArgs))
expect_identical(fit.surv[["family"]], "probit")
expect_equal(fit.surv$yhat.train, fit.probit$yhat.train)

# ---- survivalProbabilities shape, range, monotonicity (training) ----
times <- fit.probit$periods
sp <- survivalProbabilities(fit.probit, times)
n.draws <- nrow(fit.probit$yhat.train)
expect_equal(dim(sp), c(n.draws, length(times), n))
expect_true(all(sp >= 0 & sp <= 1))
sp.mean <- apply(sp, c(2L, 3L), mean)
# every posterior-mean curve is non-increasing in t
expect_true(all(apply(sp.mean, 2L, function(curve) all(diff(curve) <= 1e-9))))
# every individual draw's curve is non-increasing too (a cumulative product of
# (1 - hazard) terms in [0, 1])
expect_true(all(sp[, -1L, ] <= sp[, -length(times), ] + 1e-9))
# the default times argument is the training grid
expect_equal(survivalProbabilities(fit.probit), sp)

# ---- survivalProbabilities on newdata expands to the requested grid ----
sp.new <- survivalProbabilities(fit.probit, times, newdata = xn)
expect_equal(dim(sp.new), c(n.draws, length(times), 4L))
expect_true(all(sp.new >= 0 & sp.new <= 1))
# a subset of horizons selects the matching survival columns
sp.sub <- survivalProbabilities(fit.probit, times[c(1L, 3L)], newdata = xn)
expect_equal(sp.sub[, 1L, ], sp.new[, 1L, ])
expect_equal(sp.sub[, 2L, ], sp.new[, 3L, ])

# ---- keepTrees is required (the training design is ragged) ----
fit.notrees <- do.call(
  bart2,
  c(
    list(x, cbind(d$time, d$status), family = "hazard"),
    modifyList(fitArgs, list(keepTrees = FALSE))
  )
)
expect_false(is.null(fit.notrees[["periods"]]))
expect_error(
  survivalProbabilities(fit.notrees, times),
  "keepTrees"
)

# ---- breaks: integer count and explicit boundary vector ----
fit.k <- do.call(
  bart2,
  c(list(x, cbind(d$time, d$status), family = "hazard", breaks = 3L), fitArgs)
)
expect_true(length(fit.k$periods) <= 3L)
fit.b <- do.call(
  bart2,
  c(
    list(
      x,
      cbind(d$time, d$status),
      family = "hazard",
      breaks = c(0, 2, 4, 6)
    ),
    fitArgs
  )
)
expect_equal(fit.b$periods, c(2, 4, 6))

# ---- the N' row guard refuses an over-fine grid, naming the levers ----
expect_error(
  do.call(
    bart2,
    c(
      list(x, cbind(d$time, d$status), family = "hazard", max.rows = 50),
      fitArgs
    )
  ),
  "breaks"
)

# ---- offset / weight replication follow the binary policy ----
# a non-unit probit weight is refused after replication (the binary policy);
# an integer-count logistic weight is accepted
expect_error(
  do.call(
    bart2,
    c(
      list(
        x,
        cbind(d$time, d$status),
        family = "hazard",
        weights = runif(n) + 0.5
      ),
      fitArgs
    )
  ),
  "weight"
)
fit.wt <- do.call(
  bart2,
  c(
    list(
      x,
      cbind(d$time, d$status),
      family = "hazard.logistic",
      weights = rep(1L, n)
    ),
    fitArgs
  )
)
expect_identical(fit.wt[["family"]], "logistic")

# ---- refusals: interface, conflicting family, subset, test, non-hazard ----
# the formula interface is out of scope (matrix interface only, like aft)
expect_error(
  bart2(d$time ~ x, family = "hazard", verbose = FALSE),
  "matrix interface"
)
# a Surv response with an explicitly conflicting family errors
expect_error(
  bart2(x, surv, family = "gaussian", n.chains = 1L, verbose = FALSE),
  "aft|hazard"
)
# subset and test are refused for hazard fits
expect_error(
  dbarts(x, cbind(d$time, d$status), family = "hazard", subset = 1:10),
  "subset"
)
expect_error(
  dbarts(
    x,
    cbind(d$time, d$status),
    family = "hazard",
    test = matrix(runif(3L * p), 3L, p)
  ),
  "test"
)
# a non-hazard, non-aft fit is refused by survivalProbabilities
fit.gauss <- bart2(
  x,
  f + rnorm(n),
  n.trees = 20L,
  n.burn = 30L,
  n.samples = 30L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 7L
)
expect_error(survivalProbabilities(fit.gauss, times = 1), "aft")
# xbart and rbart_vi do not fit hazard (their family vectors are the refusal)
expect_error(xbart(x, f + rnorm(n), family = "hazard"), "should be one of")
expect_error(
  rbart_vi(
    (f + rnorm(n)) ~ x,
    family = "hazard",
    group.by = rep(1:4, length.out = n)
  ),
  "should be one of"
)

# ---- seeded recovery: fitted survival tracks the truth ----
set.seed(202L)
d.rec <- simulate("logistic")
fit.rec <- do.call(
  bart2,
  c(
    list(x, cbind(d.rec$time, d.rec$status), family = "hazard.logistic"),
    fitArgs
  )
)
# true survival curve S(k | x) = prod_{j<=k} (1 - h(j | x)) on the grid
trueHaz <- plogis(outer(baseline, f, "+")) # K x n
trueSurv <- apply(1 - trueHaz, 2L, cumprod) # K x n
estSurv <- apply(
  survivalProbabilities(fit.rec, seq_len(K)),
  c(2L, 3L),
  mean
) # K x n
# the fitted survival surface tracks the truth across (period, subject) cells
expect_true(cor(as.vector(estSurv), as.vector(trueSurv)) > 0.75)
# higher-risk subjects (larger f) survive less: the high-f half has lower mean
# fitted survival than the low-f half (a robust directional check; f is
# nonlinear, so a per-cell correlation understates the recovery)
highRisk <- f > median(f)
expect_true(mean(estSurv[, highRisk]) < mean(estSurv[, !highRisk]))
