# Grouped (random-intercept) AFT survival surface on rbart_vi. The
# engine composition GroupedResponse(AFTResponse) is exercised here through
# the R surface: a Surv / two-column response on the formula LHS with
# family = "aft", and survivalProbabilities.rbart that includes the drawn
# intercepts in S(t | x, group).

set.seed(21L)
n <- 150L
p <- 3L
x <- matrix(runif(n * p), n, p)
f <- 1.5 * sin(pi * x[, 1L]) + x[, 2L] - 0.5 * x[, 3L]
n.g <- 6L
groups <- rep(seq_len(n.g), length.out = n)
group.shift <- rnorm(n.g, 0, 0.7)[groups]
sigma.true <- 0.5
log.t <- f + group.shift + sigma.true * rnorm(n)

# ---- reduction gate: all-uncensored aft == gaussian on log T, bitwise ----
# time round-trips exp/log, so the gaussian arm is fit on log(time) - the
# identical buffer the aft path derives - and the two streams cannot diverge
# unless something consumed RNG differently
time <- exp(log.t)
log.time <- log(time)
surv.all <- structure(
  cbind(time = time, status = rep(1, n)),
  class = "Surv",
  type = "right"
)

commonArgs <- list(
  group.by = groups,
  n.samples = 60L,
  n.burn = 40L,
  n.thin = 1L,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  verbose = FALSE,
  seed = 99L
)

fit.gauss <- do.call(rbart_vi, c(list(formula = log.time ~ x), commonArgs))
fit.aft <- do.call(
  rbart_vi,
  c(list(formula = surv.all ~ x, family = "aft"), commonArgs)
)

expect_identical(fit.aft[["family"]], "aft")
expect_identical(fit.gauss$ranef, fit.aft$ranef)
expect_identical(fit.gauss$tau, fit.aft$tau)
expect_identical(fit.gauss$sigma, fit.aft$sigma)
expect_identical(fit.gauss$yhat.train, fit.aft$yhat.train)

# ---- survivalProbabilities.rbart draws (train): dims, range, monotone ----

set.seed(8L)
cens <- f + 0.4 + sigma.true * rnorm(n)
status <- as.numeric(log.t <= cens)
obs.time <- exp(ifelse(status == 1, log.t, cens))
surv <- structure(
  cbind(time = obs.time, status = status),
  class = "Surv",
  type = "right"
)

fit <- rbart_vi(
  surv ~ x,
  group.by = groups,
  family = "aft",
  n.samples = 80L,
  n.burn = 60L,
  n.thin = 1L,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  verbose = FALSE,
  seed = 3L
)
expect_identical(fit[["family"]], "aft")
expect_false(is.null(fit[["sigma"]]))

times <- c(0.5, 1, 2, 4)
sp <- survivalProbabilities(fit, times)
n.draws <- nrow(fit$yhat.train)
expect_equal(dim(sp), c(n.draws, length(times), n))
expect_true(all(sp >= 0 & sp <= 1))
sp.mean <- apply(sp, c(2L, 3L), mean)
expect_true(all(apply(sp.mean, 2L, function(curve) all(diff(curve) <= 1e-8))))
# every individual draw's curve is a normal tail, so monotone too
expect_true(all(sp[, -1L, ] <= sp[, -length(times), ] + 1e-8))

# ---- the drawn intercepts enter the curve: a high-intercept group has a
# ---- uniformly higher S(t) than a low-intercept group at the same x ----

ord <- order(fit$ranef.mean)
lowGroup <- names(fit$ranef.mean)[ord[1L]]
highGroup <- names(fit$ranef.mean)[ord[length(ord)]]
xRep <- rbind(x[1L, ], x[1L, ])
gPair <- factor(c(lowGroup, highGroup), levels = levels(fit$group.by))
sp.pair <- survivalProbabilities(fit, times, newdata = xRep, group.by = gPair)
sp.pair.mean <- apply(sp.pair, c(2L, 3L), mean) # times x 2
# the high-intercept group's curve dominates at every time; if the
# intercepts had been dropped the two curves would coincide
expect_true(all(sp.pair.mean[, 2L] >= sp.pair.mean[, 1L] - 1e-8))
expect_true(any(sp.pair.mean[, 2L] > sp.pair.mean[, 1L] + 1e-6))

# ---- newdata + group.by with an unseen group draws its intercept ----

x.new <- matrix(runif(4L * p), 4L, p)
g.new <- factor(c(1L, 2L, 3L, 99L), levels = c(levels(fit$group.by), "99"))
sp.new <- suppressWarnings(
  survivalProbabilities(fit, times, newdata = x.new, group.by = g.new)
)
expect_equal(dim(sp.new), c(n.draws, length(times), 4L))
expect_true(all(sp.new >= 0 & sp.new <= 1))

# ---- multi-chain conventions and the sigma-to-draw alignment ----

fit.chains <- rbart_vi(
  surv ~ x,
  group.by = groups,
  family = "aft",
  n.samples = 40L,
  n.burn = 40L,
  n.thin = 1L,
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 25L,
  verbose = FALSE,
  seed = 7L
)
sp.comb <- survivalProbabilities(fit.chains, times)
expect_equal(dim(sp.comb), c(80L, length(times), n))
sp.unc <- survivalProbabilities(fit.chains, times, combineChains = FALSE)
expect_equal(dim(sp.unc), c(2L, 40L, length(times), n))
# the combined result stacks the chains sample-major
expect_equal(sp.comb[1:40, , ], sp.unc[1L, , , ])
expect_equal(sp.comb[41:80, , ], sp.unc[2L, , , ])

# a pinned value is the exact normal upper tail of the ev draw (BART plus the
# drawn intercept) over the aligned per-draw sigma
ev.unc <- extract(
  fit.chains,
  type = "ev",
  sample = "train",
  combineChains = FALSE
)
sig.unc <- fit.chains$sigma
if (is.null(dim(sig.unc))) {
  sig.unc <- dbarts:::uncombineChains(as.vector(sig.unc), 2L)
}
expect_equal(
  sp.unc[2L, 17L, 3L, 5L],
  pnorm(
    (log(times[3L]) - ev.unc[2L, 17L, 5L]) / sig.unc[2L, 17L],
    lower.tail = FALSE
  )
)

# ---- recovery-under-censoring smoke: the BART component tracks the signal
# ---- and sigma stays in a sane band ----

expect_true(cor(fit$yhat.train.mean, f) > 0.6)
expect_true(mean(fit$sigma) > 0.3 && mean(fit$sigma) < 0.9)

# ---- custom-prior (callback R-loop) aft path runs and returns finite draws --

fit.custom <- rbart_vi(
  surv ~ x,
  group.by = groups,
  family = "aft",
  prior = function(x, rel.scale) dnorm(x, 0, rel.scale, log = TRUE),
  n.samples = 40L,
  n.burn = 30L,
  n.thin = 1L,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  verbose = FALSE,
  seed = 5L
)
expect_identical(fit.custom[["family"]], "aft")
expect_true(all(is.finite(fit.custom$yhat.train)))
expect_true(all(is.finite(fit.custom$sigma)))
expect_true(all(is.finite(fit.custom$tau)))
sp.custom <- survivalProbabilities(fit.custom, times)
expect_true(all(is.finite(sp.custom) & sp.custom >= 0 & sp.custom <= 1))

# ---- refusals ----

# survivalProbabilities on a non-aft rbart fit
fit.gaussian <- rbart_vi(
  log.t ~ x,
  group.by = groups,
  n.samples = 20L,
  n.burn = 20L,
  n.thin = 1L,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  verbose = FALSE,
  seed = 1L
)
expect_error(survivalProbabilities(fit.gaussian, times = 1), "aft")

# weights and subset are unsupported for aft
expect_error(
  rbart_vi(
    surv ~ x,
    group.by = groups,
    family = "aft",
    weights = runif(n) + 0.5,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "weights"
)
expect_error(
  rbart_vi(
    surv ~ x,
    group.by = groups,
    family = "aft",
    subset = seq_len(50L),
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "subset"
)

# non-positive times refused
expect_error(survivalProbabilities(fit, times = c(-1, 1)), "positive")

# ---- the public (direct-user) dbarts aft refusal stays intact: the guard
# ---- opens only for the internal rbart channel (control@bartcore.survival) --

expect_error(dbarts(log.t ~ x, family = "aft"), "matrix interface")
# a Surv left-hand side through the public formula interface is still refused
# in dbartsData, before any arithmetic trips survival's Ops.Surv guard
expect_error(
  dbarts(surv ~ x1, data = list(surv = surv.all, x1 = x[, 1L])),
  "matrix interface"
)
