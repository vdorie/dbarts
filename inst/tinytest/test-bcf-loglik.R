# extract(type = "loglik") on a BCF-shaped fit - the amplitude-coupled two-
# forest coupling, prognostic plus treatment. pointwiseLogLikelihood branches
# on $family and never sees the coupling, so what has to be pinned is that the
# location it scores at is the COMBINED one the two forests compose into:
# response.shift + a * mu + (basis %*% (b0, b1)) * tau, plus any offset. It
# reads that location off extract(type = "ev"), so every arm below rebuilds it
# from the per-forest channels instead and scores it by hand - gaussian at the
# fit's sigma (weights as precision), probit and logistic at the link of the
# combined latent location.
#
# The per-forest weight ($setForestWeights) is deliberately absent: it is a
# precision factor on one forest's leaf conditionals and enters neither the
# location nor the observation-level density (test-forest-weights.R).

n.samples <- 10L

# One BCF-shaped fit from a dbartsData carrying the two bases.
bcfLoglikFit <- function(data, family = "auto", n.chains = 1L, ...) {
  bart2(
    data,
    family = family,
    n.samples = n.samples,
    n.burn = 5L,
    n.trees = 15L,
    n.chains = n.chains,
    n.threads = 1L,
    verbose = FALSE,
    seed = 613L,
    keepSampler = TRUE,
    ...
  )
}

# The combined location of every kept draw, rebuilt by hand from the per-forest
# fits, that draw's own amplitudes and the fit's bases. extract(type = "forest")
# already carries response.scale, so response.shift and the offset are the only
# terms outside the amplitude recombination. inst/common/recombine.R does the
# same arithmetic for a single caller-env fit and folds in no offset; the local
# form takes the fit and the offset as arguments so every arm can call it.
bcfLocation <- function(fit, offset = NULL) {
  fits <- extract(fit, type = "forest")
  forestNames <- dimnames(fits)[[3L]]
  glueForest <- attr(fit$glue, "forest")
  n.obs <- dim(fits)[2L]
  location <- matrix(
    fit$fit$getCalibration(1L)[1L, "response.shift"],
    nrow(fit$glue),
    n.obs
  )
  for (k in seq_len(fit$n.forests)) {
    basis <- fit$bases[[k]]
    if (is.null(basis)) {
      basis <- matrix(1, n.obs, 1L)
    }
    g <- fit$glue[, glueForest == forestNames[k], drop = FALSE]
    location <- location + (g %*% t(basis)) * fits[,, k]
  }
  if (!is.null(offset)) {
    location <- location + rep(offset, each = nrow(location))
  }
  location
}

set.seed(613)
n <- 120L
p <- 3L
x <- matrix(runif(n * p), n, p)
z <- rbinom(n, 1L, 0.5)
bcfBases <- list(NULL, unname(model.matrix(~ factor(z) - 1L)))
mu <- 2 * sin(pi * x[, 1L]) + x[, 2L]
tau <- 1 + 2 * x[, 3L]
y <- mu + z * tau + rnorm(n, sd = 0.3)

# --- 1. gaussian: the score is dnorm at the combined location and the fit's
# own sigma draw ---
fit.g <- bcfLoglikFit(dbartsData(x, y, bases = bcfBases))
expect_identical(fit.g$n.forests, 2L)
loc.g <- bcfLocation(fit.g)
ev.g <- extract(fit.g, type = "ev")
ll.g <- extract(fit.g, type = "loglik")
expect_equal(dim(loc.g), c(n.samples, n))
expect_equal(dim(ll.g), dim(ev.g))
# the stored location IS the combined one
expect_true(max(abs(loc.g - ev.g)) < 1e-12)
hand.g <- vapply(
  seq_len(n),
  function(i) dnorm(y[i], loc.g[, i], fit.g$sigma, log = TRUE),
  numeric(n.samples)
)
expect_true(max(abs(ll.g - hand.g)) < 1e-12)

# the NEGATIVE half: it is the combined location and not the prognostic
# forest's alone, so dropping the treatment forest's contribution moves the
# score by far more than the bar above
contrib.g <- extract(fit.g, type = "forest", contribution = TRUE)
expect_equal(dim(contrib.g), c(n.samples, n, 2L))
partial.g <- loc.g - contrib.g[,, 2L]
hand.partial <- vapply(
  seq_len(n),
  function(i) dnorm(y[i], partial.g[, i], fit.g$sigma, log = TRUE),
  numeric(n.samples)
)
expect_true(max(abs(ll.g - hand.partial)) > 1)

# --- 2. gaussian with an offset and case weights: the offset rides in the
# location, the weight in the scale ---
offset <- 0.3 * x[, 1L] - 0.15
w <- rep_len(c(1, 4), n)
fit.ow <- bcfLoglikFit(
  dbartsData(x, y, bases = bcfBases, offset = offset, weights = w)
)
loc.ow <- bcfLocation(fit.ow, offset)
ev.ow <- extract(fit.ow, type = "ev")
ll.ow <- extract(fit.ow, type = "loglik")
expect_true(max(abs(loc.ow - ev.ow)) < 1e-12)
# the offset is outside the recombination: without it the rebuild is short by
# exactly the offset, on every draw
expect_true(max(abs(t(ev.ow - bcfLocation(fit.ow)) - offset)) < 1e-12)
hand.ow <- vapply(
  seq_len(n),
  function(i) {
    dnorm(y[i], loc.ow[, i], fit.ow$sigma / sqrt(w[i]), log = TRUE)
  },
  numeric(n.samples)
)
expect_true(max(abs(ll.ow - hand.ow)) < 1e-12)
# and the weight really is in the scale: the unweighted score differs on the
# w = 4 rows and agrees on the w = 1 rows
unweighted <- vapply(
  seq_len(n),
  function(i) dnorm(y[i], loc.ow[, i], fit.ow$sigma, log = TRUE),
  numeric(n.samples)
)
expect_true(max(abs(ll.ow[, w == 1] - unweighted[, w == 1])) < 1e-12)
expect_true(max(abs(ll.ow[, w == 4] - unweighted[, w == 4])) > 0.1)

# --- 3. probit: the combined location is the LATENT one, and the score is the
# bernoulli mass at its normal link ---
eta <- -0.4 + 1.5 * x[, 1L] + 0.8 * z
y.b <- rbinom(n, 1L, pnorm(eta))
fit.p <- bcfLoglikFit(dbartsData(x, y.b, bases = bcfBases), family = "probit")
expect_identical(fit.p$family, "probit")
loc.p <- bcfLocation(fit.p)
ev.p <- extract(fit.p, type = "ev")
ll.p <- extract(fit.p, type = "loglik")
expect_true(max(abs(loc.p - extract(fit.p, type = "bart"))) < 1e-12)
expect_true(max(abs(pnorm(loc.p) - ev.p)) < 1e-12)
hand.p <- vapply(
  seq_len(n),
  function(i) dbinom(y.b[i], 1L, pnorm(loc.p[, i]), log = TRUE),
  numeric(n.samples)
)
expect_true(max(abs(ll.p - hand.p)) < 1e-12)

# a binary offset shifts the same latent location, so the score follows it
fit.po <- bcfLoglikFit(
  dbartsData(x, y.b, bases = bcfBases, offset = rep(0.25, n)),
  family = "probit"
)
loc.po <- bcfLocation(fit.po, rep(0.25, n))
ll.po <- extract(fit.po, type = "loglik")
expect_true(max(abs(loc.po - extract(fit.po, type = "bart"))) < 1e-12)
hand.po <- vapply(
  seq_len(n),
  function(i) dbinom(y.b[i], 1L, pnorm(loc.po[, i]), log = TRUE),
  numeric(n.samples)
)
expect_true(max(abs(ll.po - hand.po)) < 1e-12)

# --- 4. logistic, weighted: the link is the logistic one and integer weights
# are trial counts multiplying the bernoulli score ---
w.l <- rep_len(c(1L, 3L), n)
fit.l <- bcfLoglikFit(
  dbartsData(x, y.b, bases = bcfBases, weights = w.l),
  family = "logistic"
)
expect_identical(fit.l$family, "logistic")
loc.l <- bcfLocation(fit.l)
ev.l <- extract(fit.l, type = "ev")
ll.l <- extract(fit.l, type = "loglik")
expect_true(max(abs(plogis(loc.l) - ev.l)) < 1e-12)
# the link is keyed to the family, not shared with probit
expect_true(max(abs(pnorm(loc.l) - ev.l)) > 0.01)
hand.l <- vapply(
  seq_len(n),
  function(i) w.l[i] * dbinom(y.b[i], 1L, plogis(loc.l[, i]), log = TRUE),
  numeric(n.samples)
)
expect_true(max(abs(ll.l - hand.l)) < 1e-12)

# --- 5. chains: the combined layout stacks chain blocks and each block pairs
# with its own amplitude and sigma draws, and the split layout is a reshape of
# it - the per-forest channels' own 4-d storage shape included ---
fit.mc <- bcfLoglikFit(dbartsData(x, y, bases = bcfBases), n.chains = 2L)
loc.mc <- bcfLocation(fit.mc)
ll.comb <- extract(fit.mc, type = "loglik")
ll.unc <- extract(fit.mc, type = "loglik", combineChains = FALSE)
expect_equal(dim(ll.comb), c(2L * n.samples, n))
expect_equal(dim(ll.unc), c(2L, n.samples, n))
expect_identical(ll.comb[3L, 5L], ll.unc[1L, 3L, 5L])
expect_identical(ll.comb[n.samples + 3L, 5L], ll.unc[2L, 3L, 5L])
hand.mc <- vapply(
  seq_len(n),
  function(i) dnorm(y[i], loc.mc[, i], fit.mc$sigma, log = TRUE),
  numeric(2L * n.samples)
)
expect_true(max(abs(ll.comb - hand.mc)) < 1e-12)

# the STORED layout does not matter: a fit packaged with combineChains = FALSE
# stores glue and forestFits with a leading chain margin, and extract returns
# the same numbers from either
fit.mcu <- bcfLoglikFit(
  dbartsData(x, y, bases = bcfBases),
  n.chains = 2L,
  combineChains = FALSE
)
expect_identical(length(dim(fit.mcu$forestFits)), 4L)
expect_identical(length(dim(fit.mcu$glue)), 3L)
expect_identical(extract(fit.mcu, type = "loglik"), ll.comb)
expect_identical(
  extract(fit.mcu, type = "loglik", combineChains = FALSE),
  ll.unc
)
