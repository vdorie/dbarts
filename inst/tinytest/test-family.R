# the family argument: auto dispatch, forced gaussian on 0/1 responses, and
# the public logistic family

set.seed(11)
n <- 300L
x <- matrix(runif(n * 3L), n)
f <- 3 * x[, 1L] - 1.5
y.binary <- rbinom(n, 1L, plogis(f))
y.continuous <- f + rnorm(n, 0, 0.5)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE
)

# auto preserves the existing dispatch
sampler.auto <- dbarts(y.binary ~ x, control = control)
expect_equal(sampler.auto$model@family, "probit")
expect_true(sampler.auto$control@binary)
expect_equal(
  dbarts(y.continuous ~ x, control = control)$model@family,
  "gaussian"
)

# gaussian on a 0/1 response is allowed and fits a continuous model
sampler.gauss <- dbarts(y.binary ~ x, family = "gaussian", control = control)
expect_false(sampler.gauss$control@binary)
samples.gauss <- sampler.gauss$run(50L, 50L)
expect_true(
  all(is.finite(samples.gauss$sigma)) &&
    length(unique(samples.gauss$sigma)) > 1L
)

# binary families need a 0/1 response
expect_error(
  dbarts(y.continuous ~ x, family = "probit", control = control),
  pattern = "requires a response coded 0/1"
)

# a response of one class is fitted, with the one warning saying so and no
# warning to rescale it (dec-B333)
countWarnings <- function(expr) {
  seen <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      seen[[length(seen) + 1L]] <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = seen)
}
y.ones <- rep(1, n)
res <- countWarnings(
  dbarts(y.ones ~ x, family = "probit", control = control)
)
expect_identical(length(res$warnings), 1L)
expect_true(grepl("single class", res$warnings, fixed = TRUE))
res <- countWarnings(
  dbarts(x, rep(0, n), family = "logistic", control = control)
)
expect_identical(length(res$warnings), 1L)
expect_true(grepl("single class", res$warnings, fixed = TRUE))
# the latent sampler runs at one class: every draw is finite and the fitted
# probabilities lean toward the class held
for (cls in c(0, 1)) {
  for (familyName in c("probit", "logistic")) {
    sampler <- suppressWarnings(
      dbarts(x, rep(cls, n), family = familyName, control = control)
    )
    samples <- sampler$run(20L, 50L)
    expect_true(all(is.finite(samples$train)), info = paste(cls, familyName))
    p <- if (familyName == "probit") {
      pnorm(samples$train)
    } else {
      plogis(samples$train)
    }
    expect_true(
      if (cls == 1) mean(p) > 0.5 else mean(p) < 0.5,
      info = paste(cls, familyName)
    )
  }
}
# a hazard fit's binary rows are the caller's subjects, which its refusals
# name; the binary family underneath it is never named
if (requireNamespace("survival", quietly = TRUE)) {
  expect_error(
    dbarts(x, survival::Surv(1 + rpois(n, 2), rep(0, n)), family = "hazard"),
    "family \"hazard\" needs an event; every subject is censored",
    fixed = TRUE
  )
  expect_error(
    dbarts(x, survival::Surv(rep(1, n), rep(1, n)), family = "hazard.logistic"),
    "family \"hazard.logistic\" needs a period at risk without an event",
    fixed = TRUE
  )
  expect_error(
    dbarts(
      x,
      survival::Surv(1 + rpois(n, 2), rbinom(n, 1L, 0.5)),
      family = "hazard",
      variance = TRUE
    ),
    "family \"hazard\" routes precision through its own latent channel",
    fixed = TRUE
  )
}
# a hurdle fit's zero part is a probit fit, which its refusal does not name
expect_error(
  bart(
    x,
    ifelse(y.binary == 1L, exp(y.continuous), 0),
    family = "hurdle.lognormal",
    variance = ~1,
    verbose = FALSE
  ),
  "family = \"hurdle.lognormal\" does not take a variance forest",
  fixed = TRUE
)
# every encoding of a one-class response is fitted alike, on every binary
# family and on "auto" where a categorical encoding resolves to one, with
# exactly the one warning
singleClassResponses <- list(
  logical = rep(TRUE, n),
  factor = factor(rep("a", n), levels = c("a", "b")),
  character = rep("b", n)
)
for (encoding in names(singleClassResponses)) {
  for (familyName in c("probit", "logistic", "auto")) {
    res <- countWarnings(
      dbarts(
        x,
        singleClassResponses[[encoding]],
        family = familyName,
        control = control
      )
    )
    expect_identical(
      length(res$warnings),
      1L,
      info = paste(encoding, familyName)
    )
    expect_true(
      grepl("single class", res$warnings, fixed = TRUE),
      info = paste(encoding, familyName)
    )
  }
}
# a near-constant response that is not one class keeps the warning beside
# its refusal
nearConstantWarnings <- list()
expect_error(
  withCallingHandlers(
    dbarts(x, 1e15 + runif(n, 0, 1e-3), family = "probit", control = control),
    warning = function(w) {
      nearConstantWarnings[[length(nearConstantWarnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  ),
  "requires a response coded 0/1",
  fixed = TRUE
)
expect_identical(length(nearConstantWarnings), 1L)
expect_true(grepl(
  "indistinguishable",
  conditionMessage(nearConstantWarnings[[1L]])
))
rm(
  countWarnings,
  res,
  y.ones,
  singleClassResponses,
  encoding,
  familyName,
  nearConstantWarnings
)

control.bc <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 50L)
sampler.logit <- dbarts(y.binary ~ x, family = "logistic", control = control.bc)
expect_equal(sampler.logit$model@family, "logistic")
expect_inherits(sampler.logit$model@leaf.hyperprior, "dbartsChiHyperprior")
expect_equal(sampler.logit$model@leaf.scale, pi * sqrt(3))

# fits live on the latent logit scale and recover the signal
samples.logit <- sampler.logit$run(300L, 300L)
p.hat <- rowMeans(plogis(samples.logit$train))
expect_true(cor(p.hat, plogis(f)) > 0.8)
expect_true(mean(p.hat[y.binary == 1L]) > mean(p.hat[y.binary == 0L]))

# the family survives save/load: the pointer is recreated from the stored
# state with the same response model, reproducing its trees and latents
source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)
serialized <- tempfile(fileext = ".rds")
state.logit <- sampler.logit$state
saveRDS(sampler.logit, serialized)
sampler.loaded <- readRDS(serialized)
expect_equal(sampler.loaded$model@family, "logistic")
expect_true(all(sampler.loaded$getLatents() > 0)) # omega, not probit z
sampler.loaded$storeState()
statesAgree(sampler.loaded$state, state.logit)
unlink(serialized)

# setControl does not touch the family
newControl <- sampler.logit$control
sampler.logit$setControl(newControl)
expect_equal(sampler.logit$model@family, "logistic")

# setModel cannot silently change the family either
newModel <- sampler.logit$model
newModel@family <- "gaussian"
sampler.logit$setModel(newModel)
expect_equal(sampler.logit$model@family, "logistic")

# bart forwards the argument
fit.bart <- bart(
  y.binary ~ x,
  family = "gaussian",
  n.samples = 30L,
  n.burn = 30L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_true(!is.null(fit.bart$sigma))

# probit refuses weights at the R layer through the wrapper (the bridge keeps
# the same refusal as a backstop for direct-API consumers): a weighted probit
# has no tractable latent-variable form. Logistic weights are covered in
# test-weighted-logistic.R
expect_error(
  bart(
    y.binary ~ x,
    weights = runif(n, 0.5, 1.5),
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 25L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  pattern = "probit models do not support weights"
)

# weights identically 1 are the unweighted likelihood and are treated as
# absent (SuperLearner-style callers pass obsWeights = rep(1, n)
# unconditionally); under the same creation seed the fit matches an
# unweighted one draw for draw
set.seed(7)
fit.unweighted <- bart(
  y.binary ~ x,
  n.samples = 5L,
  n.burn = 5L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
set.seed(7)
fit.unitWeights <- bart(
  y.binary ~ x,
  weights = rep(1, n),
  n.samples = 5L,
  n.burn = 5L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(fit.unitWeights$yhat.train, fit.unweighted$yhat.train)

# the wrappers record the family and transform through its link
fit.probit <- bart(
  y.binary ~ x,
  n.samples = 40L,
  n.burn = 40L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_equal(fit.probit$family, "probit")
expect_equal(extract(fit.probit, "ev"), pnorm(extract(fit.probit, "bart")))

fit.logit <- bart(
  y.binary ~ x,
  family = "logistic",
  n.samples = 40L,
  n.burn = 40L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_equal(fit.logit$family, "logistic")
latents <- extract(fit.logit, "bart")
expect_equal(extract(fit.logit, "ev"), plogis(latents))
expect_equal(
  predict(fit.logit, x, type = "ev"),
  plogis(predict(fit.logit, x, type = "bart"))
)
expect_equal(
  fitted(fit.logit),
  apply(plogis(latents), length(dim(latents)), mean)
)

# fits saved before the family element existed fall back to probit
expect_equal(dbarts:::probabilityFromLatents(0.5, list()), pnorm(0.5))
