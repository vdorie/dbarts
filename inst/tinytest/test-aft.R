# AFT log-normal survival family on the bartcore engine (src/bartcore/),
# through a two-column (time, status) response. The exact-posterior gate lives
# in benchmarks/R/aft-exact.R.

set.seed(21L)
n <- 200L
p <- 3L
x <- matrix(runif(n * p), n, p)
f <- 1.5 * sin(pi * x[, 1L]) + x[, 2L] - 0.5 * x[, 3L]
sigma.true <- 0.5
log.t <- f + sigma.true * rnorm(n)

# a seeded control makes each chain's Mersenne twister deterministic, so the
# reduction below can compare bitwise
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE,
  seed = 271L
)

# y is a log time; the sampler logs exp(y) itself, which need not return y
# bitwise, so exact comparisons read the sampler's own data@y
aftSampler <- function(y, status, weights = NULL) {
  dbarts(
    x,
    cbind(exp(y), status),
    weights = weights,
    family = "aft",
    control = control
  )
}

# ---- reduction: all-uncensored aft == gaussian on log T, bitwise ----

bc.a <- aftSampler(log.t, rep(1, n)) # every observation an event
expect_identical(bc.a$model@family, "aft")
res.a <- bc.a$run(100L, 200L)

samp.g <- dbarts(x, bc.a$data@y, control = control)
expect_identical(samp.g$model@family, "gaussian")
res.g <- samp.g$run(100L, 200L)
# the engine's own family: aft imputes a latent log-time per row, all of them
# observed here, where a gaussian engine keeps none
expect_identical(bc.a$getLatents(), bc.a$data@y)
expect_null(samp.g$getLatents())

expect_identical(res.g$train, res.a$train)
expect_identical(res.g$sigma, res.a$sigma)

# ---- recovery under censoring, and the naive downward bias corrected ----

recover <- function(censor.rate) {
  set.seed(11L)
  # censor by an independent time chosen to hit the target rate in expectation
  cens.time <- f +
    quantile(sigma.true * rnorm(2000L), 1 - censor.rate) +
    sigma.true * rnorm(n)
  status <- as.numeric(log.t <= cens.time)
  obs.log.t <- ifelse(status == 1, log.t, cens.time)

  bc <- aftSampler(obs.log.t, status)
  res <- bc$run(200L, 400L)
  fit.aft <- rowMeans(res$train)

  # ignoring the censoring underestimates the mean log-time
  samp.naive <- dbarts(x, obs.log.t, control = control)
  res.naive <- samp.naive$run(200L, 400L)
  list(
    rate = mean(status == 0),
    cor = cor(fit.aft, f),
    sigma = mean(res$sigma),
    mean.aft = mean(fit.aft),
    mean.naive = mean(rowMeans(res.naive$train)),
    lat = bc$getLatents(),
    status = status,
    obs = bc$data@y
  )
}

for (rate in c(0.2, 0.5)) {
  r <- recover(rate)
  # fit still tracks the signal
  expect_true(r$cor > 0.8)
  # sigma stays in a sane band around the truth (loose at test scale)
  expect_true(r$sigma > 0.3 && r$sigma < 0.8)
  # AFT recovers a higher mean log-time than the censoring-ignoring fit,
  # correcting its downward bias (the model extrapolates the censored tail)
  expect_true(r$mean.aft > r$mean.naive)
  # censored latents sit at or above their observed censoring time; events
  # keep their observed log event time exactly
  cens <- r$status == 0
  expect_true(all(r$lat[cens] >= r$obs[cens] - 1e-8))
  expect_equal(r$lat[!cens], r$obs[!cens])
}

# ---- setResponse under censoring redraws the latents, keeps status ----

set.seed(3L)
cens.time <- f + 0.3 + sigma.true * rnorm(n)
status <- as.numeric(log.t <= cens.time)
obs.log.t <- ifelse(status == 1, log.t, cens.time)
bc.mut <- aftSampler(obs.log.t, status)
invisible(bc.mut$run(100L, 1L))
# shift every log-time up by 1; the fit should move up with it
obs.log.t <- bc.mut$data@y
bc.mut$setResponse(obs.log.t + 1)
res.mut <- bc.mut$run(20L, 20L)
expect_equal(dim(res.mut$train), c(n, 20L))
lat.mut <- bc.mut$getLatents()
expect_true(all(lat.mut[status == 0] >= obs.log.t[status == 0] + 1 - 1e-8))
expect_equal(lat.mut[status == 1], obs.log.t[status == 1] + 1)

# ---- refusals: weights and post-creation setData on an AFT sampler ----

expect_error(
  aftSampler(log.t, rep(1, n), weights = runif(n) + 0.5),
  "weight"
)
expect_error(
  bc.a$setData(samp.g$data),
  "aft"
)

# ---- public surface: Surv / two-column ingestion, predict, and the
# ---- survivalProbabilities generic ----

set.seed(8L)
cens <- f + 0.4 + sigma.true * rnorm(n)
status.s <- as.numeric(log.t <= cens)
time.s <- exp(ifelse(status.s == 1, log.t, cens)) # observed time (not logged)

# two-column (time, status) matrix with family = "aft"
fit.2col <- bart(
  x,
  cbind(time.s, status.s),
  family = "aft",
  n.trees = 50L,
  n.burn = 100L,
  n.samples = 200L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 7L,
  keepTrees = TRUE
)
expect_identical(fit.2col[["family"]], "aft")
# aft carries sigma and returns the linear predictor E[log T | x] (no
# probability transform), so it tracks the signal on the log scale
expect_false(is.null(fit.2col[["sigma"]]))
expect_true(cor(fitted(fit.2col), f) > 0.8)

# a Surv-like object (built without importing survival) auto-dispatches to
# aft and gives an identical fit
surv <- structure(
  cbind(time = time.s, status = status.s),
  class = "Surv",
  type = "right"
)
fit.surv <- bart(
  x,
  surv,
  n.trees = 50L,
  n.burn = 100L,
  n.samples = 200L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 7L,
  keepTrees = TRUE
)
expect_identical(fit.surv[["family"]], "aft")
expect_equal(fitted(fit.2col), fitted(fit.surv))

# predict returns log-scale linear-predictor draws; median time is exp of it
x.new <- matrix(runif(5L * p), 5L, p)
pr <- predict(fit.2col, x.new)
expect_equal(ncol(pr), 5L)

# survivalProbabilities returns DRAWS (draws x times x observations); the
# posterior-mean curve is monotone decreasing in [0, 1]
times <- c(0.5, 1, 2, 4)
sp <- survivalProbabilities(fit.2col, times, newdata = x.new)
n.draws <- nrow(fit.2col$yhat.train)
expect_equal(dim(sp), c(n.draws, length(times), 5L))
expect_true(all(sp >= 0 & sp <= 1))
sp.mean <- apply(sp, c(2L, 3L), mean)
expect_true(all(apply(sp.mean, 2L, function(curve) all(diff(curve) <= 1e-8))))
# every individual draw's curve is monotone too (each is an exact normal tail)
expect_true(all(sp[, -1L, ] <= sp[, -length(times), ] + 1e-8))

# refusals through the public surface
expect_error(dbarts(x, log.t, family = "aft"), "two-column|Surv")
expect_error(
  bart(x, cbind(c(-1, time.s[-1]), status.s), family = "aft", n.chains = 1L),
  "positive"
)
# the training-fit path (no newdata) spans the training observations
sp.train <- survivalProbabilities(fit.2col, times)
expect_equal(dim(sp.train), c(n.draws, length(times), n))

# an explicitly conflicting family with a Surv response errors instead of
# silently becoming aft
expect_error(
  bart(x, surv, family = "gaussian", n.chains = 1L, verbose = FALSE),
  "aft"
)
expect_error(
  dbarts(x, surv, family = "probit"),
  "aft"
)

# a factor status: survival::Surv codes it as multi-state ("mright"), which
# is detected with a hint; a data.frame factor status hints the same way
surv.mright <- structure(
  cbind(time = time.s, status = status.s + 1),
  class = "Surv",
  type = "mright"
)
expect_error(dbarts(x, surv.mright, family = "aft"), "factor")
expect_error(
  dbarts(
    x,
    data.frame(time = time.s, status = factor(status.s)),
    family = "aft"
  ),
  "factor"
)

# the formula interface takes a Surv left-hand side directly (dec-B97): a
# plain response never wrapped in Surv() is refused for lack of one, not for
# the interface itself
surv.df <- data.frame(
  t = time.s,
  s = status.s,
  x1 = x[, 1L],
  x2 = x[, 2L],
  x3 = x[, 3L]
)
expect_error(dbarts(t ~ x1, surv.df, family = "aft"), "survival::Surv")
# a Surv-like response (built without importing survival, above) through the
# formula path with family = "auto" auto-dispatches to aft, exactly as the
# direct-response form does (Decision 2)
fit.auto.formula <- dbarts(y ~ x1, data = list(y = surv, x1 = x[, 1L]))
expect_identical(fit.auto.formula$model@family, "aft")

# as a survreg user would type them, with the real survival package: formula
# vs matrix bitwise-identical (both "auto" and explicit family = "aft"),
# subset via formula bitwise-identical to hand-subsetting
if (requireNamespace("survival", quietly = TRUE)) {
  surv.df$surv <- survival::Surv(time.s, status.s)

  fitFormula <- function(formula, data, family, subset = NULL) {
    args <- list(
      formula,
      data,
      n.trees = 50L,
      n.burn = 100L,
      n.samples = 200L,
      n.chains = 1L,
      verbose = FALSE,
      seed = 7L,
      keepTrees = TRUE
    )
    if (!missing(family)) {
      args$family <- family
    }
    if (!is.null(subset)) {
      args$subset <- subset
    }
    do.call(bart, args)
  }

  # explicit family = "aft"
  fit.formula.aft <- fitFormula(surv ~ x1 + x2 + x3, surv.df, family = "aft")
  expect_identical(fit.formula.aft[["family"]], "aft")
  expect_identical(unname(fit.formula.aft$yhat.train), fit.2col$yhat.train)

  # family = "auto" dispatches identically
  fit.formula.autoFit <- fitFormula(surv ~ x1 + x2 + x3, surv.df)
  expect_identical(fit.formula.autoFit[["family"]], "aft")
  expect_identical(
    unname(fit.formula.autoFit$yhat.train),
    fit.2col$yhat.train
  )

  # 'subset' via the formula honours the same rows a hand-subsetted matrix
  # fit would use, bitwise
  sub <- seq.int(1L, n, by = 2L)
  fit.formula.sub <- fitFormula(
    surv ~ x1 + x2 + x3,
    surv.df,
    family = "aft",
    sub
  )
  fit.matrix.sub <- bart(
    x[sub, , drop = FALSE],
    survival::Surv(time.s[sub], status.s[sub]),
    family = "aft",
    n.trees = 50L,
    n.burn = 100L,
    n.samples = 200L,
    n.chains = 1L,
    verbose = FALSE,
    seed = 7L,
    keepTrees = TRUE
  )
  expect_identical(
    unname(fit.formula.sub$yhat.train),
    fit.matrix.sub$yhat.train
  )

  # the matrix interface's own 'subset' argument (dbarts()'s aft block,
  # which subsets the status vector alongside dbartsData()'s own x/y
  # subsetting) matches the same hand-subsetted fit bitwise
  fit.matrix.subsetArg <- bart(
    x,
    survival::Surv(time.s, status.s),
    family = "aft",
    subset = sub,
    n.trees = 50L,
    n.burn = 100L,
    n.samples = 200L,
    n.chains = 1L,
    verbose = FALSE,
    seed = 7L,
    keepTrees = TRUE
  )
  expect_identical(fit.matrix.subsetArg$yhat.train, fit.matrix.sub$yhat.train)

  # a PRE-BUILT dbartsData object carrying the same Surv attributes
  # (dbartsData(Surv(...) ~ ., data), called directly) is a legitimate
  # explicit family = "aft" request too, not the unsupported-interface
  # case the matrix-interface guards otherwise refuse - identical to the
  # same object's own family = "auto" dispatch
  dataObj <- dbartsData(surv ~ x1 + x2 + x3, surv.df)
  fit.dataObj.auto <- do.call(
    bart,
    c(
      list(dataObj),
      list(
        n.trees = 50L,
        n.burn = 100L,
        n.samples = 200L,
        n.chains = 1L,
        verbose = FALSE,
        seed = 7L,
        keepTrees = TRUE
      )
    )
  )
  expect_identical(fit.dataObj.auto[["family"]], "aft")
  fit.dataObj.explicit <- do.call(
    bart,
    c(
      list(dataObj, family = "aft"),
      list(
        n.trees = 50L,
        n.burn = 100L,
        n.samples = 200L,
        n.chains = 1L,
        verbose = FALSE,
        seed = 7L,
        keepTrees = TRUE
      )
    )
  )
  expect_identical(fit.dataObj.explicit[["family"]], "aft")
  expect_identical(fit.dataObj.explicit$yhat.train, fit.dataObj.auto$yhat.train)
  # the same route with NO Surv attributes still refuses an explicit "aft"
  dataObjPlain <- dbartsData(x1 ~ x2, surv.df)
  expect_error(
    bart(dataObjPlain, family = "aft", verbose = FALSE),
    "matrix interface"
  )

  # a plain formula with family = "aft" and no Surv response
  expect_error(
    dbarts(survival::Surv(t, s) ~ x1, surv.df, family = "gaussian"),
    "aft"
  )
}

# non-aft fits are refused by the bart method
fit.gauss <- bart(
  x,
  log.t,
  n.trees = 25L,
  n.burn = 50L,
  n.samples = 50L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 7L
)
expect_error(survivalProbabilities(fit.gauss, times = 1), "aft")

# multi-chain conventions: combineChains collapses the chain margin
fit.chains <- bart(
  x,
  cbind(time.s, status.s),
  family = "aft",
  n.trees = 25L,
  n.burn = 50L,
  n.samples = 50L,
  n.chains = 2L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 7L
)
sp.comb <- survivalProbabilities(fit.chains, times)
expect_equal(dim(sp.comb), c(100L, length(times), n))
sp.unc <- survivalProbabilities(fit.chains, times, combineChains = FALSE)
expect_equal(dim(sp.unc), c(2L, 50L, length(times), n))
# the combined result is the uncombined one with chains stacked sample-major
expect_equal(sp.comb[1:50, , ], sp.unc[1L, , , ])
expect_equal(sp.comb[51:100, , ], sp.unc[2L, , , ])

# ground truth: a fit packaged uncombined (same seed, so identical draws)
# carries yhat.train and sigma with explicit chain margins; the probability
# at any (chain, sample, time, obs) is the exact normal upper tail, which
# pins the sigma-to-draw alignment
fit.chains2 <- bart(
  x,
  cbind(time.s, status.s),
  family = "aft",
  n.trees = 25L,
  n.burn = 50L,
  n.samples = 50L,
  n.chains = 2L,
  n.threads = 1L,
  combineChains = FALSE,
  verbose = FALSE,
  seed = 7L
)
sp.unc2 <- survivalProbabilities(fit.chains2, times, combineChains = FALSE)
expect_equal(sp.unc, sp.unc2)
expect_equal(
  sp.unc2[2L, 17L, 3L, 5L],
  pnorm(
    (log(times[3L]) - fit.chains2$yhat.train[2L, 17L, 5L]) /
      fit.chains2$sigma[2L, 17L],
    lower.tail = FALSE
  )
)

# ---- the censoring status is settable between draws ----

set.seed(5L)
cens.set <- f + 0.25 + sigma.true * rnorm(n)
status.set <- as.numeric(log.t <= cens.set)
obs.set <- ifelse(status.set == 1, log.t, cens.set)
all.events <- rep(1, n)
expect_true(sum(status.set == 0) > 10L)

# creation parity, bitwise: nothing advances a generator when the target
# status leaves no censored row, so a sampler created all-events and one
# created censored and then set to all events at the same response draw the
# same chain
r5.control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  seed = 271L
)
r5.aft <- function(status) {
  dbarts(x, cbind(exp(obs.set), status), family = "aft", control = r5.control)
}
s.created <- r5.aft(all.events)
s.set <- r5.aft(status.set)
s.set$setResponse(s.set$data@y, status = all.events)
# the mirror the re-creation path reads, written only after the engine accepts
expect_equal(attr(s.set$control, "bartcore.survival"), all.events)
run.created <- s.created$run(50L, 50L)
run.set <- s.set$run(50L, 50L)
expect_identical(run.created$train, run.set$train)
expect_identical(run.created$sigma, run.set$sigma)

# y and status in one call: the bounds follow the NEW response
s.joint <- r5.aft(status.set)
status.joint <- as.numeric(seq_len(n) %% 3L != 0L)
s.joint$setResponse(s.joint$data@y + 0.75, status = status.joint)
lat.joint <- s.joint$getLatents()
y.joint <- s.joint$data@y
expect_equal(lat.joint[status.joint == 1], y.joint[status.joint == 1])
expect_true(all(lat.joint[status.joint == 0] >= y.joint[status.joint == 0]))

# ---- refusals: the status is validated before anything installs ----

s.gauss <- dbarts(x, log.t, control = r5.control)
expect_error(s.gauss$setResponse(log.t, status = all.events), "non-aft")
expect_error(
  s.set$setResponse(s.set$data@y, status = all.events[-1L]),
  "length"
)
bad.value <- status.set
bad.value[3L] <- 2
expect_error(s.set$setResponse(s.set$data@y, status = bad.value), "0 .*1")
bad.na <- status.set
bad.na[3L] <- NA_real_
expect_error(s.set$setResponse(s.set$data@y, status = bad.na), "cannot be NA")
# a refusal installs nothing: the response and the mirror are the ones in force
expect_equal(attr(s.set$control, "bartcore.survival"), all.events)
expect_equal(s.set$getLatents(), s.set$data@y)
# the R5 method coerces, so an integer status reaches the engine as doubles;
# the bridge still refuses a non-real vector from a raw caller, which it would
# otherwise read as doubles
expect_silent(s.set$setResponse(s.set$data@y, status = rep(1L, n)))
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setResponse,
    s.set$getPointer(),
    s.set$data@y,
    FALSE,
    rep(1L, n)
  ),
  "numeric"
)

# ---- the state handshake: a status change after the state was stored ----

set.seed(13L)
cens.hs <- f + 0.1 + sigma.true * rnorm(n)
status.hs <- as.numeric(log.t <= cens.hs)
obs.hs <- ifelse(status.hs == 1, log.t, cens.hs)
status.hs2 <- as.numeric(log.t <= cens.hs + 0.5)
freed <- status.hs == 0 & status.hs2 == 1
newly <- status.hs == 1 & status.hs2 == 0
expect_true(any(freed))
status.hs2[which(status.hs == 1)[1L:5L]] <- 0
newly <- status.hs == 1 & status.hs2 == 0
expect_true(any(newly))

hs.control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  seed = 909L
)
s.hs <- dbarts(
  x,
  cbind(exp(obs.hs), status.hs),
  family = "aft",
  control = hs.control
)
invisible(s.hs$run(50L, 10L))
expect_inherits(s.hs$state, "bartcoreState")
# the mutators store nothing without an explicit TRUE, so the saved state is
# the one the OLD censoring structure shaped
s.hs$setResponse(s.hs$data@y, status = status.hs2)
hs.file <- tempfile(fileext = ".rds")
saveRDS(s.hs, hs.file)
s.reloaded <- readRDS(hs.file)
unlink(hs.file)
# getPointer re-creates the engine from the CURRENT status and installs the
# stored state into it
lat.hs <- s.reloaded$getLatents()
y.hs <- s.hs$data@y
# a row censored when the state was stored and an event after it keeps its
# observed log event time: the state has no business overwriting data
expect_equal(lat.hs[status.hs2 == 1], y.hs[status.hs2 == 1])
# and a row that is an event in the state and censored here comes back
# REDRAWN above its bound rather than sitting on it
expect_true(all(lat.hs[newly] > y.hs[newly]))
expect_true(all(lat.hs[status.hs2 == 0] >= y.hs[status.hs2 == 0]))
expect_silent(invisible(s.reloaded$run(0L, 1L)))

# a survival formula writes its response with Surv(); cbind() on the left is
# refused under an explicit survival family, naming Surv()
d.cbind <- data.frame(x1 = x[, 1L], time = exp(log.t), status = rep(1, n))
for (family in c("aft", "hazard")) {
  expect_error(
    bart(cbind(time, status) ~ x1, data = d.cbind, family = family),
    pattern = paste0(
      "family = \"",
      family,
      "\" takes a formula response written survival::Surv\\(time, status\\)"
    )
  )
}
rm(d.cbind, family)

# ---- an offset fit's curves carry its offset: the training rows read the
# stored channel, which carries it, and newdata the fit's offset re-evaluated
# there or the one given for it
set.seed(41L)
o.aft <- rep(c(-1, 1), n / 2L)
status.aft <- rbinom(n, 1L, 0.7)
fit.aftOff <- bart(
  x,
  cbind(exp(log.t + o.aft), status.aft),
  family = "aft",
  offset = o.aft,
  n.trees = 10L,
  n.burn = 10L,
  n.samples = 20L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
sp.aftOff <- survivalProbabilities(fit.aftOff, c(0.5, 1))
expect_equal(
  survivalProbabilities(fit.aftOff, c(0.5, 1), newdata = x),
  sp.aftOff,
  tolerance = 1e-12
)
expect_equal(
  unname(survivalProbabilities(
    fit.aftOff,
    c(0.5, 1),
    newdata = x[1:3, ],
    offset = o.aft[1:3]
  )),
  unname(sp.aftOff[,, 1:3]),
  tolerance = 1e-12
)
expect_error(
  survivalProbabilities(fit.aftOff, 1, newdata = x[1:3, ]),
  "give survivalProbabilities an 'offset' for them",
  fixed = TRUE
)
# a non-numeric offset is refused by name, not coerced: FALSE given where
# combineChains was before offset took its place is not an offset of 0
expect_error(
  survivalProbabilities(fit.aftOff, 1, x[1:3, ], FALSE),
  "'offset' must be numeric",
  fixed = TRUE
)
expect_error(
  predict(fit.aftOff, x[1:3, ], offset = TRUE),
  "'offset' must be numeric",
  fixed = TRUE
)
rm(o.aft, status.aft, fit.aftOff, sp.aftOff)

# ---- a missing time or status is a missing response, as in survreg: the
# na.action drops the row (the fit equals one on the complete rows) or fails
set.seed(43L)
time.na <- exp(log.t)
status.na <- rbinom(n, 1L, 0.7)
time.na[3L] <- NA
status.na[5L] <- NA
complete.na <- !is.na(time.na) & !is.na(status.na)
naArgs <- list(
  family = "aft",
  n.trees = 5L,
  n.burn = 5L,
  n.samples = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 5L
)
fit.naOmit <- do.call(
  bart,
  c(list(x, cbind(time.na, status.na), na.action = na.omit), naArgs)
)
fit.complete <- do.call(
  bart,
  c(list(x[complete.na, ], cbind(time.na, status.na)[complete.na, ]), naArgs)
)
expect_identical(
  unname(fit.naOmit$yhat.train),
  unname(fit.complete$yhat.train)
)
expect_identical(as.vector(unclass(fit.naOmit$na.action)), c(3L, 5L))
expect_error(
  do.call(
    bart,
    c(list(x, cbind(time.na, status.na), na.action = na.fail), naArgs)
  ),
  "missing values in object",
  fixed = TRUE
)
d.na <- data.frame(x, time = time.na, status = status.na)
fit.naFormula <- do.call(
  bart,
  c(
    list(survival::Surv(time, status) ~ ., data = d.na, na.action = na.omit),
    naArgs
  )
)
expect_identical(
  unname(fit.naFormula$yhat.train),
  unname(fit.complete$yhat.train)
)
rm(
  time.na,
  status.na,
  complete.na,
  naArgs,
  fit.naOmit,
  fit.complete,
  d.na,
  fit.naFormula
)

# ---- the matrix interface reads 'subset' once ----

# The log times, the censoring status, the predictors and the offset are
# those of one reading of 'subset', whatever its expression gives the next
# time it is read. An `id` predictor names the row each kept row is.
set.seed(31L)
onceX <- cbind(x, id = as.double(seq_len(n)))
onceTime <- exp(log.t)
onceStatus <- as.double(rbinom(n, 1L, 0.6))
onceWeights <- runif(n, 0.5, 2)
onceOffset <- rnorm(n, sd = 0.1)
subsetDraws <- list()
drawSubset <- function() {
  rows <- sample(n, 60L)
  subsetDraws[[length(subsetDraws) + 1L]] <<- rows
  rows
}
onceDoors <- list(
  dbarts = function() {
    dbarts(
      onceX,
      cbind(onceTime, onceStatus),
      family = "aft",
      subset = drawSubset(),
      offset = onceOffset,
      control = control
    )
  },
  bart = function() {
    bart(
      onceX,
      cbind(onceTime, onceStatus),
      family = "aft",
      subset = drawSubset(),
      offset = onceOffset,
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 5L,
      n.samples = 2L,
      n.burn = 0L,
      verbose = FALSE,
      samplerOnly = TRUE
    )
  },
  # a row dropped for a missing predictor leaves the status of the rows kept
  na.omit = function() {
    withMissing <- onceX
    withMissing[seq(3L, n, by = 7L), 1L] <- NA_real_
    dbarts(
      withMissing,
      cbind(onceTime, onceStatus),
      family = "aft",
      subset = drawSubset(),
      offset = onceOffset,
      na.action = na.omit,
      control = control
    )
  }
)
for (door in names(onceDoors)) {
  subsetDraws <- list()
  sampler <- onceDoors[[door]]()
  expect_identical(length(subsetDraws), 1L, info = door)
  rows <- as.integer(sampler$data@x[, "id"])
  expect_identical(
    rows,
    if (door == "na.omit") {
      setdiff(subsetDraws[[1L]], seq(3L, n, by = 7L))
    } else {
      subsetDraws[[1L]]
    },
    info = door
  )
  expect_identical(
    attr(sampler$control, "bartcore.survival"),
    onceStatus[rows],
    info = door
  )
  expect_identical(sampler$data@y, log(onceTime)[rows], info = door)
  expect_identical(sampler$data@offset, onceOffset[rows], info = door)
}
# a gaussian fit and a hazard fit read it once as well
subsetDraws <- list()
plainOnce <- dbarts(
  onceX,
  log(onceTime),
  subset = drawSubset(),
  weights = onceWeights,
  offset = onceOffset,
  control = control
)
expect_identical(length(subsetDraws), 1L)
rows <- as.integer(plainOnce$data@x[, "id"])
expect_identical(rows, subsetDraws[[1L]])
expect_identical(plainOnce$data@y, log(onceTime)[rows])
expect_identical(plainOnce$data@weights, onceWeights[rows])
expect_identical(plainOnce$data@offset, onceOffset[rows])
subsetDraws <- list()
hazardOnce <- dbarts(
  onceX,
  cbind(ceiling(onceTime), onceStatus),
  family = "hazard",
  subset = drawSubset(),
  control = control
)
expect_identical(length(subsetDraws), 1L)
expect_true(all(hazardOnce$data@x[, "id"] %in% subsetDraws[[1L]]))
rm(onceX, onceTime, onceStatus, onceWeights, onceOffset, subsetDraws)
rm(drawSubset, onceDoors, door, sampler, rows, plainOnce, hazardOnce)

# ---- a re-derived response range is that of the observed times ----

# setOffset with updateScale = TRUE read the range from the log times in
# force, which hold each chain's own draw at a censored row: two chains left
# the call under two leaf priors, neither the observed times', and the model
# recorded the first chain's. The range is the minimum and maximum of the
# observed log time less the offset, a censored row's censoring time among
# them, whatever the chains have drawn.
set.seed(19L)
rangeStatus <- rbinom(n, 1L, 0.5)
rangeOffsets <- list(
  "the offset in force" = 0.6 * (x[, 2L] - 0.5),
  "a new offset" = rep(c(-0.4, 0.2), n / 2L)
)
rangeSampler <- function() {
  sampler <- dbarts(
    x,
    cbind(exp(log.t), rangeStatus),
    offset = rangeOffsets[[1L]],
    family = "aft",
    control = dbartsControl(
      n.chains = 2L,
      n.threads = 1L,
      n.trees = 10L,
      updateState = FALSE,
      seed = 19L
    )
  )
  invisible(sampler$run(0L, 20L))
  sampler
}
# each chain's transform as the engine holds it, one row per chain
chainTransforms <- function(sampler) {
  .Call(
    dbarts:::C_dbarts_bartcore_getLeafPrior,
    sampler$getPointer(),
    0L
  )[, c("response.shift", "response.scale"), drop = FALSE]
}
expectObservedTransform <- function(sampler, offset, info) {
  bounds <- range(sampler$data@y - offset)
  transforms <- chainTransforms(sampler)
  expect_identical(transforms[1L, ], transforms[2L, ], info = info)
  expect_equal(
    unname(transforms[1L, ]),
    c(bounds[1L] + diff(bounds) / 2, diff(bounds)),
    tolerance = 1e-14,
    info = info
  )
  bounds
}
for (case in names(rangeOffsets)) {
  offset <- rangeOffsets[[case]]
  sampler <- rangeSampler()
  latents <- matrix(sampler$getLatents(), n)
  censored <- rangeStatus == 0L
  # each chain holds its own drawn times, above the observed ones
  expect_true(all(latents[censored, ] > sampler$data@y[censored]), info = case)
  expect_false(
    identical(range(latents[, 1L] - offset), range(latents[, 2L] - offset)),
    info = case
  )

  sampler$setOffset(offset, updateScale = TRUE)
  bounds <- expectObservedTransform(sampler, offset, case)
  expect_identical(
    as.vector(attr(sampler$model, "response.range")),
    bounds,
    info = case
  )
  prior <- sampler$getLeafPrior()
  expect_equal(
    c(prior$response.shift, prior$response.scale),
    c(bounds[1L] + diff(bounds) / 2, diff(bounds)),
    tolerance = 1e-14,
    info = case
  )
  # the call draws nothing and moves no drawn time
  expect_identical(matrix(sampler$getLatents(), n), latents, info = case)
  duplicate <- sampler$copy()
  expect_identical(
    chainTransforms(duplicate),
    chainTransforms(sampler),
    info = case
  )
  expect_identical(
    attr(duplicate$model, "response.range"),
    attr(sampler$model, "response.range"),
    info = case
  )
  expect_true(all(is.finite(sampler$run(0L, 5L)$train)), info = case)
  expect_true(all(is.finite(duplicate$run(0L, 5L)$train)), info = case)
}

# and through the flat entry, dbarts_sampler_setOffset
source(
  system.file("common", "capiConsumer.R", package = "dbarts"),
  local = TRUE
)
consumer <- compileCapiConsumer("aft", "the C API consumer")
if (is.null(consumer$skip)) {
  sampler <- rangeSampler()
  offset <- rangeOffsets[[2L]]
  expect_equal(
    consumer$CALL("capi_set_offset", sampler$getPointer(), offset, TRUE),
    1L
  )
  expectObservedTransform(sampler, offset, "the flat entry")
  expect_true(all(is.finite(sampler$run(0L, 5L)$train)))
}
rm(rangeStatus, rangeOffsets, rangeSampler, chainTransforms, consumer)
rm(expectObservedTransform, case, offset, sampler, latents, censored, bounds)
rm(prior, duplicate)
