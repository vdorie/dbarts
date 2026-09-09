# A variance forest under family = "aft": log-normal survival with a
# covariate-dependent dispersion, log T = f(x) + s(x) eps. The channels that
# report a residual scale read s(x) rather than the pinned sigma - the
# log-likelihood, the survival curves - each censored latent is redrawn at its
# own row's scale, and the composed fit carries the heteroscedastic refusals
# beside aft's own.

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)
source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

set.seed(24, sample.kind = "Rejection")

n <- 400L
p <- 3L
x <- matrix(
  runif(n * p),
  n,
  p,
  dimnames = list(NULL, c("x1", "x2", "x3"))
)
f <- 1.5 * x[, 2L] - 0.8 * x[, 3L]
s.true <- ifelse(x[, 1L] < 0.5, 0.25, 1.0)
log.t <- f + s.true * rnorm(n)
log.c <- f + 0.6 + 0.9 * rnorm(n)
status <- as.numeric(log.t <= log.c)
time <- exp(ifelse(status == 1, log.t, log.c))
x.test <- matrix(
  runif(20L * p),
  20L,
  p,
  dimnames = list(NULL, colnames(x))
)

# a real censoring rate, so the truncated redraw is actually exercised
expect_true(mean(status == 0) > 0.15 && mean(status == 0) < 0.6)

warnings.fit <- captureWarnings(
  fit <- bart(
    x,
    cbind(time, status),
    test = x.test,
    family = "aft",
    variance = varianceForest(n.trees = 20L),
    n.trees = 50L,
    n.burn = 300L,
    n.samples = 300L,
    n.chains = 1L,
    keepTrees = TRUE,
    verbose = FALSE,
    seed = 12L
  )
)
expect_equal(length(warnings.fit), 0L)

n.draws <- 300L
expect_identical(fit[["family"]], "aft")
expect_equal(dim(fit$s.train), c(n.draws, n))
expect_equal(dim(fit$s.test), c(n.draws, 20L))

# the surface separates the two levels of the truth
s.hat <- apply(fit$s.train, 2L, mean)
s.low <- mean(s.hat[x[, 1L] < 0.5])
s.high <- mean(s.hat[x[, 1L] >= 0.5])
expect_true(s.high > 2 * s.low)
expect_true(s.low > 0.1 && s.low < 0.6)
expect_true(s.high > 0.6 && s.high < 1.6)

# the mean surface still tracks the signal
expect_true(cor(fit$yhat.train.mean, f) > 0.8)

# sigma is the pinned constant carrying no posterior content
expect_equal(length(unique(fit$sigma)), 1L)

# ---- the log-likelihood scores at s(x_i), events and censored rows alike ----
ev <- extract(fit, type = "bart", sample = "train")
warnings.loglik <- captureWarnings(loglik <- extract(fit, type = "loglik"))
expect_equal(length(warnings.loglik), 0L)

y.rep <- rep(fit$y, each = n.draws)
loc <- as.vector(ev)
sd.rep <- as.vector(fit$s.train)
expected <- dnorm(y.rep, loc, sd.rep, log = TRUE)
censored <- rep(status, each = n.draws) == 0
expected[censored] <- pnorm(
  y.rep[censored],
  loc[censored],
  sd.rep[censored],
  lower.tail = FALSE,
  log.p = TRUE
)
expect_equal(loglik, matrix(expected, n.draws, n), tolerance = 1e-12)

# and is nowhere near the pinned scalar's answer, so the check has teeth
at.scalar <- dnorm(y.rep, loc, fit$sigma[1L], log = TRUE)
at.scalar[censored] <- pnorm(
  y.rep[censored],
  loc[censored],
  fit$sigma[1L],
  lower.tail = FALSE,
  log.p = TRUE
)
expect_true(abs(sum(loglik) - sum(at.scalar)) > 100)

# ---- survival curves divide by the surface, on training rows ----
times <- c(0.5, 1, 2)
sp <- survivalProbabilities(fit, times)
expect_equal(dim(sp), c(n.draws, length(times), n))
expect_equal(
  sp[, 2L, ],
  matrix(
    pnorm((log(times[2L]) - loc) / sd.rep, lower.tail = FALSE),
    n.draws,
    n
  ),
  tolerance = 1e-12
)
# monotone decreasing in t, on every draw
expect_true(all(sp[, -1L, ] <= sp[, -length(times), ] + 1e-8))

# ---- and at newdata, off the replayed surface predict parks on its result ----
pred <- predict(fit, x.test, type = "bart")
s.new <- attr(pred, "s")
expect_equal(dim(s.new), c(n.draws, 20L))
sp.new <- survivalProbabilities(fit, times, newdata = x.test)
expect_equal(dim(sp.new), c(n.draws, length(times), 20L))
expect_equal(
  sp.new[, 1L, ],
  matrix(
    pnorm(
      (log(times[1L]) - as.vector(pred)) / as.vector(s.new),
      lower.tail = FALSE
    ),
    n.draws,
    20L
  ),
  tolerance = 1e-12
)

# ---- without saved trees there is no surface at newdata, and it is named ----
fit.no.trees <- bart(
  x,
  cbind(time, status),
  family = "aft",
  variance = varianceForest(n.trees = 8L),
  n.trees = 20L,
  n.burn = 20L,
  n.samples = 20L,
  n.chains = 1L,
  keepSampler = TRUE,
  verbose = FALSE,
  seed = 12L
)
expect_error(
  survivalProbabilities(fit.no.trees, times, newdata = x.test),
  "replays no variance surface"
)
# the TRAINING rows are unaffected: s.train is a run channel, not a replay
expect_equal(
  dim(survivalProbabilities(fit.no.trees, times)),
  c(20L, length(times), n)
)

# ---- the composed fit carries both families' refusals ----
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  updateState = FALSE
)
sampler <- dbarts(
  x,
  cbind(time, status),
  family = "aft",
  variance = varianceForest(n.trees = 8L),
  control = control
)

# the variance forest owns the residual scale, so sigma is not settable - the
# heteroscedastic refusal, which names the forest before the family
expect_error(sampler$setSigma(2), "variance forest owns the residual scale")

# the scale leaf is calibrated once against the response transform fixed at
# creation, so a re-anchoring response or offset swap is refused
expect_error(
  sampler$setResponse(sampler$data@y, updateScale = TRUE),
  "updateScale = FALSE"
)
expect_error(
  sampler$setOffset(rep(0.1, n), updateScale = TRUE),
  "updateScale = FALSE"
)
# and the pinned swap is taken
sampler$setResponse(sampler$data@y, updateScale = FALSE)
sampler$setOffset(rep(0.1, n), updateScale = FALSE)

# aft's own refusals stand: the declined user weight channel is what frees the
# internal one, and the censoring structure is fixed at creation
expect_error(sampler$setWeights(rep(1, n)), "weight")
expect_error(sampler$setData(sampler$data), "aft")

# ---- the other latent families still refuse a variance forest ----
expect_error(
  bart(
    x,
    as.integer(log.t > median(log.t)),
    family = "probit",
    variance = TRUE,
    n.trees = 10L,
    n.samples = 10L,
    n.burn = 10L,
    n.chains = 1L,
    verbose = FALSE
  ),
  "a variance forest requires family = \"gaussian\" or \"aft\""
)

# ---- reduction: an uncensored heteroscedastic aft IS the gaussian fit ----
# every hook delegates and refreshLatents returns before drawing, so the two
# samplers must walk the same RNG stream draw for draw
control.seeded <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  updateState = FALSE,
  seed = 271L
)
samp.g <- dbarts(x, log.t, control = control.seeded, variance = TRUE)
res.g <- bartcoreRun(dbarts:::bartcoreSampler(samp.g), 100L, 100L)

samp.a <- dbarts(x, log.t, control = control.seeded, variance = TRUE)
ctrl <- samp.a$control
attr(ctrl, "bartcore.survival") <- rep(1, n) # every observation an event
samp.a$control <- ctrl
res.a <- bartcoreRun(
  dbarts:::bartcoreSampler(samp.a, family = "aft"),
  100L,
  100L
)

expect_identical(res.g$train, res.a$train)
expect_identical(res.g$variance, res.a$variance)
# the surface is not a constant, so the equality is over a forest that moved
expect_true(length(unique(res.a$variance[, 1L])) > 1L)

# ---- the exported latent replay stays scalar, and says so ----
expect_error(
  dbartsDrawLatents("aft", rep(0, 5L), rep(1, 5L), sigma = rep(0.5, 5L)),
  "single residual scale for every row"
)
