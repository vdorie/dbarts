# prior specifications as first-class objects: the dbartsPriors constructors
# and the scoped bare-name vocabulary inside the prior arguments, and the tree
# prior's plain power, base and split.probs

set.seed(42)
n <- 100L
df <- data.frame(a = runif(n), b = runif(n))
df$y <- df$a + rnorm(n, 0, 0.5)

# constructors build validated objects
prior.cgm <- dbartsPriors$cgm(power = 1.5, base = 0.9)
expect_inherits(prior.cgm, "dbartsCGMPrior")
expect_equal(prior.cgm@power, 1.5)
# NULL is absent, so it is the default
expect_equal(dbartsPriors$cgm(power = NULL)@power, 2)
expect_equal(dbartsPriors$cgm(base = NULL)@base, 0.95)
# dart() is split.probs's value and takes no depth: its first formal is 'a'
expect_error(dbartsPriors$dart(power = 3), pattern = "unused argument")
expect_equal(dbartsPriors$dart(0.3)@a, 0.3)
expect_error(dbartsPriors$cgm(power = -1), pattern = "power")
expect_error(dbartsPriors$chisq(df = -3), pattern = "df")
expect_error(dbartsPriors$normal(-2), pattern = "positive")

prior.normal <- dbartsPriors$normal(dbartsPriors$chi(1.1, 2))
expect_inherits(prior.normal@k, "dbartsChiHyperprior")
expect_equal(prior.normal@k@degreesOfFreedom, 1.1)

# string forms remain accepted for bart compatibility
expect_inherits(dbartsPriors$normal("chi(1.5, 3)")@k, "dbartsChiHyperprior")
expect_equal(dbartsPriors$normal("2.5")@k, 2.5)

# objects pass to the fitting functions and match the bare-name sugar draw
# for draw
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE
)
set.seed(1)
sampler.obj <- dbarts(
  y ~ a + b,
  df,
  power = 1.5,
  base = 0.9,
  leaf.prior = dbartsPriors$normal(3),
  control = control,
  n.trees = 25L
)
samples.obj <- sampler.obj$run(20L, 20L)
set.seed(1)
sampler.sugar <- dbarts(
  y ~ a + b,
  df,
  power = 1.5,
  base = 0.9,
  leaf.prior = normal(3),
  control = control,
  n.trees = 25L
)
samples.sugar <- sampler.sugar$run(20L, 20L)
expect_identical(samples.obj$train, samples.sugar$train)

# the vocabulary shadows the caller's environment: a masking 'normal' (as
# rstanarm exports) does not change the meaning, nor a masking 'dart' inside
# split.probs
normal <- function(...) stop("masked by another package")
dart <- function(...) stop("masked by another package")
sampler.masked <- dbarts(
  y ~ a + b,
  df,
  power = 1.5,
  base = 0.9,
  split.probs = dart(a = 0.7),
  leaf.prior = normal(3),
  control = control,
  n.trees = 25L
)
expect_inherits(sampler.masked$model@leaf.hyperprior, "dbartsFixedHyperprior")
expect_equal(sampler.masked$model@leaf.hyperprior@k, 3)
expect_equal(sampler.masked$model@tree.prior@a, 0.7)
rm(normal, dart)

# num.vars remains available inside split.probs
sampler.nv <- dbarts(
  y ~ a + b,
  df,
  split.probs = 1 / num.vars,
  control = control,
  n.trees = 25L
)
expect_equal(length(sampler.nv$model@tree.prior@splitProbabilities), 0L)

# the vocabulary is dart() plus num.vars/numvars, and nothing else: a
# caller's own 'data' or 'control' inside split.probs is the caller's
data <- c(a = 3, .default = 1)
sampler.caller <- dbarts(
  y ~ a + b,
  df,
  split.probs = data,
  control = control,
  n.trees = 25L
)
expect_equal(
  unname(sampler.caller$model@tree.prior@splitProbabilities),
  c(0.75, 0.25)
)
rm(data)

# named split probabilities resolve against the data at fit time
prior.split <- dbartsPriors$cgm(split.probs = c(a = 3, .default = 1))
sampler.split <- dbarts(
  y ~ a + b,
  df,
  split.probs = c(a = 3, .default = 1),
  control = control,
  n.trees = 25L
)
expect_equal(
  sampler.split$model@tree.prior@splitProbabilities,
  c(a = 0.75, b = 0.25)
)
expect_null(sampler.split$model@tree.prior@splitProbabilitiesSpec)

# but cannot silently enter a model without data to resolve against
expect_error(
  methods::new("dbartsModel", tree.prior = prior.split),
  pattern = "resolved against data"
)

# binary responses keep their default k hyperprior
df$y.binary <- rbinom(n, 1L, 0.5)
sampler.bin <- dbarts(y.binary ~ a + b, df, control = control, n.trees = 25L)
expect_inherits(sampler.bin$model@leaf.hyperprior, "dbartsChiHyperprior")

# DART: a Dirichlet prior over the split-variable probabilities
prior.dart <- dbartsPriors$dart(a = 0.75, alpha = 2, update.delay = 10)
expect_inherits(prior.dart, "dbartsDartPrior")
expect_inherits(prior.dart, "dbartsCGMPrior")
expect_equal(prior.dart@a, 0.75)
expect_error(
  dbartsPriors$dart(a = -1),
  "'a' must be a single positive finite number",
  fixed = TRUE
)
expect_error(
  dbartsPriors$dart(alpha = 0),
  "'alpha' must be a single positive finite number",
  fixed = TRUE
)

# a left-unset update delay resolves to half the control burn-in
control.bc <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.burn = 500L,
  updateState = TRUE
)
set.seed(7)
n.dart <- 250L
x.dart <- matrix(runif(n.dart * 10L), n.dart)
y.dart <- 6 *
  sin(pi * x.dart[, 1L]) +
  4 * (x.dart[, 2L] - 0.5)^2 +
  rnorm(n.dart, 0, 0.3)
sampler.dart <- dbarts(
  y.dart ~ x.dart,
  split.probs = dart(),
  control = control.bc,
  n.trees = 75L
)
expect_equal(sampler.dart$model@tree.prior@update.delay, 250)

# split probabilities concentrate on the signal variables
samples.dart <- sampler.dart$run(500L, 500L)
probs <- sampler.dart$state[[1L]][["dart.probabilities"]]
expect_equal(length(probs), 10L)
expect_true(sum(probs[1:2]) > 0.9)

# each kept sample records the probabilities; non-DART runs return none
expect_equal(dim(samples.dart$varprobs), c(10L, 500L))
expect_true(all(samples.dart$varprobs >= 0 & samples.dart$varprobs <= 1))
expect_equal(unname(colSums(samples.dart$varprobs)), rep(1, 500L))
expect_true(mean(colSums(samples.dart$varprobs[1:2, ])) > 0.9)
expect_null(
  dbarts(y.dart ~ x.dart, control = control.bc, n.trees = 75L)$run(
    0L,
    5L
  )$varprobs
)

# the Dirichlet machinery is fixed at creation: the sampler's own model is
# taken, and one with the ordinary tree prior is refused
expect_silent(sampler.dart$setModel(sampler.dart$model))
model.cgm <- sampler.dart$model
model.cgm@tree.prior <- methods::new("dbartsCGMPrior")
expect_error(
  sampler.dart$setModel(model.cgm),
  pattern = "cannot change the tree prior's kind"
)

# bart exposes DART through split.probs = dart() and packages varprobs
# combineChains = FALSE pinned deliberately: the shape assertions below
# expect the raw uncombined (n.chains x n.samples x n.vars) array
fit.dart <- bart(
  y.dart ~ x.dart,
  split.probs = dart(),
  n.samples = 25L,
  n.burn = 50L,
  n.trees = 25L,
  n.chains = 2L,
  n.threads = 1L,
  verbose = FALSE,
  combineChains = FALSE
)
expect_equal(dim(fit.dart$varprobs), c(2L, 25L, 10L))
expect_equal(
  unname(apply(fit.dart$varprobs, c(1L, 2L), sum)),
  matrix(1, 2L, 25L)
)
# split.probs is bart()'s own argument and warns of nothing
expect_silent(
  bart(
    y.dart ~ x.dart,
    split.probs = rep(0.1, 10L),
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )
)
expect_null(
  bart(
    y.dart ~ x.dart,
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )$varprobs
)

# dart()'s settings reach the record beside the fitting function's depth
fit.spec <- bart(
  y.dart ~ x.dart,
  split.probs = dbartsPriors$dart(a = 0.75, update.delay = 5),
  n.samples = 5L,
  n.burn = 10L,
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepSampler = TRUE
)
expect_equal(fit.spec$fit$model@tree.prior@a, 0.75)
expect_equal(fit.spec$fit$model@tree.prior@update.delay, 5)

# each scalar of a prior constructor is one number in its range, refused by
# the argument's own name: NA of either type, two values, a string and the
# range's ends, each in one line, never R's slot-class or condition error
dfKind <- "positive finite number; to hold sigma at a value, use fixed()"
refusedScalars <- list(
  list(function(v) dbartsPriors$chisq(df = v), "df", dfKind, c(0, -1, Inf)),
  list(
    function(v) dbartsPriors$chisq(3, quant = v),
    "quant",
    "number between 0 and 1",
    c(0, 1, -0.5)
  ),
  list(
    dbartsPriors$fixed,
    "value",
    "positive finite number",
    c(0, -1, Inf)
  ),
  list(dbartsPriors$chi, "df", "positive number", c(0, -1)),
  list(
    function(v) dbartsPriors$chi(1.5, scale = v),
    "scale",
    "positive number",
    c(0, -1)
  ),
  list(dbartsPriors$cgm, "power", "positive number", c(0, -1)),
  list(
    function(v) dbartsPriors$cgm(base = v),
    "base",
    "number in (0, 1)",
    c(0, 1, 1.5)
  ),
  list(
    function(v) dbartsPriors$dart(a = v),
    "a",
    "positive finite number",
    c(0, Inf)
  ),
  list(
    function(v) dbartsPriors$dart(b = v),
    "b",
    "positive finite number",
    c(0, Inf)
  ),
  list(
    function(v) dbartsPriors$dart(alpha = v),
    "alpha",
    "positive finite number",
    c(0, Inf)
  )
)
for (entry in refusedScalars) {
  message <- paste0("'", entry[[2L]], "' must be a single ", entry[[3L]])
  for (value in c(
    list(NA, NA_real_, NaN, c(1, 2), "a"),
    as.list(entry[[4L]])
  )) {
    refusal <- tryCatch(entry[[1L]](value), error = conditionMessage)
    expect_identical(refusal, message, info = deparse(value))
  }
}
expect_identical(
  tryCatch(dbartsPriors$dart(rho = c(1, 2)), error = conditionMessage),
  "'rho' must be a single positive finite number"
)
expect_identical(
  tryCatch(dbartsPriors$dart(update.delay = -1), error = conditionMessage),
  "'update.delay' must be a single non-negative number"
)
expect_identical(
  tryCatch(dbartsPriors$dart(update.alpha = NA), error = conditionMessage),
  "'update.alpha' must be TRUE/FALSE"
)
# an NA rho keeps its own words, which name NULL as the default's spelling
expect_error(dbartsPriors$dart(rho = NA), "'rho' must not be NA")
# the range's own values are taken: a huge df stands in for a held sigma, an
# infinite chi() scale is the improper limit
expect_identical(dbartsPriors$chisq(1e6, 0.9)@df, 1e6)
expect_identical(dbartsPriors$chi(1.5, Inf)@scale, Inf)
expect_identical(dbartsPriors$cgm(power = 3L)@power, 3)
# a class built directly is held to the same checks, length first
expect_error(
  new("dbartsChiSqPrior", df = c(3, 4), quantile = 0.9),
  "'df' must be a single positive finite number",
  fixed = TRUE
)
expect_error(
  new("dbartsCGMPrior", power = c(2, 3), base = 0.95),
  "'power' must be a single positive number",
  fixed = TRUE
)
expect_error(
  new("dbartsChiHyperprior", degreesOfFreedom = c(1, 2), scale = 2),
  "'df' must be a single positive number",
  fixed = TRUE
)
expect_error(new("dbartsFixedHyperprior", k = c(1, 2)), "'k' must be positive")
# and to its range: an infinite df, a quantile at 1, an infinite held value
expect_error(
  new("dbartsChiSqPrior", df = Inf, quantile = 0.9),
  "'df' must be a single positive finite number; to hold sigma at a value",
  fixed = TRUE
)
expect_error(
  new("dbartsChiSqPrior", df = 3, quantile = 1),
  "'quant' must be a single number between 0 and 1",
  fixed = TRUE
)
expect_error(
  new("dbartsFixedPrior", value = Inf),
  "'value' must be a single positive finite number",
  fixed = TRUE
)
dartInfinite <- dbartsPriors$dart()
dartInfinite@alpha <- Inf
expect_error(
  validObject(dartInfinite),
  "'alpha' must be a single positive finite number",
  fixed = TRUE
)
# a matrix of split probabilities is refused in one line naming the argument
expect_identical(
  tryCatch(
    dbartsPriors$cgm(split.probs = matrix(1:2, 1L)),
    error = conditionMessage
  ),
  "'split.probs' must be a vector, not a matrix or an array"
)
expect_error(
  bart(
    df[c("a", "b")],
    df$y,
    split.probs = matrix(1:2, 1L),
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 5L,
    n.chains = 1L,
    verbose = FALSE
  ),
  "'split.probs' must be a vector, not a matrix or an array",
  fixed = TRUE
)
# the residual prior as every door reads it: an infinite df is refused at
# once, naming fixed(), on gaussian and aft alike, and so is the 0.9-34
# shorthand; a NA reaches the same words
infiniteDf <- paste0("'df' must be a single ", dfKind)
expect_error(
  dbarts(y ~ a + b, df, family = gaussian(sigma = chisq(Inf, 0.9))),
  infiniteDf,
  fixed = TRUE
)
expect_error(
  dbartsFamilies$aft(sigma = dbartsPriors$chisq(Inf, 0.9)),
  infiniteDf,
  fixed = TRUE
)
for (sigdf in list(Inf, NA)) {
  expect_error(
    dbarts:::holdingBartNotices(bart(
      df[c("a", "b")],
      df$y,
      sigdf = sigdf,
      n.samples = 5L,
      n.burn = 5L,
      n.trees = 5L,
      n.chains = 1L,
      verbose = FALSE
    )),
    infiniteDf,
    fixed = TRUE
  )
}
# and the bridge refuses a record the constructor never saw
heldSampler <- dbarts(
  df[c("a", "b")],
  df$y,
  control = dbartsControl(n.chains = 1L, n.threads = 1L, verbose = FALSE),
  n.trees = 5L
)
infiniteModel <- heldSampler$model
infiniteModel@resid.prior@df <- Inf
expect_error(
  new("dbartsSampler", heldSampler$control, infiniteModel, heldSampler$data),
  infiniteDf,
  fixed = TRUE
)
expect_error(heldSampler$setModel(infiniteModel), infiniteDf, fixed = TRUE)
rm(
  dfKind,
  refusedScalars,
  entry,
  message,
  value,
  refusal,
  infiniteDf,
  sigdf,
  heldSampler,
  infiniteModel,
  dartInfinite
)
