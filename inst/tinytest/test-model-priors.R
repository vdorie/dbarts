# prior specifications as first-class objects: the dbartsPriors constructors
# and the scoped bare-name vocabulary inside the prior arguments

set.seed(42)
n <- 100L
df <- data.frame(a = runif(n), b = runif(n))
df$y <- df$a + rnorm(n, 0, 0.5)

# constructors build validated objects
prior.cgm <- dbartsPriors$cgm(power = 1.5, base = 0.9)
expect_inherits(prior.cgm, "dbartsCGMPrior")
expect_equal(prior.cgm@power, 1.5)
expect_error(dbartsPriors$cgm(power = -1), pattern = "power")
expect_error(dbartsPriors$chisq(df = -3), pattern = "df")
expect_error(dbartsPriors$normal(-2), pattern = "positive")

prior.normal <- dbartsPriors$normal(dbartsPriors$chi(1.1, 2))
expect_inherits(prior.normal@k, "dbartsChiHyperprior")
expect_equal(prior.normal@k@degreesOfFreedom, 1.1)

# string forms remain accepted for bart2 compatibility
expect_inherits(dbartsPriors$normal("chi(1.5, 3)")@k, "dbartsChiHyperprior")
expect_equal(dbartsPriors$normal("2.5")@k, 2.5)

# objects pass to the fitting functions and match the bare-name sugar draw
# for draw
control <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 25L,
                         updateState = FALSE)
set.seed(1)
sampler.obj <- dbarts(y ~ a + b, df, tree.prior = prior.cgm,
                      node.prior = dbartsPriors$normal(3), control = control)
samples.obj <- sampler.obj$run(20L, 20L)
set.seed(1)
sampler.sugar <- dbarts(y ~ a + b, df, tree.prior = cgm(1.5, 0.9),
                        node.prior = normal(3), control = control)
samples.sugar <- sampler.sugar$run(20L, 20L)
expect_identical(samples.obj$train, samples.sugar$train)

# the vocabulary shadows the caller's environment: a masking 'normal' (as
# rstanarm exports) does not change the meaning
normal <- function(...) stop("masked by another package")
cgm <- function(...) stop("masked by another package")
sampler.masked <- dbarts(y ~ a + b, df, tree.prior = cgm(1.5, 0.9),
                         node.prior = normal(3), control = control)
expect_inherits(sampler.masked$model@node.hyperprior, "dbartsFixedHyperprior")
expect_equal(sampler.masked$model@node.hyperprior@k, 3)
rm(normal, cgm)

# num.vars remains available inside the arguments
sampler.nv <- dbarts(y ~ a + b, df, tree.prior = cgm(2, 0.95, 1 / num.vars),
                     control = control)
expect_equal(length(sampler.nv$model@tree.prior@splitProbabilities), 0L)

# named split probabilities resolve against the data at fit time
prior.split <- dbartsPriors$cgm(split.probs = c(a = 3, .default = 1))
sampler.split <- dbarts(y ~ a + b, df, tree.prior = prior.split, control = control)
expect_equal(sampler.split$model@tree.prior@splitProbabilities,
             c(a = 0.75, b = 0.25))
expect_null(sampler.split$model@tree.prior@splitProbabilitiesSpec)

# but cannot silently enter a model without data to resolve against
expect_error(methods::new("dbartsModel", tree.prior = prior.split),
             pattern = "resolved against data")

# binary responses keep their default k hyperprior
df$y.binary <- rbinom(n, 1L, 0.5)
sampler.bin <- dbarts(y.binary ~ a + b, df, control = control)
expect_inherits(sampler.bin$model@node.hyperprior, "dbartsChiHyperprior")

# DART: a Dirichlet prior over the split-variable probabilities, bartcore
# engine only
prior.dart <- dbartsPriors$dart(a = 0.75, alpha = 2, update.delay = 10)
expect_inherits(prior.dart, "dbartsDartPrior")
expect_inherits(prior.dart, "dbartsCGMPrior")
expect_equal(prior.dart@a, 0.75)
expect_error(dbartsPriors$dart(a = -1), pattern = "'a' must be positive")
expect_error(dbartsPriors$dart(alpha = 0), pattern = "'alpha' must be positive")

expect_error(dbarts(y ~ a + b, df, tree.prior = dart(), control = control),
             pattern = "DART tree prior requires")

# a left-unset update delay resolves to half the control burn-in
control.bc <- dbartsControl(engine = "bartcore", n.chains = 1L,
                            n.threads = 1L, n.trees = 75L, n.burn = 500L,
                            updateState = TRUE)
set.seed(7)
n.dart <- 250L
x.dart <- matrix(runif(n.dart * 10L), n.dart)
y.dart <- 6 * sin(pi * x.dart[, 1L]) + 4 * (x.dart[, 2L] - 0.5)^2 +
  rnorm(n.dart, 0, 0.3)
sampler.dart <- dbarts(y.dart ~ x.dart, tree.prior = dart(),
                       control = control.bc)
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
expect_null(dbarts(y.dart ~ x.dart, control = control.bc)$run(0L, 5L)$varprobs)

# the Dirichlet machinery is fixed at creation
expect_error(sampler.dart$setModel(sampler.dart$model),
             pattern = "cannot change a DART tree prior")

# bart2 exposes DART through the dart flag and packages varprobs
fit.dart <- bart2(y.dart ~ x.dart, engine = "bartcore", dart = TRUE,
                  n.samples = 25L, n.burn = 50L, n.trees = 25L,
                  n.chains = 2L, n.threads = 1L, verbose = FALSE)
expect_equal(dim(fit.dart$varprobs), c(2L, 25L, 10L))
expect_equal(unname(apply(fit.dart$varprobs, c(1L, 2L), sum)),
             matrix(1, 2L, 25L))
expect_error(bart2(y.dart ~ x.dart, engine = "bartcore", dart = TRUE,
                   split.probs = rep(0.1, 10L)),
             pattern = "cannot be combined")
expect_error(bart2(y.dart ~ x.dart, engine = "bartcore", dart = 2),
             pattern = "must be TRUE, FALSE")
expect_null(bart2(y.dart ~ x.dart, engine = "bartcore", n.samples = 5L,
                  n.burn = 5L, n.trees = 5L, n.chains = 1L, n.threads = 1L,
                  verbose = FALSE)$varprobs)

# a full spec object overrides power/base with its own settings
fit.spec <- bart2(y.dart ~ x.dart, engine = "bartcore",
                  dart = dbartsPriors$dart(a = 0.75, update.delay = 5),
                  n.samples = 5L, n.burn = 10L, n.trees = 10L,
                  n.chains = 1L, n.threads = 1L, verbose = FALSE,
                  keepSampler = TRUE)
expect_equal(fit.spec$fit$model@tree.prior@a, 0.75)
expect_equal(fit.spec$fit$model@tree.prior@update.delay, 5)
