# A custom tau prior runs rbart_vi's Gibbs blocks in an R loop; a built-in
# prior runs them in-core (the oracle). With the custom prior set verbatim
# equal to the built-in cauchy, both draw from the same posterior, so their
# tau posteriors must agree. The response here is deliberately skewed (normal
# predictors give the Friedman f a long right tail), the case that once dumped
# the mean-vs-midpoint offset into the intercepts and inflated tau ~20-40x.

# pin the generator so suite order (an earlier file may leave sample.kind
# "Rounding") cannot shift the draws this statistical assertion depends on
oldKind <- RNGkind()
if (getRversion() >= "3.6.0") {
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
} else {
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion"))
}

n <- 800L
n.g <- 8L
p <- 5L

set.seed(11L)
x <- matrix(rnorm(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
f <- 10 * sin(pi * x[, 1L] * x[, 2L]) + 20 * (x[, 3L] - 0.5)^2 +
  10 * x[, 4L] + 5 * x[, 5L]
g <- factor(sample.int(n.g, n, replace = TRUE))
y <- f + rnorm(n.g, 0, 1)[as.integer(g)] + rnorm(n, 0, 1)
d <- data.frame(x, y = y, g = g)

# same body as rbart.priors$cauchy, under a name that misses the builtin
# lookup so rbart_vi routes through the R loop
customCauchy <- function(x, rel.scale) dcauchy(x, 0, rel.scale * 2.5, TRUE)

rbartArgs <- list(
  y ~ x1 + x2 + x3 + x4 + x5, data = d, group.by = g,
  n.trees = 50L, n.chains = 1L, n.threads = 1L,
  n.burn = 100L, n.samples = 100L, n.thin = 1L,
  keepTrees = FALSE, verbose = FALSE
)

set.seed(11L)
fit.loop <- do.call(
  dbarts::rbart_vi, c(rbartArgs, list(prior = customCauchy))
)
set.seed(11L)
fit.core <- do.call(dbarts::rbart_vi, rbartArgs)

tau.loop <- mean(fit.loop$tau)
tau.core <- mean(fit.core$tau)

# the intercepts stay near zero: tau tracks the ~1 truth, nowhere near the
# response scale the bug settled at
expect_true(tau.loop < sd(y) / 5)
# and the two posteriors agree up to Monte Carlo error (independent streams);
# the bug drove this ratio to ~20-40
expect_true(tau.loop / tau.core > 0.4 && tau.loop / tau.core < 2.5)

suppressWarnings(RNGkind(oldKind[1L], oldKind[2L], oldKind[3L]))
rm(x, f, g, y, d, customCauchy, rbartArgs, fit.loop, fit.core,
   tau.loop, tau.core, n, n.g, p, oldKind)
