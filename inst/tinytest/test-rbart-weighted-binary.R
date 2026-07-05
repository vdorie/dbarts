# rbart_vi grouped fits with binary responses and with non-uniform weights.
# The binary + weights combination is a deliberate engine refusal (a weighted
# probit was fit incorrectly by the classic engine), so the two legal halves
# are tested separately and the combination is asserted to refuse cleanly.
# rbart_vi needs group.by and defaults keepTrees = TRUE; kept tiny throughout.

set.seed(7)
n <- 200L
n.g <- 5L
x <- matrix(runif(n * 3L), n, 3L)
g <- sample(n.g, n, replace = TRUE)
b <- rnorm(n.g, 0, 0.8)
f <- 1.2 * x[, 1L] - 0.6 + b[g]

# grouped binary (probit): finite variance/ranef draws and fitted probabilities
z <- rbinom(n, 1L, pnorm(f))
fit.bin <- rbart_vi(z ~ x, group.by = g,
                    n.samples = 8L, n.burn = 5L, n.thin = 1L,
                    n.chains = 2L, n.trees = 20L, n.threads = 1L,
                    verbose = FALSE)
expect_inherits(fit.bin, "rbart")
expect_true(all(is.finite(fit.bin$tau)))
expect_true(all(is.finite(fit.bin$ranef)))
prob.bin <- fitted(fit.bin)
expect_true(all(is.finite(prob.bin)))
expect_true(all(prob.bin >= 0 & prob.bin <= 1))    # probit fits are probabilities

# grouped continuous with non-uniform weights: finite tau/ranef/fitted
y <- f + rnorm(n, 0, 0.4)
w <- runif(n, 0.5, 2)
fit.wt <- rbart_vi(y ~ x, group.by = g, weights = w,
                   n.samples = 8L, n.burn = 5L, n.thin = 1L,
                   n.chains = 2L, n.trees = 20L, n.threads = 1L,
                   verbose = FALSE)
expect_inherits(fit.wt, "rbart")
expect_true(all(is.finite(fit.wt$tau)))
expect_true(all(is.finite(fit.wt$ranef)))
expect_true(all(is.finite(fitted(fit.wt))))

# binary + weights is refused: a weighted probit is not fit
expect_error(
  rbart_vi(z ~ x, group.by = g, weights = w,
           n.samples = 8L, n.burn = 5L, n.thin = 1L,
           n.chains = 2L, n.trees = 20L, n.threads = 1L, verbose = FALSE),
  pattern = "do not support weights")

rm(fit.bin, fit.wt, prob.bin, x, g, b, f, z, y, w, n, n.g)
