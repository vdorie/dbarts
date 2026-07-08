# The custom-tau-prior path drives every kept sample through a per-sweep R
# callback fired inside one sampler run. It reproduces, draw for draw, the
# run(0, 1)-per-sample loop it replaced (verified bitwise against that loop
# before it was deleted); these pinned values are that loop's output, so a
# change that shifts the interleaving of R and sampler draws breaks them.
# n.thin = 2 exercises the thinning block logic (offset held, last sweep read).

oldKind <- RNGkind()
if (getRversion() >= "3.6.0") {
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
} else {
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion"))
}

n <- 300L
n.g <- 5L
p <- 4L

set.seed(22L)
x <- matrix(rnorm(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
f <- 5 * sin(pi * x[, 1L] * x[, 2L]) + 10 * (x[, 3L] - 0.5)^2 + 4 * x[, 4L]
g <- factor(sample.int(n.g, n, replace = TRUE))
y <- f + rnorm(n.g, 0, 1)[as.integer(g)] + rnorm(n, 0, 1)
d <- data.frame(x, y = y, g = g)

# same body as rbart.priors$cauchy, under a name that misses the builtin
# lookup so rbart_vi routes through the callback path rather than in-core
customCauchy <- function(x, rel.scale) dcauchy(x, 0, rel.scale * 2.5, TRUE)

runFit <- function() {
  set.seed(22L)
  dbarts::rbart_vi(
    y ~ x1 + x2 + x3 + x4,
    data = d,
    group.by = g,
    prior = customCauchy,
    n.trees = 30L,
    n.chains = 1L,
    n.threads = 1L,
    n.burn = 20L,
    n.samples = 25L,
    n.thin = 2L,
    keepTrees = FALSE,
    verbose = FALSE
  )
}

fit <- runFit()

# pinned draws from the deleted loop, recorded on one platform; SIMD/BLAS
# reduction order shifts the last bits elsewhere, hence the 1e-8 tolerance
tau5 <- c(
  0.21628278043302535,
  0.20047547120022122,
  0.32913812354179472,
  0.3577386605460961,
  0.91264300337492865
)
ranefRow1 <- c(
  0.015987791670784988,
  0.10784937712389432,
  0.12426941670573695,
  -0.0095328924285038429,
  -0.11437761177328708
)
yhatMeans <- c(
  11.016945998901727,
  5.1840338407404216,
  3.2384066799854492,
  19.67567655250717,
  0.16391455705261368
)

expect_equal(unname(fit$tau[1:5]), tau5, tolerance = 1e-8)
expect_equal(unname(fit$ranef[1L, 1:5]), ranefRow1, tolerance = 1e-8)
expect_equal(
  unname(colMeans(fit$yhat.train)[1:5]),
  yhatMeans,
  tolerance = 1e-8
)

# a single chain draws through R's generator, so the path reproduces exactly
fit2 <- runFit()
expect_identical(fit$tau, fit2$tau)
expect_identical(fit$ranef, fit2$ranef)
expect_identical(fit$yhat.train, fit2$yhat.train)

suppressWarnings(RNGkind(oldKind[1L], oldKind[2L], oldKind[3L]))
rm(
  x,
  f,
  g,
  y,
  d,
  customCauchy,
  runFit,
  fit,
  fit2,
  tau5,
  ranefRow1,
  yhatMeans,
  n,
  n.g,
  p,
  oldKind
)
