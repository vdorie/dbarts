# Heteroscedastic BART (docs/design/heteroscedastic.md): a `variance` selector
# adds a second forest modeling s^2(x). The fit recovers f(x) and a plausible
# s(x); a homoscedastic truth does not manufacture spurious heteroscedasticity;
# the surface is gaussian + constant-leaf only, and predict returns s(x).

set.seed(9, sample.kind = "Rejection")

# ---- a step-heteroscedastic truth: the fit recovers f(x) and s(x) ----
n <- 800L
x <- runif(n)
fTrue <- 2 * x
sTrue <- ifelse(x < 0.5, 0.3, 1.5)
y <- fTrue + sTrue * rnorm(n)

fit <- bart2(
  x,
  y,
  variance = TRUE,
  n.trees.variance = 25L,
  n.trees = 50L,
  n.samples = 400L,
  n.burn = 400L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)

# s(x) is reported, one posterior slab per training row
expect_true(!is.null(fit$s.train))
expect_equal(dim(fit$s.train), c(400L, n))

sHat <- apply(fit$s.train, 2L, mean)
fHat <- fit$yhat.train.mean

# the mean surface tracks the ramp
expect_true(cor(fHat, fTrue) > 0.9)

# s(x) is larger where the truth is noisier, and tracks the two levels within a
# factor of ~2 (the estimate attenuates with finite trees)
sLow <- mean(sHat[x < 0.5])
sHigh <- mean(sHat[x >= 0.5])
expect_true(sHigh > 2 * sLow)
expect_true(sLow > 0.15 && sLow < 0.6)
expect_true(sHigh > 0.9 && sHigh < 2.2)

# ---- predict returns s(x) alongside f(x) on new data ----
xNew <- matrix(c(0.25, 0.75), 2L, 1L)
pred <- predict(fit, xNew)
s <- attr(pred, "s")
expect_true(!is.null(s))
expect_equal(dim(s), c(400L, 2L))
sNew <- apply(s, 2L, mean)
expect_true(sNew[2L] > 2 * sNew[1L]) # x = 0.75 noisier than x = 0.25

# ---- a homoscedastic truth does not manufacture heteroscedasticity ----
set.seed(11, sample.kind = "Rejection")
xHom <- runif(n)
yHom <- 2 * xHom + 0.8 * rnorm(n)
fitHom <- bart2(
  xHom,
  yHom,
  variance = TRUE,
  n.trees.variance = 25L,
  n.trees = 50L,
  n.samples = 400L,
  n.burn = 400L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
sHom <- apply(fitHom$s.train, 2L, mean)
# no region should read as wildly more variable than another: the spread of the
# per-observation s(x) around its mean stays modest (no spurious structure)
expect_true(sd(sHom) / mean(sHom) < 0.35)
# and the recovered level is near the truth (0.8)
expect_true(mean(sHom) > 0.55 && mean(sHom) < 1.15)

# ---- a homoscedastic fit (no variance forest) carries no s channel ----
fitPlain <- bart2(
  x,
  y,
  n.trees = 50L,
  n.samples = 100L,
  n.burn = 100L,
  n.chains = 1L,
  verbose = FALSE
)
expect_null(fitPlain$s.train)

# ---- the variance forest is gaussian only ----
expect_error(
  bart2(
    x,
    as.integer(y > median(y)),
    variance = TRUE,
    family = "probit",
    n.samples = 10L,
    n.burn = 10L,
    n.chains = 1L,
    verbose = FALSE
  ),
  "variance forest requires"
)
