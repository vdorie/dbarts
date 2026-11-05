# A sweep draws its augmentation variables first: after the per-sweep hook,
# where a larger sampler makes its changes, and before the trees. So a setter
# applied inside the hook and the same setter applied between run(0, 1) calls
# give the same draws, bit for bit: the sweep draws its latents against the
# change and reads the family's working response and weights after that draw
# every sweep, which a mask installed in the hook moves. And the ordinal
# latents a run returns are drawn under the cutpoints that run records.

set.seed(20261105L)
n <- 60L
x <- matrix(runif(n * 2L), n, 2L, dimnames = list(NULL, c("x1", "x2")))
eta <- 1.5 * (x[, 1L] - 0.5)
y.binary <- as.double(eta + rnorm(n) > 0)
y.counts <- as.double(rpois(n, exp(0.5 + eta)))
y.continuous <- eta + rnorm(n, 0, 0.5)
numSweeps <- 8L
masks <- lapply(seq_len(numSweeps), function(i) as.double(runif(n) > 0.3))
offsets <- lapply(seq_len(numSweeps), function(i) rnorm(n, 0, 0.3))

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE,
  seed = 99L
)
families <- list(
  probit = list(y = y.binary, family = "probit"),
  logistic = list(y = y.binary, family = "logistic"),
  nbinom = list(y = y.counts, family = "nbinom"),
  student = list(y = y.continuous, family = dbarts:::student()),
  gaussian = list(y = y.continuous, family = "gaussian")
)
make <- function(arm) {
  dbarts::dbarts(
    x,
    arm$y,
    family = arm$family,
    n.trees = 5L,
    control = control
  )
}
between <- function(sampler, setter) {
  out <- matrix(NA_real_, n, numSweeps)
  for (i in seq_len(numSweeps)) {
    setter(sampler, i)
    out[, i] <- as.numeric(sampler$run(0L, 1L)$train)
  }
  out
}
inHook <- function(sampler, setter) {
  train <- numeric(n * numSweeps)
  .Call(
    dbarts:::C_dbarts_bartcore_runWithCallback,
    sampler$getPointer(),
    0L,
    numSweeps,
    list(train = train),
    function(i) {
      setter(sampler, i + 1L)
      FALSE
    },
    environment()
  )
  matrix(train, n, numSweeps)
}
setters <- list(
  setOffset = function(sampler, i) sampler$setOffset(offsets[[i]]),
  setActiveRows = function(sampler, i) sampler$setActiveRows(masks[[i]])
)
for (familyName in names(families)) {
  for (setterName in names(setters)) {
    arm <- families[[familyName]]
    setter <- setters[[setterName]]
    expect_identical(
      inHook(make(arm), setter),
      between(make(arm), setter),
      info = paste0(familyName, ", ", setterName, " in the hook")
    )
  }
}

# every z an ordinal run returns lies in its category's interval under the
# cutpoints that run recorded: y = k means gamma_(k-1) < z <= gamma_k
y.ordinal <- as.double(cut(eta + rnorm(n), c(-Inf, -0.4, 0.4, Inf)))
sampler <- dbarts::dbarts(
  x,
  y.ordinal,
  family = "ordinal",
  n.trees = 5L,
  control = control
)
inInterval <- logical(30L)
for (i in seq_along(inInterval)) {
  thresholds <- sampler$run(0L, 1L)$thresholds[, 1L]
  bounds <- c(-Inf, thresholds, Inf)
  z <- sampler$getLatents()
  inInterval[i] <- all(
    z > bounds[y.ordinal] & z <= bounds[y.ordinal + 1L]
  )
}
expect_true(all(inInterval))
