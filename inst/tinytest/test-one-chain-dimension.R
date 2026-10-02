# extract, predict and survivalProbabilities keep a chain dimension of
# length 1 on a one-chain fit under combineChains = FALSE (dec-A79), for
# every type and fit class; the fit's own stored fields are unchanged. Each
# uncombined shape is checked against the combineChains = TRUE one it drops
# to (identical values, one leading margin), rather than against a
# hardcoded shape alone, so a value regression fails the same test as a
# shape one.

set.seed(101)
n <- 30L
x <- matrix(rnorm(n * 2L), n, 2L)
y <- rnorm(n)
newX <- matrix(rnorm(6L * 2L), 6L, 2L)

quick <- function(y, ...) {
  suppressWarnings(bart(
    x,
    y,
    ...,
    n.samples = 8L,
    n.burn = 4L,
    n.trees = 5L,
    n.chains = 1L,
    keepTrees = TRUE,
    verbose = FALSE
  ))
}
seeded <- function(expr) {
  set.seed(5)
  expr
}
# checks that 'uncombined' is 'combined' with a leading length-1 chain
# margin, values and all
expectKeptChain <- function(uncombined, combined, info = "") {
  # a combined scalar type (sigma, k, shape) is a plain vector, with no
  # dim of its own to prepend the chain margin to
  combinedShape <- if (is.null(dim(combined))) {
    length(combined)
  } else {
    dim(combined)
  }
  expect_equal(dim(uncombined), c(1L, combinedShape), info = info)
  d <- length(dim(uncombined)) - 1L
  indices <- c(list(1L), rep(list(quote(expr = )), d))
  expect_identical(
    do.call(`[`, c(list(uncombined), indices, list(drop = TRUE))),
    combined,
    info = info
  )
}

# --- bart ---

fit <- suppressWarnings(bart(
  x,
  y,
  k = dbarts::dbartsPriors$chi(1.5, 2),
  n.samples = 8L,
  n.burn = 4L,
  n.trees = 5L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
for (type in c("ev", "bart")) {
  expectKeptChain(
    extract(fit, type, combineChains = FALSE),
    extract(fit, type),
    info = type
  )
}
expectKeptChain(
  seeded(extract(fit, "ppd", combineChains = FALSE)),
  seeded(extract(fit, "ppd")),
  info = "ppd"
)
expectKeptChain(
  extract(fit, "sigma", combineChains = FALSE),
  extract(fit, "sigma"),
  info = "sigma"
)
expectKeptChain(
  extract(fit, "k", combineChains = FALSE),
  extract(fit, "k"),
  info = "k"
)
expectKeptChain(
  extract(fit, "varcount", combineChains = FALSE),
  extract(fit, "varcount"),
  info = "varcount"
)
expectKeptChain(
  extract(fit, "loglik", combineChains = FALSE),
  extract(fit, "loglik"),
  info = "loglik"
)
expectKeptChain(
  predict(fit, newX, combineChains = FALSE),
  predict(fit, newX),
  info = "predict ev"
)
expectKeptChain(
  seeded(predict(fit, newX, "ppd", combineChains = FALSE)),
  seeded(predict(fit, newX, "ppd")),
  info = "predict ppd"
)

# fitted, residuals and summary never carry a chain margin and are
# unaffected by a one-chain fit's kept chain axis elsewhere
expect_null(dim(fitted(fit)))
expect_null(dim(residuals(fit)))
expect_true(is.numeric(summary(fit)$stats[["mean"]]))

# ci.level pools over every draw AND chain already, so it carries no chain
# margin to keep, at any chain count
band <- predict(fit, newX, ci.level = 0.9, combineChains = FALSE)
expect_equal(dim(band), c(6L, 3L))

# --- the amplitude-coupled forest arm ---

z <- rep(c(0, 1), length.out = n)
dfZ <- data.frame(y = y + z, a = x[, 1L], b = x[, 2L], z = z)
fitZ <- suppressWarnings(bart(
  y ~ a + b + z:forest(a + b),
  dfZ,
  n.samples = 8L,
  n.burn = 4L,
  n.trees = 5L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
expectKeptChain(
  extract(fitZ, "forest", combineChains = FALSE),
  extract(fitZ, "forest"),
  info = "forest extract"
)
newZ <- data.frame(a = rnorm(6L), b = rnorm(6L), z = rep(c(0, 1), 3L))
expectKeptChain(
  predict(fitZ, newZ, "forest", combineChains = FALSE),
  predict(fitZ, newZ, "forest"),
  info = "forest predict"
)
expectKeptChain(
  predict(fitZ, newZ, "bart", combineChains = FALSE),
  predict(fitZ, newZ, "bart"),
  info = "blend predict"
)

# --- multinomial: the trailing category margin rides alongside the kept
# chain one ---

category <- factor(rep(c("u", "v", "w"), length.out = n))
fitM <- quick(category, family = "multinomial")
expectKeptChain(
  extract(fitM, "ev", combineChains = FALSE),
  extract(fitM, "ev"),
  info = "multinomial ev"
)
expectKeptChain(
  extract(fitM, "varcount", combineChains = FALSE),
  extract(fitM, "varcount"),
  info = "multinomial varcount"
)
expectKeptChain(
  predict(fitM, newX, combineChains = FALSE),
  predict(fitM, newX),
  info = "multinomial predict"
)

# --- ordinal: the K - 1 threshold margin ---

fitO <- quick(factor(category, ordered = TRUE), family = "ordinal")
expectKeptChain(
  extract(fitO, "thresholds", combineChains = FALSE),
  extract(fitO, "thresholds"),
  info = "ordinal thresholds"
)
expectKeptChain(
  predict(fitO, newX, combineChains = FALSE),
  predict(fitO, newX),
  info = "ordinal predict"
)

# --- negbin: shape is a scalar-per-draw field, no observation margin ---

counts <- rpois(n, 3)
fitN <- quick(counts, family = "nbinom")
expectKeptChain(
  extract(fitN, "shape", combineChains = FALSE),
  extract(fitN, "shape"),
  info = "negbin shape"
)
expectKeptChain(
  predict(fitN, newX, combineChains = FALSE),
  predict(fitN, newX),
  info = "negbin predict"
)

# --- hurdle: sigma, k and varcount are the components' own scalar fields
# (k and varcount are lists, one entry per component), and loglik round-
# trips through hurdleParts' own recursive extract() calls ---

yPos <- ifelse(runif(n) < 0.3, 0, rlnorm(n))
fitH <- quick(yPos, family = "hurdle.lognormal")
expectKeptChain(
  extract(fitH, "sigma", combineChains = FALSE),
  extract(fitH, "sigma"),
  info = "hurdle sigma"
)
kUncombined <- extract(fitH, "k", combineChains = FALSE)
kCombined <- extract(fitH, "k")
expect_identical(names(kUncombined), names(kCombined))
for (part in names(kCombined)) {
  expectKeptChain(
    kUncombined[[part]],
    kCombined[[part]],
    info = paste("k", part)
  )
}
vcUncombined <- extract(fitH, "varcount", combineChains = FALSE)
vcCombined <- extract(fitH, "varcount")
for (part in names(vcCombined)) {
  expectKeptChain(
    vcUncombined[[part]],
    vcCombined[[part]],
    info = paste("varcount", part)
  )
}
expectKeptChain(
  extract(fitH, "loglik", combineChains = FALSE),
  extract(fitH, "loglik"),
  info = "hurdle loglik"
)
expectKeptChain(
  predict(fitH, newX, combineChains = FALSE),
  predict(fitH, newX),
  info = "hurdle predict"
)

# --- survivalProbabilities keeps the same kept chain margin ---

fitAft <- quick(cbind(abs(y) + 0.1, rbinom(n, 1L, 0.7)), family = "aft")
expectKeptChain(
  survivalProbabilities(fitAft, times = c(0.5, 1), combineChains = FALSE),
  survivalProbabilities(fitAft, times = c(0.5, 1)),
  info = "survivalProbabilities aft"
)
spNewdata <- survivalProbabilities(
  fitAft,
  times = c(0.5, 1),
  newdata = newX,
  combineChains = FALSE
)
expect_identical(dim(spNewdata), c(1L, 8L, 2L, 6L))

if (requireNamespace("survival", quietly = TRUE)) {
  status <- rbinom(n, 1L, 0.7)
  fitHaz <- quick(
    survival::Surv(sample(1:3, n, replace = TRUE), status),
    family = "hazard"
  )
  expectKeptChain(
    survivalProbabilities(fitHaz, combineChains = FALSE),
    survivalProbabilities(fitHaz),
    info = "survivalProbabilities hazard"
  )
  newHaz <- matrix(rnorm(4L * 2L), 4L, 2L)
  shNewdata <- survivalProbabilities(
    fitHaz,
    newdata = newHaz,
    combineChains = FALSE
  )
  expect_identical(dim(shNewdata)[1L:2L], c(1L, 8L))
  rm(status, fitHaz, newHaz, shNewdata)
}

rm(
  n,
  x,
  y,
  newX,
  quick,
  seeded,
  expectKeptChain,
  fit,
  type,
  band,
  z,
  dfZ,
  fitZ,
  newZ,
  category,
  fitM,
  fitO,
  counts,
  fitN,
  yPos,
  fitH,
  kUncombined,
  kCombined,
  part,
  vcUncombined,
  vcCombined,
  fitAft,
  spNewdata
)
