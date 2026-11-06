# A state holds values no chain could run from only when edited or damaged:
# setState refuses one with a sigma or k that is not positive, a response
# transform that is not finite, a generator whose index is out of range or
# whose bits are all zero, DART weights that weigh nothing where the forest
# may split, or ordinal cutpoints out of order, with the message every
# refusal gives and the sampler left as it was; a reload of one is refused
# the same way, not a crash; and a warm-start donor holding one is refused by
# its own message. The gamma memo the generator block carries is a cache and
# is recomputed on the way in.

set.seed(20261106L)
n <- 80L
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, c("x1", "x2", "x3")))
y <- x[, 1L] + 2 * x[, 2L] + rnorm(n, 0, 0.2)
df <- data.frame(x, y = y)
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  updateState = FALSE
)
make <- function(...) {
  dbarts(y ~ x1 + x2 + x3, df, n.trees = 5L, control = control, ...)
}
refused <- "state is not consistent with this sampler"

# two copies of one sampler continue it identically; the one handed an
# edited state refuses it and runs as the other does
sampler <- make(leaf.prior = normal(chi(1.5, 2)))
invisible(sampler$run(20L, 5L))
sampler$storeState()
own <- sampler$state
expect_true(is.finite(own[[1L]]$forests[[1L]]$k))
edits <- list(
  "sigma 0" = function(s) {
    s[[1L]]$sigma <- 0
    s
  },
  "sigma -1" = function(s) {
    s[[1L]]$sigma <- -1
    s
  },
  "sigma Inf" = function(s) {
    s[[1L]]$sigma <- Inf
    s
  },
  "k 0" = function(s) {
    s[[1L]]$forests[[1L]]$k <- 0
    s
  },
  "k Inf under a finite scale" = function(s) {
    s[[1L]]$forests[[1L]]$k <- Inf
    s
  },
  "fit.scale NA" = function(s) {
    s[[1L]]$fit.scale[1L] <- NA_real_
    s
  },
  "fit.scale Inf" = function(s) {
    s[[1L]]$fit.scale[2L] <- Inf
    s
  },
  "rng index -1" = function(s) {
    s[[1L]]$rng.state[1L] <- -1L
    s
  },
  "rng index -100000" = function(s) {
    s[[1L]]$rng.state[1L] <- -100000L
    s
  },
  "rng index NA" = function(s) {
    s[[1L]]$rng.state[1L] <- NA_integer_
    s
  },
  "rng index N + 1" = function(s) {
    s[[1L]]$rng.state[1L] <- 625L
    s
  },
  "rng state all zero" = function(s) {
    s[[1L]]$rng.state[2L:625L] <- 0L
    s
  }
)
for (edit in names(edits)) {
  edited <- edits[[edit]](own)
  a <- sampler$copy()
  b <- sampler$copy()
  expect_error(a$setState(edited), refused, info = edit)
  expect_error(a$setState(edited, forceUpdate = TRUE), refused, info = edit)
  expect_identical(a$run(0L, 2L), b$run(0L, 2L), info = edit)
}

# the generator's index at its ends, and the sampler's own sigma, install
for (index in c(0L, 624L)) {
  edited <- own
  edited[[1L]]$rng.state[1L] <- index
  expect_true(sampler$copy()$setState(edited), info = index)
}

# the gamma memo after the generator is recomputed from the shapes it
# records: a constant edited under its shape installs as it was stored, and
# the chain runs as an unedited one does
memo <- readBin(
  writeBin(own[[1L]]$rng.state[626L:643L], raw()),
  "double",
  n = 9L
)
expect_true(memo[1L] >= 1 && memo[4L] == memo[1L] - 0.5)
edited <- own
edited[[1L]]$rng.state[632L:633L] <- readBin(
  writeBin(NaN, raw()),
  "integer",
  n = 2L
)
a <- sampler$copy()
b <- sampler$copy()
expect_true(a$setState(edited))
a$storeState()
expect_identical(a$state[[1L]]$rng.state, own[[1L]]$rng.state)
run.a <- a$run(0L, 3L)
expect_identical(run.a, b$run(0L, 3L))
expect_true(all(is.finite(run.a$sigma)))

# a reload of an edited state is refused when the engine is re-created
reloaded <- sampler$copy()
reloaded$state <- edits[["sigma 0"]](own)
stateFile <- tempfile(fileext = ".rds")
saveRDS(reloaded, stateFile)
reloaded <- readRDS(stateFile)
unlink(stateFile)
expect_error(reloaded$run(0L, 1L), refused)

# a warm-start donor holding such a value is refused by name, in both forms
donorRefused <- paste(
  "warm-start donor holds a value this sampler cannot run from"
)
for (edit in c("sigma 0", "k Inf under a finite scale", "fit.scale Inf")) {
  recipient <- sampler$copy()
  donor <- edits[[edit]](own)
  expect_error(recipient$installTrees(donor), donorRefused, info = edit)
  expect_error(
    recipient$installTrees(donor, forceUpdate = TRUE),
    donorRefused,
    info = edit
  )
}

# k = Inf is what a draw under an infinite scale reaches, and installs there
flat <- make(leaf.prior = normal(chi(1.5, Inf)))
invisible(flat$run(5L, 1L))
flat$storeState()
infinite <- flat$state
infinite[[1L]]$forests[[1L]]$k <- Inf
expect_true(flat$copy()$setState(infinite))
expect_true(flat$copy()$installTrees(infinite))

# DART on two of three columns: a weight of 0 where the forest may split is
# refused, a column it may not split on holds anything
restricted <- make(forests = list(forest(x1 + x2)), split.probs = dart())
invisible(restricted$run(20L, 1L))
restricted$storeState()
dart.own <- restricted$state
expect_identical(dart.own[[1L]]$dart.probabilities[3L], 0)
for (weights in list(c(-5, 1, 0), c(NA, 1, 0), c(0, 0, 0), c(0, 1, 0))) {
  edited <- dart.own
  edited[[1L]]$dart.probabilities <- weights
  expect_error(
    restricted$copy()$setState(edited),
    refused,
    info = paste(weights, collapse = " ")
  )
}
for (alpha in c(0, -1, Inf)) {
  edited <- dart.own
  edited[[1L]]$dart.alpha <- alpha
  expect_error(restricted$copy()$setState(edited), refused, info = alpha)
}
edited <- dart.own
edited[[1L]]$dart.probabilities <- c(0.25, 0.25, 0.5)
expect_true(restricted$copy()$setState(edited))

# those weights, the forest's own, warm-start a forest on every column with
# its own weights started fresh
everywhere <- make(split.probs = dart())
invisible(everywhere$run(5L, 1L))
expect_true(everywhere$installTrees(dart.own))
expect_identical(everywhere$state[[1L]]$dart.probabilities, rep(1 / 3, 3L))
expect_true(all(is.finite(everywhere$run(0L, 2L)$sigma)))

# ordinal cutpoints: the first 0, each above the one before, all finite
y.ordinal <- as.double(cut(y, quantile(y, c(0, 0.3, 0.6, 1)), include.lowest = TRUE))
ordinal <- dbarts(
  x,
  y.ordinal,
  family = "ordinal",
  n.trees = 5L,
  control = control
)
invisible(ordinal$run(20L, 1L))
ordinal$storeState()
ordinal.own <- ordinal$state
gamma <- ordinal.own[[1L]]$thresholds
expect_identical(length(gamma), 2L)
expect_identical(gamma[1L], 0)
for (thresholds in list(c(0.1, gamma[2L]), c(0, 0), c(0, -gamma[2L]),
                        c(0, NA), c(0, Inf))) {
  edited <- ordinal.own
  edited[[1L]]$thresholds <- thresholds
  expect_error(
    ordinal$copy()$setState(edited),
    refused,
    info = paste(thresholds, collapse = " ")
  )
}
expect_true(ordinal$copy()$setState(ordinal.own))
