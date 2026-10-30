# The live install of one forest's tree prior, through the bridge entry the
# forest door calls: installed before any run, a forest's base and power draw
# what creating with them draws, on each of three forests and on the variance
# forest, and every refusal names its reason and leaves the draws as they
# were.

set.seed(3L)
n <- 150L
x <- data.frame(x1 = runif(n), x2 = runif(n))
z <- rbinom(n, 1L, 0.5)
dose <- runif(n, 1, 5)
y <- x$x1 + z * x$x2 + 0.2 * dose + rnorm(n, 0, 0.3)
X <- cbind(x, z = z, dose = dose)
control <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.samples = 10L,
  seed = 5L,
  updateState = FALSE
)

installTreePrior <- function(sampler, forest, base, power) {
  .Call(
    dbarts:::C_dbarts_bartcore_setForestTreePrior,
    sampler$getPointer(),
    forest,
    base,
    power
  )
}

# forest `stated` (0-based) states base 0.5 and power 3 at creation; forest 0
# takes the fitting function's tree prior
threeForests <- function(stated = NA) {
  depth <- function(i) {
    if (identical(stated, i)) list(base = 0.5, power = 3) else list()
  }
  second <- do.call(
    dbartsForests$forest,
    c(list(basis = ~ factor(z), n.trees = 5L), depth(1L))
  )
  third <- do.call(
    dbartsForests$forest,
    c(list(basis = ~dose, n.trees = 4L), depth(2L))
  )
  do.call(
    dbarts,
    c(
      list(
        X,
        y,
        forests = list(dbartsForests$forest(n.trees = 7L), second, third),
        control = control
      ),
      depth(0L)
    )
  )
}

for (forest in 0:2) {
  created <- threeForests(forest)
  installed <- threeForests()
  plain <- threeForests()
  expect_null(installTreePrior(installed, forest, 0.5, 3))
  createdDraws <- created$run(20L, 10L)
  expect_identical(installed$run(20L, 10L), createdDraws)
  expect_false(identical(plain$run(20L, 10L), createdDraws))
}

withVariance <- function(...) {
  dbarts(
    x,
    y,
    variance = dbartsForests$varianceForest(n.trees = 8L, ...),
    control = control
  )
}
created <- withVariance(base = 0.5, power = 3)
installed <- withVariance()
plain <- withVariance()
expect_null(installTreePrior(installed, -1L, 0.5, 3))
createdDraws <- created$run(20L, 10L)
expect_identical(installed$run(20L, 10L), createdDraws)
expect_false(identical(plain$run(20L, 10L), createdDraws))

refusedVariance <- withVariance()
expect_error(
  installTreePrior(refusedVariance, 1L, 0.5, 3),
  "forest index out of range",
  fixed = TRUE
)
expect_error(
  installTreePrior(refusedVariance, -1L, 1.5, 3),
  "variance base must be in (0, 1)",
  fixed = TRUE
)
expect_error(
  installTreePrior(refusedVariance, -1L, 0.5, 0),
  "variance power must be positive",
  fixed = TRUE
)

indexMessage <- "forest index out of range"
baseMessage <- "tree prior base must be in (0, 1)"
powerMessage <- "tree prior power must be positive"
lengthMessage <- "a tree prior's base and power must each be a single number"
refusals <- list(
  list(-2L, 0.5, 3, indexMessage),
  list(3L, 0.5, 3, indexMessage),
  list(NA_integer_, 0.5, 3, indexMessage),
  list(1:2, 0.5, 3, indexMessage),
  list(integer(), 0.5, 3, indexMessage),
  list(-1L, 0.5, 3, "this sampler has no variance forest"),
  list(1L, 1, 3, baseMessage),
  list(1L, 0, 3, baseMessage),
  list(1L, NA_real_, 3, baseMessage),
  list(1L, 0.5, 0, powerMessage),
  list(1L, 0.5, -1, powerMessage),
  list(1L, 0.5, NaN, powerMessage),
  list(1L, c(0.5, 0.6), 3, lengthMessage),
  list(1L, 0.5, c(3, 4), lengthMessage)
)
refused <- threeForests()
untouched <- threeForests()
for (refusal in refusals) {
  expect_error(
    installTreePrior(refused, refusal[[1L]], refusal[[2L]], refusal[[3L]]),
    refusal[[4L]],
    fixed = TRUE
  )
}
expect_identical(refused$run(20L, 10L), untouched$run(20L, 10L))
