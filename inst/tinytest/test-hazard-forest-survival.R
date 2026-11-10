# survivalProbabilities on a hazard fit with several forests: the forests
# replay separately and recombine with the fit's glue, as predict does on the
# person-period rows, and a basis is one value per subject, repeated over its
# periods as the offset is.

suppressMessages(library(survival))

set.seed(20261110L)
n <- 40L
d <- data.frame(x = rnorm(n), z = rnorm(n), w = runif(n))
d$time <- sample(1:4, n, TRUE)
d$status <- rbinom(n, 1L, 0.7)
fitArgs <- list(
  family = "hazard",
  n.samples = 12L,
  n.burn = 6L,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 8L,
  verbose = FALSE,
  keepTrees = TRUE
)
fit <- do.call(
  bart,
  c(list(Surv(time, status) ~ x + forest(basis = z), d), fitArgs)
)
K <- length(fit$periods)

# the curves of the expanded person-period rows, from predict's hazards
byPredict <- function(object, newdata, subjectBases = NULL, ...) {
  m <- nrow(newdata)
  rows <- newdata[rep(seq_len(m), times = K), , drop = FALSE]
  rows$period <- rep(seq_len(K), each = m)
  haz <- predict(
    object,
    rows,
    type = "ev",
    combineChains = FALSE,
    bases = if (!is.null(subjectBases)) rep(subjectBases, times = K),
    ...
  )
  haz <- array(haz, c(dim(haz)[1:2], m, K))
  out <- array(0, c(dim(haz)[1:2], K, m))
  for (k in seq_len(K)) {
    out[,, k, ] <- apply(1 - haz[,,, seq_len(k), drop = FALSE], 1:3, prod)
  }
  out
}

nd <- d[1:5, ]
# new subjects, the basis read from newdata
curves <- survivalProbabilities(fit, newdata = nd, combineChains = FALSE)
expect_equal(dim(curves), c(1L, 12L, K, 5L))
expect_equal(curves, byPredict(fit, nd), check.attributes = FALSE)

# a value basis, one per subject
shifted <- nd$z + 1
expect_equal(
  survivalProbabilities(
    fit,
    newdata = nd,
    bases = shifted,
    combineChains = FALSE
  ),
  byPredict(fit, nd, shifted),
  check.attributes = FALSE
)
expect_equal(
  survivalProbabilities(fit, newdata = nd, bases = list(NULL, shifted)),
  survivalProbabilities(fit, newdata = nd, bases = shifted)
)
expect_true(
  max(abs(
    survivalProbabilities(fit, newdata = nd, bases = shifted) -
      survivalProbabilities(fit, newdata = nd)
  )) >
    0
)

# a named basis is read from newdata
expect_equal(
  survivalProbabilities(fit, newdata = nd, bases = "w"),
  survivalProbabilities(fit, newdata = nd, bases = nd$w)
)

# the training subjects: their own bases, or others given
own <- survivalProbabilities(fit)
expect_equal(dim(own), c(12L, K, n))
expect_equal(
  survivalProbabilities(fit, bases = d$z),
  own,
  check.attributes = FALSE
)
expect_equal(
  survivalProbabilities(fit, bases = d$z + 1),
  survivalProbabilities(fit, newdata = d, bases = d$z + 1),
  check.attributes = FALSE
)

# a fit of two forests and no basis
plain <- do.call(
  bart,
  c(list(Surv(time, status) ~ forest(x) + forest(z), d), fitArgs)
)
expect_equal(dim(survivalProbabilities(plain)), c(12L, K, n))
expect_equal(
  survivalProbabilities(plain, newdata = nd, combineChains = FALSE),
  byPredict(plain, nd),
  check.attributes = FALSE
)

# a one-forest fit has no basis to take
one <- do.call(bart, c(list(Surv(time, status) ~ x, d), fitArgs))
expect_error(
  survivalProbabilities(one, newdata = nd, bases = nd$z),
  "only meaningful on an amplitude-coupled multi-forest fit"
)
