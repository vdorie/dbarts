# The posterior predictive draw for a weight-w logistic row is the number of
# successes among w iid bernoulli trials (binomial(w, p)), not a bernoulli
# draw multiplied by w - that older code gave a degenerate {0, w} PPD whose
# mean was coincidentally right but whose intervals/quantiles were not. See
# issue #79.

set.seed(2, sample.kind = "Rejection")
n <- 200L
x <- matrix(runif(n * 2L), n, 2L)
# a moderate slope keeps p away from 0/1, so a w=5 column's intermediate
# counts (1..4) are near-certain across hundreds of draws, not a coin flip
f <- 0.6 * x[, 1L] - 0.3
p <- plogis(f)
y <- rbinom(n, 1L, p)
w <- rep_len(c(1L, 3L, 5L), n)

fit <- bart2(
  y ~ x,
  weights = w,
  family = "logistic",
  n.samples = 500L,
  n.burn = 100L,
  n.trees = 50L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)

set.seed(11)
ppd <- extract(fit, type = "ppd")
ev <- extract(fit, type = "ev")

# 1. every column's draws stay within [0, w], and a w=5 column hits at least
# one intermediate value (not just the endpoints 0 and 5)
inRange <- vapply(
  seq_len(n),
  function(j) all(ppd[, j] >= 0 & ppd[, j] <= w[j]),
  logical(1L)
)
expect_true(all(inRange))

w5cols <- which(w == 5L)
hitsIntermediate <- vapply(
  w5cols,
  function(j) any(ppd[, j] %in% 1:4),
  logical(1L)
)
expect_true(all(hitsIntermediate))

# 3. mean sanity: colMeans(ppd) ~ w * colMeans(ev). Loose tolerance - this
# held even under the old {0, w} draw, so it only guards against a
# prob/size mixup, not the intervals the fix actually corrects.
expect_true(all(abs(colMeans(ppd) - w * colMeans(ev)) < 0.5))

# 6. reproducibility: identical seed, identical draws
set.seed(11)
ppd2 <- extract(fit, type = "ppd")
expect_equal(ppd, ppd2)

rm(inRange, w5cols, hitsIntermediate, ppd2)
rm(x, f, p, y, w, fit, ppd, ev, n)

# 2. column-alignment proof: only the observation with weight 5 may draw a
# value above 1; every weight-1 column must stay in {0, 1}. This guards the
# rep(weights, each = nrow(ev)) recycling - a transposed or mis-recycled
# size vector would leak the weight-5 support into the wrong column.
set.seed(3, sample.kind = "Rejection")
n <- 60L
x <- matrix(runif(n * 2L), n, 2L)
y <- rbinom(n, 1L, 0.5)
w <- c(5L, rep(1L, n - 1L))

fit.align <- bart2(
  y ~ x,
  weights = w,
  family = "logistic",
  n.samples = 300L,
  n.burn = 50L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
ppd.align <- extract(fit.align, type = "ppd")

expect_true(max(ppd.align[, 1L]) > 1)
expect_true(all(ppd.align[, -1L] %in% c(0, 1)))

rm(x, y, w, fit.align, ppd.align, n)

# 4. the unweighted binary path is unchanged: draws stay in {0, 1}
set.seed(4, sample.kind = "Rejection")
n <- 60L
x <- matrix(runif(n * 2L), n, 2L)
y <- rbinom(n, 1L, plogis(0.6 * x[, 1L] - 0.3))

fit.unweighted <- bart2(
  y ~ x,
  family = "logistic",
  n.samples = 100L,
  n.burn = 50L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
ppd.unweighted <- extract(fit.unweighted, type = "ppd")
expect_true(all(ppd.unweighted %in% c(0, 1)))

rm(x, y, fit.unweighted, ppd.unweighted, n)

# 5. the 3-d path (obs in the last dimension) is reachable through a public
# surface via a multi-chain fit with combineChains = FALSE; exercised once
# for shape and support, not statistical power
set.seed(5, sample.kind = "Rejection")
n <- 80L
x <- matrix(runif(n * 2L), n, 2L)
y <- rbinom(n, 1L, plogis(0.6 * x[, 1L] - 0.3))
w <- rep_len(c(1L, 3L, 5L), n)

fit.mc <- bart2(
  y ~ x,
  weights = w,
  family = "logistic",
  n.samples = 40L,
  n.burn = 20L,
  n.trees = 20L,
  n.chains = 2L,
  n.threads = 1L,
  verbose = FALSE
)
ppd.mc <- extract(fit.mc, type = "ppd", combineChains = FALSE)

expect_equal(length(dim(ppd.mc)), 3L)
expect_equal(dim(ppd.mc)[3L], n)
inRange.mc <- vapply(
  seq_len(n),
  function(j) all(ppd.mc[,, j] >= 0 & ppd.mc[,, j] <= w[j]),
  logical(1L)
)
expect_true(all(inRange.mc))

rm(x, y, w, fit.mc, ppd.mc, inRange.mc, n)
