# The posterior predictive draw for a weight-w logistic row is the number of
# successes among w iid bernoulli trials (binomial(w, p)), not a bernoulli
# draw multiplied by w - that older code gave a degenerate {0, w} PPD whose
# mean was coincidentally right but whose intervals/quantiles were not. See
# issue #79.

set.seed(2, sample.kind = "Rejection")
# Binomial counts: at each predictor point 1, 5 or 9 independent trials,
# entered as a row of successes and a row of failures, a count of 0 left out,
# so a row's count never repeats one trial. A row of weight 5 is mostly one
# side of a 4 to 5 split of nine trials.
n.points <- 200L
x.points <- matrix(runif(n.points * 2L), n.points, 2L)
# a moderate slope keeps the true p away from 0 and 1, and the fitted
# probability of each row with it; at the boundary an intermediate count is
# not near-certain, so the check below reads only the rows whose posterior
# mean probability is interior and requires that there are enough of them
p.points <- plogis(0.6 * x.points[, 1L] - 0.3)
trials <- rep_len(c(1L, 5L, 9L), n.points)
successes <- rbinom(n.points, trials, p.points)
counts <- c(successes, trials - successes)
rows <- which(counts > 0L)
point <- rep(seq_len(n.points), 2L)[rows]
x <- x.points[point, , drop = FALSE]
y <- rep(c(1L, 0L), each = n.points)[rows]
w <- counts[rows]
n <- length(rows)

fit <- bart(
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

# 1. every column's draws stay within [0, w], and an interior w=5 column hits
# at least one intermediate value (not just the endpoints 0 and 5). The old
# degenerate {0, w} draw produced NO intermediate value on any column, so
# restricting to interior columns keeps the discrimination whole.
inRange <- vapply(
  seq_len(n),
  function(j) all(ppd[, j] >= 0 & ppd[, j] <= w[j]),
  logical(1L)
)
expect_true(all(inRange))

w5cols <- which(w == 5L)
interior <- w5cols[
  colMeans(ev[, w5cols, drop = FALSE]) > 0.1 &
    colMeans(ev[, w5cols, drop = FALSE]) < 0.9
]
expect_true(length(interior) >= 30L)
hitsIntermediate <- vapply(
  interior,
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
rm(n.points, x.points, p.points, trials, successes, counts, rows, point)
rm(x, y, w, fit, ppd, ev, n)

# 2. column-alignment proof: only the observation with weight 5 may draw a
# value above 1; every weight-1 column must stay in {0, 1}. This guards the
# rep(weights, each = nrow(ev)) recycling - a transposed or mis-recycled
# size vector would leak the weight-5 support into the wrong column.
set.seed(3, sample.kind = "Rejection")
n <- 60L
x <- matrix(runif(n * 2L), n, 2L)
y <- rbinom(n, 1L, 0.5)
w <- c(5L, rep(1L, n - 1L))

fit.align <- bart(
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

fit.unweighted <- bart(
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

fit.mc <- bart(
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
