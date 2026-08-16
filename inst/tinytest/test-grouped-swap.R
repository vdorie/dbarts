# The grouped response swap. On a sampler carrying random intercepts,
# setResponse and setOffset are ACCEPTED at the pinned scale and REFUSED at
# updateScale = TRUE under a base family whose response transform is derived
# from the data (gaussian, Student-t, aft), where the effects b and their scale
# tau - held against that transform and converted by nothing - would silently
# come to mean something else. A grouped probit sampler has no such transform
# and takes TRUE as the no-op it always was. setData stays refused at every
# scale.
#
# The arithmetic of the swap itself is pinned in tests/cpp (test_model.cpp,
# "grouped response swap"): these cells gate the bridge predicate - which
# family, which conduit, which scale - and that an accepted swap reaches the
# model rather than being dropped.

set.seed(6021L)
n <- 80L
p <- 3L
x <- matrix(runif(n * p), n, p)
g <- rep_len(seq_len(4L), n)
f <- 2 * sin(pi * x[, 1L]) - x[, 2L]
b <- c(-1.2, -0.4, 0.5, 1.1)
y <- f + b[g] + rnorm(n, sd = 0.5)
# a DIFFERENT mean structure on the same design, so a swap that reached the
# model and one that did not are told apart by the fit alone
y2 <- -f + b[g] + rnorm(n, sd = 0.5)

# rbart_vi's internal control attribute is the only route to a grouped
# sampler; a fixed seed makes each chain's generator deterministic, which
# the probit no-op cell compares bitwise
groupedControl <- function(rel.scale = 1, ...) {
  ctrl <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    updateState = FALSE,
    seed = 61L,
    ...
  )
  attr(ctrl, "bartcore.groups") <- list(
    indices = as.integer(g),
    n.groups = 4L,
    prior = "cauchy",
    rel.scale = rel.scale,
    n.steps = 1L
  )
  ctrl
}

# --- gaussian: the pinned-scale swap is accepted and MOVES the fit

sampler <- dbarts(x, y, control = groupedControl(rel.scale = sd(y)))
before <- sampler$run(100L, 1L)
expect_silent(sampler$setResponse(y2))
expect_equal(sampler$data@y, y2)
after <- sampler$run(200L, 1L)
expect_true(max(abs(after$train - before$train)) > 0.5)
# and moves it TOWARD the new response, not just away from the old one
expect_true(cor(as.numeric(after$train), -f) > 0.5)
expect_true(cor(as.numeric(before$train), f) > 0.5)
# the group block is still drawing after the swap
expect_true(all(is.finite(after$ranef)))
expect_true(all(after$tau > 0))

# --- gaussian: updateScale = TRUE is refused on both response-side conduits,
# and the refusal names the two quantities it is protecting

expect_error(sampler$setResponse(y2, updateScale = TRUE), pattern = "tau")
expect_error(
  sampler$setResponse(y2, updateScale = TRUE),
  pattern = "random intercepts b"
)
expect_error(
  sampler$setOffset(rep_len(0.25, n), updateScale = TRUE),
  pattern = "tau"
)
# the same offset at the pinned scale is accepted
expect_silent(sampler$setOffset(rep_len(0.25, n)))
expect_equal(sampler$data@offset, rep_len(0.25, n))
expect_silent(sampler$setOffset(NULL))

# --- Student-t: ResponseFamily reports gaussian for it, so it is covered by
# the same predicate without a member of its own

samplerT <- dbarts(
  x,
  y,
  control = groupedControl(rel.scale = sd(y)),
  resid.dist = student(5)
)
expect_error(samplerT$setResponse(y2, updateScale = TRUE), pattern = "tau")
expect_error(
  samplerT$setOffset(rep_len(0.25, n), updateScale = TRUE),
  pattern = "tau"
)
expect_silent(samplerT$setResponse(y2))

# --- aft: a log survival time is a data-derived transform too

status <- rep_len(c(1L, 1L, 0L), n)
survTime <- exp(f + b[g] + rnorm(n, sd = 0.3))
samplerAft <- dbarts(
  x,
  cbind(survTime, status),
  control = groupedControl(rel.scale = sd(log(survTime))),
  family = "aft"
)
expect_error(
  samplerAft$setResponse(log(survTime), updateScale = TRUE),
  pattern = "tau"
)
expect_error(
  samplerAft$setOffset(rep_len(0.25, n), updateScale = TRUE),
  pattern = "tau"
)
expect_silent(samplerAft$setResponse(log(survTime) + 0.1))

# --- probit: no data-derived transform, so updateScale = TRUE is accepted and
# is exactly the no-op it always was - two identically seeded samplers, one
# swapped at TRUE and one at FALSE, draw the same chain bit for bit

yBin <- as.double(rbinom(n, 1L, pnorm(f + b[g])))
yBin2 <- as.double(rbinom(n, 1L, pnorm(-f + b[g])))
probitAnchored <- dbarts(
  x,
  yBin,
  control = groupedControl(rel.scale = 0.5),
  family = "probit"
)
probitPinned <- dbarts(
  x,
  yBin,
  control = groupedControl(rel.scale = 0.5),
  family = "probit"
)
expect_silent(probitAnchored$setResponse(yBin2, updateScale = TRUE))
probitPinned$setResponse(yBin2)
expect_identical(
  probitAnchored$run(20L, 5L)$train,
  probitPinned$run(20L, 5L)$train
)
expect_silent(probitAnchored$setOffset(rep_len(0.25, n), updateScale = TRUE))

# --- setData is refused at every scale: it may change the row count, and the
# per-observation grouping is fixed at creation

expect_error(
  sampler$setData(dbartsData(x, y2)),
  pattern = "grouped random effects"
)

# --- statistical agreement: a sampler CREATED on y2 and one created on a
# permutation of y2 and then SWAPPED to it target the same posterior. The
# permutation keeps the response transform, the leaf calibration and the tau
# prior scale identical between the two, so the only difference is which
# response the chain was conditioned on before the swap. NEVER bitwise - the
# streams have seen different data - so the comparison is of posterior means
# after a long burn.

yPermuted <- y2[sample.int(n)]
direct <- dbarts(x, y2, control = groupedControl(rel.scale = sd(y2)))
swapped <- dbarts(x, yPermuted, control = groupedControl(rel.scale = sd(y2)))
swapped$setResponse(y2)
fitsDirect <- rowMeans(direct$run(1000L, 300L)$train)
fitsSwapped <- rowMeans(swapped$run(1000L, 300L)$train)
expect_true(cor(fitsDirect, fitsSwapped) > 0.9)
expect_true(max(abs(fitsDirect - fitsSwapped)) < 0.5 * diff(range(y2)))
expect_false(identical(fitsDirect, fitsSwapped))
