# S0 pins for the latent-subset-mask arc: fixed points the active-rows
# channel must not move without disturbing. Extended at S1.

set.seed(20260812L)
n <- 200L
x <- matrix(runif(n * 2L), n, 2L, dimnames = list(NULL, c("x1", "x2")))
y.binary <- as.double(x[, 1L] + rnorm(n, 0, 0.3) > 0.5)

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  updateState = FALSE,
  rngSeed = 7L
)

# Pin: $setWeights on a probit sampler refuses post-creation, unconditionally
# on the value - even the all-ones vector that dbarts(..., weights = rep(1, n))
# accepts and normalizes away at creation is refused here
# (refuseBinaryWeightChange, R_interface_bartcore.cpp, reached from
# bartcore_setWeights before its length or value checks run).
sampler.probit <- dbarts::dbarts(
  x,
  y.binary,
  family = "probit",
  control = control
)
expect_error(
  sampler.probit$setWeights(rep(1, n)),
  "probit models do not support case weights"
)
expect_error(
  sampler.probit$setWeights(runif(n, 0.5, 1.5)),
  "probit models do not support case weights"
)

rm(sampler.probit)

# --- S1: the $setActiveRows channel ---------------------------------------
# The mask says "row i is not in the data set for this sampler this sweep": the
# row leaves every sufficient statistic, every family-level parameter update
# and its own latent draw, but keeps its leaf occupancy and its fitted value.

y <- as.double(x[, 1L] + 2 * x[, 2L] + rnorm(n, 0, 0.3))
w <- 0.5 + runif(n)
a <- as.double(seq_len(n) %% 4L != 1L)

makeSampler <- function(weights = NULL, ...) {
  dbarts::dbarts(
    x,
    y,
    weights = weights,
    control = control,
    sigma = 1,
    n.samples = 10L,
    ...
  )
}

# T1, gaussian, BITWISE: setActiveRows(a) on a sampler carrying w draws exactly
# what setWeights(w * a) draws with no mask. The comparator is w * a, never "no
# weights" - an unweighted sampler is not bitwise one carrying rep(1, n),
# because only the unweighted path takes the fused node-average roll.
masked <- makeSampler(w)
masked$setActiveRows(a)
composed <- makeSampler(w * a)
draws.masked <- masked$run(20L, 10L)
draws.composed <- composed$run(20L, 10L)
expect_identical(draws.masked$train, draws.composed$train)
expect_identical(draws.masked$sigma, draws.composed$sigma)
rm(masked, composed, draws.masked, draws.composed)

# T1, Student-t: both arms draw lambda at every row on the same build, so the
# mask annihilates the composite without moving the stream
masked.t <- makeSampler(w, resid.dist = dbarts:::student(df = 4))
masked.t$setActiveRows(a)
composed.t <- makeSampler(w * a, resid.dist = dbarts:::student(df = 4))
draws.masked.t <- masked.t$run(20L, 10L)
draws.composed.t <- composed.t$run(20L, 10L)
expect_identical(draws.masked.t$train, draws.composed.t$train)
expect_identical(draws.masked.t$sigma, draws.composed.t$sigma)
rm(masked.t, composed.t, draws.masked.t, draws.composed.t)

# T2(b), the all-ones normalizer, on BOTH surfaces this slice ships - the R5
# method and the dbarts::: bridge entry. The normalizer is in the ENGINE, so
# neither surface can be the one that is right. (The flat entry is the third
# surface; it lands with the dbarts.h reshape.)
#
# The comparator is deliberately an UNWEIGHTED sampler: gaussian-unweighted,
# probit and ordinal are the three configurations still on the fused
# node-average path, which any non-null precision pointer leaves. A weighted
# arm would pass on the values alone and prove nothing about the normalizer.
plain <- makeSampler()
draws.plain <- plain$run(20L, 10L)

r5.ones <- makeSampler()
r5.ones$setActiveRows(rep(1, n))
expect_identical(r5.ones$run(20L, 10L)$train, draws.plain$train)

# the bridge entry drives its own engine handle (bartcoreSampler builds a
# fresh sampler), so its arm runs through that handle end to end
bc.plain <- dbarts:::bartcoreSampler(makeSampler())
bc.ones <- dbarts:::bartcoreSampler(makeSampler())
dbarts:::bartcoreSetActiveRows(bc.ones, rep(1, n))
expect_identical(
  dbarts:::bartcoreRun(bc.ones, 20L, 10L)$train,
  dbarts:::bartcoreRun(bc.plain, 20L, 10L)$train
)

# a mask returning to all ones CLEARS, restoring the pre-mask pointer BY
# IDENTITY - so the fused path comes back - and NULL clears the same way
returned <- makeSampler()
returned$setActiveRows(a)
returned$setActiveRows(rep(1, n))
expect_identical(returned$run(20L, 10L)$train, draws.plain$train)

cleared <- makeSampler()
cleared$setActiveRows(a)
cleared$setActiveRows(NULL)
expect_identical(cleared$run(20L, 10L)$train, draws.plain$train)
rm(r5.ones, bc.plain, bc.ones, returned, cleared)

# A refusal must install nothing AND perturb nothing: a fractional value, an
# NA, and a wrong length each leave the sampler exactly where it was. Unweighted
# again, so a partial application would show up by leaving the fused path.
for (bad in list(replace(a, 2L, 0.5), replace(a, 2L, NA_real_), a[-1L])) {
  refused <- makeSampler()
  expect_error(refused$setActiveRows(bad))
  expect_identical(refused$run(20L, 10L)$train, draws.plain$train)
}
rm(refused, bad)

# The bridge enforces the same values with no R validation in front of it,
# which is the point of putting the scan in the engine
bc.bad <- dbarts:::bartcoreSampler(makeSampler())
expect_error(
  dbarts:::bartcoreSetActiveRows(bc.bad, replace(a, 2L, 0.5)),
  "exactly 0 or 1"
)
expect_identical(
  dbarts:::bartcoreRun(bc.bad, 20L, 10L)$train,
  draws.plain$train
)
rm(bc.bad)

# The two setters are absolute and INDEPENDENT: either order leaves the same
# composite, and setResponse/setOffset do not disturb the mask.
w2 <- 0.5 + runif(n)
order.aw <- makeSampler(w)
order.aw$setActiveRows(a)
order.aw$setWeights(w2)
order.wa <- makeSampler(w)
order.wa$setWeights(w2)
order.wa$setActiveRows(a)
expect_identical(order.aw$run(20L, 10L)$train, order.wa$run(20L, 10L)$train)

survives <- makeSampler(w)
survives$setActiveRows(a)
survives$setOffset(rep(0.25, n))
kept <- makeSampler(w * a)
kept$setOffset(rep(0.25, n))
expect_identical(survives$run(20L, 10L)$train, kept$run(20L, 10L)$train)

# setData CLEARS: the mask is length-n and setData may change n. Compared
# against a twin that took the same setData with no mask installed, since
# setData cold-initializes state a never-swapped sampler does not have.
newData <- dbarts::dbartsData(x, y, weights = w)
dropped <- makeSampler(w)
dropped$setActiveRows(a)
dropped$setData(newData)
swapped <- makeSampler(w)
swapped$setData(newData)
expect_identical(dropped$run(20L, 10L)$train, swapped$run(20L, 10L)$train)
rm(order.aw, order.wa, survives, kept, dropped, swapped, newData, w2)

# F6: an all-zeros mask is ACCEPTED and runs - the natural answer for a host
# whose stratum has emptied. Every forest sits at its prior and every row still
# receives a fit, which is what makes the channel worth more than compaction.
empty <- makeSampler(w)
empty$setActiveRows(rep(0, n))
draws.empty <- empty$run(20L, 10L)
expect_true(all(is.finite(draws.empty$train)))
expect_true(all(is.finite(draws.empty$sigma) & draws.empty$sigma > 0))
expect_equal(dim(draws.empty$train), c(n, 10L))
rm(empty, draws.empty)

# T2(c), probit: substituting arbitrary labels at the INACTIVE rows leaves
# every active row's recorded draw bitwise. Because the truncated-normal
# primitive is a rejection sampler, this fails outright if an inactive row's
# latent is drawn and discarded rather than skipped.
probit.a <- dbarts::dbarts(x, y.binary, family = "probit", control = control)
probit.a$setActiveRows(a)
probit.b <- dbarts::dbarts(
  x,
  ifelse(a == 0, 1 - y.binary, y.binary),
  family = "probit",
  control = control
)
probit.b$setActiveRows(a)
train.a <- probit.a$run(20L, 10L)$train
train.b <- probit.b$run(20L, 10L)$train
expect_identical(train.a[a == 1, ], train.b[a == 1, ])
# An INACTIVE row's own fit is not claimed bitwise and is not: the running
# residual reconstructs totalFits as working - residual, and the working
# response there is the row's own stale latent, so the reported value carries
# reconstruction rounding. The trees, the leaf draws and every active row are
# untouched.
expect_true(max(abs(train.a[a == 0, ] - train.b[a == 0, ])) < 1e-12)
rm(probit.a, probit.b, train.a, train.b)

# The refusal for a family v1 does not build names that family and the reason
expect_error(
  dbarts::dbarts(
    x,
    y.binary,
    family = "logistic",
    control = control
  )$setActiveRows(a),
  "not implemented for logistic"
)

# Decorators compose with no edit of their own, which is a claim and so a pin.
# Heteroscedastic: the mean weights divide the composed weight by s^2(x), so a
# masked row is zero there too, and the variance forest carries the same
# composed weight.
hetero <- dbarts::dbarts(
  x,
  y,
  variance = TRUE,
  n.trees.variance = 10L,
  control = control,
  sigma = 1,
  n.samples = 10L
)
hetero$setActiveRows(a)
draws.hetero <- hetero$run(20L, 10L)
expect_true(all(is.finite(draws.hetero$train)))
rm(hetero, draws.hetero)

# A GP leaf under a mask: the pin is STRUCTURAL, never an equality - the clean
# and positive-subset paths consume different numbers of variates per leaf, so
# installing a mask moves the stream at every leaf by construction.
gp <- dbarts::dbarts(
  x,
  y,
  control = dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    updateState = FALSE,
    rngSeed = 7L
  ),
  sigma = 1,
  n.samples = 10L,
  node.prior = dbarts:::gp("x1", max.leaf.size = 100L)
)
gp$setActiveRows(a)
draws.gp <- gp$run(20L, 10L)
expect_true(all(is.finite(draws.gp$train)))
expect_equal(dim(draws.gp$train), c(n, 10L))
# and it must DETECT a broken routing rather than merely running: with the
# zero-weight branch bypassed an inactive row carries infinite noise variance,
# every containing leaf's marginal degenerates, no tree ever splits and the fit
# collapses to the response mean.
expect_true(sum(draws.gp$varcount) > 0L)
expect_true(
  sqrt(mean((rowMeans(draws.gp$train)[a == 1] - y[a == 1])^2)) <
    0.6 * sd(y[a == 1])
)
rm(gp, draws.gp)

rm(
  n,
  x,
  y,
  y.binary,
  w,
  a,
  control,
  makeSampler,
  plain,
  draws.plain
)
