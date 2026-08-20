# Fixed points the active-row mask channel must not move without disturbing,
# plus the channel's own semantics pinned one assertion at a time below.

set.seed(20260812L)
n <- 200L
x <- matrix(runif(n * 2L), n, 2L, dimnames = list(NULL, c("x1", "x2")))
y.binary <- as.double(x[, 1L] + rnorm(n, 0, 0.3) > 0.5)

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  updateState = FALSE,
  seed = 7L
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

# --- the $setActiveRows channel ---------------------------------------
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

# Gaussian, bitwise: setActiveRows(a) on a sampler carrying w draws exactly
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

# Student-t, bitwise: both arms draw lambda at every row on the same build,
# so the mask annihilates the composite without moving the stream
masked.t <- makeSampler(w, resid.dist = dbarts:::student(df = 4))
masked.t$setActiveRows(a)
composed.t <- makeSampler(w * a, resid.dist = dbarts:::student(df = 4))
draws.masked.t <- masked.t$run(20L, 10L)
draws.composed.t <- composed.t$run(20L, 10L)
expect_identical(draws.masked.t$train, draws.composed.t$train)
expect_identical(draws.masked.t$sigma, draws.composed.t$sigma)
rm(masked.t, composed.t, draws.masked.t, draws.composed.t)

# Bayesian causal forest, bitwise: a BCF sampler's two forests share one
# gaussian response, which composes the mask exactly as the single-forest
# gaussian arm above, and the mask reaches both forests through
# composeForestWeights with no code of its own - so setActiveRows(a) on a
# BCF sampler carrying w draws exactly what setWeights(w * a) draws with no
# mask, on train and sigma both.
z.bcf <- as.double(seq_len(n) %% 2L)
masked.bcf <- dbarts::dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z.bcf))),
  weights = w,
  control = control,
  sigma = 1
)
masked.bcf$setActiveRows(a)
composed.bcf <- dbarts::dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z.bcf))),
  weights = w * a,
  control = control,
  sigma = 1
)
draws.masked.bcf <- masked.bcf$run(20L, 10L)
draws.composed.bcf <- composed.bcf$run(20L, 10L)
expect_identical(draws.masked.bcf$train, draws.composed.bcf$train)
expect_identical(draws.masked.bcf$sigma, draws.composed.bcf$sigma)
rm(masked.bcf, composed.bcf, draws.masked.bcf, draws.composed.bcf, z.bcf)

# The all-ones normalizer, on BOTH R-facing surfaces - the R-level method and the
# dbarts::: bridge entry. The normalizer is in the ENGINE, so neither surface
# can be the one that is right. (The flat entry is the third surface; it
# lands with the dbarts.h reshape.)
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
  expect_error(
    refused$setActiveRows(bad),
    "'active' (must be all 0 or 1|cannot be NA|must have the same length as 'y')"
  )
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

# An all-zeros mask is ACCEPTED and runs - the natural answer for a host
# whose stratum has emptied. Every forest sits at its prior and every row still
# receives a fit, which is what makes the channel worth more than compaction.
empty <- makeSampler(w)
empty$setActiveRows(rep(0, n))
draws.empty <- empty$run(20L, 10L)
expect_true(all(is.finite(draws.empty$train)))
expect_true(all(is.finite(draws.empty$sigma) & draws.empty$sigma > 0))
expect_equal(dim(draws.empty$train), c(n, 10L))
rm(empty, draws.empty)

# Probit: substituting arbitrary labels at the INACTIVE rows leaves
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

# Decorators compose with no edit of their own, which is a claim and so a pin.
# Heteroscedastic: the mean weights divide the composed weight by s^2(x), so a
# masked row is zero there too, and the variance forest carries the same
# composed weight. Pinned BITWISE against setWeights(w * a) rather than for
# finiteness, which a mask that never reached the variance forest's sufficient
# statistics would also satisfy.
heteroSampler <- function(weights) {
  dbarts::dbarts(
    x,
    y,
    weights = weights,
    variance = dbarts::varianceForest(n.trees = 10L),
    control = control,
    sigma = 1,
    n.samples = 10L
  )
}
hetero <- heteroSampler(w)
hetero$setActiveRows(a)
hetero.composed <- heteroSampler(w * a)
draws.hetero <- hetero$run(20L, 10L)
draws.hetero.composed <- hetero.composed$run(20L, 10L)
expect_true(all(is.finite(draws.hetero$train)))
expect_identical(draws.hetero$train, draws.hetero.composed$train)
expect_identical(draws.hetero$varcount, draws.hetero.composed$varcount)
rm(hetero, hetero.composed, draws.hetero, draws.hetero.composed, heteroSampler)

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
    seed = 7L
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

# --- logistic, nbinom and aft -----------------------------------------
# Each family's oracle substitutes arbitrary in-support responses at the
# INACTIVE rows and checks that every ACTIVE row's recorded draw stays
# bitwise, at the SAMPLER level. Because every one of these latents comes
# from a rejection sampler, an arm fails outright if an inactive row's draw
# is taken and discarded rather than skipped. An inactive row's own reported
# fit is not claimed bitwise, for the reason the probit arm above states.

# Ordinal, at the sampler level beside the kernel-level coverage in
# tests/cpp: the cutpoint pass is the surface the kernels cover one at a time,
# and both its proposal (a count tally) and its target (a category-wise
# likelihood) read the mask, so a mask reaching the latents but not the
# cutpoints parts here.
levels.ord <- c("lo", "mid", "hi")
z.ord <- 2 * (x[, 1L] - 0.5) + rnorm(n)
codes.ord <- 1L + (z.ord > 0) + (z.ord > 0.8)
codes.other <- codes.ord
codes.other[a == 0] <- 1L + (codes.ord[a == 0] %% 3L)
ordinalSampler <- function(codes) {
  sampler <- dbarts::dbarts(
    x,
    ordered(levels.ord[codes], levels = levels.ord),
    family = "ordinal",
    control = control
  )
  sampler$setActiveRows(a)
  sampler
}

# logistic: the channel is not redundant with the count weights - a zero count
# is refused at creation, no count change is accepted after it, and a
# zero-count row would still consume one Polya-Gamma variate. The inactive
# rows' TRIAL COUNTS are substituted alongside their labels, since the count is
# what the Polya-Gamma draw itself reads (the label enters only the working
# response), so the arm sees the draw and not just the statistic.
counts.binary <- rep(1, n)
counts.binary.other <- counts.binary
counts.binary.other[a == 0] <- 3
logisticSampler <- function(y.logistic, trials = counts.binary, mask = a) {
  sampler <- dbarts::dbarts(
    x,
    y.logistic,
    weights = trials,
    family = "logistic",
    control = control
  )
  if (!is.null(mask)) {
    sampler$setActiveRows(mask)
  }
  sampler
}

# nbinom: the substituted counts move the dispersion kernel too, so this arm
# also fails unless the count histogram behind L_k is rebuilt over the ACTIVE
# rows at every mask change.
counts <- as.double(rpois(n, 3))
counts.other <- counts
counts.other[a == 0] <- counts[a == 0] + 5
nbinomSampler <- function(y.count, mask = a) {
  sampler <- dbarts::dbarts(x, y.count, family = "nbinom", control = control)
  if (!is.null(mask)) {
    sampler$setActiveRows(mask)
  }
  sampler
}

for (arms in list(
  list(ordinalSampler(codes.ord), ordinalSampler(codes.other)),
  list(
    logisticSampler(y.binary),
    logisticSampler(ifelse(a == 0, 1 - y.binary, y.binary), counts.binary.other)
  ),
  list(nbinomSampler(counts), nbinomSampler(counts.other))
)) {
  train.a <- arms[[1L]]$run(20L, 10L)$train
  train.b <- arms[[2L]]$run(20L, 10L)$train
  expect_identical(train.a[a == 1, ], train.b[a == 1, ])
  expect_true(max(abs(train.a[a == 0, ] - train.b[a == 0, ])) < 1e-12)
}
rm(arms, train.a, train.b)

# aft: the censored rows' log-time redraw is the latent. The substituted times
# stay strictly inside the observed range and leave the extremes alone, because
# the response transform is the FULL-data one by design and a new extreme would
# move both arms' scale rather than just one row.
status <- as.double(seq_len(n) %% 3L != 1L)
aftHandle <- function(logTime, mask = a) {
  sampler <- dbarts::dbarts(
    x,
    logTime,
    control = control,
    sigma = 1,
    n.samples = 10L
  )
  ctrl <- sampler$control
  attr(ctrl, "bartcore.survival") <- status
  sampler$control <- ctrl
  handle <- dbarts:::bartcoreSampler(sampler, family = "aft")
  if (!is.null(mask)) {
    dbarts:::bartcoreSetActiveRows(handle, mask)
  }
  handle
}
y.other <- y
substitutable <- a == 0
substitutable[c(which.min(y), which.max(y))] <- FALSE
y.other[substitutable] <- mean(y)
aft.train.a <- dbarts:::bartcoreRun(aftHandle(y), 20L, 10L)$train
aft.train.b <- dbarts:::bartcoreRun(aftHandle(y.other), 20L, 10L)$train
expect_identical(aft.train.a[a == 1, ], aft.train.b[a == 1, ])
expect_true(max(abs(aft.train.a[a == 0, ] - aft.train.b[a == 0, ])) < 1e-12)

# The all-ones normalizer for logistic, nbinom and aft: the engine's
# normalizer clears an all-ones mask, and each serves its pre-mask precision
# pointer by identity when it does - a new branch per family, so a new arm
# per family.
expect_identical(
  logisticSampler(y.binary, mask = rep(1, n))$run(20L, 10L)$train,
  logisticSampler(y.binary, mask = NULL)$run(20L, 10L)$train
)
expect_identical(
  nbinomSampler(counts, rep(1, n))$run(20L, 10L)$train,
  nbinomSampler(counts, NULL)$run(20L, 10L)$train
)
expect_identical(
  dbarts:::bartcoreRun(aftHandle(y, rep(1, n)), 20L, 10L)$train,
  dbarts:::bartcoreRun(aftHandle(y, NULL), 20L, 10L)$train
)

# --- multinomial, GLOBAL only -----------------------------------------
# The mask lands on the softmax COUPLING, not on the response, which holds no
# precisions of its own: an inactive row's K interleaved Polya-Gamma draws are
# skipped and its composed precision is zero in every category. The R-level object
# never builds a multinomial sampler, so every arm here drives the dbarts:::
# handle, which is the surface this family's mask ships on.
K.mn <- 3L
labels.mn <- codes.ord - 1L
counts.mn <- matrix(0L, n, K.mn)
counts.mn[cbind(seq_len(n), labels.mn + 1L)] <- 1L + seq_len(n) %% 3L
# The substituted arm moves the successes AND the trials at the inactive rows:
# PG(n_i, psi) sums n_i variates, so a mask that reached the working response
# but not the draw count would part here, exactly as it does for logistic
# above.
counts.mn.other <- counts.mn
counts.mn.other[a == 0, ] <- counts.mn[a == 0, c(2L, 3L, 1L)] + 2L

multinomialHandle <- function(counts, mask = a) {
  handle <- dbarts:::bartcoreMultinomialCountSampler(
    dbarts::dbarts(x, as.double(labels.mn), control = control),
    counts,
    K = K.mn
  )
  if (!is.null(mask)) {
    dbarts:::bartcoreSetActiveRows(handle, mask)
  }
  handle
}

# Substituting the inactive rows' counts leaves every ACTIVE row's recorded
# softmax bitwise. Polya-Gamma is a rejection sampler, so this fails outright
# if an inactive row's K draws are taken and discarded, not skipped.
train.mn.a <- dbarts:::bartcoreRun(multinomialHandle(counts.mn), 20L, 10L)$train
train.mn.b <- dbarts:::bartcoreRun(
  multinomialHandle(counts.mn.other),
  20L,
  10L
)$train
expect_identical(train.mn.a[a == 1, , ], train.mn.b[a == 1, , ])
expect_true(max(abs(train.mn.a[a == 0, , ] - train.mn.b[a == 0, , ])) < 1e-12)

# The all-ones normalizer clears here too, and the coupling serves its
# unmasked precisions when it does
expect_identical(
  dbarts:::bartcoreRun(multinomialHandle(counts.mn, rep(1, n)), 20L, 10L)$train,
  dbarts:::bartcoreRun(multinomialHandle(counts.mn, NULL), 20L, 10L)$train
)

# an all-zeros mask runs: every category forest sits at its prior and every row
# still gets its K reported probabilities, which stay a simplex
train.mn.empty <- dbarts:::bartcoreRun(
  multinomialHandle(counts.mn, rep(0, n)),
  20L,
  10L
)$train
expect_true(all(is.finite(train.mn.empty)))
expect_true(max(abs(apply(train.mn.empty, c(1L, 3L), sum) - 1)) < 1e-12)

# the engine's value scan is under this surface too
expect_error(
  dbarts:::bartcoreSetActiveRows(
    multinomialHandle(counts.mn, NULL),
    replace(a, 2L, 0.5)
  ),
  "exactly 0 or 1"
)

# PER-FOREST masking is refused permanently and on model grounds, not for want
# of an implementation: the per-forest weight entry is the only per-forest,
# per-observation channel a caller can reach, and a softmax sampler refuses it
# with that reason.
expect_error(
  dbarts:::bartcoreSetForestWeights(multinomialHandle(counts.mn, NULL), 1L, a),
  "applies to every category"
)

rm(
  K.mn,
  labels.mn,
  counts.mn,
  counts.mn.other,
  multinomialHandle,
  train.mn.a,
  train.mn.b,
  train.mn.empty,
  n,
  x,
  y,
  y.binary,
  w,
  a,
  control,
  makeSampler,
  plain,
  draws.plain,
  levels.ord,
  z.ord,
  codes.ord,
  codes.other,
  ordinalSampler,
  counts.binary,
  counts.binary.other,
  logisticSampler,
  counts,
  counts.other,
  nbinomSampler,
  status,
  aftHandle,
  y.other,
  substitutable,
  aft.train.a,
  aft.train.b
)
