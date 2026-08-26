source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that predict fails if sampler not saved
bartFit <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = 20,
  nskip = 5,
  ntree = 5L,
  verbose = FALSE
)
# the pattern matters: bartFit is built with a namespace-qualified call, whose
# call[[1L]] is a `::` call rather than a bare symbol. An unanchored
# expect_error passed even while the guard died with "the condition has length
# > 1" instead of naming 'keeptrees'.
expect_error(
  predict(bartFit, testData$x, n.threads = 1L),
  pattern = "keeptrees"
)
expect_error(extract(bartFit, type = "trees"), pattern = "keeptrees")
# the exact refusal stem is shared by predict, extract(type = "trees") and
# plotTree, naming only the spelling this fit's own surface used
expect_error(
  predict(bartFit, testData$x),
  "requires the fit's saved trees; refit with keeptrees = TRUE",
  fixed = TRUE
)
expect_error(
  extract(bartFit, type = "trees"),
  "requires the fit's saved trees; refit with keeptrees = TRUE",
  fixed = TRUE
)
expect_error(
  plotTree(bartFit),
  "requires the fit's saved trees; refit with keeptrees = TRUE",
  fixed = TRUE
)

bart2Fit <- dbarts::bart2(
  testData$x,
  testData$y,
  n.samples = 20L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = FALSE,
  keepTrainingFits = FALSE
)
expect_error(
  predict(bart2Fit, testData$x, n.threads = 1L),
  pattern = "keepTrees"
)
expect_error(extract(bart2Fit, type = "trees"), pattern = "keepTrees")
expect_error(extract(bart2Fit, sample = "train"), pattern = "keepTrainingFits")
expect_error(
  predict(bart2Fit, testData$x),
  "requires the fit's saved trees; refit with keepTrees = TRUE",
  fixed = TRUE
)
expect_error(
  extract(bart2Fit, type = "trees"),
  "requires the fit's saved trees; refit with keepTrees = TRUE",
  fixed = TRUE
)
expect_error(
  plotTree(bart2Fit),
  "requires the fit's saved trees; refit with keepTrees = TRUE",
  fixed = TRUE
)

rbartFitNoTrees <- dbarts::rbart_vi(
  testData$y ~ testData$x,
  group.by = rep(1:2, length.out = nrow(testData$x)),
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = FALSE
)
expect_error(
  predict(
    rbartFitNoTrees,
    testData$x,
    group.by = rep(1:2, length.out = nrow(testData$x))
  ),
  "requires the fit's saved trees; refit with keepTrees = TRUE",
  fixed = TRUE
)
expect_error(
  extract(rbartFitNoTrees, type = "trees"),
  "requires the fit's saved trees; refit with keepTrees = TRUE",
  fixed = TRUE
)
expect_error(
  plotTree(rbartFitNoTrees),
  "requires the fit's saved trees; refit with keepTrees = TRUE",
  fixed = TRUE
)

rm(bartFit, bart2Fit, rbartFitNoTrees)

# predict.bart/predict.rbart did not call refuseUnusedGenericArgs: a caller
# typing the sibling family's offset formal name ('offset.test', used by
# predict.bartNegbin) instead of this fit's own 'offset' had it silently
# vanish into '...' instead of applied
bart2FitKT <- dbarts::bart2(
  testData$x,
  testData$y,
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_error(
  predict(bart2FitKT, testData$x, offset.test = rep(0.1, nrow(testData$x))),
  "'offset.test' is not used by predict on a bart fit",
  fixed = TRUE
)
# the derived surface vocabulary: a name that is a formal on a sibling
# predict method and foreign here is refused rather than silently dropped
expect_error(
  predict(bart2FitKT, testData$x, sample = "train"),
  "'sample' is not used by predict on a bart fit: the fit's stored train and test channels are extract's 'sample'",
  fixed = TRUE
)
expect_error(
  predict(bart2FitKT, testData$x, contribution = TRUE),
  "'contribution' is not used by predict on a bart fit: the per-observation contribution decomposition belongs to extract(type = \"forest\")",
  fixed = TRUE
)
expect_error(
  predict(bart2FitKT, testData$x, group.by = 1L),
  "'group.by' is not used by predict on a bart fit: 'group.by' is the grouped (rbart_vi) fit's own predict argument",
  fixed = TRUE
)
rm(bart2FitKT)

groupBy <- rep(1:2, length.out = nrow(testData$x))
rbartFitKT <- dbarts::rbart_vi(
  testData$y ~ testData$x,
  group.by = groupBy,
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_error(
  predict(
    rbartFitKT,
    testData$x,
    group.by = groupBy,
    offset.test = rep(0.1, nrow(testData$x))
  ),
  "'offset.test' is not used by predict on a rbart fit",
  fixed = TRUE
)
# group.by follows '...' now: a positional third argument binds to 'type'
# instead, and a missing one - named or positional - is refused by name
expect_error(
  predict(rbartFitKT, testData$x, groupBy),
  "'group.by' must be given by name",
  fixed = TRUE
)
expect_error(
  predict(rbartFitKT, testData$x, type = "ev"),
  "'group.by' must be given by name",
  fixed = TRUE
)
# the pre-1.0 'value' shim is deleted; the old spelling is refused by name
# rather than accepted and folded onto 'type'
expect_error(
  predict(rbartFitKT, testData$x, group.by = groupBy, value = "ev"),
  "'value' is not used by predict on a rbart fit: predict's channel argument is named 'type'",
  fixed = TRUE
)
# the derived surface vocabulary: a name foreign to predict.rbart is refused,
# 'forest'/'contribution' by the same single-forest reason ordinal/negbin/
# hurdle already carry
expect_error(
  predict(rbartFitKT, testData$x, group.by = groupBy, sample = "train"),
  "'sample' is not used by predict on a rbart fit: the fit's stored train and test channels are extract's 'sample'",
  fixed = TRUE
)
expect_error(
  predict(rbartFitKT, testData$x, group.by = groupBy, bases = list()),
  "'bases' is not used by predict on a rbart fit: only an amplitude-coupled multi-forest fit takes bases at the predicted rows",
  fixed = TRUE
)
expect_error(
  predict(rbartFitKT, testData$x, group.by = groupBy, forest = 1L),
  "'forest' is not used by predict on a rbart fit: this selects among an amplitude-coupled fit's co-fit forests; this fit has a single forest",
  fixed = TRUE
)
expect_error(
  predict(rbartFitKT, testData$x, group.by = groupBy, contribution = TRUE),
  "'contribution' is not used by predict on a rbart fit: this selects among an amplitude-coupled fit's co-fit forests; this fit has a single forest",
  fixed = TRUE
)
# 'group.by' is rbart's own formal, so it is NOT foreign here - the derived
# refusal above must not fire on the one method that declares it
expect_true(is.numeric(predict(
  rbartFitKT,
  testData$x,
  group.by = groupBy,
  type = "ev"
)))
rm(rbartFitKT, groupBy)

# extend the offset.test refusal to the four own-class fits, and add the
# 'weights' and offset-channel refusals the same signature reshape gained
n <- 40L
xSmall <- matrix(rnorm(n * 2L), n, 2L)

multinomialFitKT <- dbarts::bart2(
  xSmall,
  factor(sample(letters[1:3], n, replace = TRUE)),
  family = "multinomial",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  predict(multinomialFitKT, xSmall, offset.test = 0),
  "'offset.test' is not used by predict on a bartMultinomial fit",
  fixed = TRUE
)
expect_error(
  predict(multinomialFitKT, xSmall, weights = rep(1, n)),
  "'weights' is not used by predict on a bartMultinomial fit",
  fixed = TRUE
)
expect_error(
  predict(multinomialFitKT, xSmall),
  "requires the fit's saved trees; refit with keepTrees = TRUE",
  fixed = TRUE
)
expect_error(
  predict(multinomialFitKT, xSmall, sample = "train"),
  "'sample' is not used by predict on a bartMultinomial fit: the fit's stored train and test channels are extract's 'sample'",
  fixed = TRUE
)
expect_error(
  predict(multinomialFitKT, xSmall, bases = list()),
  "'bases' is not used by predict on a bartMultinomial fit: only an amplitude-coupled multi-forest fit takes bases at the predicted rows",
  fixed = TRUE
)
expect_error(
  predict(multinomialFitKT, xSmall, group.by = 1L),
  "'group.by' is not used by predict on a bartMultinomial fit: 'group.by' is the grouped (rbart_vi) fit's own predict argument",
  fixed = TRUE
)
rm(multinomialFitKT)

ordinalFitKT <- dbarts::bart2(
  xSmall,
  ordered(
    sample(c("lo", "mid", "hi"), n, replace = TRUE),
    levels = c("lo", "mid", "hi")
  ),
  family = "ordinal",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  predict(ordinalFitKT, xSmall, offset.test = 0),
  "this fit has no out-of-sample offset channel",
  fixed = TRUE
)
expect_error(
  predict(ordinalFitKT, xSmall, offset = 0),
  "this fit has no out-of-sample offset channel",
  fixed = TRUE
)
expect_error(
  predict(ordinalFitKT, xSmall, weights = rep(1, n)),
  "'weights' is not used by predict on a bartOrdinal fit",
  fixed = TRUE
)
expect_error(
  predict(ordinalFitKT, xSmall),
  "requires the fit's saved trees; refit with keepTrees = TRUE",
  fixed = TRUE
)
expect_error(
  predict(ordinalFitKT, xSmall, sample = "train"),
  "'sample' is not used by predict on a bartOrdinal fit: the fit's stored train and test channels are extract's 'sample'",
  fixed = TRUE
)
expect_error(
  predict(ordinalFitKT, xSmall, bases = list()),
  "'bases' is not used by predict on a bartOrdinal fit: only an amplitude-coupled multi-forest fit takes bases at the predicted rows",
  fixed = TRUE
)
expect_error(
  predict(ordinalFitKT, xSmall, group.by = 1L),
  "'group.by' is not used by predict on a bartOrdinal fit: 'group.by' is the grouped (rbart_vi) fit's own predict argument",
  fixed = TRUE
)
rm(ordinalFitKT)

negbinFitKT <- dbarts::bart2(
  xSmall,
  rpois(n, 3),
  family = "nbinom",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  predict(negbinFitKT, xSmall, offset.test = rep(0, n)),
  "'offset.test' is not used by predict on a bartNegbin fit",
  fixed = TRUE
)
expect_error(
  predict(negbinFitKT, xSmall, weights = rep(1, n)),
  "'weights' is not used by predict on a bartNegbin fit",
  fixed = TRUE
)
expect_error(
  predict(negbinFitKT, xSmall),
  "requires the fit's saved trees; refit with keepTrees = TRUE",
  fixed = TRUE
)
expect_error(
  predict(negbinFitKT, xSmall, sample = "train"),
  "'sample' is not used by predict on a bartNegbin fit: the fit's stored train and test channels are extract's 'sample'",
  fixed = TRUE
)
expect_error(
  predict(negbinFitKT, xSmall, bases = list()),
  "'bases' is not used by predict on a bartNegbin fit: only an amplitude-coupled multi-forest fit takes bases at the predicted rows",
  fixed = TRUE
)
expect_error(
  predict(negbinFitKT, xSmall, group.by = 1L),
  "'group.by' is not used by predict on a bartNegbin fit: 'group.by' is the grouped (rbart_vi) fit's own predict argument",
  fixed = TRUE
)
rm(negbinFitKT)

hurdleFitKT <- dbarts::bart2(
  xSmall,
  ifelse(runif(n) < 0.5, 0, rlnorm(n)),
  family = "hurdle.lognormal",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  predict(hurdleFitKT, xSmall, offset.test = rep(0, n)),
  "this fit has no out-of-sample offset channel",
  fixed = TRUE
)
expect_error(
  predict(hurdleFitKT, xSmall, offset = rep(0, n)),
  "this fit has no out-of-sample offset channel",
  fixed = TRUE
)
expect_error(
  predict(hurdleFitKT, xSmall, weights = rep(1, n)),
  "'weights' is not used by predict on a bartHurdle fit",
  fixed = TRUE
)
expect_error(
  predict(hurdleFitKT, xSmall),
  "requires the fit's saved trees; refit with keepTrees = TRUE",
  fixed = TRUE
)
expect_error(
  predict(hurdleFitKT, xSmall, sample = "train"),
  "'sample' is not used by predict on a bartHurdle fit: the fit's stored train and test channels are extract's 'sample'",
  fixed = TRUE
)
expect_error(
  predict(hurdleFitKT, xSmall, bases = list()),
  "'bases' is not used by predict on a bartHurdle fit: only an amplitude-coupled multi-forest fit takes bases at the predicted rows",
  fixed = TRUE
)
expect_error(
  predict(hurdleFitKT, xSmall, group.by = 1L),
  "'group.by' is not used by predict on a bartHurdle fit: 'group.by' is the grouped (rbart_vi) fit's own predict argument",
  fixed = TRUE
)
rm(hurdleFitKT, n, xSmall)

# --- extract arm: ci.level/newdata/n.threads refused on extract.bart and
# --- extract.rbart, BELOW the type == "trees" branch, and a representative
# --- foreign name on each own-class extract ---
bartFitKT <- dbarts::bart2(
  testData$x,
  testData$y,
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_error(
  extract(bartFitKT, type = "ev", ci.level = 0.9),
  "'ci.level' is not used by extract on a bart fit: extract returns the draws that fitted() and predict() take a band over",
  fixed = TRUE
)
expect_error(
  extract(bartFitKT, type = "ev", newdata = testData$x),
  "'newdata' is not used by extract on a bart fit: predict(object, newdata) is the read at new rows",
  fixed = TRUE
)
expect_error(
  extract(bartFitKT, type = "ev", n.threads = 1L),
  "'n.threads' is not used by extract on a bart fit: extract reads stored channels and replays nothing",
  fixed = TRUE
)
# the acceptance half: extract(type = "trees", newdata = ) still forwards to
# getTrees below the new refusal call, so 'newdata' keeps working there - the
# root of every tree holds every row, so the largest 'n' is the supplied
# row count rather than the training one
treesAtNewdata <- extract(
  bartFitKT,
  type = "trees",
  newdata = testData$x[1:5, , drop = FALSE]
)
expect_equal(max(treesAtNewdata$n), 5L)
rm(bartFitKT)

groupByKT <- rep(1:2, length.out = nrow(testData$x))
rbartFitKT2 <- dbarts::rbart_vi(
  testData$y ~ testData$x,
  group.by = groupByKT,
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_error(
  extract(rbartFitKT2, type = "ev", ci.level = 0.9),
  "'ci.level' is not used by extract on a rbart fit: extract returns the draws that fitted() and predict() take a band over",
  fixed = TRUE
)
expect_error(
  extract(rbartFitKT2, type = "ev", newdata = testData$x),
  "'newdata' is not used by extract on a rbart fit: predict(object, newdata) is the read at new rows",
  fixed = TRUE
)
expect_error(
  extract(rbartFitKT2, type = "ev", n.threads = 1L),
  "'n.threads' is not used by extract on a rbart fit: extract reads stored channels and replays nothing",
  fixed = TRUE
)
rm(rbartFitKT2, groupByKT)

n2 <- 40L
xSmall2 <- matrix(rnorm(n2 * 2L), n2, 2L)
multinomialFitKT2 <- dbarts::bart2(
  xSmall2,
  factor(sample(letters[1:3], n2, replace = TRUE)),
  family = "multinomial",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  extract(multinomialFitKT2, type = "ev", ci.level = 0.9),
  "'ci.level' is not used by extract on a bartMultinomial fit: extract returns the draws that fitted() and predict() take a band over",
  fixed = TRUE
)
rm(multinomialFitKT2)

ordinalFitKT2 <- dbarts::bart2(
  xSmall2,
  ordered(
    sample(c("lo", "mid", "hi"), n2, replace = TRUE),
    levels = c("lo", "mid", "hi")
  ),
  family = "ordinal",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  extract(ordinalFitKT2, type = "ev", ci.level = 0.9),
  "'ci.level' is not used by extract on a bartOrdinal fit: extract returns the draws that fitted() and predict() take a band over",
  fixed = TRUE
)
rm(ordinalFitKT2)

negbinFitKT2 <- dbarts::bart2(
  xSmall2,
  rpois(n2, 3),
  family = "nbinom",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  extract(negbinFitKT2, type = "ev", ci.level = 0.9),
  "'ci.level' is not used by extract on a bartNegbin fit: extract returns the draws that fitted() and predict() take a band over",
  fixed = TRUE
)
rm(negbinFitKT2)

hurdleFitKT2 <- dbarts::bart2(
  xSmall2,
  ifelse(runif(n2) < 0.5, 0, rlnorm(n2)),
  family = "hurdle.lognormal",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  extract(hurdleFitKT2, type = "ev", ci.level = 0.9),
  "'ci.level' is not used by extract on a bartHurdle fit: extract returns the draws that fitted() and predict() take a band over",
  fixed = TRUE
)
rm(hurdleFitKT2, n2, xSmall2)

# extract.dbartsSampler's own catch-all: its '...' forwards nowhere, so any
# supplied name is refused rather than only a listed few
samplerKT <- dbarts::dbarts(
  testData$y ~ testData$x,
  control = dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 5L
  )
)
expect_error(
  extract(samplerKT, sample = "train"),
  "'sample' is not used by extract on a dbartsSampler: this method returns the sampler's coded predictor matrix",
  fixed = TRUE
)
expect_error(
  extract(samplerKT, combineChains = FALSE),
  "'combineChains' is not used by extract on a dbartsSampler: this method returns the sampler's coded predictor matrix",
  fixed = TRUE
)
rm(samplerKT)

# --- fitted/residuals arm: combineChains refused on all six fitted,
# --- ci.level refused on all six residuals (section 7), sample refused on
# --- the three residuals methods that did not already refuse it, and
# --- n.threads on one fitted and one residuals ---
n3 <- 40L
xSmall3 <- matrix(rnorm(n3 * 2L), n3, 2L)
bartFitPlain <- dbarts::bart2(
  testData$x,
  testData$y,
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  fitted(bartFitPlain, combineChains = FALSE),
  "'combineChains' is not used by fitted on a bart fit: the per-chain draws are extract(object, combineChains = FALSE)",
  fixed = TRUE
)
expect_error(
  fitted(bartFitPlain, n.threads = 1L),
  "'n.threads' is not used by fitted on a bart fit: fitted summarizes stored channels and replays nothing",
  fixed = TRUE
)
expect_error(
  residuals(bartFitPlain, ci.level = 0.9),
  "'ci.level' is not used by residuals on a bart fit: residuals are the observed response minus the posterior-mean fit",
  fixed = TRUE
)
expect_error(
  residuals(bartFitPlain, n.threads = 1L),
  "'n.threads' is not used by residuals on a bart fit: residuals summarize stored channels and replay nothing",
  fixed = TRUE
)
rm(bartFitPlain)

groupByPlain <- rep(1:2, length.out = nrow(testData$x))
rbartFitPlain <- dbarts::rbart_vi(
  testData$y ~ testData$x,
  group.by = groupByPlain,
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = FALSE
)
expect_error(
  fitted(rbartFitPlain, combineChains = FALSE),
  "'combineChains' is not used by fitted on a rbart fit: the per-chain draws are extract(object, combineChains = FALSE)",
  fixed = TRUE
)
expect_error(
  residuals(rbartFitPlain, ci.level = 0.9),
  "'ci.level' is not used by residuals on a rbart fit: residuals are the observed response minus the posterior-mean fit",
  fixed = TRUE
)
rm(rbartFitPlain, groupByPlain)

multinomialFitPlain <- dbarts::bart2(
  xSmall3,
  factor(sample(letters[1:3], n3, replace = TRUE)),
  family = "multinomial",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  fitted(multinomialFitPlain, combineChains = FALSE),
  "'combineChains' is not used by fitted on a bartMultinomial fit: the per-chain draws are extract(object, combineChains = FALSE)",
  fixed = TRUE
)
expect_error(
  residuals(multinomialFitPlain, ci.level = 0.9),
  "'ci.level' is not used by residuals on a bartMultinomial fit: residuals are the observed response minus the posterior-mean fit",
  fixed = TRUE
)
expect_error(
  residuals(multinomialFitPlain, sample = "train"),
  "'sample' is not used by residuals; residuals are always against the training response",
  fixed = TRUE
)
rm(multinomialFitPlain)

ordinalFitPlain <- dbarts::bart2(
  xSmall3,
  ordered(
    sample(c("lo", "mid", "hi"), n3, replace = TRUE),
    levels = c("lo", "mid", "hi")
  ),
  family = "ordinal",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  fitted(ordinalFitPlain, combineChains = FALSE),
  "'combineChains' is not used by fitted on a bartOrdinal fit: the per-chain draws are extract(object, combineChains = FALSE)",
  fixed = TRUE
)
expect_error(
  residuals(ordinalFitPlain, ci.level = 0.9),
  "'ci.level' is not used by residuals on a bartOrdinal fit: residuals are the observed response minus the posterior-mean fit",
  fixed = TRUE
)
expect_error(
  residuals(ordinalFitPlain, sample = "train"),
  "'sample' is not used by residuals; residuals are always against the training response",
  fixed = TRUE
)
rm(ordinalFitPlain)

negbinFitPlain <- dbarts::bart2(
  xSmall3,
  rpois(n3, 3),
  family = "nbinom",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  fitted(negbinFitPlain, combineChains = FALSE),
  "'combineChains' is not used by fitted on a bartNegbin fit: the per-chain draws are extract(object, combineChains = FALSE)",
  fixed = TRUE
)
expect_error(
  residuals(negbinFitPlain, ci.level = 0.9),
  "'ci.level' is not used by residuals on a bartNegbin fit: residuals are the observed response minus the posterior-mean fit",
  fixed = TRUE
)
expect_error(
  residuals(negbinFitPlain, sample = "train"),
  "'sample' is not used by residuals; residuals are always against the training response",
  fixed = TRUE
)
rm(negbinFitPlain)

hurdleFitPlain <- dbarts::bart2(
  xSmall3,
  ifelse(runif(n3) < 0.5, 0, rlnorm(n3)),
  family = "hurdle.lognormal",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  fitted(hurdleFitPlain, combineChains = FALSE),
  "'combineChains' is not used by fitted on a bartHurdle fit: the per-chain draws are extract(object, combineChains = FALSE)",
  fixed = TRUE
)
expect_error(
  residuals(hurdleFitPlain, ci.level = 0.9),
  "'ci.level' is not used by residuals on a bartHurdle fit: residuals are the observed response minus the posterior-mean fit",
  fixed = TRUE
)
rm(hurdleFitPlain, n3, xSmall3)

rm(testData)
