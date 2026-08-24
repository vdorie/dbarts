# The PUBLIC multinomial surface (docs/design/multinomial.md):
# dbarts(family = "multinomial") builds a real, runnable
# K-forest softmax sampler; the response and both category offsets ride the
# data object, so a re-created or reloaded sampler carries them; and every
# channel the softmax gives no meaning to is refused by a message naming the
# capability rather than a C entry point.
#
# Fixtures use genuinely varying predictor columns throughout: a constant
# column is NOT inert on this engine (it shifts the tree prior), so it can
# never stand in for "a column that does not split".

set.seed(303)
n <- 120L
p <- 4L
K <- 3L
x <- matrix(rnorm(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
nTest <- 9L
x.test <- matrix(rnorm(nTest * p), nTest, p)
colnames(x.test) <- colnames(x)

eta <- cbind(0.9 * x[, 1L], -0.6 * x[, 2L], 0.4 * x[, 3L])
prob <- exp(eta) / rowSums(exp(eta))
labels <- factor(
  apply(prob, 1L, function(row) sample.int(K, 1L, prob = row)),
  levels = seq_len(K),
  labels = c("low", "mid", "high")
)
oneHot <- matrix(0L, n, K, dimnames = list(NULL, levels(labels)))
oneHot[cbind(seq_len(n), as.integer(labels))] <- 1L

control <- dbartsControl(
  n.chains = 1L,
  n.trees = 25L,
  n.samples = 8L,
  n.burn = 0L,
  updateState = FALSE,
  verbose = FALSE
)

# --- creation from a factor response ----------------------------------------
# a K >= 3 factor is the single-trial special case: one-hot counts, every trial
# 1, and the levels carried on the count matrix's own column names
sampler <- dbarts(
  x,
  labels,
  test = x.test,
  family = "multinomial",
  control = control
)
expect_inherits(sampler, "dbartsSampler")
expect_identical(sampler$model@family, "multinomial")
expect_false(sampler$control@binary)
expect_identical(sampler$data@counts, oneHot)
expect_identical(colnames(sampler$data@counts), levels(labels))
# G1a: y is the TRIALS vector, so every length(data@y) reader stays meaningful
expect_identical(sampler$data@y, rep(1.0, n))
expect_identical(length(sampler$data@y), n)
expect_null(sampler$data@offset.category)
expect_null(sampler$data@offset.category.test)

# --- family gate 1: the degeneracy warning skips a count-carrying object -----
# a single-trial response has an identically constant trials vector, which the
# precision-degeneracy ratio reads as a degenerate response on every such fit
degenerateWarnings <- character(0L)
withCallingHandlers(
  dbartsData(x, counts = oneHot),
  warning = function(w) {
    degenerateWarnings <<- c(degenerateWarnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
expect_identical(length(degenerateWarnings), 0L)
# and the check still bites on the same vector WITHOUT the counts, which is
# what makes the bypass a bypass rather than a deletion
plainWarnings <- character(0L)
withCallingHandlers(
  dbartsData(x, rep(1.0, n)),
  warning = function(w) {
    plainWarnings <<- c(plainWarnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
expect_identical(length(plainWarnings), 1L)
expect_true(grepl("indistinguishable", plainWarnings[1L], fixed = TRUE))

# --- family gate 2: no degenerate starting sigma ------------------------------
# the plain fit of the same constant response installs sigma ~ 1e-16; the
# multinomial one is a fixed-unit-scale family and installs none at all
expect_true(is.na(sampler$data@sigma))
expect_false(is.na(suppressWarnings(dbarts(x, rep(1.0, n)))$data@sigma))

# --- it runs, and reports K softmax probabilities ---------------------------
samples <- sampler$run(6L, 8L)
expect_identical(dim(samples$train), c(n, K, 8L))
expect_equal(apply(samples$train[,, 1L], 1L, sum), rep(1.0, n))
expect_identical(dim(samples$test), c(nTest, K, 8L))
expect_true(all(samples$train > 0 & samples$train < 1))

# --- a grouped-count response ------------------------------------------------
set.seed(11)
grouped <- matrix(rpois(n * K, 2), n, K)
grouped[rowSums(grouped) == 0L, 1L] <- 1L
groupedSampler <- dbarts(x, grouped, family = "multinomial", control = control)
expect_identical(groupedSampler$data@y, as.double(rowSums(grouped)))
expect_identical(storage.mode(groupedSampler$data@counts), "integer")
expect_identical(dim(groupedSampler$run(4L, 8L)$train), c(n, K, 8L))

# --- dbartsSpec reaches the same family through the data object -------------
# a counts-carrying data object IS the declaration, so "auto" resolves to it
specData <- dbartsData(x, counts = oneHot)
spec <- dbartsSpec(specData, control = control)
expect_identical(spec$family, "multinomial")
expect_identical(
  dbartsSpec(specData, control = control, family = "multinomial")$family,
  "multinomial"
)

# --- the response slots are constrained --------------------------------------
expect_error(
  dbartsData(x, counts = oneHot[, 1L, drop = FALSE]),
  "at least two categories"
)
expect_error(dbartsData(x, counts = oneHot - 1L), "non-negative")
expect_error(dbartsData(x, counts = oneHot + 0.5), "whole numbers")
expect_error(dbartsData(x, counts = matrix(0L, n, K)), "at least one trial")
expect_error(
  dbartsData(x, counts = oneHot[-1L, , drop = FALSE]),
  "same number of rows"
)
expect_error(
  dbartsData(x, counts = oneHot, offset.category = matrix(0.0, n, K + 1L)),
  "dimensions of 'counts'"
)
expect_error(
  dbartsData(x, counts = oneHot, offset.category = matrix(Inf, n, K)),
  "finite"
)
expect_error(
  dbartsData(x, offset.category = matrix(0.0, n, K)),
  "supply one, or use 'offset'"
)
expect_error(
  dbartsData(x, counts = oneHot, offset.category.test = matrix(0.0, nTest, K)),
  "'x.test' is null"
)
# a counts-carrying object fits only the family whose response it is
expect_error(
  dbartsSpec(specData, control = control, family = "gaussian"),
  "only family \"multinomial\" fits"
)
expect_error(
  dbarts(
    x,
    labels[seq_len(n)],
    family = "multinomial",
    control = control,
    weights = runif(n, 0.5, 1.5)
  ),
  "row-wise replication"
)
# and the token needs a response to name
expect_error(
  dbartsSpec(
    dbartsData(x, rnorm(n)),
    control = control,
    family = "multinomial"
  ),
  "needs an n x K count-matrix response"
)
# a formula LHS carries no count matrix
formulaFrame <- data.frame(y = labels, x1 = x[, 1L], x2 = x[, 2L])
expect_error(
  dbarts(y ~ x1 + x2, formulaFrame, family = "multinomial", control = control),
  "matrix interface"
)

# --- compositions the K-forest factory does not build are named -------------
expect_error(
  dbarts(
    x,
    labels,
    family = "multinomial",
    control = control,
    monotone = c(1, 0, 0, 0)
  ),
  "does not support 'monotone'"
)
expect_error(
  dbarts(
    x,
    labels,
    family = "multinomial",
    control = control,
    tree.prior = dart
  ),
  "does not support a DART tree prior"
)
expect_error(
  dbarts(
    x,
    labels,
    family = "multinomial",
    control = control,
    node.prior = normal(chi(1.25))
  ),
  "does not support a 'k' hyperprior"
)
expect_error(
  dbarts(x, labels, family = "multinomial", control = control, variance = TRUE),
  "requires family = \"gaussian\""
)

# --- the refusal matrix ------------------------------------------------------
# every refusal names the CAPABILITY and, where one exists, the channel that
# serves the caller instead - never a C entry point, which R cannot call
expect_error(sampler$setResponse(rep(1.0, n)), "\\$setCounts")
expect_error(sampler$setResponse(rep(1.0, n)), "n x K count matrix")
expect_error(sampler$setOffset(rep(0.1, n)), "\\$setCategoryOffset")
expect_error(sampler$setOffset(rep(0.1, n)), "null direction")
expect_error(sampler$setWeights(rep(1.0, n)), "row-wise replication")
expect_error(sampler$setSigma(1.0), "no residual scale")
expect_error(sampler$setData(sampler$data), "fix their data at creation")
expect_error(sampler$setModel(sampler$model), "calibrated at creation")
expect_error(
  sampler$setCalibration(prior.scale = 1.0),
  "softmax calibration map"
)
expect_error(sampler$setForestWeights(1L, rep(1.0, n)), "log-sum-exp")
expect_error(sampler$setForestBasis(1L, rep(1.0, n)), "carry no amplitudes")
expect_error(sampler$getFitsWithoutOffset(), "softmax probabilities")
# the four channels whose refusal is R-canonical - raised before the .Call, so
# the caller is told which R method serves them - name no C entry point. The
# rest are the bridge's own backstops and keep its caller-name prefix, which is
# the shared-guard convention (docs/design/error-style.md R7/R8), not a leak.
multinomialRefusals <- vapply(
  list(
    function() sampler$setResponse(rep(1.0, n)),
    function() sampler$setOffset(rep(0.1, n)),
    function() sampler$setWeights(rep(1.0, n)),
    function() sampler$setSigma(1.0)
  ),
  function(call) {
    tryCatch(
      {
        call()
        ""
      },
      error = conditionMessage
    )
  },
  ""
)
expect_false(any(grepl("bartcore_", multinomialRefusals, fixed = TRUE)))
# a written decline, not an omission: the PG omegas are interleaved
# one-vs-rest augmentation variables refreshed against a moving margin, so
# nothing marginally meaningful is read here
expect_null(sampler$getLatents())

# --- the channels that stay open --------------------------------------------
expect_silent(sampler$setPredictor(x))
expect_silent(sampler$setPredictor(x[, 2L] + 0.05, 2L))
expect_silent(sampler$setCutPoints(sort(unique(x[, 3L]))[c(10L, 40L, 70L)], 3L))
expect_silent(sampler$setTestPredictor(x.test))
expect_silent(sampler$setActiveRows(rep(1.0, n)))
expect_silent(sampler$setActiveRows(NULL))
# a flat test offset is refused for the train-side reason, naming the matrix
expect_error(sampler$setTestOffset(rep(0.1, nTest)), "category test offset")

# --- the three response channels --------------------------------------------
set.seed(41)
swapped <- matrix(0L, n, K, dimnames = list(NULL, levels(labels)))
swapped[cbind(seq_len(n), sample.int(K, n, replace = TRUE))] <- 1L
sampler$setCounts(swapped)
expect_identical(sampler$data@counts, swapped)
expect_identical(sampler$data@y, as.double(rowSums(swapped)))
expect_identical(dim(sampler$run(0L, 4L)$train), c(n, K, 4L))
# n and K are fixed at creation
expect_error(sampler$setCounts(swapped[, seq_len(2L)]), "3 categories")
expect_error(
  sampler$setCounts(swapped[-1L, , drop = FALSE]),
  "same number of rows"
)
# a rejected matrix leaves the sampler untouched
expect_error(sampler$setCounts(swapped - 1L), "non-negative")
expect_identical(sampler$data@counts, swapped)

categoryOffset <- matrix(rnorm(n * K, 0, 0.3), n, K)
sampler$setCategoryOffset(categoryOffset)
expect_identical(sampler$data@offset.category, categoryOffset)
categoryTestOffset <- matrix(rnorm(nTest * K, 0, 0.3), nTest, K)
sampler$setCategoryTestOffset(categoryTestOffset)
expect_identical(sampler$data@offset.category.test, categoryTestOffset)
expect_error(sampler$setCategoryOffset(matrix(0.0, n, K + 1L)), "3 categories")
expect_error(sampler$setCategoryOffset(matrix(NA_real_, n, K)), "finite")
sampler$setCategoryOffset(NULL)
expect_null(sampler$data@offset.category)
sampler$setCategoryTestOffset(NULL)
expect_null(sampler$data@offset.category.test)

# the three channels are the softmax's own, and are refused by name elsewhere
plain <- dbarts(x, rnorm(n), control = control)
expect_error(plain$setCounts(oneHot), "carries no count response")
expect_error(
  plain$setCategoryOffset(matrix(0.0, n, K)),
  "carries no count response"
)
expect_error(
  plain$setCategoryTestOffset(matrix(0.0, nTest, K)),
  "carries no count response"
)

# --- the K-matrix predict arm ------------------------------------------------
predictions <- sampler$predict(x.test)
expect_identical(dim(predictions)[1:2], c(nTest, K))
expect_equal(apply(predictions[,, 1L], 1L, sum), rep(1.0, nTest))
# an all-zero per-category matrix is the offset-free surface, exactly
expect_identical(sampler$predict(x.test, matrix(0.0, nTest, K)), predictions)
# a scalar broadcasts across the whole matrix, which the softmax absorbs
expect_equal(sampler$predict(x.test, 0.25), predictions)
# a flat vector is the softmax's null direction and stays refused
expect_error(
  sampler$predict(x.test, rep(0.25, nTest)),
  "must be a per-category matrix"
)
expect_error(sampler$predict(x.test, rep(0.25, nTest)), "3 categories")
expect_error(
  sampler$predict(x.test, matrix(0.0, nTest, K + 1L)),
  "3 categories"
)

# --- the single-location entry points refuse a count-carrying data object ----
# every channel bart2() and xbart() report - the chain/sample reshaping, the
# sigma vector, the losses - is written against ONE location per observation,
# and a K-location fit would be reshaped and scored without a word. The refusal
# is at the entry, not at a family branch: family = "auto" resolves counts to
# multinomial during spec resolution, so no branch keyed on the token sees it.
countsData <- dbartsData(x, counts = oneHot)
expect_error(
  bart2(
    countsData,
    n.samples = 4L,
    n.burn = 2L,
    n.chains = 1L,
    n.trees = 10L,
    verbose = FALSE
  ),
  "does not support a data object carrying an n x K count matrix"
)
expect_error(
  bart2(
    countsData,
    n.samples = 4L,
    n.burn = 2L,
    n.chains = 1L,
    n.trees = 10L,
    verbose = FALSE
  ),
  "dbarts\\(family = \"multinomial\"\\)"
)
# and under an explicit token, which reaches a different refusal but must also
# never fit
expect_error(
  bart2(
    countsData,
    family = "multinomial",
    n.samples = 4L,
    n.burn = 2L,
    n.chains = 1L,
    n.trees = 10L,
    verbose = FALSE
  ),
  "count matrix"
)
expect_error(
  xbart(
    countsData,
    n.samples = 4L,
    n.burn = c(2L, 1L, 1L),
    n.reps = 1L,
    n.trees = 10L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "does not support a data object carrying an n x K count matrix"
)
# bart(), the frozen shim, reaches the same dispatch through dbartsData's
# passthrough branch and is gated at its own entry too
expect_error(
  bart(countsData, ndpost = 4L, nskip = 2L, ntree = 10L, verbose = FALSE),
  "does not support a data object carrying an n x K count matrix"
)
# rbart_vi needs no gate of its own: grouped random effects and the softmax's
# own augmentation are not shown to interleave, so spec resolution refuses the
# composition by name before any packaging sees a K-location fit
expect_error(
  rbart_vi(
    countsData,
    group.by = rep_len(1L:3L, n),
    n.samples = 20L,
    n.burn = 10L,
    n.chains = 1L,
    n.trees = 10L,
    verbose = FALSE
  ),
  "grouped random effects"
)
# an ordinary data object is untouched by the gate
expect_silent(invisible(bart2(
  x,
  rnorm(n),
  n.samples = 4L,
  n.burn = 2L,
  n.chains = 1L,
  n.trees = 10L,
  verbose = FALSE
)))

# --- 'counts' is a response, and a formula already names one -----------------
# refused at BOTH layers: dbarts() at its own entry, and dbartsData(), which
# owns response ingestion. The model frame has already applied 'subset' by the
# time a formula response exists, so a count matrix at the caller's full row
# count could not be aligned to the kept rows the way 'weights' is
countsFrame <- data.frame(y = rnorm(n), x1 = x[, 1L], x2 = x[, 2L])
expect_error(
  dbartsData(y ~ x1 + x2, countsFrame, counts = oneHot),
  "cannot be given with a formula"
)
expect_error(
  dbartsData(
    y ~ x1 + x2,
    countsFrame,
    counts = oneHot,
    offset.category = matrix(0.0, n, K)
  ),
  "cannot be given with a formula"
)
# the matrix paths carry 'subset' through to the counts, as documented
subsetRows <- seq_len(40L)
subsetData <- dbartsData(x, counts = oneHot, subset = subsetRows)
expect_identical(nrow(subsetData@counts), 40L)
expect_identical(subsetData@counts, oneHot[subsetRows, , drop = FALSE])
expect_identical(subsetData@y, as.double(rowSums(oneHot[subsetRows, ])))
subsetOffset <- matrix(rnorm(n * K), n, K)
subsetOffsetData <- dbartsData(
  x,
  counts = oneHot,
  subset = subsetRows,
  offset.category = subsetOffset
)
expect_identical(
  subsetOffsetData@offset.category,
  subsetOffset[subsetRows, , drop = FALSE]
)

# --- a lone NA on the FLAT predict path still reads as "no offset" -----------
# the multinomial arm takes the per-category matrix, and adding it must not
# move what the single-location path does with the scalar NA its older
# spelling defaulted to
oneRow <- x[1L, , drop = FALSE]
expect_identical(plain$predict(oneRow, NA_real_), plain$predict(oneRow))
# the collapse is the one-row case only, then as now: a scalar NA against more
# rows is recycled to a full NA vector before the test can see it, and stays
# one - restoring the collapse must not quietly widen it either
expect_true(anyNA(plain$predict(x.test, NA_real_)))

# --- the migration guard: an object serialized before the slots existed ------
# a dbartsData written under the pre-S2 class definition carries no 'counts',
# 'offset.category' or 'offset.category.test' attribute at all, which is what
# removing them here reproduces; every internal read must answer NULL rather
# than raising "no slot of name"
oldData <- plain$data
attributes(oldData)$counts <- NULL
attributes(oldData)$offset.category <- NULL
attributes(oldData)$offset.category.test <- NULL
expect_false(methods::.hasSlot(oldData, "counts"))
expect_null(dbarts:::dataCounts(oldData))
expect_null(dbarts:::dataSlotOrNULL(oldData, "offset.category"))
expect_null(dbarts:::dataSlotOrNULL(oldData, "offset.category.test"))
oldSampler <- dbarts:::dbartsSampler$new(control, plain$model, oldData)
expect_false(dbarts:::samplerCarriesCounts(oldSampler))
expect_identical(dim(oldSampler$run(4L, 8L)$train), c(n, 8L))
oldSampler$setResponse(rnorm(n))
expect_identical(length(oldSampler$data@y), n)
expect_error(oldSampler$setCounts(oneHot), "carries no count response")

# --- state round trip: STRUCTURAL, not bitwise -------------------------------
# omega is a per-sweep latent redrawn against whatever margins the restored
# forests present, so the restored chain reconstructs the MODEL - the trees and
# their leaf parameters - rather than the draw stream (multinomial.md)
roundTrip <- dbarts(
  x,
  labels,
  test = x.test,
  family = "multinomial",
  control = control
)
invisible(roundTrip$run(6L, 8L))
roundTrip$setCategoryOffset(categoryOffset)
roundTrip$setCategoryTestOffset(categoryTestOffset)
roundTrip$storeState()
savedState <- roundTrip$state
savedPredictions <- roundTrip$predict(x.test, matrix(0.0, nTest, K))

serialized <- tempfile(fileext = ".rds")
saveRDS(roundTrip, serialized)
restored <- readRDS(serialized)
# the three slots are DATA and ride the object, so re-creation reads them where
# creation did - there is no reapply step to forget
expect_identical(restored$data@counts, roundTrip$data@counts)
expect_identical(restored$data@offset.category, categoryOffset)
expect_identical(restored$data@offset.category.test, categoryTestOffset)
expect_identical(restored$data@y, roundTrip$data@y)
# and the re-created engine replays the same saved forests
expect_identical(
  restored$predict(x.test, matrix(0.0, nTest, K)),
  savedPredictions
)
restored$storeState()
source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)
statesAgree(restored$state, savedState)
# a resident category offset makes an unnamed predict a refusal rather than a
# silent report of the offset-free surface
expect_error(restored$predict(x.test), "cannot be inferred")
unlink(serialized)

# --- copy() carries the response ---------------------------------------------
duplicate <- roundTrip$copy()
expect_identical(duplicate$data@counts, roundTrip$data@counts)
expect_identical(duplicate$data@offset.category, categoryOffset)
expect_error(duplicate$setResponse(rep(1.0, n)), "n x K count matrix")
expect_identical(dim(duplicate$run(0L, 4L)$train), c(n, K, 4L))

rm(sampler, groupedSampler, roundTrip, restored, duplicate, plain, oldSampler)
