# A warm start from a donor's SAVED sample installs that sample's own scale
# surface, not the donor's final one: the mean and variance saved buffers are
# index-aligned, so sample k names one (mean forest, scale surface) pair
# (docs/design/heteroscedastic.md). There is no behavioral readout of a LIVE
# variance surface - predict replays the saved slots, which a warm start leaves
# untouched - so the pin is on the state's variance blocks, sliced out of the
# donor's own saved block by hand rather than by the code under test.

set.seed(37, sample.kind = "Rejection")

n <- 300L
xHet <- cbind(runif(n), runif(n))
sHet <- ifelse(xHet[, 1L] < 0.5, 0.3, 1.3)
yHet <- 2 * xHet[, 2L] + sHet * rnorm(n)

nSamples <- 4L
nVarianceTrees <- 5L
hetControl <- function(keepTrees) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    n.samples = nSamples,
    keepTrees = keepTrees,
    updateState = FALSE
  )
}
hetSampler <- function(keepTrees, x = xHet) {
  dbarts(
    x,
    yHet,
    control = hetControl(keepTrees),
    variance = varianceForest(n.trees = nVarianceTrees)
  )
}

donor <- hetSampler(TRUE)
donorSamples <- donor$run(60L, nSamples)
donor$storeState()
donorState <- donor$state

# the flat blocks are tree-major and the saved buffer is slot-major, so slot
# `slot`'s trees are entries slot * nVarianceTrees + 1:nVarianceTrees; `sizes`
# gives each tree's node count and `values` holds eight bytes per node
savedSlice <- function(state, slot, numTrees) {
  sizes <- state[[1L]][["variance.saved.sizes"]]
  trees <- slot * numTrees + seq_len(numTrees)
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  nodes <- unlist(lapply(trees, function(t) starts[t]:ends[t]))
  bytes <- unlist(lapply(nodes, function(i) ((i - 1L) * 8L + 1L):(i * 8L)))
  list(
    vars = state[[1L]][["variance.saved.vars"]][nodes],
    values = state[[1L]][["variance.saved.values"]][bytes],
    sizes = sizes[trees],
    flags = state[[1L]][["variance.saved.flags"]][nodes]
  )
}
liveVariance <- function(state) {
  list(
    vars = state[[1L]][["variance.vars"]],
    values = state[[1L]][["variance.values"]],
    sizes = state[[1L]][["variance.sizes"]],
    flags = state[[1L]][["variance.flags"]]
  )
}

## ---- slot k's mean forest arrives with slot k's own scale surface ----
# sample 3 of 4, so the pin separates the named slot from the donor's final
# sweep, whose saved trees are its live ones
k <- 3L
expected <- savedSlice(donorState, k - 1L, nVarianceTrees)
# non-vacuity: the named sample's surface is not the donor's current one, and
# the state does carry the saved channel at all
expect_false(is.null(donorState[[1L]][["variance.saved.vars"]]))
expect_false(identical(expected, liveVariance(donorState)))

dest <- hetSampler(FALSE)
dest$installTrees(donor, samples = k)
expect_identical(liveVariance(dest$state), expected)

# the two halves come from ONE sample: the mean forest reproduces the donor's
# fit at the same k
expect_equal(as.vector(dest$predict(xHet)), donorSamples$train[, k])

## ---- a live-sourced start is unchanged ----
liveDonor <- hetSampler(FALSE)
invisible(liveDonor$run(60L, nSamples))
liveDonor$storeState()
liveDest <- hetSampler(FALSE)
liveDest$installTrees(liveDonor)
expect_identical(liveVariance(liveDest$state), liveVariance(liveDonor$state))

## ---- the pooled-categorical mask channel rides the same slice ----
# a variance tree splitting on a >63-level column keeps its rule's words in a
# side channel, so the mask block must be sliced with the trees. Pinned at slot
# 0, where the expected value is a PREFIX of the saved mask block: a general
# slot would force this test to re-derive the engine's per-tree mask-word count
# (pooled split nodes x words per column), duplicating engine logic in an
# oracle.
set.seed(23, sample.kind = "Rejection")
nWide <- 400L
wideData <- data.frame(
  g = factor(sample(80L, nWide, replace = TRUE)),
  z = runif(nWide)
)
# the SCALE rides the wide factor, so the variance trees split on it
wideSd <- ifelse(as.integer(wideData$g) < 40L, 0.25, 1.4)
wideData$y <- wideData$z + wideSd * rnorm(nWide)
nWideSamples <- 12L
wideSampler <- function(keepTrees) {
  dbarts(
    y ~ g + z,
    wideData,
    control = dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 10L,
      n.samples = nWideSamples,
      keepTrees = keepTrees,
      updateState = FALSE
    ),
    variance = varianceForest(n.trees = nVarianceTrees)
  )
}
wideDonor <- wideSampler(TRUE)
invisible(wideDonor$run(60L, nWideSamples))
wideDonor$storeState()
wideDest <- wideSampler(FALSE)
wideDest$installTrees(wideDonor, samples = 1L)
wideMasks <- wideDest$state[[1L]][["variance.masks"]]
# non-vacuity: some variance tree does split on the wide column, and slot 0's
# words are not the donor's live ones
expect_true(length(wideMasks) > 0L)
expect_false(identical(wideMasks, wideDonor$state[[1L]][["variance.masks"]]))
expect_identical(
  wideMasks,
  wideDonor$state[[1L]][["variance.saved.masks"]][seq_along(wideMasks)]
)

## ---- refusals, each leaving the destination usable ----
# a donor state carrying no saved variance surface at all: the mean forests name
# a slot the scale surface cannot supply, so the install is refused rather than
# falling back to the live surface
strippedState <- donorState
for (block in c("vars", "values", "sizes", "flags")) {
  strippedState[[1L]][[paste0("variance.saved.", block)]] <- NULL
}
expect_error(
  dest$installTrees(strippedState, samples = 2L),
  "saved variance buffer"
)

# a state whose LIVE variance block is shorter than the saved buffer's stride:
# a bound-only check would accept it at a high slot and slice ACROSS two
# sweeps' trees, handing the install a correctly sized mixture
splicedState <- donorState
liveSizes <- splicedState[[1L]][["variance.sizes"]]
keptTrees <- seq_len(nVarianceTrees - 2L)
liveEnds <- cumsum(liveSizes)
keptNodes <- seq_len(liveEnds[length(keptTrees)])
splicedState[[1L]][["variance.sizes"]] <- liveSizes[keptTrees]
splicedState[[1L]][["variance.vars"]] <-
  splicedState[[1L]][["variance.vars"]][keptNodes]
splicedState[[1L]][["variance.values"]] <-
  splicedState[[1L]][["variance.values"]][seq_len(8L * length(keptNodes))]
splicedState[[1L]][["variance.flags"]] <-
  splicedState[[1L]][["variance.flags"]][keptNodes]
expect_error(
  dest$installTrees(splicedState, samples = nSamples),
  "saved variance buffer"
)

# the destination is still the sampler it was before either refusal
expect_identical(liveVariance(dest$state), expected)
expect_true(all(is.finite(dest$run(0L, 3L)$sigma)))

## ---- a saved slot the destination's rows no longer support ----
# saved variance trees are exempt from the occupancy pass that live ones take,
# so a donor that kept sweeps and then had its rows moved can hold a slot whose
# region is empty against the destination. That is refused by the variance
# install, not silently installed: an unoccupied scale leaf reports a scale the
# data never supported.
strandDest <- hetSampler(FALSE)
strandDest$setPredictor(rep(0.97, n), 1L, forceUpdate = TRUE)
expect_error(
  strandDest$installTrees(donor, samples = 1L),
  "variance surface is incompatible"
)
expect_true(all(is.finite(strandDest$run(0L, 3L)$sigma)))

## ---- setState is held to the column mask installTrees is ----
# A `variance = <subset>` forest may split only on its own columns:
# splitVariableLogProbability prices a rule against an availability menu that
# drops the rest, so a tree carrying a forbidden split is mis-scored for as
# long as it lives. setState used to admit one where installTrees refused it;
# the two entries now run the one predicate over every forest a state carries
# and report the one message.
maskRefusal <- paste0(
  "warm-start donor holds a tree that splits on a variable outside this ",
  "forest's allowed column set; the donor's fit is incompatible with the ",
  "column restriction (a forest's own column subset or a restricted ",
  "variance forest) in force here"
)
refusalOf <- function(expr) {
  tryCatch(
    {
      expr
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
}
restrictedSampler <- function() {
  dbarts(
    xHet,
    yHet,
    control = hetControl(FALSE),
    variance = varianceForest(vars = 2L, n.trees = nVarianceTrees)
  )
}
maskDonor <- hetSampler(FALSE)
invisible(maskDonor$run(60L, 0L))
maskDonor$storeState()
maskRecipient <- restrictedSampler()
invisible(maskRecipient$run(60L, 0L))
maskRecipient$storeState()
# non-vacuity: the donor really does split on the column the recipient forbids,
# and the recipient's own variance trees never do
expect_true(any(maskDonor$state[[1L]][["variance.vars"]] == 1L))
expect_false(any(maskRecipient$state[[1L]][["variance.vars"]] == 1L))

maskBefore <- maskRecipient$state
expect_identical(
  refusalOf(maskRecipient$setState(maskDonor$state)),
  maskRefusal
)
# the same donor at the other entry, in the same words: the two cannot disagree
expect_identical(refusalOf(maskRecipient$installTrees(maskDonor)), maskRefusal)
# transactional: both refusals validate before they mutate, so the recipient is
# the sampler it was, and still runs
maskRecipient$storeState()
expect_identical(maskRecipient$state, maskBefore)
expect_true(all(is.finite(maskRecipient$run(0L, 3L)$variance)))

# ---- and a compliant state still round trips ----
# the recipient's own state is inside the mask, so the gate is a no-op and the
# restore reproduces it exactly
maskRecipient$storeState()
maskSaved <- maskRecipient$state
maskRestored <- restrictedSampler()
maskRestored$setState(maskSaved)
maskRestored$storeState()
expect_identical(maskRestored$state, maskSaved)
maskRestoredRun <- maskRestored$run(0L, 3L)
expect_true(all(is.finite(maskRestoredRun$variance)))
expect_true(all(maskRestoredRun$variance > 0))
