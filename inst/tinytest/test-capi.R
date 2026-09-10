# The flat C API (inst/include/dbarts/dbarts.h): compiles a small consumer
# against the installed headers, resolves every entry point through
# R_GetCCallable as a LinkingTo package would, and drives the
# conditional-sampling workout stan4bart performs. Skips wherever the
# consumer cannot be compiled.

source(
  system.file("common", "capiConsumer.R", package = "dbarts"),
  local = TRUE
)
consumer <- compileCapiConsumer("capi", "the C API consumer")
if (!is.null(consumer$skip)) {
  exit_file(consumer$skip)
}
consumerSource <- consumer$consumerSource
includeDir <- consumer$includeDir
dll <- consumer$dll
CALL <- consumer$CALL

# the signature token, resolved BOTH ways: through the stubs and through the
# raw R_GetCCallable canary, which is the un-stubbed per-symbol path a consumer
# that declines DBARTS_USE_STUBS still relies on. The two must agree, and the
# installed library must carry the token this consumer compiled against - the
# only runtime signal that separates a stale consumer binary from a fresh one
# while the version constants stay put.
hashes <- CALL("capi_hash")
expect_true(hashes$raw.agrees)
expect_true(hashes$matches.header)
# and the token really moved on each re-signing: these are the literals the
# pre-reshape header and the post-reshape header baked, so a token blind to
# either change would still read one of them
expect_false(identical(hashes$text, "0x1a911c00bb26dcd7"))
expect_false(identical(hashes$text, "0xcd88efcd67de55d7"))
# the layout-blind token, back when the fold covered the entry-point signatures
# alone: a token that still could not see a struct's layout, an ABI enumerator
# or the callback's parameters would read this one
expect_false(identical(hashes$text, "0x85bd1ef04beb3848"))
# the token the header baked before dbarts_sampler_predict took a thread
# count: a signature that gained a parameter must move it
expect_false(identical(hashes$text, "0x6c9776ae1197e8f5"))
# the token before the freeze slice: the family enum, the four get renames,
# the widened printEvery and printTrees' useLiveTrees all moved it
expect_false(identical(hashes$text, "0x66d33f1613892406"))
# the token before the shape freeze: seven void setters became int, getTrees
# and printTrees took the forest index second, and setForestBasis renamed its
# basis parameter, so a token blind to a return type, an argument order or a
# parameter name would still read this one
expect_false(identical(hashes$text, "0x0939c0224353505b"))
# the token before dbarts_drawLatents' ordinal threshold parameter was renamed
# off the split-candidate grid's spelling: a parameter rename is an ABI
# acknowledgment here, so a token blind to one would still read this
expect_false(identical(hashes$text, "0x5a32aa4cd3872d55"))
# the token from the first half of that rename, when only the threshold vector
# itself had moved and its length still spelled numCutpoints: a token blind to
# one parameter of a signature would still read this
expect_false(identical(hashes$text, "0xc4a2d83f6050bb1f"))
# the token before dbarts_column_type gained its ordered-factor enumerator: an
# appended enumerator moves it, which is what makes the enumerator list the
# place a new column kind is declared
expect_false(identical(hashes$text, "0xb6c0e97dc0688991"))
# the token before dbarts_predictor_source gained its code channel: an
# appended field moves it, since the fold carries each struct's size
expect_false(identical(hashes$text, "0xe14b499a84f501d2"))
# and again when the code channel gained the width that bounds it: a second
# appended field moves the same size the first did
expect_false(identical(hashes$text, "0x37288e7c56449b34"))
# and again when dbarts_results LOST two fields: a pre-1.0-0 removal shifts
# every field below it and shrinks the struct, both of which the fold sees
expect_false(identical(hashes$text, "0xca7b56a64c812b8d"))
# the token the header baked while it still declared the four SEXP entries and
# the twenty other entries this release trims: half the surface leaving is what
# the signature half of the fold sees, and the callback's parameter list and
# dbarts_forest_calibration's layout leaving the fold is what the rest sees
expect_false(identical(hashes$text, "0x616ffcda8c947777"))
# the token before the per-draw callback: one entry added to the surface and
# one struct added to the layout fold, so a token blind to either half would
# still read this
expect_false(identical(hashes$text, "0xab4909b71853c7df"))
# and it does NOT move for doc text outside what it folds, which the token
# cannot see
expect_identical(hashes$text, "0x6380bf095d5cae3f")

# the two version components did NOT move: no version of this API has shipped,
# so whatever they read at the first release becomes the initial contract, and
# the hash above is what acknowledges a pre-release change
versions <- CALL("capi_versions")
expect_equal(versions, c(1L, 0L))

# the enum round trip: this consumer's compiled-in dbarts_family numbering
# agrees with the installed header's, in header (not alphabetical) order
familyConstants <- CALL("capi_family_constants")
expect_equal(unname(familyConstants), 0:8)

set.seed(0)
n <- 150L
p <- 4L
x <- matrix(runif(n * p), n, p)
y <- 2 + x[, 1L] + 0.5 * sin(4 * x[, 2L]) + rnorm(n, 0, 0.4)
x.test <- matrix(runif(20L * p), 20L, p)

nSamples <- 7L
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  updateState = FALSE,
  seed = 99L
)
# THE HANDLE: this header has no creation entry. A consumer builds the
# sampler in R and reads the handle out of the object's own external pointer,
# which is what every ptr below is - so the R object has to stay reachable
# for as long as the C side holds one, and a fresh sampler is a fresh
# dbarts() call rather than a second creation off one specification.
spec <- dbarts(x, y, test = x.test, control = control)

# queries and a run into caller-owned buffers
ptr1 <- spec$getPointer()
dims <- CALL("capi_dims", ptr1)
expect_equal(dims, c(n, p, 20L, 1L, 25L, 0L, 0L, 0L))
expect_equal(CALL("capi_sampler_family", ptr1), familyConstants[["gaussian"]])

CALL("capi_sample_trees_from_prior", ptr1)
r1 <- CALL("capi_run", ptr1, 5L, nSamples, TRUE, TRUE)
expect_equal(length(r1$sigma), as.integer(nSamples))
expect_true(all(is.finite(r1$sigma)) && all(r1$sigma > 0))
expect_equal(length(r1$train), n * nSamples)
expect_true(all(is.finite(r1$train)))
expect_equal(length(r1$test), 20L * nSamples)
expect_equal(length(r1$varcount), p * nSamples)
expect_true(all(colSums(matrix(r1$varcount, p)) > 0))
# training fits are on the response scale
expect_true(abs(mean(r1$train) - mean(y)) < 1)

# a seeded single-chain sampler reproduces bitwise across creations
spec2 <- dbarts(x, y, test = x.test, control = control)
ptr2 <- spec2$getPointer()
CALL("capi_sample_trees_from_prior", ptr2)
r2 <- CALL("capi_run", ptr2, 5L, nSamples, TRUE, TRUE)
expect_identical(r1, r2)

# the write guard on the size-first results struct: a caller whose struct
# predates a field (structSize pinned below it) is never written past its
# declared size, even for fields the sampler produces (varcount, test) - a
# size-blind write would crash on the poisoned pointers the canary installs
# past the boundary
specGuard <- dbarts(x, y, test = x.test, control = control)
ptrGuard <- specGuard$getPointer()
CALL("capi_sample_trees_from_prior", ptrGuard)
expect_true(CALL("capi_run_guard", ptrGuard, 5L, nSamples))

# the zero-structSize guard: a caller that forgets to set results.structSize
# (leaves it 0) is rejected outright, not silently given an all-skip no-op that
# leaves its buffers uninitialized - the flat-API footgun that garbaged a
# consumer's Gibbs loop
specZero <- dbarts(x, y, test = x.test, control = control)
ptrZero <- specZero$getPointer()
CALL("capi_sample_trees_from_prior", ptrZero)
expect_error(
  CALL("capi_run_zero_structsize", ptrZero, 5L, nSamples),
  "structSize"
)

# the per-observation log-likelihood channel: for a gaussian sampler it must
# equal dnorm(y, train, sigma, log = TRUE) recomputed on the same draws (train
# is n x nSamples, sigma constant within a draw), pairing observation-fastest
specLL <- dbarts(x, y, test = x.test, control = control)
ptrLL <- specLL$getPointer()
CALL("capi_sample_trees_from_prior", ptrLL)
rLL <- CALL("capi_run_loglik", ptrLL, 5L, nSamples)
expect_equal(length(rLL$loglik), n * nSamples)
expectedLL <- dnorm(
  rep(y, times = nSamples),
  rLL$train,
  rep(rLL$sigma, each = n),
  log = TRUE
)
expect_equal(rLL$loglik, expectedLL, tolerance = 1e-12)

# THE PER-DRAW CALLBACK. A registered C observer fires once per RECORDED draw
# per chain: five burn-in sweeps produce no call, so the count is nSamples and
# not nSamples + 5, and the draw indices run 0, 1, ... within the chain. The
# consumer's status word is how it reports a disagreement - a callback cannot
# raise - so a nonzero sum here is a malformed draw, not a failed comparison.
specDraw <- dbarts(x, y, test = x.test, control = control)
ptrDraw <- specDraw$getPointer()
CALL("capi_draw_reset", -1L)
CALL("capi_set_draw_callback", ptrDraw, TRUE)
rDraw <- CALL("capi_run", ptrDraw, 5L, nSamples, TRUE, TRUE)
drawReport <- CALL("capi_draw_report")
expect_equal(sum(drawReport$calls), as.integer(nSamples))
expect_equal(drawReport$calls[1L], as.integer(nSamples))
expect_equal(drawReport$last.draw.index[1L], nSamples - 1L)
expect_equal(sum(drawReport$status), 0L)
# the library fills structSize with its own sizeof, which for a consumer built
# against this installed header is the size it compiled
expect_equal(drawReport$struct.size, drawReport$expected.struct.size)
expect_equal(drawReport$num.observations, n)
expect_equal(drawReport$num.predictors, p)
expect_equal(drawReport$num.reported.locations, 1L)
# and the draw the callback saw IS the draw the run recorded
expect_equal(drawReport$last.sigma[1L], rDraw$sigma[nSamples])

# a null function clears: the next run fires nothing
CALL("capi_set_draw_callback", ptrDraw, FALSE)
CALL("capi_draw_reset", -1L)
invisible(CALL("capi_run", ptrDraw, 0L, 3L, TRUE, FALSE))
expect_equal(sum(CALL("capi_draw_report")$calls), 0L)

# and a second registration takes effect on the run after it
CALL("capi_set_draw_callback", ptrDraw, TRUE)
CALL("capi_draw_reset", -1L)
invisible(CALL("capi_run", ptrDraw, 0L, 4L, TRUE, FALSE))
expect_equal(sum(CALL("capi_draw_report")$calls), 4L)

# a nonzero return ABORTS the run at that draw, and the entry still returns
# normally: no condition is raised, so the caller learns of the stop from its
# own context and discards the results. This sampler is inconsistent with them
# afterwards and is not reused.
CALL("capi_draw_reset", 2L)
invisible(CALL("capi_run", ptrDraw, 0L, 6L, TRUE, FALSE))
expect_equal(sum(CALL("capi_draw_report")$calls), 3L)
rm(ptrDraw, specDraw)
invisible(gc(FALSE))

# two chains on two threads: each chain counts its own draws, and each chain
# writes only its own slot - the header's (chainIndex, drawIndex) discipline,
# which is what makes a callback safe with no engine lock
controlChains <- dbartsControl(
  n.chains = 2L,
  n.threads = 2L,
  n.trees = 25L,
  updateState = FALSE,
  seed = 7L
)
specChains <- dbarts(x, y, control = controlChains)
ptrChains <- specChains$getPointer()
CALL("capi_draw_reset", -1L)
CALL("capi_set_draw_callback", ptrChains, TRUE)
invisible(CALL("capi_run", ptrChains, 2L, 4L, TRUE, FALSE))
reportChains <- CALL("capi_draw_report")
expect_equal(reportChains$calls[1:2], c(4L, 4L))
expect_equal(sum(reportChains$calls), 8L)
expect_equal(reportChains$last.draw.index[1:2], c(3L, 3L))
expect_equal(sum(reportChains$status), 0L)
rm(ptrChains, specChains)
invisible(gc(FALSE))

# THE R ROUTE. bart()'s own 'callback' argument builds a per-run hook from
# the SAME two pointers, independent of dbarts_sampler_setDrawCallback above:
# a hook registered through that C entry does not fire here, and one
# registered through 'callback' does not fire through capi_run's own
# 'callback' slot either - the two channels never see each other's draws.
callbackFn <- CALL("capi_draw_function")
callbackContext <- CALL("capi_draw_context")

# exactly n.samples calls at bart()'s own n.burn/n.samples defaults (500,
# 500): the burn-in regression. A hook that saw the burn phase too would
# report 1000 calls, not 500 - the split bart()'s runWithBurnIn makes
# installs the callback on the kept-sample run alone. Kept small everywhere
# BUT n.burn/n.samples so 1000 total sweeps stays fast.
set.seed(101)
nBurnIn <- 16L
xBurnIn <- matrix(runif(nBurnIn * 2L), nBurnIn, 2L)
yBurnIn <- xBurnIn[, 1L] + rnorm(nBurnIn, 0, 0.2)
CALL("capi_draw_reset", -1L)
fitBurnIn <- bart(
  xBurnIn,
  yBurnIn,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  callback = list(fn = callbackFn, context = callbackContext),
  verbose = FALSE
)
reportBurnIn <- CALL("capi_draw_report")
expect_equal(sum(reportBurnIn$calls), 500L)
expect_equal(reportBurnIn$calls[1L], 500L)
expect_equal(reportBurnIn$last.draw.index[1L], 499L)
expect_equal(sum(reportBurnIn$status), 0L)
# a callback fit's automatic keepFits = FALSE means no train channel at all
expect_null(fitBurnIn$yhat.train)

# per-chain counts under n.threads > 1, driven the same way through bart()
nMulti <- 20L
xMulti <- matrix(runif(nMulti * 2L), nMulti, 2L)
yMulti <- xMulti[, 1L] + rnorm(nMulti, 0, 0.2)
CALL("capi_draw_reset", -1L)
fitMulti <- bart(
  xMulti,
  yMulti,
  n.chains = 2L,
  n.threads = 2L,
  n.trees = 5L,
  n.burn = 3L,
  n.samples = 6L,
  callback = list(fn = callbackFn, context = callbackContext),
  verbose = FALSE
)
reportMulti <- CALL("capi_draw_report")
expect_equal(reportMulti$calls[1:2], c(6L, 6L))
expect_equal(sum(reportMulti$calls), 12L)
expect_equal(reportMulti$last.draw.index[1:2], c(5L, 5L))
expect_equal(sum(reportMulti$status), 0L)

# a stop-flag run: the callback aborts partway through the kept-sample run,
# and bart() surfaces that as an error distinct from an interrupt, exactly as
# capi_run does above
nStop <- 20L
xStop <- matrix(runif(nStop * 2L), nStop, 2L)
yStop <- xStop[, 1L] + rnorm(nStop, 0, 0.2)
CALL("capi_draw_reset", 2L)
expect_error(
  bart(
    xStop,
    yStop,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.burn = 0L,
    n.samples = 6L,
    callback = list(fn = callbackFn, context = callbackContext),
    verbose = FALSE
  ),
  "callback"
)
expect_equal(sum(CALL("capi_draw_report")$calls), 3L)
CALL("capi_draw_reset", -1L)

# null buffers skip quantities
r3 <- CALL("capi_run", ptr2, 0L, 2L, FALSE, FALSE)
expect_null(r3$train)
expect_null(r3$test)
expect_equal(length(r3$sigma), 2L)

# the Gibbs conditioning surface: an external offset enters the training
# fits, and with a fixed residual prior an externally set sigma is held
offset <- rep(100, n)
CALL("capi_set_offset", ptr2, offset, TRUE)
rOffset <- CALL("capi_run", ptr2, 30L, 3L, TRUE, FALSE)
expect_true(abs(mean(rOffset$train) - mean(y)) < 3)

specFixed <- dbarts(x, y, resid.prior = fixed(1), control = control)
ptrFixed <- specFixed$getPointer()
expect_equal(CALL("capi_set_sigma", ptrFixed, 0.37), 1L)
rFixed <- CALL("capi_run", ptrFixed, 0L, 3L, FALSE, FALSE)
expect_equal(unique(rFixed$sigma), 0.37)

# the remaining conditioning hooks: a replacement response moves the fits,
# gaussian weights install, and leaf parameters redraw from the prior. Each
# answers 1 - the accepting half of the capability channel, without which a 0
# somewhere below would not discriminate
yShifted <- y + 10
expect_equal(CALL("capi_set_response", ptr2, yShifted, TRUE), 1L)
rShifted <- CALL("capi_run", ptr2, 10L, 3L, TRUE, FALSE)
expect_true(abs(mean(rShifted$train) - mean(yShifted)) < 3)
CALL("capi_set_response", ptr2, y, TRUE)

# thinning, thread, and verbosity controls apply to subsequent runs
CALL("capi_set_run_controls", ptr2, 1L, 2L, FALSE)
rThinned <- CALL("capi_run", ptr2, 0L, 2L, FALSE, FALSE)
expect_equal(length(rThinned$sigma), 2L)
CALL("capi_set_run_controls", ptr2, 1L, 1L, FALSE)

# printEvery divides the iteration count in the print condition, so 0 is a
# division by zero in the engine rather than "never print"
expect_error(
  CALL("capi_set_verbose", ptr2, TRUE, 0L),
  "dbarts_sampler_setVerbose: printEvery must be at least 1"
)
CALL("capi_set_verbose", ptr2, FALSE, 1L)

# a live-tree dump goes through the R console without touching state
printed <- capture.output(CALL("capi_print_trees", ptr2, FALSE, 0L))
expect_true(is.character(printed))

# useLiveTrees forces the live branch on a sampler that has tree storage on
# but nothing recorded yet: FALSE still hits the empty-store refusal there,
# TRUE bypasses it and prints
specLiveTrees <- dbarts(x, y, test = x.test, control = control)
ptrLiveTrees <- specLiveTrees$getPointer()
CALL("capi_set_tree_storage", ptrLiveTrees, TRUE, 2L)
expect_true(
  length(capture.output(CALL("capi_print_trees", ptrLiveTrees, TRUE, 0L))) > 0L
)
expect_error(
  CALL("capi_print_trees", ptrLiveTrees, FALSE, 0L),
  "holds no recorded draws"
)
rm(ptrLiveTrees, specLiveTrees)
invisible(gc(FALSE))

# latents: absent for gaussian, sign-locked to the response for probit
expect_null(CALL("capi_get_latents", ptr1))

yBinary <- rbinom(n, 1L, pnorm(scale(y)))
specBinary <- dbarts(x, yBinary, control = control)
ptrBinary <- specBinary$getPointer()
invisible(CALL("capi_run", ptrBinary, 3L, 2L, FALSE, FALSE))
latents <- CALL("capi_get_latents", ptrBinary)
expect_equal(length(latents), n)
expect_true(all(latents[yBinary == 1] > 0))
expect_equal(
  CALL("capi_sampler_family", ptrBinary),
  familyConstants[["probit"]]
)
expect_true(all(latents[yBinary == 0] <= 0))

# probit log-likelihood: the recorded train fits are the latent location eta,
# so the channel must be the stable log dbinom(y, 1, pnorm(eta)) - log Phi(eta)
# for a success, log Phi(-eta) for a failure - and it agrees with the R twin's
# dbinom(y, 1, pnorm(eta)) form wherever that stays finite
rBinLL <- CALL("capi_run_loglik", ptrBinary, 3L, 2L)
etaBin <- rBinLL$train
yRep <- rep(as.double(yBinary), times = 2L)
stableLL <- ifelse(
  yRep == 1,
  pnorm(etaBin, log.p = TRUE),
  pnorm(etaBin, lower.tail = FALSE, log.p = TRUE)
)
expect_equal(rBinLL$loglik, stableLL, tolerance = 1e-12)
twinLL <- dbinom(yRep, 1L, pnorm(etaBin), log = TRUE)
finite <- is.finite(twinLL)
expect_equal(rBinLL$loglik[finite], twinLL[finite], tolerance = 1e-9)

# dbarts_sampler_setSigma carries the same family rule as the R bridge entry:
# a probit sampler's sigma is pinned at 1 by the model definition, so a write
# would persist (no redraw corrects it) and rescale every leaf posterior
# precision. The flat entry answers this in the capability channel rather than
# by unwinding - no sigma would have worked - while the R5 twin still raises.
# ptrFixed above is the permitted case - gaussian with resid.prior = fixed(),
# the outer-Gibbs conduit - and answers 1 there
expect_equal(CALL("capi_set_sigma", ptrBinary, 0.5), 0L)

# the nbinom magnitude cap at the flat funnel, which has no R layer ahead of it
# to state the rule: the dispersion grid's count histogram is sized from the
# largest count, so the bound is an allocation bound, and creation and the one
# conduit that swaps y here both carry it. The counts are built without drawing,
# so the stream stays where the seeds above put it
yCount <- as.double(seq_len(n) %% 7L)
specCount <- dbarts(x, yCount, family = "nbinom", control = control)
ptrCount <- specCount$getPointer()
capRefusal <- "counts no larger than 1000000"
expect_error(
  CALL("capi_set_response", ptrCount, replace(yCount, 1L, 1000001), FALSE),
  capRefusal
)

# the dispersion channel from C: the results slot appended to dbarts_results.
# The slot is NA-poisoned before the run, so an unfilled channel cannot pass
# for a filled one
disp <- CALL("capi_run_dispersion", ptrCount, 2L, 3L)
expect_true(disp$present)
expect_equal(length(disp$recorded), 3L)
expect_true(all(is.finite(disp$recorded)))
expect_true(all(disp$recorded > 0))
# a caller whose structSize predates the field is never written past, on the
# one family that HAS a dispersion to write
expect_true(disp$guarded)
expect_equal(CALL("capi_sampler_family", ptrCount), familyConstants[["nbinom"]])

rm(ptrCount, specCount, capRefusal)
invisible(gc(FALSE))

# tree storage, prediction, and the state round trip stan4bart's
# predict-after-reload uses
CALL("capi_set_tree_storage", ptr1, TRUE, nSamples)
r4 <- CALL("capi_run", ptr1, 0L, nSamples, FALSE, TRUE)
expect_equal(CALL("capi_dims", ptr1)[6L], as.integer(nSamples))

pred1 <- CALL("capi_predict", ptr1, x.test, NULL)
expect_equal(length(pred1), 20L * nSamples)
expect_equal(pred1, as.double(r4$test))

predOffset <- CALL("capi_predict", ptr1, x.test, rep(2, 20L))
expect_equal(predOffset, pred1 + 2)

# the per-call thread count: it does not persist and cannot move a value, so
# every count answers with pred1 bit for bit. 0 resolves to the sampler's own
# count, which capi_predict itself already passes.
for (nThreads in c(0L, 1L, 2L)) {
  expect_identical(CALL("capi_predict_threads", ptr1, x.test, nThreads), pred1)
}
rm(nThreads)

# recorded draws oldest first, so the second run's draws are the last of them
r5 <- CALL("capi_run", ptr1, 0L, 2L, FALSE, TRUE)
pred5 <- CALL("capi_predict", ptr1, x.test, NULL)
expect_equal(length(pred5), 20L * nSamples)
expect_equal(tail(pred5, 20L * 2L), as.double(r5$test))

# ---------------------------------------------------------------------------
# The self-describing predictor source. One struct carries the four predictor
# entries' arguments, so a C consumer hands the sampler compressed-column
# storage for prediction and test data without densifying it, and every
# argument declares its own width and its own CSC column count.
# ---------------------------------------------------------------------------

# a CSC column over a dense vector: the rows differing from the implicit value
# are stored, and every other row reads the implicit one
cscColumn <- function(values, implicit) {
  stored <- which(values != implicit)
  list(i = as.integer(stored - 1L), x = as.double(values[stored]))
}

# an R-built dbarts_predictor_source, member for member, so a malformed
# argument is as easy to hand the entries as a well-formed one
makeSource <- function(
  numRows,
  numColumns,
  dense = NULL,
  cscColumns = NULL,
  map = NULL,
  types = NULL,
  counts = NULL,
  refs = NULL,
  numCscColumns = NULL,
  codes = NULL,
  numCodeColumns = NULL
) {
  spec <- list(
    numRows = as.integer(numRows),
    numColumns = as.integer(numColumns)
  )
  if (!is.null(dense)) {
    spec$dense <- as.double(dense)
  }
  if (!is.null(cscColumns)) {
    sizes <- vapply(cscColumns, function(column) length(column$i), integer(1L))
    spec$cscColumnPointers <- as.integer(c(0L, cumsum(sizes)))
    spec$cscRowIndices <- as.integer(unlist(lapply(cscColumns, `[[`, "i")))
    spec$cscValues <- as.double(unlist(lapply(cscColumns, `[[`, "x")))
  }
  if (!is.null(numCscColumns)) {
    spec$numCscColumns <- as.integer(numCscColumns)
  }
  if (!is.null(map)) {
    spec$columnSources <- as.integer(map)
  }
  if (!is.null(types)) {
    spec$columnTypes <- as.integer(types)
  }
  if (!is.null(counts)) {
    spec$categoryCounts <- as.integer(counts)
  }
  if (!is.null(refs)) {
    spec$referenceCodes <- as.integer(refs)
  }
  if (!is.null(codes)) {
    spec$denseCodes <- as.integer(codes)
  }
  if (!is.null(numCodeColumns)) {
    spec$numDenseCodeColumns <- as.integer(numCodeColumns)
  }
  spec
}

set.seed(909L)
nSrc <- 120L
levelsSrc <- c("a", "b", "c")
labelsSrc <- sample(levelsSrc, nSrc, replace = TRUE)
frameSrc <- data.frame(x1 = rnorm(nSrc))
frameSrc$f <- factor(labelsSrc, levels = levelsSrc)
ySrc <- 2 * match(labelsSrc, levelsSrc) + rnorm(nSrc, 0, 0.2)
controlSrc <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  updateState = FALSE,
  seed = 17L
)
specSrc <- dbarts(frameSrc, ySrc, control = controlSrc)
ptrSrc <- specSrc$getPointer()
CALL("capi_set_tree_storage", ptrSrc, TRUE, 4L)
rSrc <- CALL("capi_run", ptrSrc, 20L, 4L, FALSE, FALSE)
# the categorical column is split on, so the implicit-value legs below are not
# vacuous: a wrong reference has somewhere to show up
expect_true(sum(matrix(rSrc$varcount, 2L)[2L, ]) > 0L)

codesSrc <- as.matrix(specSrc$data@x)
storage.mode(codesSrc) <- "double"

nTestSrc <- 30L
x1TestSrc <- rnorm(nTestSrc)
x1TestSrc[c(2L, 5L, 9L, 14L, 21L, 27L)] <- 0 # implicit rows of an ordinal CSC
codesTestSrc <- as.double(sample(0:2, nTestSrc, replace = TRUE))
codesTestSrc[c(1L, 4L, 8L)] <- 1 # implicit rows under reference "b"
denseTestSrc <- cbind(x1TestSrc, codesTestSrc)
predDenseSrc <- CALL("capi_predict", ptrSrc, denseTestSrc, NULL)

# the sparse-predict oracle: a CSC source predicts bitwise identically to the
# dense source holding the materialized same values, on an ordinal CSC column,
# a categorical CSC column with a NONZERO reference, an all-CSC design, and a
# column whose every row is implicit
srcOrdinal <- makeSource(
  nTestSrc,
  2L,
  dense = codesTestSrc,
  cscColumns = list(cscColumn(x1TestSrc, 0)),
  map = c(-1L, 0L),
  types = c(0L, 1L)
)
expect_identical(CALL("capi_predict_source", ptrSrc, srcOrdinal), predDenseSrc)

srcCategorical <- makeSource(
  nTestSrc,
  2L,
  dense = x1TestSrc,
  cscColumns = list(cscColumn(codesTestSrc, 1)),
  map = c(0L, -1L),
  types = c(0L, 1L),
  counts = c(0L, 3L),
  refs = c(-1L, 1L)
)
expect_identical(
  CALL("capi_predict_source", ptrSrc, srcCategorical),
  predDenseSrc
)

srcBoth <- makeSource(
  nTestSrc,
  2L,
  cscColumns = list(cscColumn(x1TestSrc, 0), cscColumn(codesTestSrc, 1)),
  map = c(-1L, -2L),
  types = c(0L, 1L),
  counts = c(0L, 3L),
  refs = c(-1L, 1L)
)
expect_identical(CALL("capi_predict_source", ptrSrc, srcBoth), predDenseSrc)

codesAllImplicit <- rep(1, nTestSrc)
predAllImplicit <- CALL(
  "capi_predict",
  ptrSrc,
  cbind(x1TestSrc, codesAllImplicit),
  NULL
)
srcAllImplicit <- makeSource(
  nTestSrc,
  2L,
  dense = x1TestSrc,
  cscColumns = list(cscColumn(codesAllImplicit, 1)),
  map = c(0L, -1L),
  types = c(0L, 1L),
  counts = c(0L, 3L),
  refs = c(-1L, 1L)
)
expect_identical(
  CALL("capi_predict_source", ptrSrc, srcAllImplicit),
  predAllImplicit
)

# the code channel: a caller whose factor column is already integer codes
# hands it over as integers, and the answer is the dense one bitwise. Both
# blocks are indexed by the SAME columnSources[j], within the channel the
# sampler's kind for that column selects - so under the identity map the code
# block is indexed by predictor position and its column 0 goes unread
codesIntSrc <- as.integer(codesTestSrc)
srcCodedIdentity <- makeSource(
  nTestSrc,
  2L,
  dense = c(x1TestSrc, rep(0, nTestSrc)),
  codes = c(rep(0L, nTestSrc), codesIntSrc)
)
expect_identical(
  CALL("capi_predict_source", ptrSrc, srcCodedIdentity),
  predDenseSrc
)

# and an explicit map is what makes the packed form legal: each column names
# its index WITHIN its own channel, so both blocks are one column wide
srcCodedPacked <- makeSource(
  nTestSrc,
  2L,
  dense = x1TestSrc,
  codes = codesIntSrc,
  map = c(0L, 0L)
)
expect_identical(
  CALL("capi_predict_source", ptrSrc, srcCodedPacked),
  predDenseSrc
)

# the code channel composes with CSC storage rather than densifying it: the
# ordinal column stays sparse beside the coded factor column
srcCodedCsc <- makeSource(
  nTestSrc,
  2L,
  codes = codesIntSrc,
  cscColumns = list(cscColumn(x1TestSrc, 0)),
  map = c(-1L, 0L),
  types = c(0L, 1L)
)
expect_identical(
  CALL("capi_predict_source", ptrSrc, srcCodedCsc),
  predDenseSrc
)

# a null code channel is the double path every caller written before the
# field takes, unchanged
srcNoCodes <- makeSource(nTestSrc, 2L, dense = denseTestSrc)
expect_identical(CALL("capi_predict_source", ptrSrc, srcNoCodes), predDenseSrc)

# a code channel does not supply a NUMERIC dense column's values, so a source
# that leaves denseValues out is refused rather than read through a null
srcCodesOnly <- makeSource(
  nTestSrc,
  2L,
  codes = c(rep(0L, nTestSrc), codesIntSrc)
)
expect_error(
  CALL("capi_predict_source", ptrSrc, srcCodesOnly),
  "names no denseValues"
)

# and the code block declares its own extent, so an index past it is refused
# rather than read past the caller's array: the identity map with a PACKED
# one-column code block is exactly that mistake
srcCodesNarrow <- makeSource(
  nTestSrc,
  2L,
  dense = c(x1TestSrc, rep(0, nTestSrc)),
  codes = codesIntSrc
)
expect_error(
  CALL("capi_predict_source", ptrSrc, srcCodesNarrow),
  "names code column 1, but the source declares 1"
)

# the channel's own missing marker is R's integer NA, and it reaches the
# entries as the NA every double-typed rule tests for - here the refusal a
# test NA takes against a training column that carried none
codesWithNA <- codesIntSrc
codesWithNA[3L] <- NA_integer_
srcCodedNA <- makeSource(
  nTestSrc,
  2L,
  dense = x1TestSrc,
  codes = codesWithNA,
  map = c(0L, 0L)
)
expect_error(
  CALL("capi_predict_source", ptrSrc, srcCodedNA),
  "missing values but the training column had none"
)

# the negative half of that oracle: force the categorical column's implicit
# rows to read the storage's zero rather than the declared reference and the
# answers must part company
srcWrongReference <- makeSource(
  nTestSrc,
  2L,
  dense = x1TestSrc,
  cscColumns = list(cscColumn(codesTestSrc, 1)),
  map = c(0L, -1L),
  types = c(0L, 1L),
  counts = c(0L, 3L),
  refs = c(-1L, 0L)
)
expect_false(identical(
  CALL("capi_predict_source", ptrSrc, srcWrongReference),
  predDenseSrc
))

# ... while an UNDECLARED reference on the same ordinal column is accepted and
# reads exactly as a source declaring none at all: "< 0" is the absence a code
# type with no sentinel cannot express
srcUndeclared <- makeSource(
  nTestSrc,
  2L,
  dense = codesTestSrc,
  cscColumns = list(cscColumn(x1TestSrc, 0)),
  map = c(-1L, 0L),
  types = c(0L, 1L),
  refs = c(-1L, -1L)
)
expect_identical(
  CALL("capi_predict_source", ptrSrc, srcUndeclared),
  CALL("capi_predict_source", ptrSrc, srcOrdinal)
)
expect_identical(
  CALL("capi_predict_source", ptrSrc, srcUndeclared),
  predDenseSrc
)

# the self-description matrix: an argument that does not describe itself is
# refused rather than read to the sampler's own width. The dense entries used
# to infer p from the sampler and consume whatever lay past a narrow caller's
# matrix - a measured one-unit swing on the response scale, silently
wideTest <- makeSource(nTestSrc, 3L, dense = cbind(denseTestSrc, x1TestSrc))
expect_error(CALL("capi_predict_source", ptrSrc, wideTest), "columns")
# DBARTS_COLUMN_ORDERED_FACTOR (2) is a declared column type a source may
# carry: the entrance accepts it and answers exactly as the same source
# declaring DBARTS_COLUMN_ORDINAL does, since a source's declared types are
# validated and then discarded - the STORE's types are what route.
expect_identical(
  CALL(
    "capi_predict_source",
    ptrSrc,
    makeSource(nTestSrc, 2L, dense = denseTestSrc, types = c(2L, 1L))
  ),
  predDenseSrc
)
# a column type outside {ordinal, categorical, ordered factor}, a declared
# level count or a reference code past the engine's category limit, and a
# source naming a CSC column it does not carry (without the bound, the ~v
# decode reads past the caller's own pointer array)
expect_error(
  CALL(
    "capi_predict_source",
    ptrSrc,
    makeSource(nTestSrc, 2L, dense = denseTestSrc, types = c(0L, 7L))
  ),
  "DBARTS_COLUMN_ORDERED_FACTOR"
)
expect_error(
  CALL(
    "capi_predict_source",
    ptrSrc,
    makeSource(nTestSrc, 2L, dense = denseTestSrc, counts = c(0L, 70000L))
  ),
  "categoryCounts"
)
expect_error(
  CALL(
    "capi_predict_source",
    ptrSrc,
    makeSource(nTestSrc, 2L, dense = denseTestSrc, refs = c(-1L, 70000L))
  ),
  "referenceCodes"
)
expect_error(
  CALL(
    "capi_predict_source",
    ptrSrc,
    makeSource(
      nTestSrc,
      2L,
      dense = x1TestSrc,
      cscColumns = list(cscColumn(codesTestSrc, 1)),
      map = c(0L, -2L),
      numCscColumns = 1L,
      types = c(0L, 1L),
      counts = c(0L, 3L),
      refs = c(-1L, 1L)
    )
  ),
  "CSC column"
)

# the read guard on the caller-filled struct, the input-side twin of the
# results write guard: a caller whose struct predates the typing channel pins
# structSize below it, and those members - pointed at an unmapped page here -
# must never be read
expect_identical(
  CALL("capi_predict_truncated", ptrSrc, denseTestSrc),
  predDenseSrc
)

# forest addressing on a single-forest sampler: forest 0 is exactly today's
# behavior, and an index past the last forest is an error on both entries that
# take one rather than a read past the last forest (the engine's printers
# index their forest unchecked, so the bridge's check is the only guard)
expect_equal(CALL("capi_num_trees", ptrSrc, 0L), 20L)
expect_error(CALL("capi_num_trees", ptrSrc, 1L), "forest index out of range")
expect_error(
  CALL("capi_print_trees", ptrSrc, FALSE, 1L),
  "forest index out of range"
)
# both print branches take the index: ptr2 has no tree storage (the live
# trees print), ptrSrc has it (the saved ones)
expect_true(
  length(capture.output(CALL("capi_print_trees", ptr2, FALSE, 0L))) > 0L
)
expect_true(
  length(capture.output(CALL("capi_print_trees", ptrSrc, FALSE, 0L))) > 0L
)

# useLiveTrees overrides the store on a sampler that HAS one: forcing the
# live branch prints a different (unsaved, most-recent) tree than the saved
# read the same call without it takes
printSrcLive <- capture.output(CALL("capi_print_trees", ptrSrc, TRUE, 0L))
printSrcSaved <- capture.output(CALL("capi_print_trees", ptrSrc, FALSE, 0L))
expect_true(length(printSrcLive) > 0L)
expect_false(identical(printSrcLive, printSrcSaved))


# the Student-t df channel from C: the results slot appended to
# dbarts_results after the dispersion one, on a sampler whose error law is
# selected by the model's resid.dist rather than by the family string. The
# slot is NA-poisoned before the run, so an unfilled channel cannot pass for
# a filled one
specT <- dbarts(x, y, resid.dist = student(df = 5), control = control)
ptrT <- specT$getPointer()
dfT <- CALL("capi_run_residual_df", ptrT, 2L, 3L)
expect_true(dfT$present)
expect_equal(length(dfT$recorded), 3L)
# a FIXED df repeats the value the sampler was created with, every draw
expect_equal(dfT$recorded, rep(5, 3L))
# a caller whose structSize predates the field is never written past, on the
# one error law that HAS a df to write
expect_true(dfT$guarded)
# and it is ABSENT off the error law: a gaussian sampler handed the buffer
# leaves the poisoned slot exactly as it found it
dfG <- CALL("capi_run_residual_df", ptr1, 2L, 3L)
expect_true(all(is.na(dfG$recorded)))
# a Student-t residual sampler's family IS gaussian: resid.dist selects the
# error law, not a family of its own
expect_equal(CALL("capi_sampler_family", ptrT), familyConstants[["gaussian"]])
rm(specT, ptrT, dfT, dfG)
invisible(gc(FALSE))

# ---------------------------------------------------------------------------
# COPY-ON-SET. The value setters copy into a buffer the sampler owns, so the
# caller's array is free the moment the call returns. The probe hands the
# entry a buffer the consumer owns, OVERWRITES that buffer before returning,
# and reads a channel that goes back to the installed vector every draw: the
# per-observation log-likelihood is dnorm(y_i, ., .) with the offset in the
# mean, so a setter that RETAINED the pointer scores the clobbered values.
# (The fits alone would not discriminate: gaussian bakes y into its working
# scale at set time, so the leak shows only where the vector is re-read.)
# ---------------------------------------------------------------------------
specCopyA <- dbarts(x, y, control = control)
specCopyB <- dbarts(x, y, control = control)
yConditioned <- y + 3
yClobber <- y - 99
expect_equal(
  CALL(
    "capi_set_response_clobber",
    specCopyA$getPointer(),
    yConditioned,
    yClobber
  ),
  1L
)
expect_equal(
  CALL("capi_set_response", specCopyB$getPointer(), yConditioned, FALSE),
  1L
)
expect_identical(
  CALL("capi_run_loglik", specCopyA$getPointer(), 5L, 3L),
  CALL("capi_run_loglik", specCopyB$getPointer(), 5L, 3L)
)
# and the clobbered values are not merely equal to the conditioned ones: a
# sampler actually set from them scores differently, so the assertion above
# has somewhere to fail
specCopyC <- dbarts(x, y, control = control)
expect_equal(
  CALL("capi_set_response", specCopyC$getPointer(), yClobber, FALSE),
  1L
)
expect_false(identical(
  CALL("capi_run_loglik", specCopyC$getPointer(), 5L, 3L),
  CALL("capi_run_loglik", specCopyB$getPointer(), 5L, 3L)
))
# the offset conduit carries the same contract, and the same channel reads it
specOffA <- dbarts(x, y, control = control)
specOffB <- dbarts(x, y, control = control)
offsetConditioned <- rep(2, n)
expect_equal(
  CALL(
    "capi_set_offset_clobber",
    specOffA$getPointer(),
    offsetConditioned,
    rep(-50, n)
  ),
  1L
)
expect_equal(
  CALL("capi_set_offset", specOffB$getPointer(), offsetConditioned, FALSE),
  1L
)
expect_identical(
  CALL("capi_run_loglik", specOffA$getPointer(), 5L, 3L),
  CALL("capi_run_loglik", specOffB$getPointer(), 5L, 3L)
)
rm(specCopyA, specCopyB, specCopyC, specOffA, specOffB)
invisible(gc(FALSE))

# ---------------------------------------------------------------------------
# THE HANDLE ACROSS A RESTORE. The handle is the address in the R object's
# external pointer, so it is valid only until that object replaces its
# pointer - which is exactly what a save/load round trip makes it do. A
# consumer re-reads the handle afterwards and drives the re-created engine.
# ---------------------------------------------------------------------------
specRestore <- dbarts(x, y, control = control)
handleBefore <- specRestore$getPointer()
invisible(CALL("capi_run", handleBefore, 5L, 2L, FALSE, FALSE))
specRestore$storeState()
revived <- unserialize(serialize(specRestore, NULL))
handleAfter <- revived$getPointer()
# a deserialized external pointer carries no address, so getPointer built a
# new engine from the stored state and the old handle is not it
expect_false(identical(handleBefore, handleAfter))
rRestored <- CALL("capi_run", handleAfter, 0L, 2L, TRUE, FALSE)
expect_equal(length(rRestored$train), n * 2L)
expect_true(all(is.finite(rRestored$train)))
expect_equal(CALL("capi_dims", handleAfter)[1L], n)
rm(specRestore, revived, handleBefore, handleAfter, rRestored)
invisible(gc(FALSE))

# ---------------------------------------------------------------------------
# DESTROY. The engine goes early, the R object is left in its dead-pointer
# state, and its own methods re-create from a stored state or refuse for want
# of one. A second destroy is the one call a destroyed handle still takes.
# ---------------------------------------------------------------------------
specDestroy <- dbarts(x, y, control = control)
invisible(CALL("capi_run", specDestroy$getPointer(), 5L, 2L, FALSE, FALSE))
specDestroy$storeState()
handleDestroyed <- specDestroy$getPointer()
CALL("capi_destroy", handleDestroyed)
CALL("capi_destroy", handleDestroyed) # a no-op, not a double free
# the R object's pointer now reads DEAD, which is the branch its own methods
# take to re-create; without this the re-creation below would prove nothing,
# since a live pointer is returned as it stands
expect_false(.Call(dbarts:::C_dbarts_bartcore_isValidPointer, handleDestroyed))
handleRecreated <- specDestroy$getPointer()
expect_false(identical(handleDestroyed, handleRecreated))
expect_true(all(is.finite(
  CALL("capi_run", handleRecreated, 0L, 2L, TRUE, FALSE)$train
)))

# and without a stored state the R object refuses rather than re-creating
specNoState <- dbarts(x, y, control = control)
CALL("capi_destroy", specNoState$getPointer())
expect_error(
  specNoState$getPointer(),
  "cannot be re-created without a stored state"
)
rm(specDestroy, specNoState, handleDestroyed, handleRecreated)
invisible(gc(FALSE))

# The handshake, D3: every stub enforces major-equality plus a minor floor by
# default on its first resolution, and checks the exact-ABI hash too only when
# the consumer opts in with DBARTS_REQUIRE_EXACT_ABI before including the
# header. Five arms, each its own temp dir and object name (the loader caches
# by path) and the same compile pattern as the consumer above, so each skips
# wherever that one skips: (a) a wrong hash ALONE now passes silently - the
# half proving the default gate does not look at it; (b) the same wrong hash
# plus the opt-in macro raises the ABI mismatch; (c) and (d) a wrong major or
# minor raise the version mismatch regardless of the opt-in; (e) the opt-in
# macro alone, hash untouched, loads and calls clean - the configuration both
# consumers ship. Arm (a) probes through capi_versions rather than
# capi_dims: the version accessors are allocation-free stubs, so the
# handshake still runs on their first resolution without touching a sampler.
# The other arms drive capi_dims on the handle spec already holds - this
# consumer owns no sampler, so nothing an arm resolves outlives its
# dyn.unload.
compileHandshakeConsumer <- function(label, extraFlags) {
  dir <- tempfile(paste0("capi-", label))
  dir.create(dir)
  file.copy(consumerSource, file.path(dir, paste0(label, ".c")))
  writeLines(
    sprintf('PKG_CPPFLAGS = -I"%s" %s', includeDir, extraFlags),
    file.path(dir, "Makevars")
  )
  owd <- setwd(dir)
  output <- tryCatch(
    suppressWarnings(system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "SHLIB", paste0(label, ".c")),
      stdout = TRUE,
      stderr = TRUE
    )),
    error = function(e) e
  )
  setwd(owd)
  lib <- file.path(dir, paste0(label, .Platform$dynlib.ext))
  if (!file.exists(lib)) {
    if (nzchar(Sys.getenv("CI", ""))) {
      stop(
        "could not compile the ",
        label,
        " C API consumer under CI:\n",
        paste(output, collapse = "\n")
      )
    }
    return(NULL)
  }
  lib
}

libWrongHashAlone <- compileHandshakeConsumer(
  "stale-hash-alone",
  "-DDBARTS_C_API_HASH=0x0123456789abcdefULL"
)
if (!is.null(libWrongHashAlone)) {
  dllWrongHashAlone <- dyn.load(libWrongHashAlone)
  versionsWrongHashAlone <- .Call(
    getNativeSymbolInfo("capi_versions", dllWrongHashAlone)
  )
  expect_equal(versionsWrongHashAlone, c(1L, 0L))
  dyn.unload(libWrongHashAlone)
}

libWrongHashExact <- compileHandshakeConsumer(
  "stale-hash-exact",
  "-DDBARTS_C_API_HASH=0x0123456789abcdefULL -DDBARTS_REQUIRE_EXACT_ABI"
)
if (!is.null(libWrongHashExact)) {
  dllWrongHashExact <- dyn.load(libWrongHashExact)
  expect_error(
    .Call(getNativeSymbolInfo("capi_dims", dllWrongHashExact), ptr1),
    "dbarts C ABI mismatch"
  )
  dyn.unload(libWrongHashExact)
}

libWrongMajor <- compileHandshakeConsumer(
  "stale-major",
  "-DDBARTS_C_API_MAJOR=99"
)
if (!is.null(libWrongMajor)) {
  dllWrongMajor <- dyn.load(libWrongMajor)
  expect_error(
    .Call(getNativeSymbolInfo("capi_dims", dllWrongMajor), ptr1),
    "dbarts C API version mismatch"
  )
  dyn.unload(libWrongMajor)
}

libWrongMinor <- compileHandshakeConsumer(
  "stale-minor",
  "-DDBARTS_C_API_MINOR=99"
)
if (!is.null(libWrongMinor)) {
  dllWrongMinor <- dyn.load(libWrongMinor)
  expect_error(
    .Call(getNativeSymbolInfo("capi_dims", dllWrongMinor), ptr1),
    "dbarts C API version mismatch"
  )
  dyn.unload(libWrongMinor)
}

libCorrectExact <- compileHandshakeConsumer(
  "correct-exact",
  "-DDBARTS_REQUIRE_EXACT_ABI"
)
if (!is.null(libCorrectExact)) {
  dllCorrectExact <- dyn.load(libCorrectExact)
  dimsCorrectExact <-
    .Call(getNativeSymbolInfo("capi_dims", dllCorrectExact), ptr1)
  expect_equal(dimsCorrectExact[1L], n)
  dyn.unload(libCorrectExact)
}
