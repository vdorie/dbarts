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
# the token before the multinomial run and predict: two accessors appended to
# the surface and the draw struct's numReportedLocations renamed, so a token
# blind to an appended entry or a field name would still read this
expect_false(identical(hashes$text, "0x6380bf095d5cae3f"))
# the token before the nbinom shape fields: two struct fields renamed, which only
# the by-name layout fold sees
expect_false(identical(hashes$text, "0xa7415a6f1bcc93c3"))
# the token before the swap entries' status took names: a fourth enumerator
# list joined the fold with no signature moving, so a token that folded only
# the lists it already knew would still read this
expect_false(identical(hashes$text, "0xa33182bf349aa60d"))
# the token before the sweep hook: one entry, a parameter, two enum lists
# (the statuses and the draw channels) with DBARTS_DRAW_ALL's text, the chain
# view's layout and operations, both callback types' text and the draw's
# latents joined the fold
expect_false(identical(hashes$text, "0xa84b67742536fbc4"))
# the tokens of this header with one parameter retyped at the same width: the
# sweep hook's isBurnIn unsigned, the view's setOffset taking a non-const
# offset, the draw callback a non-const draw. A token blind to a callback's
# or an operation's parameter list would read one of these
expect_false(identical(hashes$text, "0x0337e88b95749642"))
expect_false(identical(hashes$text, "0x1945d8d3f58ded0d"))
expect_false(identical(hashes$text, "0x42ac4b080ae4b8dd"))
# and this header's with DBARTS_DRAW_ALL naming only today's ten channels
expect_false(identical(hashes$text, "0xd2db09fcb6caecbf"))
# and it does NOT move for doc text outside what it folds, which the token
# cannot see
expect_identical(hashes$text, "0x7220b30db0d94fd4")

# the two version components did NOT move: no version of this API has shipped,
# so whatever they read at the first release becomes the initial contract, and
# the hash above is what acknowledges a pre-release change
versions <- CALL("capi_versions")
expect_equal(versions, c(1L, 0L))

# the enum round trip: this consumer's compiled-in dbarts_family numbering
# agrees with the installed header's, in header (not alphabetical) order
familyConstants <- CALL("capi_family_constants")
expect_equal(unname(familyConstants), 0:8)
# and the swap entries' status: none, done, scale kept
expect_identical(CALL("capi_swap_constants"), 0:2)

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
  updateState = FALSE,
  seed = 99L
)
# THE HANDLE: this header has no creation entry. A consumer builds the
# sampler in R and reads the handle out of the object's own external pointer,
# which is what every ptr below is - so the R object has to stay reachable
# for as long as the C side holds one, and a fresh sampler is a fresh
# dbarts() call rather than a second creation off one specification.
spec <- dbarts(x, y, test = x.test, control = control, n.trees = 25L)

# queries and a run into caller-owned buffers
ptr1 <- spec$getPointer()
dims <- CALL("capi_dims", ptr1)
expect_equal(dims, c(n, p, 20L, 1L, 25L, 0L, 0L, 0L, 1L, 1L))
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
spec2 <- dbarts(x, y, test = x.test, control = control, n.trees = 25L)
ptr2 <- spec2$getPointer()
CALL("capi_sample_trees_from_prior", ptr2)
r2 <- CALL("capi_run", ptr2, 5L, nSamples, TRUE, TRUE)
expect_identical(r1, r2)

# the write guard on the size-first results struct: a caller whose struct
# predates a field (structSize pinned below it) is never written past its
# declared size, even for fields the sampler produces (varcount, test) - a
# size-blind write would crash on the poisoned pointers the canary installs
# past the boundary
specGuard <- dbarts(x, y, test = x.test, control = control, n.trees = 25L)
ptrGuard <- specGuard$getPointer()
CALL("capi_sample_trees_from_prior", ptrGuard)
expect_true(CALL("capi_run_guard", ptrGuard, 5L, nSamples))

# the zero-structSize guard: a caller that forgets to set results.structSize
# (leaves it 0) is rejected outright, not silently given an all-skip no-op that
# leaves its buffers uninitialized - the flat-API footgun that garbaged a
# consumer's Gibbs loop
specZero <- dbarts(x, y, test = x.test, control = control, n.trees = 25L)
ptrZero <- specZero$getPointer()
CALL("capi_sample_trees_from_prior", ptrZero)
expect_error(
  CALL("capi_run_zero_structsize", ptrZero, 5L, nSamples),
  "structSize"
)

# the per-observation log-likelihood channel: for a gaussian sampler it must
# equal dnorm(y, train, sigma, log = TRUE) recomputed on the same draws (train
# is n x nSamples, sigma constant within a draw), pairing observation-fastest
specLL <- dbarts(x, y, test = x.test, control = control, n.trees = 25L)
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
specDraw <- dbarts(x, y, test = x.test, control = control, n.trees = 25L)
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
expect_equal(drawReport$num.fitted.values.per.observation, 1L)
expect_equal(drawReport$num.variable.count.forests, 1L)
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
# own context and discards the results.
CALL("capi_draw_reset", 2L)
stopRun <- tryCatch(
  {
    CALL("capi_run", ptrDraw, 0L, 6L, TRUE, FALSE)
    "returned"
  },
  interrupt = function(cond) "an interrupt"
)
expect_equal(stopRun, "returned")
expect_equal(sum(CALL("capi_draw_report")$calls), 3L)
# and a negative one, the callback's failure, stops it the same way
CALL("capi_draw_reset", 2L)
CALL("capi_draw_stop_status", -1L)
stopRun <- tryCatch(
  {
    CALL("capi_run", ptrDraw, 0L, 6L, TRUE, FALSE)
    "returned"
  },
  interrupt = function(cond) "an interrupt"
)
expect_equal(stopRun, "returned")
expect_equal(sum(CALL("capi_draw_report")$calls), 3L)

# the run drops the saved trees it wrote: of a full store of ten draws, a
# run stopped at its draw 5 overwrote six, and the sampler predicts from the
# four it left, oldest first
CALL("capi_set_tree_storage", ptrDraw, TRUE, 10L)
CALL("capi_draw_reset", -1L)
invisible(CALL("capi_run", ptrDraw, 0L, 10L, FALSE, FALSE))
predFull <- CALL("capi_predict", ptrDraw, x.test, NULL)
expect_equal(length(predFull), 20L * 10L)
CALL("capi_draw_reset", 5L)
invisible(CALL("capi_run", ptrDraw, 0L, 10L, FALSE, FALSE))
expect_equal(sum(CALL("capi_draw_report")$calls), 6L)
expect_identical(
  CALL("capi_predict", ptrDraw, x.test, NULL),
  predFull[(6L * 20L + 1L):(10L * 20L)]
)
CALL("capi_set_draw_callback", ptrDraw, FALSE)
rm(predFull, ptrDraw, specDraw)
invisible(gc(FALSE))

# A RAISING callback, on a fresh sampler because a run that stops mid-stream
# is not reused. dbarts_draw_callback admits an R error from an INLINE run -
# one chain here, so every call reaches this thread - and dbarts_sampler_run
# protects itself against the jump: the run stops at the raising draw, the
# entry does NOT return, and the condition reaches R carrying the callback's
# own message.
specRaise <- dbarts(x, y, control = control, n.trees = 25L)
ptrRaise <- specRaise$getPointer()

# registered but idle first, so raising is the only difference between the two
# runs below
CALL("capi_set_raising_callback", ptrRaise, -1L)
invisible(CALL("capi_run_raising", ptrRaise, 0L, 4L))
idleReport <- CALL("capi_raise_report")
expect_equal(idleReport$calls, 4L)
expect_true(idleReport$returned)

CALL("capi_set_raising_callback", ptrRaise, 2L)
expect_error(
  CALL("capi_run_raising", ptrRaise, 0L, 6L),
  "callback refused draw 2"
)
raiseReport <- CALL("capi_raise_report")
# three draws arrived and the sixth never did, and the entry never returned:
# the jump really came out of the middle of the run rather than after it
expect_equal(raiseReport$calls, 3L)
expect_false(raiseReport$returned)

# and the handle survives it. The callback is cleared and the SAME sampler
# runs again - the "a later run or a destroy is safe" half of the contract -
# and the cleared hook fires nothing.
CALL("capi_set_draw_callback", ptrRaise, FALSE)
CALL("capi_draw_reset", -1L)
afterRaise <- CALL("capi_run", ptrRaise, 0L, 3L, TRUE, FALSE)
expect_equal(length(afterRaise$sigma), 3L)
expect_true(all(is.finite(afterRaise$sigma)) && all(afterRaise$sigma > 0))
expect_equal(sum(CALL("capi_draw_report")$calls), 0L)

# an explicit destroy after the raise, which is the other way out
CALL("capi_destroy", ptrRaise)
rm(ptrRaise, specRaise)
invisible(gc(FALSE))

# and the same from a MULTI-CHAIN inline run: two chains on one thread, so
# every call still reaches this thread, and the jump comes out of the middle of
# the engine's chain loop rather than its first iteration - where the per-chain
# results vector and the cancel closure are live. Chains run in order, so the
# refused draw is chain 0's.
controlRaiseChains <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  updateState = FALSE,
  seed = 11L
)
specRaiseChains <- dbarts(x, y, control = controlRaiseChains, n.trees = 25L)
ptrRaiseChains <- specRaiseChains$getPointer()
CALL("capi_set_raising_callback", ptrRaiseChains, 1L)
expect_error(
  CALL("capi_run_raising", ptrRaiseChains, 1L, 5L),
  "callback refused draw 1"
)
chainRaiseReport <- CALL("capi_raise_report")
expect_equal(chainRaiseReport$calls, 2L)
expect_false(chainRaiseReport$returned)

# the session is intact afterwards, which is the check standing in for a
# memory-check run here: a FRESH sampler builds and runs ...
specAfterRaise <- dbarts(x, y, control = control, n.trees = 25L)
ptrAfterRaise <- specAfterRaise$getPointer()
afterFresh <- CALL("capi_run", ptrAfterRaise, 2L, 3L, TRUE, FALSE)
expect_true(all(is.finite(afterFresh$sigma)) && all(afterFresh$sigma > 0))
rm(ptrAfterRaise, specAfterRaise)

# ... and the sampler the raise came out of destroys cleanly
CALL("capi_destroy", ptrRaiseChains)
rm(ptrRaiseChains, specRaiseChains, controlRaiseChains)
invisible(gc(FALSE))

# two chains on two threads: each chain counts its own draws, and each chain
# writes only its own slot - the header's (chainIndex, drawIndex) discipline,
# which is what makes a callback safe with no engine lock
controlChains <- dbartsControl(
  n.chains = 2L,
  n.threads = 2L,
  updateState = FALSE,
  seed = 7L
)
specChains <- dbarts(x, y, control = controlChains, n.trees = 25L)
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
stopped <- tryCatch(
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
  interrupt = function(cond) "an interrupt",
  error = conditionMessage
)
expect_true(is.character(stopped) && grepl("callback", stopped))
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

# a named leaf-prior sd keeps its meaning in response units across a
# re-anchor made through the flat entry, which never returns to R: the engine
# holds the sd and restates it as k against the data's new scale. A held sd's
# spread stays the sd; a prior on the sd has its chi scale restated, and the
# k drawn so far is left to lag until its next draw. The reads are the
# engine's own calibration, since the R object is not told of the move.
engineLeafPrior <- function(sampler) {
  .Call(dbarts:::C_dbarts_bartcore_getLeafPrior, sampler$getPointer(), 0L)
}
# an offset that moves the range of the response less the offset
rangeOffset <- 4 * x[, 3L] - 1
namedControl <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  updateState = FALSE,
  seed = 99L
)
heldSd <- dbarts(
  x,
  y,
  control = namedControl,
  leaf.prior = normal(sd = 0.6),
  n.trees = 25L
)
heldBefore <- engineLeafPrior(heldSd)
expect_true(max(abs(heldBefore[, "prior.sd"] / 0.6 - 1)) < 1e-12)
expect_equal(
  CALL("capi_set_offset", heldSd$getPointer(), rangeOffset, TRUE),
  1L
)
heldAfter <- engineLeafPrior(heldSd)
expect_true(all(
  abs(heldAfter[, "prior.scale"] / heldBefore[, "prior.scale"] - 1) > 0.05
))
expect_true(max(abs(heldAfter[, "prior.sd"] / 0.6 - 1)) < 1e-12)
expect_identical(heldAfter[, "k"], heldAfter[, "prior.scale"] / 0.6)
expect_identical(heldAfter[, "named.sd"], c(0.6, 0.6))

drawnSd <- dbarts(
  x,
  y,
  control = namedControl,
  leaf.prior = normal(sd = invchi(3, 0.6)),
  n.trees = 25L
)
invisible(CALL("capi_run", drawnSd$getPointer(), 20L, 1L, FALSE, FALSE))
drawnBefore <- engineLeafPrior(drawnSd)
expect_identical(
  drawnBefore[, "k.prior.scale"],
  drawnBefore[, "prior.scale"] / 0.6
)
# the chains have left their common start
expect_true(drawnBefore[1L, "k"] != drawnBefore[2L, "k"])
expect_equal(
  CALL("capi_set_offset", drawnSd$getPointer(), rangeOffset, TRUE),
  1L
)
drawnAfter <- engineLeafPrior(drawnSd)
expect_true(all(
  abs(drawnAfter[, "prior.scale"] / drawnBefore[, "prior.scale"] - 1) > 0.05
))
expect_identical(
  drawnAfter[, "k.prior.scale"],
  drawnAfter[, "prior.scale"] / 0.6
)
expect_identical(drawnAfter[, "k"], drawnBefore[, "k"])
# and without the scale update nothing is restated
expect_equal(CALL("capi_set_offset", drawnSd$getPointer(), offset, FALSE), 1L)
expect_identical(engineLeafPrior(drawnSd), drawnAfter)

# the swap status. A re-derivation with a row out of the likelihood and the
# rows in at fewer than two distinct values of the response less the offset
# has too little to read: the entry keeps the scale, answers 2 and raises
# nothing, and the sampler is where the same swap without updateScale leaves
# it, in its calibration and in the sweeps that follow. Each sampler is
# created with every row in and loses its rows afterwards, so creation has
# nothing to say; the vectors handed over span another range than creation's,
# so an entry that re-derived from every row would move the calibration.
keepRows <- 1:3
keepWeights <- replace(numeric(n), keepRows, 1)
keepOffset <- replace(rangeOffset, keepRows, 0.25)
keepResponse <- replace(3 * y + 10, keepRows, 7)
keepSampler <- function() {
  sampler <- dbarts(
    x,
    replace(y, keepRows, 1.5),
    control = namedControl,
    n.trees = 25L
  )
  sampler$setWeights(keepWeights)
  sampler
}
countConditions <- function(expr, class = "condition") {
  count <- 0L
  value <- withCallingHandlers(expr, condition = function(cond) {
    if (inherits(cond, class)) {
      count <<- count + 1L
    }
    if (inherits(cond, "warning")) invokeRestart("muffleWarning")
  })
  list(value = value, count = count)
}
kept <- keepSampler()
twin <- keepSampler()
keptBefore <- engineLeafPrior(kept)
flatKeep <- countConditions(c(
  CALL("capi_set_offset", kept$getPointer(), keepOffset, TRUE),
  CALL("capi_set_response", kept$getPointer(), keepResponse, TRUE)
))
expect_identical(flatKeep$value, c(2L, 2L))
expect_identical(flatKeep$count, 0L)
expect_identical(
  c(
    CALL("capi_set_offset", twin$getPointer(), keepOffset, FALSE),
    CALL("capi_set_response", twin$getPointer(), keepResponse, FALSE)
  ),
  c(1L, 1L)
)
expect_identical(engineLeafPrior(kept), keptBefore)
expect_identical(engineLeafPrior(kept), engineLeafPrior(twin))
expect_identical(
  CALL("capi_run", kept$getPointer(), 0L, 5L, TRUE, FALSE),
  CALL("capi_run", twin$getPointer(), 0L, 5L, TRUE, FALSE)
)
# the flat keeps spent nothing: the R method's first keep on the same sampler
# still warns, once, and its second does not
methodKeep <- countConditions(
  kept$setOffset(keepOffset, updateScale = TRUE),
  "dbartsScaleFallbackWarning"
)
expect_identical(methodKeep$count, 1L)
methodKeep <- countConditions(
  kept$setOffset(keepOffset, updateScale = TRUE),
  "warning"
)
expect_identical(methodKeep$count, 0L)
# and a flat keep after it is still a 2, as a swap that does not ask is a 1
expect_identical(
  CALL("capi_set_offset", kept$getPointer(), keepOffset, TRUE),
  2L
)
expect_identical(
  CALL("capi_set_offset", kept$getPointer(), keepOffset, FALSE),
  1L
)
# two rows in at two values are enough: 1, and the transform is theirs
twoRows <- dbarts(x, y, control = namedControl, n.trees = 25L)
twoRows$setWeights(replace(numeric(n), 4:5, 1))
expect_identical(CALL("capi_set_response", twoRows$getPointer(), y, TRUE), 1L)
expect_identical(
  engineLeafPrior(twoRows)[, "response.scale"],
  rep(max(y[4:5]) - min(y[4:5]), 2L)
)
# a constant response with every row in takes its window and answers 1
constantRows <- dbarts(x, y, control = namedControl, n.trees = 25L)
expect_identical(
  CALL("capi_set_response", constantRows$getPointer(), rep(3, n), TRUE),
  1L
)
expect_identical(
  engineLeafPrior(constantRows)[1L, c("response.scale", "response.shift")],
  c(response.scale = 1, response.shift = 3)
)
# a count sampler under a mask keeping rows of one count keeps its center
# counts written out, not drawn: nothing here moves the session's stream
keepCounts <- replace(as.double(seq_len(n) %% 5L), keepRows, 2)
countRows <- dbarts(
  x,
  keepCounts,
  family = "nbinom",
  control = namedControl,
  n.trees = 25L
)
countBefore <- engineLeafPrior(countRows)
countRows$setActiveRows(keepWeights)
expect_identical(
  CALL("capi_set_offset", countRows$getPointer(), rangeOffset / 10, TRUE),
  2L
)
expect_identical(
  engineLeafPrior(countRows)[, "response.shift"],
  countBefore[, "response.shift"]
)
# a probit sampler takes nothing from its response, so it has nothing to keep
probitRows <- dbarts(x, as.double(y > 2), control = namedControl, n.trees = 25L)
probitRows$setActiveRows(keepWeights)
expect_identical(
  CALL("capi_set_offset", probitRows$getPointer(), rangeOffset / 10, TRUE),
  1L
)
rm(kept, twin, twoRows, constantRows, countRows, probitRows, countConditions)

specFixed <- dbarts(
  x,
  y,
  family = gaussian(sigma = fixed(1)),
  control = control,
  n.trees = 25L
)
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
specLiveTrees <- dbarts(x, y, test = x.test, control = control, n.trees = 25L)
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
specBinary <- dbarts(x, yBinary, control = control, n.trees = 25L)
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
# ptrFixed above is the permitted case - gaussian with a fixed sigma prior,
# the outer-Gibbs conduit - and answers 1 there
expect_equal(CALL("capi_set_sigma", ptrBinary, 0.5), 0L)

# the nbinom magnitude cap at the flat funnel, which has no R layer ahead of it
# to state the rule: the shape grid's count histogram is sized from the
# largest count, so the bound is an allocation bound, and creation and the one
# conduit that swaps y here both carry it. The counts are built without drawing,
# so the stream stays where the seeds above put it
yCount <- as.double(seq_len(n) %% 7L)
specCount <- dbarts(
  x,
  yCount,
  family = "nbinom",
  control = control,
  n.trees = 25L
)
ptrCount <- specCount$getPointer()
capRefusal <- "counts no larger than 1000000"
expect_error(
  CALL("capi_set_response", ptrCount, replace(yCount, 1L, 1000001), FALSE),
  capRefusal
)

# the shape channel from C: the results slot appended to dbarts_results.
# The slot is NA-poisoned before the run, so an unfilled channel cannot pass
# for a filled one
disp <- CALL("capi_run_shape", ptrCount, 2L, 3L)
expect_true(disp$present)
expect_equal(length(disp$recorded), 3L)
expect_true(all(is.finite(disp$recorded)))
expect_true(all(disp$recorded > 0))
# a caller whose structSize predates the field is never written past, on the
# one family that HAS a shape to write
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
  updateState = FALSE,
  seed = 17L
)
specSrc <- dbarts(frameSrc, ySrc, control = controlSrc, n.trees = 20L)
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
# dbarts_results after the shape one, on a sampler whose error
# distribution is selected by the model's resid.df attribute rather than
# by the family string. The slot is NA-poisoned before the run, so an
# unfilled channel cannot pass for a filled one
specT <- dbarts(
  x,
  y,
  family = student(df = 5),
  control = control,
  n.trees = 25L
)
ptrT <- specT$getPointer()
dfT <- CALL("capi_run_residual_df", ptrT, 2L, 3L)
expect_true(dfT$present)
expect_equal(length(dfT$recorded), 3L)
# a FIXED df repeats the value the sampler was created with, every draw
expect_equal(dfT$recorded, rep(5, 3L))
# a caller whose structSize predates the field is never written past, on the
# one error distribution that HAS a df to write
expect_true(dfT$guarded)
# and it is ABSENT off the error distribution: a gaussian sampler handed
# the buffer leaves the poisoned slot exactly as it found it
dfG <- CALL("capi_run_residual_df", ptr1, 2L, 3L)
expect_true(all(is.na(dfG$recorded)))
# a Student-t residual sampler's family IS gaussian: family = student()
# selects the error distribution, not a family of its own
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
specCopyA <- dbarts(x, y, control = control, n.trees = 25L)
specCopyB <- dbarts(x, y, control = control, n.trees = 25L)
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
specCopyC <- dbarts(x, y, control = control, n.trees = 25L)
expect_equal(
  CALL("capi_set_response", specCopyC$getPointer(), yClobber, FALSE),
  1L
)
expect_false(identical(
  CALL("capi_run_loglik", specCopyC$getPointer(), 5L, 3L),
  CALL("capi_run_loglik", specCopyB$getPointer(), 5L, 3L)
))
# the offset conduit carries the same contract, and the same channel reads it
specOffA <- dbarts(x, y, control = control, n.trees = 25L)
specOffB <- dbarts(x, y, control = control, n.trees = 25L)
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
specRestore <- dbarts(x, y, control = control, n.trees = 25L)
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
specDestroy <- dbarts(x, y, control = control, n.trees = 25L)
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
specNoState <- dbarts(x, y, control = control, n.trees = 25L)
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
  compileCapiSource(
    dir,
    paste0(label, ".c"),
    paste0("the ", label, " C API consumer")
  )
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

# the flat conduits refuse a non-finite value, as creation does, and the
# sampler runs on unharmed
specFinite <- dbarts(x, y, control = control, n.trees = 25L)
ptrFinite <- specFinite$getPointer()
expect_error(
  CALL("capi_set_response", ptrFinite, replace(y, 1L, Inf), FALSE),
  "dbarts_sampler_setResponse: response contains non-finite values"
)
expect_error(
  CALL("capi_set_offset", ptrFinite, replace(numeric(n), 1L, NaN), FALSE),
  "dbarts_sampler_setOffset: offset contains non-finite values"
)
expect_error(
  CALL("capi_set_sigma", ptrFinite, Inf),
  "dbarts_sampler_setSigma: sigma must be finite and positive"
)
finiteRun <- specFinite$run(0L, 2L)
expect_true(all(is.finite(finiteRun$sigma)))
rm(specFinite, ptrFinite, finiteRun)

# a negative-binomial row of count w draws w (y + r) Polya-Gamma units a
# sweep, so the flat response conduit holds a counted sampler to the cap
# $setResponse states, r at the drawn shape's largest grid point (50); the
# counts reach the sampler through the bridge entry, which the fitting
# functions do not yet open
countY <- as.double(rpois(n, 3))
nbinomCounts <- dbarts(
  x,
  countY,
  family = "nbinom",
  control = control,
  n.trees = 5L
)
.Call(
  dbarts:::C_dbarts_bartcore_setWeights,
  nbinomCounts$getPointer(),
  replace(rep(1, n), 1L, 2),
  FALSE,
  FALSE,
  NULL
)
expect_error(
  CALL(
    "capi_set_response",
    nbinomCounts$getPointer(),
    replace(countY, 1L, 499975),
    FALSE
  ),
  paste0(
    "dbarts_sampler_setResponse: family \"nbinom\" requires 'weights' times ",
    "(y + r) no larger than 1000000"
  ),
  fixed = TRUE
)
# the bound itself is taken: 2 (499950 + 50) is the cap
expect_equal(
  CALL(
    "capi_set_response",
    nbinomCounts$getPointer(),
    replace(countY, 1L, 499950),
    FALSE
  ),
  1L
)
rm(countY, nbinomCounts)

# a binomial sampler (the counts slot under logistic) refuses a flat response
# on the flat conduit in $setResponse's words
binomialControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 2L,
  n.burn = 0L,
  updateState = FALSE,
  seed = 5L
)
binaryFit <- dbarts(
  x,
  as.double(y > median(y)),
  family = "logistic",
  control = binomialControl,
  n.trees = 5L
)
binomialData <- binaryFit$data
binomialData@counts <- cbind(rep(1L, n), rep(2L, n))
binomialData@y <- rep(3, n)
ptrBinomial <- .Call(
  dbarts:::C_dbarts_bartcore_create,
  binomialControl,
  binaryFit$model,
  binomialData,
  "logistic",
  NULL,
  NULL
)
expect_error(
  CALL("capi_set_response", ptrBinomial, rep(1, n), FALSE),
  paste0(
    "dbarts_sampler_setResponse: its response is the n x 2 count pair of ",
    "successes and failures, which a flat vector cannot express; pass the ",
    "pair to the R sampler's $setResponse, or make a new sampler"
  ),
  fixed = TRUE
)
rm(binaryFit, binomialData, ptrBinomial)

# a sampler with gp leaves that holds saved draws refuses a re-derived range
# on both flat conduits, by an R error naming the entry and before anything is
# read; holding the range answers 1 and leaves the saved draws as they were
specKernel <- dbarts(
  x,
  y,
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 3L,
    n.burn = 10L,
    keepTrees = TRUE,
    updateState = FALSE,
    seed = 99L
  ),
  leaf.prior = gp(columns = 2L),
  n.trees = 8L
)
invisible(specKernel$run())
ptrKernel <- specKernel$getPointer()
kernelDraws <- specKernel$predict(x.test)
refusedKernel <- function(entry) {
  paste0(
    entry,
    ": 'updateScale' cannot be TRUE for a sampler with gp leaves that holds ",
    "saved draws: re-derive the scale during burn-in, before draws are saved ",
    "(burn-in sweeps, or keepTrees off until the scale is fixed), or make a ",
    "new sampler"
  )
}
expect_error(
  CALL("capi_set_response", ptrKernel, 3 * y + 10, TRUE),
  refusedKernel("dbarts_sampler_setResponse"),
  fixed = TRUE
)
expect_error(
  CALL("capi_set_offset", ptrKernel, rep(-5, n), TRUE),
  refusedKernel("dbarts_sampler_setOffset"),
  fixed = TRUE
)
expect_error(
  CALL("capi_set_offset", ptrKernel, NULL, TRUE),
  refusedKernel("dbarts_sampler_setOffset"),
  fixed = TRUE
)
expect_identical(specKernel$predict(x.test), kernelDraws)
# a call that is served answers 1; one wrongly refused answers its message
# here, so that it is counted and the file runs on
served <- function(...) tryCatch(CALL(...), error = conditionMessage)
expect_equal(served("capi_set_response", ptrKernel, y + 1, FALSE), 1L)
expect_identical(specKernel$predict(x.test), kernelDraws)
expect_equal(served("capi_set_offset", ptrKernel, rep(1, n), FALSE), 1L)
expect_identical(specKernel$predict(x.test), kernelDraws)
expect_equal(served("capi_set_offset", ptrKernel, NULL, FALSE), 1L)
expect_identical(specKernel$predict(x.test), kernelDraws)
# a family with no response scale has none to re-derive, so a probit, logistic
# or ordinal sampler with gp leaves and saved draws takes both calls with
# updateScale true (dec-B481), the saved draws as they were
yBinary <- as.double(y > median(y))
yOrdinal <- as.double(cut(y, quantile(y, 0:3 / 3), include.lowest = TRUE))
for (kind in c("probit", "logistic", "ordinal")) {
  responseOf <- if (kind == "ordinal") yOrdinal else yBinary
  specFree <- dbarts(
    x,
    responseOf,
    family = kind,
    control = dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.samples = 3L,
      n.burn = 10L,
      keepTrees = TRUE,
      updateState = FALSE,
      seed = 99L
    ),
    leaf.prior = gp(columns = 2L),
    n.trees = 8L
  )
  invisible(specFree$run())
  ptrFree <- specFree$getPointer()
  freeDraws <- specFree$predict(x.test)
  expect_equal(
    served("capi_set_response", ptrFree, rev(responseOf), TRUE),
    1L,
    info = kind
  )
  expect_equal(
    served("capi_set_offset", ptrFree, rep(-5, n), TRUE),
    1L,
    info = kind
  )
  expect_equal(served("capi_set_offset", ptrFree, NULL, TRUE), 1L, info = kind)
  expect_identical(specFree$predict(x.test), freeDraws, info = kind)
}
rm(yBinary, yOrdinal, responseOf, specFree, ptrFree, freeDraws, kind)
# with nothing saved the same entries take the re-derived range
specBare <- dbarts(
  x,
  y,
  control = control,
  leaf.prior = gp(columns = 2L),
  n.trees = 25L
)
ptrBare <- specBare$getPointer()
expect_equal(served("capi_set_response", ptrBare, 3 * y + 10, TRUE), 1L)
expect_equal(served("capi_set_offset", ptrBare, rep(-5, n), TRUE), 1L)
rm(specKernel, ptrKernel, kernelDraws, refusedKernel, specBare, ptrBare)
rm(served)

# a store too large to hold raises from the flat setter rather than aborting,
# and leaves the store at its previous capacity: the allocation's exception
# leaves the entry through the guard every flat entry is registered under,
# which names the entry and carries the exception's message
specStorage <- dbarts(x, y, control = control, n.trees = 25L)
ptrStorage <- specStorage$getPointer()
CALL("capi_set_tree_storage", ptrStorage, TRUE, 3L)
expect_error(
  CALL("capi_set_tree_storage", ptrStorage, TRUE, -1L),
  "^dbarts_sampler_setTreeStorage: .+"
)
invisible(specStorage$run(0L, 3L))
expect_identical(dim(specStorage$predict(x))[2L], 3L)
rm(specStorage, ptrStorage)

# a store switched on through the flat API survives a save and load: the
# state the R object stores carries it, and the re-created sampler takes its
# capacity, so predict and run work as on the live one
storageControl <- control
storageControl@keepTrees <- FALSE
specFlatStore <- dbarts(x, y, control = storageControl, n.trees = 25L)
CALL("capi_set_tree_storage", specFlatStore$getPointer(), TRUE, 4L)
invisible(specFlatStore$run(3L, 4L))
specFlatStore$storeState()
livePredict <- specFlatStore$predict(x)
storeFile <- tempfile(fileext = ".rds")
saveRDS(specFlatStore, storeFile)
reloaded <- readRDS(storeFile)
unlink(storeFile)
expect_identical(reloaded$predict(x), livePredict)
expect_silent(invisible(reloaded$run(0L, 1L)))
rm(specFlatStore, reloaded, storageControl, livePredict, storeFile)
invisible(gc(FALSE))

# ---- THE FITTED-VALUE AND SPLIT-COUNT WIDTHS. Every buffer below is sized by
# the documented layout from an F and a V computed in R - F the count matrix's
# column count (1 off multinomial), V the sampler's own forest count - and
# carries a canaried tail K bodies long, so a write past the documented size
# shows as a tail that moved. Every channel the flat struct carries is then
# compared bitwise against what the R route returns on the same seed.

tailFactor <- 3L
widths <- function(sampler, numFitted) {
  c(numFitted, dbarts:::bartcoreNumForests(sampler$getPointer()))
}
channelOk <- function(channel) {
  isTRUE(channel$tail.intact)
}

set.seed(77)
nW <- 40L
xW <- matrix(runif(nW * 3L), nW, 3L)
zW <- rbinom(nW, 1L, 0.5)
wW <- runif(nW)
yW <- xW[, 1L] + zW * (1 + xW[, 2L]) + wW * xW[, 3L] + rnorm(nW, sd = 0.2)
xTestW <- xW[1:6, , drop = FALSE]
KW <- 3L
countsW <- matrix(rpois(nW * KW, 1.2), nW, KW)
countsW[rowSums(countsW) == 0L, 1L] <- 1L
storage.mode(countsW) <- "integer"
widthControl <- function(keepTrees = TRUE) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 4L,
    keepTrees = keepTrees,
    updateState = FALSE,
    seed = 21L
  )
}
gaussianW <- function(keepTrees = TRUE) {
  dbarts(xW, yW, test = xTestW, control = widthControl(keepTrees), n.trees = 8L)
}
multinomialW <- function(keepTrees = TRUE) {
  dbarts(
    dbartsData(xW, counts = countsW, test = xTestW),
    family = "multinomial",
    control = widthControl(keepTrees),
    n.trees = 8L
  )
}
bcfW <- function() {
  dbarts(
    xW,
    yW,
    forests = list(forest(), forest(basis = ~ factor(zW))),
    control = widthControl(),
    n.trees = 8L
  )
}
threeW <- function() {
  dbarts(
    xW,
    yW,
    forests = list(
      forest(),
      forest(basis = ~ factor(zW)),
      forest(basis = ~wW)
    ),
    control = widthControl(),
    n.trees = 8L
  )
}

# the accessors, each against its R-side count: F is K, 1, 1, 1 and V is K, 1,
# 1, 2 on multinomial, gaussian, probit and BCF
probitW <- dbarts(
  xW,
  as.numeric(yW > 1),
  family = "probit",
  control = widthControl(),
  n.trees = 8L
)
for (case in list(
  list(sampler = multinomialW(), F = KW),
  list(sampler = gaussianW(), F = 1L),
  list(sampler = probitW, F = 1L),
  list(sampler = bcfW(), F = 1L)
)) {
  expected <- widths(case$sampler, case$F)
  dimsW <- CALL("capi_dims", case$sampler$getPointer())
  expect_equal(dimsW[9L], expected[1L])
  expect_equal(dimsW[10L], expected[2L])
}
expect_equal(CALL("capi_dims", multinomialW()$getPointer())[9:10], c(KW, KW))
expect_equal(CALL("capi_dims", bcfW()$getPointer())[9:10], c(1L, 2L))

# the gaussian parity check comes first: a flat run and an R run of two
# identically seeded samplers agree on every channel both carry
flatG <- gaussianW()
rG <- gaussianW()
widthG <- widths(flatG, 1L)
outG <- CALL(
  "capi_run_canaried",
  flatG$getPointer(),
  3L,
  4L,
  widthG[1L],
  widthG[2L],
  tailFactor
)
runG <- rG$run(3L, 4L)
expect_identical(outG$sigma$body, as.vector(runG$sigma))
expect_identical(outG$train$body, as.vector(runG$train))
expect_identical(outG$test$body, as.vector(runG$test))
expect_identical(outG$varcount$body, as.vector(runG$varcount))
expect_true(all(is.finite(outG$logLikelihood$body)))
for (channel in c("k", "varprobs", "shape", "residualDf")) {
  expect_true(outG[[channel]]$body.untouched, info = channel)
}
for (channel in names(outG)) {
  expect_true(channelOk(outG[[channel]]), info = channel)
}

# multinomial: K fitted values per observation and K split-count sets, the
# softmax probabilities with every row of K summing to 1, sigma the pinned 1,
# and the log-likelihood the engine's quiet NaN in every word
flatM <- multinomialW()
rM <- multinomialW()
widthM <- widths(flatM, KW)
outM <- CALL(
  "capi_run_canaried",
  flatM$getPointer(),
  3L,
  4L,
  widthM[1L],
  widthM[2L],
  tailFactor
)
runM <- rM$run(3L, 4L)
expect_true(all(outM$sigma$body == 1))
expect_identical(outM$train$body, as.vector(runM$train))
expect_identical(outM$test$body, as.vector(runM$test))
expect_true(all(abs(apply(runM$train, c(1L, 3L, 4L), sum) - 1) < 1e-12))
expect_identical(outM$varcount$body, as.vector(runM$varcount))
expect_true(outM$logLikelihood$all.quiet.nan)
for (channel in c("k", "varprobs", "shape", "residualDf")) {
  expect_true(outM[[channel]]$body.untouched, info = channel)
}
for (channel in names(outM)) {
  expect_true(channelOk(outM[[channel]]), info = channel)
}

# BCF and a three-forest sampler: one split-count set per mean forest, the
# prognostic forest's first
flatB <- bcfW()
widthB <- widths(flatB, 1L)
outB <- CALL(
  "capi_run_canaried",
  flatB$getPointer(),
  3L,
  4L,
  widthB[1L],
  widthB[2L],
  tailFactor
)
expect_identical(outB$varcount$body, as.vector(bcfW()$run(3L, 4L)$varcount))
expect_true(channelOk(outB$varcount))
flatT <- threeW()
widthT <- widths(flatT, 1L)
expect_equal(widthT[2L], 3L)
expect_equal(CALL("capi_dims", flatT$getPointer())[10L], 3L)
outT <- CALL(
  "capi_run_canaried",
  flatT$getPointer(),
  3L,
  4L,
  widthT[1L],
  widthT[2L],
  tailFactor
)
runT <- threeW()$run(3L, 4L)
expect_identical(outT$varcount$body, as.vector(runT$varcount))
expect_identical(dim(runT$varcount), c(3L, 3L, 4L, 2L))
expect_true(channelOk(outT$varcount))

# the callback sees F and V equal to the accessors on every multi-forest model
for (build in list(multinomialW, bcfW, threeW)) {
  sampler <- build()
  ptrW <- sampler$getPointer()
  dimsW <- CALL("capi_dims", ptrW)
  CALL("capi_draw_reset", -1L)
  CALL("capi_set_draw_callback", ptrW, TRUE)
  invisible(CALL(
    "capi_run_canaried",
    ptrW,
    0L,
    2L,
    dimsW[9L],
    dimsW[10L],
    tailFactor
  ))
  CALL("capi_set_draw_callback", ptrW, FALSE)
  report <- CALL("capi_draw_report")
  expect_equal(report$num.fitted.values.per.observation, dimsW[9L])
  expect_equal(report$num.variable.count.forests, dimsW[10L])
}

# predict, against $predict on the same sampler, with and without tree storage
for (keepTrees in c(TRUE, FALSE)) {
  sampler <- gaussianW(keepTrees)
  invisible(sampler$run(3L, 4L))
  flatPredict <- CALL(
    "capi_predict_canaried",
    sampler$getPointer(),
    xTestW,
    NULL,
    1L,
    tailFactor
  )
  expect_equal(flatPredict$status, 1L)
  expect_identical(flatPredict$out$body, as.vector(sampler$predict(xTestW)))
  expect_true(channelOk(flatPredict$out))

  sampler <- multinomialW(keepTrees)
  invisible(sampler$run(3L, 4L))
  flatPredict <- CALL(
    "capi_predict_canaried",
    sampler$getPointer(),
    xTestW,
    NULL,
    KW,
    tailFactor
  )
  expect_equal(flatPredict$status, 1L)
  expect_identical(flatPredict$out$body, as.vector(sampler$predict(xTestW)))
  expect_true(channelOk(flatPredict$out))
}

# the offset arms. A gaussian offset still adds after the fit
sampler <- gaussianW()
invisible(sampler$run(3L, 4L))
offsetG <- seq(-1, 1, length.out = nrow(xTestW))
flatPredict <- CALL(
  "capi_predict_canaried",
  sampler$getPointer(),
  xTestW,
  offsetG,
  1L,
  tailFactor
)
expect_identical(
  flatPredict$out$body,
  as.vector(sampler$predict(xTestW, offset = offsetG))
)

# a multinomial one is rows x K, entering before the softmax, bitwise R's
sampler <- multinomialW()
invisible(sampler$run(3L, 4L))
offsetM <- cbind(0.5 * xTestW[, 1L], -0.3, 0.2 * xTestW[, 2L])
flatPredict <- CALL(
  "capi_predict_canaried",
  sampler$getPointer(),
  xTestW,
  offsetM,
  KW,
  tailFactor
)
expect_equal(flatPredict$status, 1L)
expect_identical(
  flatPredict$out$body,
  as.vector(sampler$predict(xTestW, offset = offsetM))
)
expect_true(channelOk(flatPredict$out))
expect_false(identical(
  flatPredict$out$body,
  as.vector(sampler$predict(xTestW))
))
# a non-finite entry raises, as R's does
expect_error(
  CALL(
    "capi_predict_canaried",
    sampler$getPointer(),
    xTestW,
    replace(offsetM, 2L, Inf),
    KW,
    tailFactor
  ),
  "dbarts_sampler_predict: offset requires every category offset entry"
)
# and the conduits a softmax has no use for answer 0
expect_equal(CALL("capi_set_response", sampler$getPointer(), yW, FALSE), 0L)
expect_equal(
  CALL("capi_set_offset", sampler$getPointer(), numeric(nW), FALSE),
  0L
)

# a sampler given a category offset from R - a train one, or only a test one -
# refuses a predict that names none, and an all-zero matrix asks for the
# offset-free surface exactly as R's does
trainOffsetW <- cbind(0.4 * xW[, 1L], -0.2, 0.3 * xW[, 3L])
testOffsetW <- cbind(0.1, 0.5 * xTestW[, 2L], -0.4)
zeroOffsetW <- matrix(0, nrow(xTestW), KW)
withTrain <- multinomialW()
withTrain$setCategoryOffset(trainOffsetW, updateState = FALSE)
withTest <- multinomialW()
withTest$setCategoryTestOffset(testOffsetW)
for (sampler in list(withTrain, withTest)) {
  invisible(sampler$run(3L, 4L))
  expect_error(
    CALL(
      "capi_predict_canaried",
      sampler$getPointer(),
      xTestW,
      NULL,
      KW,
      tailFactor
    ),
    "dbarts_sampler_predict: this sampler carries an n x K category offset"
  )
  expect_error(sampler$predict(xTestW), "carries an n x K category offset")
  flatPredict <- CALL(
    "capi_predict_canaried",
    sampler$getPointer(),
    xTestW,
    zeroOffsetW,
    KW,
    tailFactor
  )
  expect_identical(
    flatPredict$out$body,
    as.vector(sampler$predict(xTestW, offset = zeroOffsetW))
  )
}
rm(
  flatG,
  rG,
  flatM,
  rM,
  flatB,
  flatT,
  runT,
  probitW,
  sampler,
  withTrain,
  withTest,
  outG,
  outM,
  outB,
  outT,
  runG,
  runM,
  flatPredict
)
invisible(gc(FALSE))

# EVERY CHANNEL THE FIT CARRIES reaches a registered callback, whether or not
# the run stores it: the recording callback notes the channels each call hands
# over and copies each one at draw 2 of chain 2. Each arm runs 3 draws on 2
# chains and 2 threads, with results NULL and with every dbarts_results field,
# and through the R run at keepFits TRUE and FALSE; under one seed each copy
# equals the array of the route that stores it (the flat route: train, test,
# varcount, varprobs, the log-likelihood; the R route: the thresholds, the
# variance pair, forestFits and glue). So storing or not moves no draw either.
set.seed(31)
nC <- 60L
pC <- 3L
nTestC <- 8L
xC <- matrix(runif(nC * pC), nC, pC)
xTestC <- matrix(runif(nTestC * pC), nTestC, pC)
zC <- rbinom(nC, 1L, 0.5)
yC <- xC[, 1L] + zC * (1 + xC[, 2L]) + (0.2 + xC[, 3L]) * rnorm(nC, sd = 0.3)
ordinalC <- cut(xC[, 1L] + rnorm(nC, sd = 0.2), 3L, labels = FALSE)
timeC <- exp(yC / 2)
statusC <- rbinom(nC, 1L, 0.7)
channelControl <- function(keepFits = TRUE) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 2L,
    updateState = FALSE,
    seed = 17L,
    keepFits = keepFits
  )
}
channelArms <- list(
  dart = list(
    build = function(keepFits = TRUE) {
      dbarts(
        xC,
        yC,
        test = xTestC,
        split.probs = dart(),
        control = channelControl(keepFits),
        n.trees = 8L
      )
    },
    channels = c("train", "test", "varprobs", "logLikelihood", "varcount")
  ),
  ordinal = list(
    build = function(keepFits = TRUE) {
      dbarts(
        xC,
        ordinalC,
        test = xTestC,
        family = "ordinal",
        control = channelControl(keepFits),
        n.trees = 8L
      )
    },
    channels = c(
      "train",
      "test",
      "logLikelihood",
      "thresholds",
      "varcount",
      "latents"
    )
  ),
  variance = list(
    build = function(keepFits = TRUE) {
      dbarts(
        xC,
        yC,
        test = xTestC,
        variance = varianceForest(n.trees = 4L),
        control = channelControl(keepFits),
        n.trees = 8L
      )
    },
    channels = c(
      "train",
      "test",
      "variance",
      "varianceTest",
      "logLikelihood",
      "varcount"
    )
  ),
  aft = list(
    build = function(keepFits = TRUE) {
      dbarts(
        xC,
        cbind(timeC, statusC),
        test = xTestC,
        family = "aft",
        variance = varianceForest(n.trees = 4L),
        control = channelControl(keepFits),
        n.trees = 8L
      )
    },
    channels = c(
      "train",
      "test",
      "variance",
      "varianceTest",
      "logLikelihood",
      "varcount",
      "latents"
    )
  ),
  bcf = list(
    build = function(keepFits = TRUE) {
      dbarts(
        xC,
        yC,
        forests = list(forest(), forest(basis = ~ factor(zC))),
        control = channelControl(keepFits),
        n.trees = 8L
      )
    },
    channels = c("train", "forestFits", "glue", "varcount")
  )
)
recordingCallback <- list(
  fn = CALL("capi_recording_function"),
  context = CALL("capi_recording_context")
)
# chain 2's draw 2 of an array whose leading extent is perDraw, draws and then
# chains varying slowest
drawOf <- function(values, perDraw) {
  slot <- 1L * 3L + 1L
  as.vector(values)[slot * perDraw + seq_len(perDraw)]
}
flatRecording <- function(sampler, every, channels = NA_integer_) {
  CALL("capi_recording_reset", 1L, 1L)
  CALL("capi_set_recording_callback", sampler$getPointer(), TRUE, channels)
  stored <- CALL("capi_run_every", sampler$getPointer(), 5L, 3L, every)
  CALL("capi_set_recording_callback", sampler$getPointer(), FALSE, NA_integer_)
  list(stored = stored, report = CALL("capi_recording_report", 2L, 3L))
}
rRecording <- function(sampler) {
  CALL("capi_recording_reset", 1L, 1L)
  stored <- sampler$run(5L, 3L, callback = recordingCallback)
  list(stored = stored, report = CALL("capi_recording_report", 2L, 3L))
}
for (armName in names(channelArms)) {
  arm <- channelArms[[armName]]
  flatNone <- flatRecording(arm$build(), FALSE)
  flatEvery <- flatRecording(arm$build(), TRUE)
  rKept <- rRecording(arm$build(TRUE))
  rDropped <- rRecording(arm$build(FALSE))
  reports <- list(
    flat.null = flatNone$report,
    flat.every = flatEvery$report,
    r.keepFits = rKept$report,
    r.noFits = rDropped$report
  )
  for (route in names(reports)) {
    report <- reports[[route]]
    info <- paste(armName, route)
    expect_identical(report$calls, c(3L, 3L), info = info)
    expect_true(report$steady, info = info)
    for (chain in 1:2) {
      expect_identical(report$channels[[chain]], arm$channels, info = info)
    }
    # every route hands over the same draw
    expect_identical(report$copies, flatNone$report$copies, info = info)
    expect_identical(report$sigma, flatNone$report$sigma, info = info)
  }
  copies <- flatNone$report$copies
  flat <- flatEvery$stored
  expect_identical(copies$train, drawOf(flat$train, nC), info = armName)
  expect_identical(
    copies$varcount,
    as.double(drawOf(flat$varcount, length(copies$varcount))),
    info = armName
  )
  if ("test" %in% arm$channels) {
    expect_identical(copies$test, drawOf(flat$test, nTestC), info = armName)
  }
  if ("varprobs" %in% arm$channels) {
    expect_identical(copies$varprobs, drawOf(flat$varprobs, pC), info = armName)
  }
  if ("logLikelihood" %in% arm$channels) {
    expect_identical(
      copies$logLikelihood,
      drawOf(flat$logLikelihood, nC),
      info = armName
    )
  }
  stored <- rKept$stored
  if ("thresholds" %in% arm$channels) {
    expect_identical(
      copies$thresholds,
      drawOf(stored$thresholds, flatNone$report$num.ordinal.thresholds),
      info = armName
    )
  }
  if ("variance" %in% arm$channels) {
    expect_identical(copies$variance, drawOf(stored$variance, nC))
    expect_identical(copies$varianceTest, drawOf(stored$varianceTest, nTestC))
  }
  if ("forestFits" %in% arm$channels) {
    expect_identical(copies$forestFits, drawOf(stored$forestFits, 2L * nC))
    expect_identical(copies$glue, drawOf(stored$glue, 3L))
  }
  # a variance forest has no single residual scale: sigma is NaN in every
  # draw the callback sees and every draw the flat run stores; elsewhere the
  # callback's sigma is the stored one
  if ("variance" %in% arm$channels) {
    expect_true(all(is.nan(flatNone$report$sigma)), info = armName)
    expect_true(all(is.nan(flat$sigma)), info = armName)
    expect_null(stored$sigma, info = armName)
  } else {
    expect_identical(as.vector(flatNone$report$sigma), flat$sigma)
    expect_identical(as.vector(flatNone$report$sigma), as.vector(stored$sigma))
  }
}

# A callback gets the channels it names that the fit carries, the ones the run
# stores, and no other. Each flag alone, with results NULL, hands over that
# channel where the fit carries it and nothing else, the draw it hands over
# the one every channel's registration does; with every dbarts_results array
# stored and no channel named, the ones those arrays hold.
channelFlags <- CALL("capi_draw_channel_flags")
for (armName in names(channelArms)) {
  arm <- channelArms[[armName]]
  everyCopy <- flatRecording(arm$build(), FALSE)$report$copies
  for (flag in names(channelFlags)) {
    named <- if (flag == "variance") c("variance", "varianceTest") else flag
    report <- flatRecording(arm$build(), FALSE, channelFlags[[flag]])$report
    info <- paste(armName, flag)
    expected <- intersect(arm$channels, named)
    for (chain in 1:2) {
      expect_identical(report$channels[[chain]], expected, info = info)
    }
    expect_identical(report$copies[expected], everyCopy[expected], info = info)
  }
  report <- flatRecording(arm$build(), TRUE, 0L)$report
  for (chain in 1:2) {
    expect_identical(
      report$channels[[chain]],
      intersect(
        arm$channels,
        c("train", "test", "varprobs", "logLikelihood", "varcount")
      ),
      info = armName
    )
  }
}
rm(everyCopy, named, expected)

# the scalar sigma is NaN keyed on the variance forest, not on a sigma held
# fixed: a probit sampler stores and hands over the pinned 1, a gaussian one at
# sigma = fixed(1) the 1 it holds
for (pinned in list(
  probit = function() {
    dbarts(
      xC,
      as.numeric(yC > 1),
      family = "probit",
      control = channelControl(),
      n.trees = 8L
    )
  },
  fixed = function() {
    dbarts(
      xC,
      yC,
      family = gaussian(sigma = fixed(1)),
      control = channelControl(),
      n.trees = 8L
    )
  }
)) {
  recorded <- flatRecording(pinned(), TRUE)
  expect_true(all(recorded$stored$sigma == 1))
  expect_true(all(recorded$report$sigma == 1))
}
rm(
  flatNone,
  flatEvery,
  rKept,
  rDropped,
  reports,
  report,
  copies,
  flat,
  stored,
  recorded
)
CALL("capi_recording_reset", 0L, 0L)
invisible(gc(FALSE))

# kIsSampled follows the model as it now stands, through setLeafPrior and
# setModel alike: a fixed-k sampler given a k hyperprior reports 1 and writes
# results.k, its draws varying, and the R run returns k; given a fixed k again
# it reports 0, leaves results.k at the caller's sentinel and the R run returns
# no k
kControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE,
  seed = 23L
)
kFixed <- function() dbarts(xC, yC, control = kControl, n.trees = 8L)
kDrawn <- function() {
  dbarts(
    xC,
    yC,
    leaf.prior = normal(k = chi(1.5, 2)),
    control = kControl,
    n.trees = 8L
  )
}
expectKDrawn <- function(sampler, drawn, info) {
  expect_equal(CALL("capi_dims", sampler$getPointer())[7L], as.integer(drawn))
  flatK <- CALL("capi_run_every", sampler$getPointer(), 0L, 4L, TRUE)$k
  if (drawn) {
    expect_true(all(is.finite(flatK) & flatK > 0), info = info)
    expect_true(length(unique(flatK)) > 1L, info = info)
  } else {
    expect_identical(flatK, rep(-1, 4L), info = info)
  }
  expect_identical(!is.null(sampler$run(0L, 2L)$k), drawn, info = info)
}
sK <- kFixed()
expectKDrawn(sK, FALSE, "created with k fixed")
sK$setLeafPrior(normal(k = chi(1.5, 2)))
expectKDrawn(sK, TRUE, "setLeafPrior to a k hyperprior")
sK$setLeafPrior(normal(k = 2))
expectKDrawn(sK, FALSE, "setLeafPrior back to k = 2")
sK <- kFixed()
fixedModel <- sK$model
sK$setModel(kDrawn()$model)
expectKDrawn(sK, TRUE, "setModel to a k hyperprior")
sK$setModel(fixedModel)
expectKDrawn(sK, FALSE, "setModel back to k fixed")
sK <- kDrawn()
expectKDrawn(sK, TRUE, "created with a k hyperprior")
rm(sK, fixedModel)
invisible(gc(FALSE))

# ---- the sweep hook: a compiled toy Gibbs step through the chain view ----
# Each chain's next offset and sigma are a fixed function of its own last
# draw (consumer.c, capi_gibbsStep), handed back before each sweep through
# the chain view; chunked, the same step goes through the sampler-wide
# setters between one-sweep runs.
set.seed(5L)
nH <- 80L
xH <- matrix(runif(nH * 3L), nH, 3L)
yH <- 1 + xH[, 1L] + sin(3 * xH[, 2L]) + rnorm(nH, 0, 0.3)
hookSampler <- function(chains, threads, binary = FALSE) {
  dbarts(
    xH,
    if (binary) as.numeric(yH > 1.5) else yH,
    control = dbartsControl(
      n.chains = chains,
      n.threads = threads,
      updateState = FALSE,
      seed = 7L
    ),
    n.trees = 10L
  )
}
gibbsRun <- function(
  sampler,
  numSweeps,
  chunked,
  setsSigma = TRUE,
  stop = c(-1L, -1L, 1L)
) {
  ptr <- sampler$getPointer()
  CALL("capi_gibbs_reset", ptr, setsSigma, stop[1L], stop[2L], stop[3L])
  status <- CALL("capi_set_gibbs_hooks", ptr, if (chunked) 1L else 2L)
  result <- CALL("capi_gibbs_run", ptr, numSweeps, chunked)
  CALL("capi_set_gibbs_hooks", ptr, 0L)
  c(result, status = status, CALL("capi_gibbs_report"))
}

# at one chain the two routes are the same model, bit for bit, on a sampler
# that takes sigma and on one that pins it
for (binary in c(FALSE, TRUE)) {
  hooked <- gibbsRun(hookSampler(1L, 1L, binary), 6L, FALSE, !binary)
  chunked <- gibbsRun(hookSampler(1L, 1L, binary), 6L, TRUE, !binary)
  plain <- CALL(
    "capi_run",
    hookSampler(1L, 1L, binary)$getPointer(),
    0L,
    6L,
    TRUE,
    FALSE
  )
  expect_identical(hooked$status, 1L, info = paste("binary", binary))
  expect_false(hooked$unexpected, info = paste("binary", binary))
  expect_identical(hooked$train, chunked$train, info = paste("binary", binary))
  expect_identical(hooked$sigma, chunked$sigma, info = paste("binary", binary))
  # and the step moves the draws it is checked by
  expect_false(identical(hooked$train, plain$train))
}

# four chains, each on its own offset, the same at one thread and two
fourInline <- gibbsRun(hookSampler(4L, 1L), 6L, FALSE)
fourThreaded <- gibbsRun(hookSampler(4L, 2L), 6L, FALSE)
expect_false(fourThreaded$unexpected)
expect_identical(fourInline$train, fourThreaded$train)
expect_identical(fourInline$sigma, fourThreaded$sigma)

# the view's refusals touch nothing and raise nothing: a non-finite offset,
# a sigma of 0, and on a probit sampler any sigma, which it pins
gaussianProbe <- CALL("capi_view_probe", hookSampler(2L, 2L)$getPointer(), 3L)
expect_identical(gaussianProbe$status, c(1L, -2L, -2L, 0L))
probitProbe <- CALL(
  "capi_view_probe",
  hookSampler(2L, 2L, TRUE)$getPointer(),
  3L
)
expect_identical(probitProbe$status, c(1L, -2L, -1L, -1L))
expect_identical(
  probitProbe$train,
  CALL(
    "capi_run",
    hookSampler(2L, 2L, TRUE)$getPointer(),
    0L,
    3L,
    TRUE,
    FALSE
  )$train
)
# a softmax, whose offset no chain can write, registers nothing
softmaxProbe <- CALL("capi_view_probe", multinomialW()$getPointer(), 2L)
expect_identical(softmaxProbe$status, c(0L, NA, NA, NA))

# a hook's nonzero return stops the run, which returns normally: a positive
# value, a deliberate stop, and a negative one, the hook's failure
for (stopStatus in c(1L, -1L)) {
  stopped <- gibbsRun(
    hookSampler(2L, 2L),
    6L,
    FALSE,
    stop = c(1L, 3L, stopStatus)
  )
  expect_true(stopped$stopped, info = stopStatus)
}

# an interrupt in a loop of one-sweep hooked runs on two threads raises the
# interrupt condition, and a fresh sampler then draws what one did before it
setInterruptAfter <- function(polls) {
  .Call(
    dbarts:::C_dbarts_bartcore_setMonotoneCountHooks,
    NA_real_,
    FALSE,
    as.integer(polls)
  )
}
reference <- gibbsRun(hookSampler(2L, 2L), 6L, FALSE)
interrupted <- local({
  ptr <- hookSampler(2L, 2L)$getPointer()
  CALL("capi_gibbs_reset", ptr, TRUE, -1L, -1L, 1L)
  CALL("capi_set_gibbs_hooks", ptr, 2L)
  on.exit({
    setInterruptAfter(0L)
    CALL("capi_set_gibbs_hooks", ptr, 0L)
  })
  setInterruptAfter(3L)
  tryCatch(
    {
      for (i in 1:5) {
        CALL("capi_gibbs_run", ptr, 1L, FALSE)
      }
      "not interrupted"
    },
    interrupt = function(cond) class(cond)
  )
})
expect_identical(interrupted, c("interrupt", "condition"))
again <- gibbsRun(hookSampler(2L, 2L), 6L, FALSE)
expect_identical(again$train, reference$train)
expect_identical(again$sigma, reference$sigma)

# while chains hold their own offsets a re-derivation of the transform is
# refused on every route; a sampler-wide offset, or a clear, puts every chain
# back on one offset and lifts it
refusal <- "sweep hook wrote"
heldOffsets <- function() {
  sampler <- hookSampler(2L, 2L)
  ptr <- sampler$getPointer()
  CALL("capi_gibbs_reset", ptr, TRUE, -1L, -1L, 1L)
  CALL("capi_set_gibbs_hooks", ptr, 2L)
  CALL("capi_gibbs_run", ptr, 3L, FALSE)
  sampler
}
held <- heldOffsets()
ptr <- held$getPointer()
expect_error(CALL("capi_set_response", ptr, yH, TRUE), refusal)
expect_error(held$setResponse(yH, updateScale = TRUE), refusal)
expect_error(
  held$setWeights(c(0.5, rep(1, nH - 1L)), updateScale = TRUE),
  refusal
)
expect_error(
  held$setActiveRows(c(0, rep(1, nH - 1L)), updateScale = TRUE),
  refusal
)
expect_equal(CALL("capi_set_response", ptr, yH, FALSE), 1L)
expect_equal(CALL("capi_set_offset", ptr, rep(0.1, nH), FALSE), 1L)
expect_equal(CALL("capi_set_response", ptr, yH, TRUE), 1L)
CALL("capi_set_gibbs_hooks", ptr, 0L)
held <- heldOffsets()
ptr <- held$getPointer()
CALL("capi_set_gibbs_hooks", ptr, 0L)
expect_equal(CALL("capi_set_response", ptr, yH, TRUE), 1L)
rm(held, ptr)

# a C++ consumer's exceptions, from its sweep hook and its draw callback,
# inline and on a worker, come back to R as errors carrying their what(), and
# the sampler runs on
throwingDir <- tempfile("capi-throwing")
dir.create(throwingDir)
file.copy(
  system.file("tinytest", "capi", "throwing.cpp", package = "dbarts"),
  file.path(throwingDir, "throwing.cpp")
)
writeLines(
  sprintf('PKG_CPPFLAGS = -I"%s"', includeDir),
  file.path(throwingDir, "Makevars")
)
throwingLib <- compileCapiSource(
  throwingDir,
  "throwing.cpp",
  "the throwing C++ consumer"
)
throwingDll <- dyn.load(throwingLib)
THROW <- function(name, ...) {
  .Call(getNativeSymbolInfo(name, throwingDll), ...)
}
for (threads in 1:2) {
  for (which in 1:2) {
    info <- paste("threads", threads, "callback", which)
    sampler <- hookSampler(2L, threads)
    ptr <- sampler$getPointer()
    THROW("throwing_register", ptr, which)
    expect_error(
      THROW("throwing_run", ptr, 4L),
      if (which == 1L) "the sweep hook threw" else "the draw callback threw",
      info = info
    )
    THROW("throwing_register", ptr, 0L)
    expect_true(CALL("capi_run_plain", ptr, 0L, 2L), info = info)
  }
}
dyn.unload(throwingLib)
rm(sampler, ptr, throwingDll)
invisible(gc(FALSE))
