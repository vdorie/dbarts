# The flat C API (inst/include/dbarts/dbarts.h): compiles a small consumer
# against the installed headers, resolves every entry point through
# R_GetCCallable as a LinkingTo package would, and drives the
# conditional-sampling workout stan4bart performs. Skips wherever the
# consumer cannot be compiled.

consumerSource <- system.file(
  "tinytest",
  "capi",
  "consumer.c",
  package = "dbarts"
)
if (consumerSource == "") {
  exit_file("consumer source not installed")
}

buildDir <- tempfile("capi")
dir.create(buildDir)
file.copy(consumerSource, file.path(buildDir, "consumer.c"))

includeDir <- system.file("include", package = "dbarts")
headerPath <- file.path(includeDir, "dbarts", "dbarts.h")
if (!nzchar(includeDir) || !file.exists(headerPath)) {
  msg <- paste0("dbarts.h not found under includeDir '", includeDir, "'")
  if (nzchar(Sys.getenv("CI", ""))) stop(msg) else exit_file(msg)
}

# system2's env= is not reliably passed through to the child process on
# Windows; a Makevars in the build dir is the portable channel for
# PKG_CPPFLAGS across all platforms including Rtools.
writeLines(
  sprintf('PKG_CPPFLAGS = -I"%s"', includeDir),
  file.path(buildDir, "Makevars")
)
owd <- setwd(buildDir)
compileOutput <- tryCatch(
  suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "SHLIB", "consumer.c"),
    stdout = TRUE,
    stderr = TRUE
  )),
  error = function(e) e
)
setwd(owd)

sharedLib <- file.path(buildDir, paste0("consumer", .Platform$dynlib.ext))
if (!file.exists(sharedLib)) {
  if (nzchar(Sys.getenv("CI", ""))) {
    stop(
      "could not compile the C API consumer under CI:\n",
      paste(compileOutput, collapse = "\n")
    )
  }
  exit_file("could not compile the C API consumer")
}

dll <- dyn.load(sharedLib)
CALL <- function(name, ...) .Call(getNativeSymbolInfo(name, dll), ...)

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
# and it does NOT move for doc text outside DBARTS_C_API_LIST, which the token
# cannot see: this is the literal baked for the three entry points appended
# after the reshape - the dispersion getter and the two augmentation forms
expect_identical(hashes$text, "0x85bd1ef04beb3848")

# the two version components did NOT move: no version of this API has shipped,
# so whatever they read at the first release becomes the initial contract, and
# the hash above is what acknowledges a pre-release change
versions <- CALL("capi_versions")
expect_equal(versions, c(1L, 0L))

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
spec <- dbarts(x, y, test = x.test, control = control)

# creation, queries, and a run into caller-owned buffers
ptr1 <- CALL("capi_create", spec$control, spec$model, spec$data, "")
dims <- CALL("capi_dims", ptr1)
expect_equal(dims, c(n, p, 20L, 1L, 25L, 0L, 0L, 0L))

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
ptr2 <- CALL("capi_create", spec$control, spec$model, spec$data, "")
CALL("capi_sample_trees_from_prior", ptr2)
r2 <- CALL("capi_run", ptr2, 5L, nSamples, TRUE, TRUE)
expect_identical(r1, r2)

# the write guard on the size-first results struct: a caller whose struct
# predates a field (structSize pinned below it) is never written past its
# declared size, even for fields the sampler produces (varcount, test) - a
# size-blind write would crash on the poisoned pointers the canary installs
# past the boundary
ptrGuard <- CALL("capi_create", spec$control, spec$model, spec$data, "")
CALL("capi_sample_trees_from_prior", ptrGuard)
expect_true(CALL("capi_run_guard", ptrGuard, 5L, nSamples))

# the zero-structSize guard: a caller that forgets to set results.structSize
# (leaves it 0) is rejected outright, not silently given an all-skip no-op that
# leaves its buffers uninitialized - the flat-API footgun that garbaged a
# consumer's Gibbs loop
ptrZero <- CALL("capi_create", spec$control, spec$model, spec$data, "")
CALL("capi_sample_trees_from_prior", ptrZero)
expect_error(
  CALL("capi_run_zero_structsize", ptrZero, 5L, nSamples),
  "structSize"
)

# the per-observation log-likelihood channel: for a gaussian sampler it must
# equal dnorm(y, train, sigma, log = TRUE) recomputed on the same draws (train
# is n x nSamples, sigma constant within a draw), pairing observation-fastest
ptrLL <- CALL("capi_create", spec$control, spec$model, spec$data, "")
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
ptrFixed <- CALL(
  "capi_create",
  specFixed$control,
  specFixed$model,
  specFixed$data,
  ""
)
rm(specFixed)
invisible(gc(FALSE)) # creation preserves what the sampler borrows
CALL("capi_set_sigma", ptrFixed, 0.37)
rFixed <- CALL("capi_run", ptrFixed, 0L, 3L, FALSE, FALSE)
expect_equal(unique(rFixed$sigma), 0.37)

# transactional predictor mutation: the replacement column is drawn from a
# fixed local seed so it stays compatible with ptr2's current tree state (a
# fully random column can strand a split and be rejected); the stream is
# restored so later draws in this file are unaffected
xNew <- x
savedSeed <- .Random.seed
set.seed(1L)
xNew[, 3L] <- runif(n)
assign(".Random.seed", savedSeed, envir = globalenv())
rm(savedSeed)
expect_true(CALL("capi_set_predictor", ptr2, xNew))
expect_true(CALL("capi_update_predictor", ptr2, matrix(x[, 3L], n), 2L))

# dbarts_sampler_setPredictor and dbarts_sampler_updatePredictor take a
# heteroscedastic sampler's transactional predictor mutation (forceUpdate = 0)
# and return the engine's verdict, matching the R bridge entries: the
# two-phase revalidation reaches the variance forest, so the entry either
# re-routes every variance tree or rolls the whole change back
nVar <- 60L
xVar <- matrix(runif(nVar * 2L), nVar, 2L)
yVar <- rnorm(nVar)
controlVar <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 8L,
  updateState = FALSE,
  seed = 5L
)
specVar <- dbarts(
  xVar,
  yVar,
  control = controlVar,
  variance = varianceForest(n.trees = 4L)
)
ptrVar <- CALL("capi_create", specVar$control, specVar$model, specVar$data, "")
CALL("capi_sample_trees_from_prior", ptrVar)
xVarNew <- xVar
xVarNew[, 1L] <- runif(nVar)

# the transactional entries answer with the engine's verdict where they used to
# refuse: an accept re-routes every variance tree with the mean ones
expect_true(CALL("capi_set_predictor", ptrVar, xVarNew))
expect_true(is.logical(CALL(
  "capi_update_predictor",
  ptrVar,
  matrix(xVar[, 2L], nVar),
  1L
)))
# and the DECLINE arm: a two-level replacement column empties leaves, so the
# transaction rolls back rather than installing (against the existing grid,
# since refreshing the cuts from a collapsed column is refused for its cut
# count first)
expect_false(CALL(
  "capi_update_predictor_fixed_cuts",
  ptrVar,
  matrix(ifelse(seq_len(nVar) %% 2L == 0L, 0.25, 0.75), nVar),
  0L
))
# ... while the forced flavors stay what they were
expect_true(CALL("capi_set_predictor_forced", ptrVar, xVarNew))
expect_true(CALL(
  "capi_update_predictor_forced",
  ptrVar,
  matrix(xVar[, 2L], nVar),
  1L
))

rm(nVar, xVar, yVar, controlVar, specVar, ptrVar, xVarNew)
invisible(gc(FALSE))

# test-data replacement and removal
CALL("capi_set_test_predictors", ptr2, NULL)
expect_equal(CALL("capi_dims", ptr2)[3L], 0L)
CALL("capi_set_test_predictors", ptr2, x.test)
expect_equal(CALL("capi_dims", ptr2)[3L], 20L)

# the remaining conditioning hooks: a replacement response moves the fits,
# gaussian weights install, and leaf parameters redraw from the prior
yShifted <- y + 10
CALL("capi_set_response", ptr2, yShifted, TRUE)
rShifted <- CALL("capi_run", ptr2, 10L, 3L, TRUE, FALSE)
expect_true(abs(mean(rShifted$train) - mean(yShifted)) < 3)
CALL("capi_set_response", ptr2, y, TRUE)

# the sampler borrows what it is handed, so the vector must stay live
weights <- rep(c(0.5, 1.5), length.out = n)
CALL("capi_set_weights", ptr2, weights)
rWeighted <- CALL("capi_run", ptr2, 10L, 3L, TRUE, FALSE)
expect_true(all(is.finite(rWeighted$train)))

CALL("capi_sample_node_parameters_from_prior", ptr2)
rPriorParams <- CALL("capi_run", ptr2, 0L, 2L, TRUE, FALSE)
expect_true(all(is.finite(rPriorParams$train)))

# thinning, thread, and verbosity controls apply to subsequent runs
CALL("capi_set_run_controls", ptr2, 1L, 2L, FALSE)
rThinned <- CALL("capi_run", ptr2, 0L, 2L, FALSE, FALSE)
expect_equal(length(rThinned$sigma), 2L)
CALL("capi_set_run_controls", ptr2, 1L, 1L, FALSE)

# a live-tree dump goes through the R console without touching state
printed <- capture.output(CALL("capi_print_trees", ptr2, 0L))
expect_true(is.character(printed))

# a test offset adds to the recorded test fits without entering the trees:
# saved-tree replays bit-match the recorded fits (asserted above), so with
# an offset installed the two must differ by exactly the offset
ptrC <- CALL("capi_create", spec$control, spec$model, spec$data, "")
CALL("capi_set_tree_storage", ptrC, TRUE, 2L)
testOffset <- rep(1.5, 20L) # borrowed: must outlive the run below
CALL("capi_set_test_offset", ptrC, testOffset)
rOffsetTest <- CALL("capi_run", ptrC, 3L, 2L, FALSE, TRUE)
predNoOffset <- CALL("capi_predict", ptrC, x.test, NULL)
expect_equal(as.double(rOffsetTest$test), predNoOffset + 1.5)
CALL("capi_set_test_offset", ptrC, NULL)
rm(ptrC, testOffset, weights)
invisible(gc(FALSE))

# latents: absent for gaussian, sign-locked to the response for probit
expect_null(CALL("capi_get_latents", ptr1))

yBinary <- rbinom(n, 1L, pnorm(scale(y)))
specBinary <- dbarts(x, yBinary, control = control)
ptrBinary <- CALL(
  "capi_create",
  specBinary$control,
  specBinary$model,
  specBinary$data,
  "probit"
)
invisible(CALL("capi_run", ptrBinary, 3L, 2L, FALSE, FALSE))
latents <- CALL("capi_get_latents", ptrBinary)
expect_equal(length(latents), n)
expect_true(all(latents[yBinary == 1] > 0))
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
# precision. ptrFixed above is the permitted case - gaussian with
# resid.prior = fixed(), the outer-Gibbs conduit
expect_error(
  CALL("capi_set_sigma", ptrBinary, 0.5),
  "response family fixes the residual standard deviation"
)

# dbarts_sampler_setWeights carries the same family rule too, and used to carry
# none: every family but gaussian refuses a weight change, so this entry
# installed a vector a probit/ordinal/aft/nbinom sampler never reads and took a
# logistic one's negative weight into a division by zero. ptr2 above is the
# permitted case - a gaussian sampler taking weights mid-chain
expect_error(
  CALL("capi_set_weights", ptrBinary, rep(1, n)),
  "probit models do not support case weights"
)

# the nbinom magnitude cap at the flat funnel, which has no R layer ahead of it
# to state the rule: the dispersion grid's count histogram is sized from the
# largest count, so the bound is an allocation bound, and creation and the one
# conduit that swaps y here both carry it. The counts are built without drawing,
# so the stream stays where the seeds above put it
yCount <- as.double(seq_len(n) %% 7L)
specCount <- dbarts(x, yCount, family = "nbinom", control = control)
dataOverBound <- specCount$data
dataOverBound@y <- replace(yCount, 1L, 1000001)
capRefusal <- "counts no larger than 1000000"
expect_error(
  CALL(
    "capi_create",
    specCount$control,
    specCount$model,
    dataOverBound,
    "nbinom"
  ),
  capRefusal
)
ptrCount <- CALL(
  "capi_create",
  specCount$control,
  specCount$model,
  specCount$data,
  "nbinom"
)
expect_error(
  CALL("capi_set_response", ptrCount, dataOverBound@y, FALSE),
  capRefusal
)
# and the weight refusal under a second family, whose message names it
expect_error(
  CALL("capi_set_weights", ptrCount, rep(1, n)),
  "nbinom \\(count\\) models do not support case weights"
)

# the dispersion channel from C: the results slot appended to dbarts_results
# and the mid-sweep getter. The slot is NA-poisoned before the run, so an
# unfilled channel cannot pass for a filled one
disp <- CALL("capi_run_dispersion", ptrCount, 2L, 3L)
expect_true(disp$present)
expect_equal(length(disp$recorded), 3L)
expect_true(all(is.finite(disp$recorded)))
expect_true(all(disp$recorded > 0))
# a caller whose structSize predates the field is never written past, on the
# one family that HAS a dispersion to write
expect_true(disp$guarded)
# the getter reads the same scalar mid-sweep, without serializing state
live <- CALL("capi_dispersion", ptrCount)
expect_equal(length(live), 1L)
expect_equal(live, disp$recorded[3L])
# and it is ABSENT off the family: the entry answers 0, reported here as NULL
expect_null(CALL("capi_dispersion", ptr1))

rm(ptrCount, specCount, dataOverBound)
invisible(gc(FALSE))

# the wrapped augmentation entries, the flat twins of dbartsDrawLatents and
# dbartsWorkingResponse: one implementation behind two surfaces, drawing from
# R's own stream, so under one seed they agree bitwise. The stream is restored
# after the block so later draws in this file are unaffected
drawFlat <- function(family, fit, y, offset = NULL, r = NULL, cuts = NULL) {
  CALL("capi_draw_latents", family, fit, y, NULL, offset, 1, r, cuts, NULL)
}
workFlat <- function(family, latent, y, offset = NULL, r = NULL) {
  CALL("capi_working_response", family, latent, y, NULL, offset, r)
}

augFit <- as.vector(0.4 * scale(y))
augOffset <- rep_len(c(-0.3, 0.2), n)
yDouble <- as.double(yBinary)
augSeed <- .Random.seed
set.seed(11)
flatLatent <- drawFlat("probit", augFit, yDouble, offset = augOffset)
set.seed(11)
rLatent <- dbartsDrawLatents("probit", augFit, yDouble, offset = augOffset)
expect_equal(flatLatent, as.vector(rLatent))

# the two arguments the flat forms carry that the R surface derives for itself:
# a scalar its family's law requires, and the cut points with their own length
yOrdinal <- 1 + (seq_len(n) %% 3L)
cuts <- c(-0.5, 0.7)
set.seed(14)
flatCount <- drawFlat("nbinom", augFit, yCount, r = 3)
flatOrdinal <- drawFlat("ordinal", augFit, yOrdinal, cuts = cuts)
set.seed(14)
rCount <- dbartsDrawLatents("nbinom", augFit, yCount, dispersion = 3)
rOrdinal <- dbartsDrawLatents("ordinal", augFit, yOrdinal, cutpoints = cuts)
expect_equal(flatCount, as.vector(rCount))
expect_equal(flatOrdinal, as.vector(rOrdinal))

# the offset convention is load-bearing and the flat form carries it: 'fit' is
# the location WITHOUT the offset, so drawing at (fit, offset) must reproduce
# drawing at (fit + offset, none) - an entry that dropped the argument would
# agree with neither
set.seed(12)
withOffset <- drawFlat("probit", augFit, yDouble, offset = augOffset)
set.seed(12)
folded <- drawFlat("probit", augFit + augOffset, yDouble)
expect_identical(withOffset, folded)

# the working response, deterministic given the latent: the flat form matches
# the R helper on a PRECISION family, where the quantity is not the latent
set.seed(13)
omega <- as.vector(dbartsDrawLatents(
  "logistic",
  augFit,
  yDouble,
  offset = augOffset
))
expect_equal(
  workFlat("logistic", omega, yDouble, offset = augOffset),
  as.vector(dbartsWorkingResponse(
    "logistic",
    omega,
    yDouble,
    offset = augOffset
  ))
)

# and the flat forms carry the guard rail rather than exporting the sharp edge:
# a law's required parameter has no default, a Polya-Gamma working response
# divides by its latent, and the support rule is the shared one
expect_error(drawFlat("ordinal", augFit, yOrdinal), "requires cut points")
expect_error(workFlat("logistic", rep(0, n), yDouble), "must be positive")
expect_error(drawFlat("probit", augFit, rep(2, n)), "must be coded 0 or 1")
assign(".Random.seed", augSeed, envir = globalenv())
rm(augFit, augOffset, yDouble, yOrdinal, cuts, omega, augSeed)
rm(drawFlat, workFlat)

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

state <- CALL("capi_store_state", ptr1)
expect_inherits(state, "bartcoreState")

ptr3 <- CALL("capi_create", spec$control, spec$model, spec$data, "")
CALL("capi_set_tree_storage", ptr3, TRUE, nSamples)
CALL("capi_set_state", ptr3, state)
pred3 <- CALL("capi_predict", ptr3, x.test, NULL)
expect_identical(pred3, pred1)

# tree introspection: saved trees carry a sample column, live trees do not
trees <- CALL("capi_get_trees", ptr1, FALSE, NULL, 0L)
expect_inherits(trees, "data.frame")
expect_equal(names(trees), c("sample", "tree", "n", "var", "value"))
expect_equal(sort(unique(trees$tree)), 1:25)
expect_true(all(trees$var[trees$var < 0] == -1))
expect_true(all(trees$var <= p))

treesLive <- CALL("capi_get_trees", ptr1, TRUE, NULL, 0L)
expect_equal(names(treesLive), c("tree", "n", "var", "value"))
# every tree's root row sees all observations
expect_true(all(treesLive$n[!duplicated(treesLive$tree)] == n))

# extract-trees-consistency: the all-samples table must equal the per-sample
# gather stan4bart concatenates - in particular the replayed n counts. A stale
# consumer that dropped the getTrees trainingData parameter read an
# indeterminate replay source and returned inconsistent n; the restored 8-arg
# ABI replays the retained creation spec identically for either call shape
perSample <- do.call(
  rbind,
  lapply(
    seq_len(nSamples),
    function(i) CALL("capi_get_trees", ptr1, FALSE, as.integer(i), 0L)
  )
)
row.names(perSample) <- row.names(trees)
expect_equal(trees, perSample)

# a per-sweep callback that sets sigma reproduces the setSigma + run(0, 1)
# loop bitwise (fixed residual prior, so a set sigma is held and recorded)
sigmas <- runif(6L, 0.2, 0.9)
specFixed2 <- dbarts(x, y, resid.prior = fixed(1), control = control)
ptrCbA <- CALL(
  "capi_create",
  specFixed2$control,
  specFixed2$model,
  specFixed2$data,
  ""
)
ptrCbB <- CALL(
  "capi_create",
  specFixed2$control,
  specFixed2$model,
  specFixed2$data,
  ""
)
rCbA <- CALL("capi_run_with_callback", ptrCbA, 0L, 6L, sigmas, -1L)
expect_equal(rCbA$sigma, sigmas)

sigmaB <- numeric(0)
trainB <- numeric(0)
varcountB <- integer(0)
for (i in seq_len(6L)) {
  CALL("capi_set_sigma", ptrCbB, sigmas[i])
  rb <- CALL("capi_run", ptrCbB, 0L, 1L, TRUE, FALSE)
  sigmaB <- c(sigmaB, rb$sigma)
  trainB <- c(trainB, rb$train)
  varcountB <- c(varcountB, rb$varcount)
}
expect_identical(rCbA$sigma, sigmaB)
expect_identical(rCbA$train, trainB)
expect_identical(rCbA$varcount, varcountB)

# returning 0 halts the run early: the callback stops on its third invocation
ptrStop <- CALL(
  "capi_create",
  specFixed2$control,
  specFixed2$model,
  specFixed2$data,
  ""
)
rStop <- CALL("capi_run_with_callback", ptrStop, 0L, 10L, NULL, 3L)
expect_equal(rStop$count, 3L)

# registration is refused while chains would run on worker threads
controlMT <- dbartsControl(
  n.chains = 2L,
  n.threads = 2L,
  n.trees = 25L,
  updateState = FALSE,
  seed = 99L
)
specMT <- dbarts(x, y, control = controlMT)
ptrMT <- CALL("capi_create", specMT$control, specMT$model, specMT$data, "")
expect_error(CALL("capi_run_with_callback", ptrMT, 0L, 2L, NULL, -1L))

# tau/groupEffects on a grouped random-intercept fit, configured through the
# internal control attribute rbart_vi's in-core path sets
numGroups <- 3L
controlG <- spec$control
attr(controlG, "bartcore.groups") <- list(
  indices = as.integer(rep_len(seq_len(numGroups), n)),
  n.groups = numGroups,
  prior = "cauchy",
  rel.scale = sd(y),
  n.steps = 1L
)
ptrG <- CALL("capi_create", controlG, spec$model, spec$data, "")
rG <- CALL("capi_run_grouped", ptrG, 5L, 4L, numGroups)
expect_equal(length(rG$tau), 4L)
expect_true(all(is.finite(rG$tau)) && all(rG$tau > 0))
expect_equal(length(rG$ranef), numGroups * 4L)
expect_true(all(is.finite(rG$ranef)))

# the predictor entries take a design carrying a CSC-backed CATEGORICAL
# column: the engine keys the replacement's nonzero pattern on
# {i : code != refCode}, so a whole-column swap through the flat API lands
# the same codes a dense factor of those values would
if (requireNamespace("Matrix", quietly = TRUE)) {
  set.seed(808L)
  nSparse <- 120L
  levelsSparse <- c("a", "b", "c")
  labelsSparse <- sample(levelsSparse, nSparse, replace = TRUE)
  frameSparse <- data.frame(x1 = rnorm(nSparse))
  frameSparse$f <- dbarts::sparseFactor(
    labelsSparse,
    levels = levelsSparse,
    reference = "b"
  )
  ySparse <- rnorm(nSparse) + match(labelsSparse, levelsSparse)
  specSparse <- dbarts(
    frameSparse,
    ySparse,
    control = dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 10L,
      updateState = FALSE,
      seed = 11L
    )
  )
  ptrSparse <- CALL(
    "capi_create",
    specSparse$control,
    specSparse$model,
    specSparse$data,
    ""
  )
  invisible(CALL("capi_run", ptrSparse, 5L, 3L, FALSE, FALSE))
  codesSparse <- as.matrix(specSparse$data@x)
  storage.mode(codesSparse) <- "double"
  newCodesSparse <- codesSparse
  newCodesSparse[, 2L] <- (codesSparse[, 2L] + 1) %% 3
  expect_true(CALL(
    "capi_update_predictor",
    ptrSparse,
    matrix(newCodesSparse[, 2L], nSparse),
    1L
  ))
  rSparse <- CALL("capi_run", ptrSparse, 0L, 3L, FALSE, FALSE)
  expect_true(all(is.finite(rSparse$sigma)))
  expect_true(CALL("capi_set_predictor", ptrSparse, newCodesSparse))
  rm(ptrSparse)
}

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
  numCscColumns = NULL
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
ptrSrc <- CALL("capi_create", specSrc$control, specSrc$model, specSrc$data, "")
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

# rule 5 at the flat funnel: a reference level is the code a column's IMPLICIT
# rows take, which only a categorical column has, so declaring one against a
# store-ORDINAL column is refused at EVERY entry that takes a source - and 0,
# a legal code, is refused exactly as any other declared value is
for (declared in c(2L, 0L)) {
  trainReference <- makeSource(
    nSrc,
    2L,
    dense = codesSrc[, 2L],
    cscColumns = list(cscColumn(codesSrc[, 1L], 0)),
    map = c(-1L, 0L),
    types = c(0L, 1L),
    refs = c(declared, -1L)
  )
  testReference <- makeSource(
    nTestSrc,
    2L,
    dense = codesTestSrc,
    cscColumns = list(cscColumn(x1TestSrc, 0)),
    map = c(-1L, 0L),
    types = c(0L, 1L),
    refs = c(declared, -1L)
  )
  reason <- "reference level only for a categorical"
  expect_error(
    CALL("capi_set_predictor_source", ptrSrc, trainReference),
    reason
  )
  expect_error(
    CALL("capi_update_predictor_source", ptrSrc, trainReference, c(0L, 1L)),
    reason
  )
  expect_error(
    CALL("capi_set_test_predictors_source", ptrSrc, testReference),
    reason
  )
  expect_error(CALL("capi_predict_source", ptrSrc, testReference), reason)
}

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
wideTrain <- makeSource(nSrc, 3L, dense = cbind(codesSrc, codesSrc[, 1L]))
expect_error(CALL("capi_predict_source", ptrSrc, wideTest), "columns")
expect_error(
  CALL("capi_set_test_predictors_source", ptrSrc, wideTest),
  "columns"
)
expect_error(CALL("capi_set_predictor_source", ptrSrc, wideTrain), "columns")
expect_error(
  CALL("capi_update_predictor_source", ptrSrc, wideTrain, c(0L, 1L)),
  "columns"
)
# a mutation source must also declare the sampler's own row count
expect_error(
  CALL(
    "capi_set_predictor_source",
    ptrSrc,
    makeSource(nSrc - 1L, 2L, dense = codesSrc[-1L, ])
  ),
  "rows"
)
# a column type outside {ordinal, categorical}, a declared level count or a
# reference code past the engine's category limit, and a source naming a CSC
# column it does not carry (without the bound, the ~v decode reads past the
# caller's own pointer array)
expect_error(
  CALL(
    "capi_predict_source",
    ptrSrc,
    makeSource(nTestSrc, 2L, dense = denseTestSrc, types = c(0L, 7L))
  ),
  "columnTypes"
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
# behavior, and an index past the last forest is an error on all three tree
# queries rather than a read past the last forest (the engine's printers index
# their forest unchecked, so the bridge's check is the only guard)
expect_equal(CALL("capi_num_forests", ptrSrc), 1L)
expect_equal(CALL("capi_num_trees", ptrSrc, 0L), 20L)
expect_error(CALL("capi_num_trees", ptrSrc, 1L), "forest index out of range")
expect_error(
  CALL("capi_get_trees", ptrSrc, FALSE, NULL, 1L),
  "forest index out of range"
)
expect_error(
  CALL("capi_print_trees", ptrSrc, 1L),
  "forest index out of range"
)
# both print branches take the index: ptr2 has no tree storage (the live
# trees print), ptrSrc has it (the saved ones)
expect_true(length(capture.output(CALL("capi_print_trees", ptr2, 0L))) > 0L)
expect_true(length(capture.output(CALL("capi_print_trees", ptrSrc, 0L))) > 0L)

# the honesty falsifier: a CSC MUTATION source is materialized exactly as the
# R bridge materializes its own, so it draws the same accept/reject verdict as
# its dense equivalent and the sampler that took it draws bitwise identically.
# Handing the engine a non-dense view instead would return unsupportedSource,
# which the 1/0 mapping would report as an ordinary rollback
ptrMutDense <- CALL(
  "capi_create",
  specSrc$control,
  specSrc$model,
  specSrc$data,
  ""
)
ptrMutCsc <- CALL(
  "capi_create",
  specSrc$control,
  specSrc$model,
  specSrc$data,
  ""
)
newCodesSrc <- codesSrc
newCodesSrc[, 2L] <- (codesSrc[, 2L] + 1) %% 3
verdictDense <- CALL("capi_set_predictor", ptrMutDense, newCodesSrc)
verdictCsc <- CALL(
  "capi_set_predictor_source",
  ptrMutCsc,
  makeSource(
    nSrc,
    2L,
    dense = newCodesSrc[, 1L],
    cscColumns = list(cscColumn(newCodesSrc[, 2L], 1)),
    map = c(0L, -1L),
    types = c(0L, 1L),
    counts = c(0L, 3L),
    refs = c(-1L, 1L)
  )
)
expect_identical(verdictDense, verdictCsc)
expect_true(verdictDense)
expect_identical(
  CALL("capi_run", ptrMutDense, 5L, 3L, TRUE, FALSE),
  CALL("capi_run", ptrMutCsc, 5L, 3L, TRUE, FALSE)
)
rm(ptrMutDense, ptrMutCsc)
invisible(gc(FALSE))

# the leaf-covariate refusal at both flat test entries: a linear leaf reads its
# covariate's contiguous raw values, which CSC storage does not serve. The
# refusal must leave the installed test store INTACT - a discarded refusal
# would report the PREVIOUS rows as the new test set
set.seed(505L)
nLeaf <- 90L
xLeaf <- matrix(rnorm(nLeaf * 2L), nLeaf, 2L)
yLeaf <- 2 * xLeaf[, 1L] + rnorm(nLeaf, 0, 0.3)
controlLeaf <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  updateState = FALSE,
  seed = 23L
)
specLeaf <- dbarts(xLeaf, yLeaf, node.prior = linear(1L), control = controlLeaf)
xTestLeaf <- matrix(rnorm(12L * 2L), 12L, 2L)
sparseLeafSource <- makeSource(
  12L,
  2L,
  dense = xTestLeaf[, 2L],
  cscColumns = list(cscColumn(xTestLeaf[, 1L], 0)),
  map = c(-1L, 0L)
)
ptrLeafA <- CALL(
  "capi_create",
  specLeaf$control,
  specLeaf$model,
  specLeaf$data,
  ""
)
ptrLeafB <- CALL(
  "capi_create",
  specLeaf$control,
  specLeaf$model,
  specLeaf$data,
  ""
)
CALL("capi_set_test_predictors", ptrLeafA, xTestLeaf)
CALL("capi_set_test_predictors", ptrLeafB, xTestLeaf)
reasonLeaf <- "leaf covariate column cannot be a sparse test column"
expect_error(
  CALL("capi_set_test_predictors_source", ptrLeafB, sparseLeafSource),
  reasonLeaf
)
expect_error(
  CALL("capi_predict_source", ptrLeafB, sparseLeafSource),
  reasonLeaf
)
expect_equal(CALL("capi_dims", ptrLeafB)[3L], 12L)
expect_identical(
  CALL("capi_run", ptrLeafA, 10L, 2L, FALSE, TRUE)$test,
  CALL("capi_run", ptrLeafB, 10L, 2L, FALSE, TRUE)$test
)
rm(ptrLeafA, ptrLeafB)
invisible(gc(FALSE))

# the calibration surface from C: the reader fills only the members the caller
# both carries (by structSize) and points somewhere, and the writer answers a
# capability with a return value while a malformed value raises
calibration <- CALL("capi_forest_calibration", ptrSrc, 0L, 0L, FALSE)
expect_equal(calibration$accepted, 1L)
expect_true(all(calibration$prior.scale > 0))
expect_equal(calibration$prior.sd, calibration$prior.scale / calibration$k)
expect_true(all(calibration$response.scale > 0))
expect_equal(calibration$k.has.hyperprior, 0L)
expect_equal(calibration$leaf.model, 0L) # DBARTS_LEAF_CONSTANT
# the five calibration-map fields are filled, and on a single-forest sampler
# what they carry is NaN - the same "not map-derived" signal the R matrix gives
mapFields <- c(
  "amplitude.prior.variance",
  "amplitude.prior.scale",
  "node.scale.factor",
  "node.scale.divisor",
  "basis.row.norm"
)
for (field in mapFields) {
  expect_true(all(is.nan(calibration[[field]])), info = field)
}
# an index past the last forest touches nothing and answers 0
expect_equal(
  CALL("capi_forest_calibration", ptrSrc, 1L, 0L, FALSE)$accepted,
  0L
)
# an omitting caller (structSize below leafModel, that member poisoned) and a
# null member are both skipped, and everything else still fills
calibrationPartial <- CALL("capi_forest_calibration", ptrSrc, 0L, 1L, TRUE)
expect_equal(calibrationPartial$accepted, 1L)
expect_equal(calibrationPartial$leaf.model, -1L)
expect_equal(calibrationPartial$k, -1)
expect_identical(calibrationPartial$prior.scale, calibration$prior.scale)
# a PRE-APPEND caller: structSize pinned at the calibration-map boundary, the
# five appended members poisoned. It still reads all EIGHT original fields, and
# its five buffers keep the sentinel they went in with - which is what
# DBARTS_HAS_FIELD promises a consumer compiled against the older struct, and
# what lets sizeof(dbarts_forest_calibration) move without a source break.
calibrationPreAppend <- CALL("capi_forest_calibration", ptrSrc, 0L, 2L, FALSE)
expect_equal(calibrationPreAppend$accepted, 1L)
expect_identical(calibrationPreAppend$prior.scale, calibration$prior.scale)
expect_identical(calibrationPreAppend$prior.sd, calibration$prior.sd)
expect_identical(calibrationPreAppend$prior.mean, calibration$prior.mean)
expect_identical(calibrationPreAppend$k, calibration$k)
expect_identical(
  calibrationPreAppend$response.scale,
  calibration$response.scale
)
expect_identical(
  calibrationPreAppend$response.shift,
  calibration$response.shift
)
expect_identical(
  calibrationPreAppend$k.has.hyperprior,
  calibration$k.has.hyperprior
)
expect_identical(calibrationPreAppend$leaf.model, calibration$leaf.model)
for (field in mapFields) {
  expect_equal(
    calibrationPreAppend[[field]],
    rep_len(-1, length(calibration[[field]])),
    info = field
  )
}
expect_error(
  CALL("capi_forest_calibration_zero_structsize", ptrSrc),
  "structSize"
)
expect_equal(CALL("capi_set_forest_prior_scale", ptrSrc, 0L, 2.5), 1L)
expect_equal(
  CALL("capi_forest_calibration", ptrSrc, 0L, 0L, FALSE)$prior.scale,
  2.5,
  tolerance = 1e-12
)
# a read-then-write is inert, an index past the last forest refuses without
# raising, and a malformed scale raises
expect_equal(CALL("capi_set_forest_prior_scale", ptrSrc, 0L, 2.5), 1L)
expect_equal(CALL("capi_set_forest_prior_scale", ptrSrc, 1L, 2.5), 0L)
expect_error(
  CALL("capi_set_forest_prior_scale", ptrSrc, 0L, -1),
  "positive finite"
)
expect_error(
  CALL("capi_set_forest_prior_scale", ptrSrc, 0L, NaN),
  "positive finite"
)
rm(ptrSrc)
invisible(gc(FALSE))

# the active-row mask from C: a masked row leaves the data set for the sweep,
# an all-ones mask installs nothing, a fractional value is refused (0, not a
# raise), and NULL clears. The probe is a capability, never a family switch,
# so a probit sampler takes one too, and a genuine mask must move ITS draws
# exactly as it does the gaussian ones
ptrMaskA <- CALL("capi_create", spec$control, spec$model, spec$data, "")
ptrMaskB <- CALL("capi_create", spec$control, spec$model, spec$data, "")
ptrMaskC <- CALL("capi_create", spec$control, spec$model, spec$data, "")
expect_equal(
  CALL("capi_set_active_rows", ptrMaskA, as.double(seq_len(n) %% 2L)),
  1L
)
rMaskA <- CALL("capi_run", ptrMaskA, 5L, 3L, TRUE, FALSE)
rMaskB <- CALL("capi_run", ptrMaskB, 5L, 3L, TRUE, FALSE)
expect_false(identical(rMaskA$train, rMaskB$train))
expect_equal(CALL("capi_set_active_rows", ptrMaskC, rep(1, n)), 1L)
expect_identical(CALL("capi_run", ptrMaskC, 5L, 3L, TRUE, FALSE), rMaskB)
expect_equal(CALL("capi_set_active_rows", ptrMaskA, rep(0.5, n)), 0L)
expect_equal(CALL("capi_set_active_rows", ptrMaskA, NULL), 1L)
ptrProbitMaskA <- CALL(
  "capi_create",
  specBinary$control,
  specBinary$model,
  specBinary$data,
  "probit"
)
ptrProbitMaskB <- CALL(
  "capi_create",
  specBinary$control,
  specBinary$model,
  specBinary$data,
  "probit"
)
expect_equal(
  CALL("capi_set_active_rows", ptrProbitMaskA, as.double(seq_len(n) %% 2L)),
  1L
)
rProbitMaskA <- CALL("capi_run", ptrProbitMaskA, 5L, 3L, TRUE, FALSE)
rProbitMaskB <- CALL("capi_run", ptrProbitMaskB, 5L, 3L, TRUE, FALSE)
expect_false(identical(rProbitMaskA$train, rProbitMaskB$train))
rm(ptrMaskA, ptrMaskB, ptrMaskC, ptrProbitMaskA, ptrProbitMaskB)
invisible(gc(FALSE))

# the two-forest (BCF) surface, docs/design/bcf.md: a treatment vector on the
# data object and the treatment forest's configuration on the control make a
# Bayesian causal forest through this same creation entry point. The whole
# mutation and reporting surface is driven from the consumer - the acceptances,
# the refusals, and the reason each refusal names - so what is checked is the
# flat API's own agreement with the R bridge rather than an R restatement of it.
# The response-side legs also pin the one branch of that shared guard this
# surface cannot reach: a coupling whose response is its own count matrix is
# refused with a message naming the counts channel, conditioned on the
# capability and not on the forest count, so a BCF sampler - which opts into
# the response conduit - must still be ACCEPTED at updateScale = FALSE and must
# still refuse with the generic wording at TRUE
set.seed(21L)
nBCF <- 200L
xBCF <- matrix(runif(nBCF * 3L), nBCF, 3L)
zBCF <- as.double(rbinom(nBCF, 1L, 0.5))
yBCF <- 2 * xBCF[, 1L] + zBCF * (1 + xBCF[, 2L]) + rnorm(nBCF, 0, 0.2)
specBCF <- dbartsSpec(
  dbartsData(xBCF, yBCF),
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.trees = 25L,
    n.samples = 3L,
    updateState = FALSE,
    seed = 41L
  ),
  forests = list(forest(), forest(basis = ~ factor(zBCF), n.trees = 15L))
)
ptrBCF <- CALL("capi_create", specBCF$control, specBCF$model, specBCF$data, "")
expect_equal(CALL("capi_dims", ptrBCF)[5L], 25L)
invisible(CALL("capi_run", ptrBCF, 5L, 3L, FALSE, FALSE))

# borrowed until replaced, so these must outlive the sampler below
offsetBCF <- rep(0, nBCF)
weightsBCF <- rep(1, nBCF)
bcfLegs <- CALL(
  "capi_bcf_surface",
  ptrBCF,
  yBCF,
  offsetBCF,
  weightsBCF,
  zBCF,
  xBCF[1:5, , drop = FALSE]
)
expect_equal(length(bcfLegs), 19L)
for (leg in names(bcfLegs)) {
  expect_true(bcfLegs[[leg]], info = leg)
}

# the accepted mutations leave the sampler running
rBCF <- CALL("capi_run", ptrBCF, 0L, 3L, TRUE, FALSE)
expect_true(all(is.finite(rBCF$train)))

# the varcount contract on a two-forest sampler, from the consumer's side.
# dbarts_results declares no forest count, so the engine writes exactly the
# numPredictors x numSamples x numChains slab the consumer sized - the widened
# per-forest channel belongs to the R run bridge, which declares one - and the
# bytes are the PROGNOSTIC forest's. The length alone cannot falsify that (the
# consumer computes it from p, samples and chains), so the value is pinned
# against the same model's live per-forest read at the same seed: a fresh flat
# sampler and a fresh R-level sampler over the identical (control, model, data),
# each run 5 burn-in and 3 kept draws, must agree on forest 1's counts
pBCF <- ncol(xBCF)
expect_equal(length(rBCF$varcount), pBCF * 3L * 2L)

ptrBCFCounts <- CALL(
  "capi_create",
  specBCF$control,
  specBCF$model,
  specBCF$data,
  ""
)
rBCFCounts <- CALL("capi_run", ptrBCFCounts, 5L, 3L, FALSE, FALSE)
r5BCF <- new("dbartsSampler", specBCF$control, specBCF$model, specBCF$data)
r5Counts <- r5BCF$run(5L, 3L)
flatCounts <- matrix(rBCFCounts$varcount, pBCF)
for (chain in seq_len(2L)) {
  expect_equal(
    flatCounts[, chain * 3L],
    as.integer(r5BCF$getForestVariableCounts(1L)[, chain])
  )
}
# and the whole flat slab is the R-level channel's forest-1 slab, so what the flat
# caller gets is a projection of the widened channel rather than a third thing
expect_equal(
  as.integer(flatCounts),
  as.integer(r5Counts$varcount[, 1L, , ])
)

# forest addressing on a two-forest sampler, which is what the index is FOR:
# the tree queries name the forest they read instead of silently answering for
# forest 0, and each reads its own forest's tree count
expect_equal(CALL("capi_num_forests", ptrBCF), 2L)
expect_equal(CALL("capi_num_trees", ptrBCF, 0L), 25L)
expect_equal(CALL("capi_num_trees", ptrBCF, 1L), 15L)
expect_error(CALL("capi_num_trees", ptrBCF, 2L), "forest index out of range")
treesMu <- CALL("capi_get_trees", ptrBCF, TRUE, NULL, 0L)
treesTau <- CALL("capi_get_trees", ptrBCF, TRUE, NULL, 1L)
expect_equal(sort(unique(treesMu$tree)), 1:25)
expect_equal(sort(unique(treesTau$tree)), 1:15)
expect_error(
  CALL("capi_get_trees", ptrBCF, TRUE, NULL, 2L),
  "forest index out of range"
)
# the flat forwarding must actually carry the forest index through to the
# engine's printer rather than silently answering for forest 0: the live
# (unstored) branch here, exercised through ptrBCF, which never enabled tree
# storage, and the saved branch below, on a holder that did. Comparing
# lengths would not falsify a forest-0 pin (both forests have live trees to
# print), so this is a content pin - forest 0 (mu, all 3 predictors) and
# forest 1 (tau, the z-factor basis) print different splits
printBCFLive0 <- capture.output(CALL("capi_print_trees", ptrBCF, 0L))
printBCFLive1 <- capture.output(CALL("capi_print_trees", ptrBCF, 1L))
expect_true(length(printBCFLive1) > 0L)
expect_false(identical(printBCFLive0, printBCFLive1))
expect_error(
  CALL("capi_print_trees", ptrBCF, 2L),
  "forest index out of range"
)

# the same forwarding, through the saved-tree branch this time: a BCF holder
# with tree storage enabled and a run behind it takes printSavedTree instead
# of printTree (bartcore/sampler.hpp printTrees dispatches on keepTrees), a
# separate code path the live-branch pin above cannot reach
ptrBCFSaved <- CALL(
  "capi_create",
  specBCF$control,
  specBCF$model,
  specBCF$data,
  ""
)
CALL("capi_set_tree_storage", ptrBCFSaved, TRUE, 2L)
invisible(CALL("capi_run", ptrBCFSaved, 5L, 2L, FALSE, FALSE))
expect_true(CALL("capi_dims", ptrBCFSaved)[6L] > 0L)
printBCFSaved0 <- capture.output(CALL("capi_print_trees", ptrBCFSaved, 0L))
printBCFSaved1 <- capture.output(CALL("capi_print_trees", ptrBCFSaved, 1L))
expect_true(length(printBCFSaved0) > 0L)
expect_true(length(printBCFSaved1) > 0L)
expect_false(identical(printBCFSaved0, printBCFSaved1))
rm(ptrBCFSaved)
invisible(gc(FALSE))

# the per-forest precision weight from C: 1 = accepted, 0 = refused - the
# shipped polarity, which two landed plans recorded inverted - and an accepted
# non-degenerate weight really installs, while an all-ones one is bitwise inert
# and a refused call leaves the sampler where it was. The weights are BORROWED,
# so these vectors outlive the samplers below
weightsForest <- rep(c(0.25, 1.75), length.out = nBCF)
onesForest <- rep(1, nBCF)
ptrW1 <- CALL("capi_create", specBCF$control, specBCF$model, specBCF$data, "")
ptrW2 <- CALL("capi_create", specBCF$control, specBCF$model, specBCF$data, "")
ptrW3 <- CALL("capi_create", specBCF$control, specBCF$model, specBCF$data, "")
expect_equal(CALL("capi_set_forest_weights", ptrW1, 1L, weightsForest), 1L)
expect_equal(CALL("capi_set_forest_weights", ptrW3, 1L, onesForest), 1L)
rW1 <- CALL("capi_run", ptrW1, 5L, 3L, TRUE, FALSE)
rW2 <- CALL("capi_run", ptrW2, 5L, 3L, TRUE, FALSE)
rW3 <- CALL("capi_run", ptrW3, 5L, 3L, TRUE, FALSE)
expect_false(identical(rW1$train, rW2$train))
expect_identical(rW3, rW2)
# the refusals: a forest past the last one, and a sampler whose coupling
# admits no such weight at all
expect_equal(CALL("capi_set_forest_weights", ptrW2, 2L, weightsForest), 0L)
expect_equal(CALL("capi_set_forest_weights", ptr1, 0L, rep(1, n)), 0L)
# and the refused call moved nothing: the next draws still agree with the
# sampler that was never asked
expect_identical(
  CALL("capi_run", ptrW2, 0L, 3L, TRUE, FALSE),
  CALL("capi_run", ptrW3, 0L, 3L, TRUE, FALSE)
)

# the amplitude pair from C, ragged by construction: the basis forest carries
# two amplitudes and the intercept forest one, and the values agree with the
# R-level reader on the same draws
amplitudes0 <- CALL("capi_forest_amplitudes", ptrW2, 0L)
amplitudes1 <- CALL("capi_forest_amplitudes", ptrW2, 1L)
expect_equal(amplitudes0$count, 1L)
expect_equal(amplitudes1$count, 2L)
expect_equal(amplitudes0$accepted, 1L)
expect_equal(length(amplitudes1$values), 2L * 2L) # two amplitudes, two chains
expect_true(all(is.finite(c(amplitudes0$values, amplitudes1$values))))
expect_error(
  CALL("capi_forest_amplitudes", ptrW2, 2L),
  "forest index out of range"
)

# destruction runs through the finalizer
rm(ptr1, ptr2, ptr3, ptrBinary, ptrFixed, ptrLL)
rm(ptrCbA, ptrCbB, ptrStop, ptrMT, ptrG, ptrBCF, ptrW1, ptrW2, ptrW3)
invisible(gc(FALSE))
rm(offsetBCF, weightsBCF, yBCF, zBCF, weightsForest, onesForest)
