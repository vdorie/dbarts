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
owd <- setwd(buildDir)
compileOutput <- tryCatch(
  suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "SHLIB", "consumer.c"),
    env = paste0("PKG_CPPFLAGS=-I", shQuote(includeDir)),
    stdout = TRUE,
    stderr = TRUE
  )),
  error = function(e) e
)
setwd(owd)

sharedLib <- file.path(buildDir, paste0("consumer", .Platform$dynlib.ext))
if (!file.exists(sharedLib)) {
  exit_file("could not compile the C API consumer")
}

dll <- dyn.load(sharedLib)
CALL <- function(name, ...) .Call(getNativeSymbolInfo(name, dll), ...)

expect_equal(CALL("capi_version"), 1L)

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
  rngSeed = 99L
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

# test-data replacement and removal
CALL("capi_set_test_predictors", ptr2, NULL)
expect_equal(CALL("capi_dims", ptr2)[3L], 0L)
CALL("capi_set_test_predictors", ptr2, x.test)
expect_equal(CALL("capi_dims", ptr2)[3L], 20L)

# the remaining conditioning hooks: a replacement response moves the fits,
# gaussian weights install, and leaf parameters redraw from the prior
yShifted <- y + 10
CALL("capi_set_response", ptr2, yShifted)
rShifted <- CALL("capi_run", ptr2, 10L, 3L, TRUE, FALSE)
expect_true(abs(mean(rShifted$train) - mean(yShifted)) < 3)
CALL("capi_set_response", ptr2, y)

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
printed <- capture.output(CALL("capi_print_trees", ptr2))
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
trees <- CALL("capi_get_trees", ptr1, FALSE)
expect_inherits(trees, "data.frame")
expect_equal(names(trees), c("sample", "tree", "n", "var", "value"))
expect_equal(sort(unique(trees$tree)), 1:25)
expect_true(all(trees$var[trees$var < 0] == -1))
expect_true(all(trees$var <= p))

treesLive <- CALL("capi_get_trees", ptr1, TRUE)
expect_equal(names(treesLive), c("tree", "n", "var", "value"))
# every tree's root row sees all observations
expect_true(all(treesLive$n[!duplicated(treesLive$tree)] == n))

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
  rngSeed = 99L
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

# destruction runs through the finalizer
rm(ptr1, ptr2, ptr3, ptrBinary, ptrFixed)
rm(ptrCbA, ptrCbB, ptrStop, ptrMT, ptrG)
invisible(gc(FALSE))
expect_true(TRUE)
