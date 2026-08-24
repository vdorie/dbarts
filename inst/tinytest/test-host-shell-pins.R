# The dbartsSampler host-shell guard (docs/design/multinomial-mutation-arc.md
# section 1): today's census of hostFor-guarded vs unguarded methods, the
# placeholder host's degenerate predict, its K-dependent resolved family, and
# the save/reload failure the three K-forest/host-shell families
# (multinomial, ordinal, nbinom) share. These assertions exist to MOVE - a
# later change that gives one of these families a real public surface should
# invert or delete the corresponding one here, not leave it standing
# unnoticed.

# --- the census: every own method of dbartsSampler, classified by which
# guard (if either) its body calls first. Read the SOURCE rather than probe
# each of the 22 mutators by call, since most need contrived arguments this
# file has no stake in constructing; a spot check below confirms the guard
# still fires at runtime, not just in text.
gen <- getRefClass("dbartsSampler")
ownMethods <- gen$def@refMethods
infrastructure <- c(
  "initialize",
  "adoptPointer",
  "refuseHostMutation",
  "refuseHostRead",
  "reapplyForestWeights",
  "getPointer",
  "show"
)
# envRefClass's own inherited machinery, not part of dbartsSampler's surface
inherited <- c(
  ".objectPackage",
  ".objectParent",
  "callSuper",
  "copy#envRefClass",
  "export",
  "field",
  "getClass",
  "getRefClass",
  "import",
  "initFields",
  "show#envRefClass",
  "trace",
  "untrace",
  "usingMethods"
)
own <- setdiff(names(ownMethods), inherited)

classifyGuard <- function(name) {
  text <- paste(deparse(body(ownMethods[[name]])), collapse = "\n")
  if (grepl("refuseHostMutation(", text, fixed = TRUE)) {
    return("mutation")
  }
  if (grepl("refuseHostRead(", text, fixed = TRUE)) {
    return("read")
  }
  "unguarded"
}
guard <- vapply(own, classifyGuard, character(1L))

hostMutationMethods <- sort(own[guard == "mutation"])
hostReadMethods <- sort(own[guard == "read"])
substantiveUnguardedMethods <- sort(setdiff(
  own[guard == "unguarded"],
  infrastructure
))

expect_equal(length(own), 48L)
expect_equal(length(hostMutationMethods), 25L)
expect_equal(length(hostReadMethods), 2L)
expect_equal(length(substantiveUnguardedMethods), 14L)
expect_equal(length(infrastructure), 7L)

expect_identical(
  hostMutationMethods,
  sort(c(
    "run",
    "sampleTreesFromPrior",
    "sampleNodeParametersFromPrior",
    "growFromRoot",
    "setControl",
    "setModel",
    "setData",
    "setResponse",
    "setOffset",
    "setWeights",
    "setCounts",
    "setCategoryOffset",
    "setCategoryTestOffset",
    "setActiveRows",
    "setForestWeights",
    "setForestBasis",
    "setSigma",
    "setPredictor",
    "setCutPoints",
    "setTestPredictor",
    "setTestPredictorAndOffset",
    "setTestOffset",
    "setCalibration",
    "setState",
    "installTrees"
  ))
)
expect_identical(
  hostReadMethods,
  sort(c("getDispersion", "getFitsWithoutOffset"))
)
expect_identical(
  substantiveUnguardedMethods,
  sort(c(
    "copy",
    "predict",
    "predictForests",
    "getLatents",
    "getSigmas",
    "getSumsOfSquaredResiduals",
    "getForestFits",
    "getForestAmplitudes",
    "getForestVariableCounts",
    "getCalibration",
    "storeState",
    "printTrees",
    "getTrees",
    "plotTree"
  ))
)

# --- spot check: the guard fires at runtime, on an actual host shell,
# consistent with the source-level classification above ---
set.seed(94011)
nSpot <- 60L
xSpot <- matrix(runif(nSpot * 2L), nSpot, 2L)
ySpot <- factor(
  c("a", "b", "c")[1L + (xSpot[, 1L] > 0.5) + (xSpot[, 1L] > 0.8)],
  levels = c("a", "b", "c")
)
fitSpot <- bart2(
  xSpot,
  ySpot,
  family = "multinomial",
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = 2L,
  n.samples = 2L,
  verbose = FALSE,
  keepTrees = TRUE
)
hostRefusalSpot <- "host sampler of a bart2\\(family = \"multinomial\"\\) fit"
expect_error(fitSpot$fit$setResponse(rnorm(nSpot)), hostRefusalSpot)
expect_error(fitSpot$fit$getDispersion(), hostRefusalSpot)
# unguarded: answers (however uselessly), never the host-shell refusal
predSpot <- fitSpot$fit$predict(xSpot)
expect_true(is.numeric(predSpot))

# --- fixtures for the remaining pins: one multinomial K = 3, one K = 2, one
# ordinal, one nbinom, each with a host shell retained via keepTrees ---
n.trees <- 15L
n.burn <- 15L
n.samples <- 15L

set.seed(9401)
n3 <- 120L
p3 <- 3L
x3 <- matrix(runif(n3 * p3), n3, p3)
eta3 <- cbind(
  2 * (x3[, 1L] - 0.5),
  x3[, 2L] - x3[, 3L],
  1.2 * (x3[, 1L] - x3[, 3L])
)
probs3 <- exp(eta3) / rowSums(exp(eta3))
labels3 <- vapply(
  seq_len(n3),
  function(i) sample.int(3L, 1L, prob = probs3[i, ]) - 1L,
  integer(1L)
)
y3 <- factor(c("a", "b", "c")[labels3 + 1L], levels = c("a", "b", "c"))
set.seed(9402)
fit3 <- bart2(
  x3,
  y3,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE,
  keepTrees = TRUE
)

set.seed(9403)
n2 <- 100L
p2 <- 2L
x2 <- matrix(runif(n2 * p2), n2, p2)
labels2 <- rbinom(n2, 1L, plogis(2 * (x2[, 1L] - 0.5) + x2[, 2L]))
y2 <- factor(c("no", "yes")[labels2 + 1L], levels = c("no", "yes"))
set.seed(9404)
fit2 <- bart2(
  x2,
  y2,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE,
  keepTrees = TRUE
)

set.seed(99)
nOrd <- 80L
xOrd <- matrix(runif(nOrd * 3L), nOrd, 3L)
zOrd <- 2 * (xOrd[, 1L] - 0.5) + rnorm(nOrd)
codesOrd <- 1L + (zOrd > 0) + (zOrd > 0.8)
yOrd <- ordered(c("lo", "mid", "hi")[codesOrd], levels = c("lo", "mid", "hi"))
fitOrd <- bart2(
  xOrd,
  yOrd,
  family = "ordinal",
  n.trees = 10L,
  n.chains = 1L,
  n.burn = 10L,
  n.samples = 10L,
  verbose = FALSE,
  keepTrees = TRUE
)

set.seed(199)
nNb <- 80L
xNb <- matrix(runif(nNb * 3L), nNb, 3L)
etaNb <- 0.9 * (xNb[, 1L] - 0.5)
yNb <- rnbinom(nNb, size = 5L, mu = exp(etaNb))
fitNb <- bart2(
  xNb,
  yNb,
  family = "nbinom",
  n.trees = 10L,
  n.chains = 1L,
  n.burn = 10L,
  n.samples = 10L,
  verbose = FALSE,
  keepTrees = TRUE
)

# --- $copy() no longer launders the host-shell guard: hostFor transfers,
# so the copy refuses exactly as the original does ---
expect_true(length(fit3$fit$hostFor) == 1L)
dupe3 <- fit3$fit$copy()
expect_identical(dupe3$hostFor, fit3$fit$hostFor)
expect_error(dupe3$setResponse(rnorm(n3)), hostRefusalSpot)
expect_error(dupe3$run(0L, 1L), hostRefusalSpot)

# --- ordinal's $fit is no longer a host shell (pointer adoption): hostFor
# stays empty, on the original and its copy alike, and a mutation through
# either succeeds ---
expect_true(length(fitOrd$fit$hostFor) == 0L)
dupeOrd <- fitOrd$fit$copy(shallow = TRUE)
expect_identical(dupeOrd$hostFor, fitOrd$fit$hostFor)
expect_silent(dupeOrd$setData(dbartsData(xOrd, as.double(codesOrd))))

# --- the K = 3 host is created and never run: predict answers a single,
# constant value with no error - a plausible-looking number, not a warning ---
predHost3 <- fit3$fit$predict(x3)
expect_true(is.numeric(predHost3))
expect_equal(length(unique(as.vector(predHost3))), 1L)

# --- the resolved host family is K-dependent, and only reachable under an
# explicit family = "multinomial" token: K = 2 is probit, K >= 3 is
# gaussian, and family = "auto" on a 2-level factor never reaches either -
# it announces probit directly and returns class "bart" ---
expect_equal(fit3$fit$model@family, "gaussian")
expect_false(fit3$fit$control@binary)
expect_equal(fit2$fit$model@family, "probit")
expect_true(fit2$fit$control@binary)

autoWarnings <- 0L
fitAuto <- withCallingHandlers(
  bart2(
    x2,
    y2,
    family = "auto",
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    n.burn = 2L,
    n.samples = 2L,
    verbose = FALSE
  ),
  message = function(m) {
    autoWarnings <<- autoWarnings + 1L
    invokeRestart("muffleMessage")
  }
)
expect_equal(autoWarnings, 1L)
expect_identical(class(fitAuto), "bart")
expect_false(inherits(fitAuto, "bartMultinomial"))

# --- save/reload: multinomial's $fit is still a host shell, so predict
# still fails the same way and the documented bart escape
# ($fit$storeState(), man/bart.Rd) still does not help - predict reaches
# through $bc, which storeState never touches ---
pinSaveReloadFailure <- function(fit, xtest) {
  fit$fit$storeState()
  tempFile <- tempfile()
  saveRDS(fit, tempFile)
  reloaded <- readRDS(tempFile)
  unlink(tempFile)
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_error(predict(reloaded, xtest), "NULL external pointer")
  reloaded$fit$storeState()
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_error(predict(reloaded, xtest), "NULL external pointer")
}
pinSaveReloadFailure(fit3, x3[1:5, ])

# --- ordinal/nbinom's $fit is the engine that ran (pointer adoption), so
# predict now routes through it and getPointer's re-creation branch rebuilds
# the engine from the stored state: save/reload/predict round-trips ---
checkSaveReloadPredicts <- function(fit, xtest) {
  before <- predict(fit, xtest)
  fit$fit$storeState()
  tempFile <- tempfile()
  saveRDS(fit, tempFile)
  reloaded <- readRDS(tempFile)
  unlink(tempFile)
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_equal(predict(reloaded, xtest), before)
}
checkSaveReloadPredicts(fitOrd, xOrd[1:5, ])
checkSaveReloadPredicts(fitNb, xNb[1:5, ])
