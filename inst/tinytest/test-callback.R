# The R surface of the per-draw callback: 'callback' on bart(), dbarts() and
# the sampler's run method, the paired 'keepFits' control slot, its automatic
# default, the four family refusals, and what a keepFits = FALSE fit costs
# the returned object. No compiled callback is registered here - a fake but
# genuinely-typed external pointer stands in wherever a shape check is what's
# under test, and every run in this file either never fires a callback
# (samplerOnly, or no callback at all) or runs with keepFits = FALSE and NO
# callback, so nothing here ever dereferences a bogus function address; the
# compiled counting consumer that actually fires one lives in test-capi.R.

set.seed(7, sample.kind = "Rejection")
n <- 24L
x <- matrix(runif(n * 2L), n, 2L)
y <- x[, 1L] + rnorm(n, 0, 0.2)

# a real, but non-function, external pointer: fine wherever the callback is
# validated-but-never-invoked (typeof() is all 'callback' checks)
fakePtr <- dbarts(x, y)$getPointer()
expect_equal(typeof(fakePtr), "externalptr")
validCallback <- list(fn = fakePtr, context = NULL)

# ---- argument validation, on all three surfaces ----

expect_error(
  bart(x, y, callback = "not a list", samplerOnly = TRUE),
  "'callback'"
)
expect_error(
  bart(x, y, callback = list(context = NULL), samplerOnly = TRUE),
  "'callback'"
)
expect_error(
  bart(x, y, callback = list(fn = 1, context = NULL), samplerOnly = TRUE),
  "'callback\\$fn'"
)
expect_error(
  bart(
    x,
    y,
    callback = list(fn = fakePtr, context = "nope"),
    samplerOnly = TRUE
  ),
  "'callback\\$context'"
)
# NULL fn inside a supplied list is refused too - "no callback" is spelled
# callback = NULL at the top, not a list with a null fn
expect_error(
  bart(x, y, callback = list(fn = NULL, context = NULL), samplerOnly = TRUE),
  "'callback\\$fn'"
)

expect_error(
  dbarts(x, y, callback = list(fn = 1, context = NULL)),
  "'callback\\$fn'"
)

samplerForRun <- dbarts(
  x,
  y,
  control = dbartsControl(n.chains = 1L, n.threads = 1L)
)
expect_error(
  samplerForRun$run(0L, 1L, callback = list(fn = 1, context = NULL)),
  "'callback\\$fn'"
)

# ---- the control slot: validity ----

expect_true(methods::new("dbartsControl")@keepFits)
expect_error(dbartsControl(keepFits = NA), "'keepFits'")
expect_error(dbartsControl(keepFits = c(TRUE, FALSE)), "'keepFits'")
expect_true(dbartsControl(keepFits = TRUE)@keepFits)
expect_false(dbartsControl(keepFits = FALSE)@keepFits)

# ---- the automatic default: callback supplied, keepFits unnamed -> FALSE;
# an explicit keepFits always wins; no callback -> TRUE regardless ----

specNoCallback <- bart(x, y, samplerOnly = TRUE)
expect_true(specNoCallback$control@keepFits)

specAuto <- bart(x, y, callback = validCallback, samplerOnly = TRUE)
expect_false(specAuto$control@keepFits)

specExplicitTrue <- bart(
  x,
  y,
  callback = validCallback,
  keepFits = TRUE,
  samplerOnly = TRUE
)
expect_true(specExplicitTrue$control@keepFits)

specExplicitFalse <- bart(x, y, keepFits = FALSE, samplerOnly = TRUE)
expect_false(specExplicitFalse$control@keepFits)

# dbartsControl()'s own default is unaffected by anything - it has no
# 'callback' to react to
expect_true(dbartsControl()@keepFits)

# ---- the four family refusals: automatic and explicit alike ----

y.multi <- factor(sample(letters[1:3], n, replace = TRUE))
y.ordinal <- factor(sample(1:3, n, replace = TRUE), ordered = TRUE)
y.count <- rpois(n, 4)
y.hurdle <- c(rep(0, n %/% 2L), abs(rnorm(n - n %/% 2L)) + 0.1)

expect_error(
  bart(
    x,
    y.multi,
    family = "multinomial",
    keepFits = FALSE,
    samplerOnly = TRUE
  ),
  "keepFits"
)
expect_error(
  bart(x, y.ordinal, family = "ordinal", keepFits = FALSE, samplerOnly = TRUE),
  "keepFits"
)
expect_error(
  bart(x, y.count, family = "nbinom", keepFits = FALSE, samplerOnly = TRUE),
  "keepFits"
)
expect_error(
  bart(x, y.hurdle, family = "hurdle.lognormal", keepFits = FALSE),
  "keepFits"
)
# the automatic path, not just the explicit one, is named - a callback that
# nobody asked to disable keepFits for still gets refused rather than
# silently breaking the multinomial packager
expect_error(
  bart(
    x,
    y.multi,
    family = "multinomial",
    callback = validCallback,
    samplerOnly = TRUE
  ),
  "callback"
)

# ---- the returned object under keepFits = FALSE, no callback: a REAL run,
# small enough to be fast and safe (nothing here ever calls a callback) ----

fitDropped <- bart(
  x,
  y,
  keepFits = FALSE,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.burn = 2L,
  n.samples = 4L,
  verbose = FALSE
)
expect_null(fitDropped$yhat.train)
expect_null(fitDropped$yhat.train.mean)
expect_false("yhat.train" %in% names(fitDropped))
expect_false("yhat.train.mean" %in% names(fitDropped))

expect_error(plot(fitDropped), "keepFits")
expect_error(extract(fitDropped, sample = "train"), "keepFits")
expect_error(fitted(fitDropped), "keepFits")
expect_error(residuals(fitDropped), "keepFits")

# a keepFits = FALSE fit with test data also loses yhat.test/yhat.test.mean,
# unlike keepTrainingFits = FALSE, which never touched the test channel
xTest <- matrix(runif(6L * 2L), 6L, 2L)
fitDroppedWithTest <- bart(
  x,
  y,
  test = xTest,
  keepFits = FALSE,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.burn = 2L,
  n.samples = 4L,
  verbose = FALSE
)
expect_null(fitDroppedWithTest$yhat.test)
expect_false("yhat.test" %in% names(fitDroppedWithTest))
expect_error(
  extract(fitDroppedWithTest, type = "ev", sample = "test"),
  "keepFits"
)

# keepTrainingFits = FALSE alone (no keepFits involved) still gives its
# longstanding message, unaffected by any of the above
fitTrainingFitsOnly <- bart(
  x,
  y,
  keepTrainingFits = FALSE,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.burn = 2L,
  n.samples = 4L,
  verbose = FALSE
)
expect_null(fitTrainingFitsOnly$yhat.train)
expect_error(plot(fitTrainingFitsOnly), "keepTrainingFits")

# heteroscedastic + keepFits = FALSE: predict(type = "ppd") cannot tell this
# fit apart from a homoscedastic one once s.train is dropped and keepTrees
# is FALSE (its live replay carries no variance surface either) - refused
# by name rather than silently sampling without s(x). hasVariance is the
# fix: it survives keepFits = FALSE where s.train does not.
fitHetero <- bart(
  x,
  y,
  variance = TRUE,
  keepFits = FALSE,
  keepSampler = TRUE,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.burn = 2L,
  n.samples = 4L,
  verbose = FALSE
)
expect_true(fitHetero$hasVariance)
expect_null(fitHetero$s.train)
expect_error(
  predict(fitHetero, x, type = "ppd"),
  "keepFits"
)
