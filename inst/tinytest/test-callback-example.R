# The running-mean recipe the component vignette carries
# (vignette("dbarts-as-a-component"), recipe 7), in the plain-C form
# inst/tinytest/capi/consumer.c holds beside the counting consumer
# test-capi.R drives - the same source file, compiled the same way, only the
# entry points used differ.
#
# This is the test that the RECIPE is correct, not merely that it compiles:
# a callback-accumulated running mean of the training channel must equal the
# same seeded fit's own yhat.train.mean, with keepFits = TRUE forced so that
# channel comes back for the comparison (a callback fit defaults it FALSE).

source(
  system.file("common", "capiConsumer.R", package = "dbarts"),
  local = TRUE
)
consumer <- compileCapiConsumer("capi-mean", "the callback-example consumer")
if (!is.null(consumer$skip)) exit_file(consumer$skip)
CALL <- consumer$CALL

meanFn <- CALL("capi_mean_function")

n <- 40L
nChains <- 2L
set.seed(303, sample.kind = "Rejection")
x <- matrix(runif(n * 2L), n, 2L)
y <- x[, 1L] + rnorm(n, 0, 0.2)

# the ONLY per-observation allocation, per the vignette: REAL(acc) is what
# capi_mean_context_new takes, so acc must not be reassigned before the run
acc <- numeric(n * nChains)
ctx <- CALL("capi_mean_context_new", acc, n, nChains)

fit <- bart(
  x,
  y,
  n.chains = nChains,
  n.threads = nChains,
  n.trees = 20L,
  n.burn = 20L,
  n.samples = 40L,
  keepFits = TRUE, # explicit: overrides the callback's automatic FALSE
  callback = list(fn = meanFn, context = ctx),
  verbose = FALSE
)

expect_equal(CALL("capi_mean_status"), 0L)
# The tolerance covers summation order and nothing else: an incremental mean
# per chain then rowMeans, against yhat.train.mean's own reduction over the
# whole n x S x C array. Measured, that gap is 1.1e-16 absolute on values of
# order 1, so 1e-12 keeps four orders of headroom and still fails a mistaken
# slice offset, a doubled chain or an off-by-one draw count.
means <- rowMeans(matrix(acc, n, nChains))
expect_equal(means, fit$yhat.train.mean, tolerance = 1e-12)
