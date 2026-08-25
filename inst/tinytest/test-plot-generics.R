# smoke coverage for the plotting generics: plot.bart on a gaussian fit, on a
# binary bart2 fit, and plot.rbart on a grouped fit. Rendered to a null device;
# each plot is silent (no output, message, or warning) at these small sizes.

set.seed(5)
n <- 120L
x <- matrix(runif(n * 3L), n, 3L)
y.cont <- 2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.3)
z.bin <- rbinom(n, 1L, pnorm(1.5 * x[, 1L] - 0.5))
g <- sample(4L, n, replace = TRUE)

# gaussian bart
fit.bart <- bart(
  x,
  y.cont,
  ntree = 10L,
  nskip = 20L,
  ndpost = 20L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE,
  keeptrees = TRUE
)
pdf(NULL)
expect_silent(plot(fit.bart))
dev.off()

# binary bart2
fit.bart2 <- bart2(
  x,
  z.bin,
  n.trees = 10L,
  n.burn = 20L,
  n.samples = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
pdf(NULL)
expect_silent(plot(fit.bart2))
dev.off()

# grouped rbart_vi
fit.rbart <- rbart_vi(
  y.cont ~ x,
  group.by = g,
  n.samples = 8L,
  n.burn = 5L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 10L,
  n.threads = 1L,
  verbose = FALSE
)
pdf(NULL)
expect_silent(plot(fit.rbart))
dev.off()

# plotTree now dispatches at the fit level (previously reachable only through
# $fit$plotTree). Kept bart, bart2, and rbart fits plot a single tree, the
# sampler dispatches directly, and a fit without kept trees errors.
pt.bart2 <- bart2(
  x,
  z.bin,
  n.trees = 10L,
  n.burn = 20L,
  n.samples = 20L,
  n.chains = 2L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
pt.rbart <- rbart_vi(
  y.cont ~ x,
  group.by = g,
  n.samples = 8L,
  n.burn = 5L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 10L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
pdf(NULL)
expect_silent(plotTree(fit.bart, treeNum = 1L))
expect_silent(plotTree(fit.bart$fit, treeNum = 2L))
expect_silent(plotTree(pt.bart2, treeNum = 1L, chainNum = 2L, sampleNum = 5L))
expect_silent(plotTree(pt.rbart, treeNum = 1L))
dev.off()

expect_error(
  plotTree(fit.bart2, treeNum = 1L),
  pattern = "requires the trees to be kept"
)
expect_error(
  plotTree(pt.rbart, treeNum = 1L, chainNum = 5L),
  pattern = "'chainNum' must be a single chain index"
)

# the print methods summarize a fit to the console, and keep doing so with
# keepCall = FALSE (previously just "NULL()")
# the synopsis is the last block of the printout, one labelled line per fact,
# and each label names the value under it: a fit whose control survived
# reports its own n.trees and n.burn, and a fit's family is its own
expect_equal(
  tail(capture.output(print(fit.bart)), 5L),
  c(
    "family: gaussian",
    "n.chains: 1",
    "n.trees: 10",
    "n.burn: 20",
    "kept draws (per chain): 20"
  )
)
expect_equal(
  tail(capture.output(print(fit.rbart)), 5L),
  c(
    "family: gaussian",
    "n.chains: 1",
    "n.trees: 10",
    "n.burn: 5",
    "kept draws (per chain): 8"
  )
)
# a fit with no kept trees keeps no control, so the two control-borne lines
# drop out - and the family line is still the fit's own, not gaussian
expect_equal(
  tail(capture.output(print(fit.bart2)), 3L),
  c("family: probit", "n.chains: 1", "kept draws (per chain): 20")
)
fit.noCall <- bart2(
  x,
  y.cont,
  n.trees = 10L,
  n.burn = 5L,
  n.samples = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepCall = FALSE
)
noCallOutput <- capture.output(print(fit.noCall))
expect_true(!any(grepl("NULL()", noCallOutput, fixed = TRUE)))
expect_true(any(grepl("^family:", noCallOutput)))
rm(fit.noCall, noCallOutput)

# keepTrainingFits = FALSE: plot/fitted/residuals must stop early and name
# the control flag, instead of dying inside apply() on a NULL yhat.train
# (plot) or silently returning NA (fitted/residuals)
fit.noTrainFits <- bart2(
  x,
  y.cont,
  n.trees = 10L,
  n.burn = 5L,
  n.samples = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrainingFits = FALSE
)
expect_error(plot(fit.noTrainFits), pattern = "keepTrainingFits")
expect_error(fitted(fit.noTrainFits), pattern = "keepTrainingFits")
expect_error(residuals(fit.noTrainFits), pattern = "keepTrainingFits")

fit.rbart.noTrainFits <- rbart_vi(
  y.cont ~ x,
  group.by = g,
  n.samples = 8L,
  n.burn = 5L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 10L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrainingFits = FALSE
)
expect_error(plot(fit.rbart.noTrainFits), pattern = "keepTrainingFits")
expect_error(fitted(fit.rbart.noTrainFits), pattern = "keepTrainingFits")
expect_error(residuals(fit.rbart.noTrainFits), pattern = "keepTrainingFits")

rm(fit.noTrainFits, fit.rbart.noTrainFits)

# residuals is always against the training response; a caller-supplied
# 'sample' collided with the fixed sample = "train" residuals.bart/.rbart
# forward to fitted, raising a raw 'formal argument "sample" matched by
# multiple actual arguments' - refused by name instead
residualsSampleReason <- paste0(
  "'sample' is not used by residuals; residuals are always against the ",
  "training response"
)
expect_error(
  residuals(fit.bart, sample = "train"),
  residualsSampleReason,
  fixed = TRUE
)
expect_error(
  residuals(fit.rbart, sample = "train"),
  residualsSampleReason,
  fixed = TRUE
)
rm(residualsSampleReason)

rm(fit.bart, fit.bart2, fit.rbart, pt.bart2, pt.rbart)


# multi-chain fits carry a matrix-shaped 'sigma', so the trace panel draws one
# bridged line per chain rather than a single scatter
fit.chains <- bart(
  x,
  y.cont,
  ntree = 10L,
  nskip = 10L,
  ndpost = 20L,
  nchain = 2L,
  nthread = 1L,
  verbose = FALSE
)
pdf(NULL)
expect_silent(plot(fit.chains))
dev.off()

fit.rbart.chains <- rbart_vi(
  y.cont ~ x,
  group.by = g,
  n.samples = 20L,
  n.burn = 5L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 10L,
  n.threads = 1L,
  verbose = FALSE
)
pdf(NULL)
expect_silent(plot(fit.rbart.chains))
dev.off()

# a binary grouped fit has no residual scale, so plot.rbart takes its
# probability-interval branch instead of the E(Y | x) one
fit.rbart.bin <- rbart_vi(
  z.bin ~ x,
  group.by = g,
  n.samples = 20L,
  n.burn = 5L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 10L,
  n.threads = 1L,
  verbose = FALSE
)
pdf(NULL)
expect_silent(plot(fit.rbart.bin))
dev.off()

rm(fit.chains, fit.rbart.chains, fit.rbart.bin)


# a namespace-qualified call stores the `dbarts::bart` call in call[[1L]],
# which as.character() splits into c("::", "dbarts", "bart"); the guards that
# name the missing argument used to compare against that and die with "the
# condition has length > 1" instead
fit.qualified <- dbarts::bart(
  x,
  y.cont,
  ntree = 10L,
  nskip = 5L,
  ndpost = 5L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE,
  keeptrainfits = FALSE
)
expect_error(plot(fit.qualified), pattern = "keeptrainfits")

fit2.qualified <- dbarts::bart2(
  x,
  y.cont,
  n.trees = 10L,
  n.burn = 5L,
  n.samples = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrainingFits = FALSE
)
expect_error(plot(fit2.qualified), pattern = "keepTrainingFits")

rm(fit.qualified, fit2.qualified)

rm(x, y.cont, z.bin, g, n)
