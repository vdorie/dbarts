# smoke coverage for the plotting generics: plot.bart on a gaussian fit, on a
# binary bart2 fit, and plot.rbart on a grouped fit. Rendered to a null device;
# each plot is silent (no output, message, or warning) at these small sizes.

set.seed(5)
n <- 120L
x <- matrix(runif(n * 3L), n, 3L)
y.cont <- 2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.3)
z.bin  <- rbinom(n, 1L, pnorm(1.5 * x[, 1L] - 0.5))
g      <- sample(4L, n, replace = TRUE)

# gaussian bart
fit.bart <- bart(x, y.cont, ntree = 10L, nskip = 20L, ndpost = 20L,
                 nchain = 1L, nthread = 1L, verbose = FALSE, keeptrees = TRUE)
pdf(NULL)
expect_silent(plot(fit.bart))
dev.off()

# binary bart2
fit.bart2 <- bart2(x, z.bin, n.trees = 10L, n.burn = 20L, n.samples = 20L,
                   n.chains = 1L, n.threads = 1L, verbose = FALSE)
pdf(NULL)
expect_silent(plot(fit.bart2))
dev.off()

# grouped rbart_vi
fit.rbart <- rbart_vi(y.cont ~ x, group.by = g, n.samples = 8L, n.burn = 5L,
                      n.thin = 1L, n.chains = 1L, n.trees = 10L, n.threads = 1L,
                      verbose = FALSE)
pdf(NULL)
expect_silent(plot(fit.rbart))
dev.off()

# plotTree now dispatches at the fit level (previously reachable only through
# $fit$plotTree). Kept bart, bart2, and rbart fits plot a single tree, the
# sampler dispatches directly, and a fit without kept trees errors.
pt.bart2 <- bart2(x, z.bin, n.trees = 10L, n.burn = 20L, n.samples = 20L,
                  n.chains = 2L, n.threads = 1L, verbose = FALSE,
                  keepTrees = TRUE)
pt.rbart <- rbart_vi(y.cont ~ x, group.by = g, n.samples = 8L, n.burn = 5L,
                     n.thin = 1L, n.chains = 1L, n.trees = 10L, n.threads = 1L,
                     verbose = FALSE, keepTrees = TRUE)
pdf(NULL)
expect_silent(plotTree(fit.bart, treeNum = 1L))
expect_silent(plotTree(fit.bart$fit, treeNum = 2L))
expect_silent(plotTree(pt.bart2, treeNum = 1L, chainNum = 2L, sampleNum = 5L))
expect_silent(plotTree(pt.rbart, treeNum = 1L))
dev.off()

expect_error(plotTree(fit.bart2, treeNum = 1L),
             pattern = "requires the trees to be kept")
expect_error(plotTree(pt.rbart, treeNum = 1L, chainNum = 5L),
             pattern = "chainNum must be a single chain index")

rm(fit.bart, fit.bart2, fit.rbart, pt.bart2, pt.rbart, x, y.cont, z.bin, g, n)
