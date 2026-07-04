library(dbarts, quietly = TRUE)

set.seed(99)
n <- 250L
x1 <- runif(n)
x2 <- runif(n, -1, 1)
x3 <- runif(n)
g  <- factor(sample(letters[1:3], n, replace = TRUE))
mu <- ifelse(x1 > 0.5, x2, 0)
y  <- mu + rnorm(n, 0, 0.2)
df <- data.frame(x1, x2, x3, y)

# designation validation happens when the prior resolves against the data
expect_error(dbarts(y ~ x1 + x2 + x3, df, node.prior = linear("zz")),
             pattern = "unrecognized column")
expect_error(dbarts(y ~ x1 + x2 + x3, df, node.prior = linear(10)),
             pattern = "out of range")
expect_error(dbarts(y ~ x1 + x2 + x3, df, node.prior = linear(c(2, 2))),
             pattern = "duplicate")
expect_error(dbarts(y ~ x1 + x2 + g, data.frame(df, g),
                    node.prior = linear("g")),
             pattern = "must be continuous")
# an unresolved designation cannot enter a model object directly
expect_error(new("dbartsModel", node.prior = dbarts:::linear("x2")),
             pattern = "resolved against data")

# fitting recovers the varying slope structure
control <- dbartsControl(n.trees = 20L, n.chains = 1L, n.samples = 40L,
                         n.burn = 150L, keepTrees = TRUE,
                         updateState = FALSE)
set.seed(0)
sampler <- dbarts(y ~ x1 + x2 + x3, df, test = df[1:5, c("x1", "x2", "x3")],
                  node.prior = linear("x2"), control = control)
samples <- sampler$run()
fits <- rowMeans(samples$train)
expect_true(sum((fits - mu)^2) < 0.2 * sum((mean(y) - mu)^2))

# recorded test fits match a saved-tree replay of the same rows
predictions <- sampler$predict(as.matrix(df[1:5, c("x1", "x2", "x3")]))
expect_equal(predictions, samples$test, tolerance = 1e-10)

# getTrees reports one slope column per covariate, NA on internal nodes
trees <- sampler$getTrees(treeNums = 1:2, sampleNums = 1L)
expect_true("beta.x2" %in% names(trees))
expect_true(all(is.na(trees$beta.x2[trees$var > 0])))
expect_true(all(!is.na(trees$beta.x2[trees$var == -1])))

# plotTree labels linear leaves with their coefficients
pdf(NULL)
expect_silent(sampler$plotTree(1L, sampleNum = 1L))
dev.off()

# the leaf covariate designation is fixed at creation: a replacement model
# with a constant node prior is refused
model.const <- sampler$model
model.const@node.prior <- dbarts:::normal(2)
expect_error(sampler$setModel(model.const),
             pattern = "fixed when a sampler is created")

# state serialization carries the slope arrays: a restored sampler
# continues bitwise identically
control.state <- dbartsControl(n.chains = 2L, n.threads = 1L, n.trees = 10L,
                               n.samples = 5L, updateState = FALSE)
sampler.state <- dbarts(y ~ x1 + x2 + x3, df, node.prior = linear("x2"),
                        control = control.state)
invisible(sampler.state$run(30L, 2L))
sampler.state$storeState()
expect_true("tree.params" %in% names(sampler.state$state[[1L]]))
sampler.restored <- dbarts(y ~ x1 + x2 + x3, df, node.prior = linear("x2"),
                           control = control.state)
sampler.restored$setState(sampler.state$state)
expect_identical(sampler.state$run(0L, 3L), sampler.restored$run(0L, 3L))

# the mutable-data surface stays live under linear leaves
set.seed(1)
sampler.mut <- dbarts(y ~ x1 + x2 + x3, df, node.prior = linear("x2"),
                      control = control)
invisible(sampler.mut$run(50L, 5L))
x2.new <- df$x2 * 1.1
expect_silent(sampler.mut$setPredictor(x2.new, "x2", forceUpdate = TRUE))
more <- sampler.mut$run(0L, 5L)
expect_true(all(is.finite(more$train)))

# a probit response composes with linear leaves
z <- rbinom(n, 1L, pnorm(mu / 0.5))
df.binary <- data.frame(x1, x2, x3, z)
set.seed(2)
sampler.binary <- dbarts(z ~ x1 + x2 + x3, df.binary,
                         node.prior = linear("x2"), control = control)
samples.binary <- sampler.binary$run(100L, 20L)
expect_true(all(is.finite(samples.binary$train)))

# linear leaves ride the data-handle views: a full-rows view matches the
# raw-data path bitwise, standardizing with the parent's constants
control.view <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 15L,
                              updateState = FALSE)
sampler.view <- dbarts(y ~ x1 + x2 + x3, df, node.prior = linear("x2"),
                       control = control.view)
handle <- dbarts:::bartcoreDataHandle(sampler.view$control, sampler.view$data)
view <- dbarts:::bartcoreSamplerFromHandle(handle, sampler.view$control,
                                           sampler.view$model,
                                           sampler.view$data,
                                           trainRows = seq_len(n))
full <- dbarts:::bartcoreSampler(sampler.view)
set.seed(7)
samples.view <- dbarts:::bartcoreRun(view, 40L, 20L)
set.seed(7)
samples.full <- dbarts:::bartcoreRun(full, 40L, 20L)
expect_identical(samples.view$sigma, samples.full$sigma)
expect_identical(samples.view$train, samples.full$train)

# a proper fold serves its held-out rows through the gathered covariates
testRows <- seq(1L, n, by = 4L)
set.seed(11)
fold <- dbarts:::bartcoreSamplerFromHandle(handle, sampler.view$control,
                                           sampler.view$model,
                                           sampler.view$data,
                                           setdiff(seq_len(n), testRows),
                                           testRows)
samples.fold <- dbarts:::bartcoreRun(fold, 150L, 100L)
expect_true(all(is.finite(samples.fold$test)))
expect_true(cor(rowMeans(samples.fold$test), mu[testRows]) > 0.7)

# views still refuse raw-predictor mutation under linear leaves
expect_error(dbarts:::bartcoreSetPredictor(fold, sampler.view$data@x),
             pattern = "views hold none")

# xbart accepts a linear node prior, with its k standing in for a missing
# k argument and the k grid overriding per cell
xbart.linear <- xbart(y ~ x1 + x2 + x3, df, node.prior = linear("x2", k = 3),
                      n.samples = 60L, n.burn = c(60L, 30L, 0L), n.reps = 2L,
                      n.trees = 15L, n.threads = 1L, seed = 1L)
expect_true(all(is.finite(xbart.linear)))
xbart.grid <- xbart(y ~ x1 + x2 + x3, df, node.prior = linear("x2"),
                    k = c(1, 4), n.samples = 60L, n.burn = c(60L, 30L, 0L),
                    n.reps = 2L, n.trees = 15L, n.threads = 1L, seed = 1L,
                    drop = FALSE)
expect_equal(dim(xbart.grid)[3L], 2L)
expect_error(xbart(y ~ x1 + x2 + g, data.frame(df, g),
                   node.prior = linear("g"), n.threads = 1L),
             pattern = "must be continuous")

rm(sampler, sampler.state, sampler.restored, sampler.mut, sampler.binary,
   samples, samples.binary, more, trees, predictions, fits, control,
   control.state, df, df.binary, x1, x2, x3, g, y, z, mu, x2.new, n,
   control.view, sampler.view, handle, view, full, samples.view,
   samples.full, testRows, fold, samples.fold, xbart.linear, xbart.grid)
