source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

df <- with(testData, data.frame(x, y))
df$X10 <- as.factor(paste0("C", 1 + round(4 * df$X10, 0)))

fitCall <- quote(dbarts::bart(
  testData$x,
  testData$y,
  ndpost = 1,
  nskip = 0,
  ntree = 1,
  keeptrees = TRUE,
  verbose = FALSE
))

# test that works defaults, no column names
bartFit <- eval(fitCall)
expect_equal(
  length(bartFit$fit$model@tree.prior@splitProbabilities),
  0L
)


fitCall$splitprobs <- quote(1 / numvars)

bartFit <- eval(fitCall)
expect_equal(
  length(bartFit$fit$model@tree.prior@splitProbabilities),
  0L
)

fitCall$splitprobs <- quote(1)

bartFit <- eval(fitCall)
expect_equal(
  length(bartFit$fit$model@tree.prior@splitProbabilities),
  0L
)

rm(bartFit)


# test that works specific values, no column names
probs <- c(2, rep.int(1, ncol(testData$x) - 1L))

fitCall$splitprobs <- quote(probs)
bartFit <- eval(fitCall)

expect_equal(bartFit$fit$model@tree.prior@splitProbabilities, probs / sum(probs))

probs <- probs[-1L]
expect_error(
  eval(fitCall),
  "length of input \\(9\\) does not equal number of columns in model matrix \\(10\\)"
)

probs <- c(-1, probs)
expect_error(
  eval(fitCall),
  "'splitProbabilities' must form a simplex"
)

probs[1L] <- NA_real_
expect_error(
  eval(fitCall),
  "missing values for columns 1"
)

rm(bartFit, probs)

rm(fitCall)


fitCall <- quote(dbarts::bart(
  y ~ .,
  df,
  ndpost = 3,
  nskip = 0,
  ntree = 2,
  verbose = FALSE,
  keeptrees = TRUE
))

# test that works with column names
fitCall$splitprobs <- quote(c(X4 = 2, X10 = 1.5, .default = 1))
bartFit <- eval(fitCall)

split.probs <- bartFit$fit$model@tree.prior@splitProbabilities

expect_equal(length(split.probs), ncol(df) - 2L + nlevels(df$X10))
expect_equal(split.probs[["X4"]], 2 * split.probs[["X1"]])
x10_values <- startsWith(names(split.probs), "X10.")
expect_equal(sum(x10_values), nlevels(df$X10))
expect_true(
  all(split.probs[["X4"]] ==  2 * split.probs[x10_values] / 1.5)
)
default_values <- !x10_values
default_values[names(split.probs) == "X4"] <- FALSE
expect_true(
  all(split.probs[default_values] == split.probs[default_values][1L])
)


fitCall$splitprobs <- quote(c(X4 = -1, X10 = 1.5, .default = 1))
expect_error(
  eval(fitCall),
  "'splitProbabilities' must form a simplex"
)

fitCall$splitprobs <- quote(c(X4 = NA_real_, X10 = 1.5, .default = 1))
expect_error(
  eval(fitCall),
  "missing values for columns 'X4'"
)

rm(default_values, x10_values, split.probs, bartFit)

rm(fitCall)


# test that split probabilities sample from prior
set.seed(0L)
n.trees <- 200L
control <- dbarts::dbartsControl(
  n.burn = 0L, n.samples = 1L, n.thin = 1L,
  n.trees = n.trees, keepTrees = FALSE,
  n.chains = 1L, n.threads = 1L,
  updateState = FALSE, verbose = FALSE
)
sampler <- dbarts::dbarts(
  y ~ ., df, control = control,
  tree.prior = cgm(split.probs = c(X4 = 2, .default = 1))
)
sampler$sampleTreesFromPrior()

trees <- sampler$getTrees()
treeTable <- table(trees$var)
treeTable <- treeTable[setdiff(names(treeTable), "-1")]
names(treeTable) <- colnames(sampler$data@x)
expect_true(
  abs(treeTable[["X4"]] - 2 * mean(treeTable[setdiff(names(treeTable), "X4")]))
    / n.trees < 0.05
)

rm(treeTable, trees, sampler, control, n.trees)


# test that split probabilities sample from posterior

# X6 is uncorrelated
set.seed(0L)
n.trees <- 5L
control <- dbarts::dbartsControl(
  n.burn = 1000L, n.samples = 100L, n.thin = 1L,
  n.trees = n.trees, keepTrees = FALSE,
  n.chains = 1L, n.threads = 1L,
  updateState = FALSE, verbose = FALSE
)
sampler <- dbarts::dbarts(
  y ~ ., df, control = control,
  tree.prior = cgm(split.probs = c(X6 = 2, .default = 1))
)
samples <- sampler$run(200L, 100L)


varcounts <- apply(samples$varcount, 1L, mean)
names(varcounts) <- colnames(sampler$data@x)

expect_true(all(varcounts[["X6"]] <= varcounts[paste0("X", 1:5)]))
expect_true(all(varcounts[["X6"]] >= varcounts[setdiff(names(varcounts), paste0("X", 1:6))] - 0.02))

rm(varcounts, samples, sampler, control, n.trees)

rm(testData)

