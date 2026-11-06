# A run that stops early - an interrupt, a draw callback's nonzero return, a
# sweep closure's TRUE, an error - drops the saved trees it wrote and keeps
# the earlier runs' draws it did not overwrite: predict and getTrees serve
# them, oldest first, and a sampler stored and reloaded after the stop serves
# the same.

set.seed(20261106L)
n <- 50L
numDraws <- 20L
x <- matrix(runif(n * 2L), n, 2L, dimnames = list(NULL, c("x1", "x2")))
y <- x[, 1L] + rnorm(n, 0, 0.2)
x.test <- x[1L:5L, ]
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = numDraws,
  keepTrees = TRUE,
  updateState = FALSE
)
# the sweep-closure entry, the closure stopping the run before sweep stopAt
runStoppedAt <- function(sampler, numBurnIn, numSamples, stopAt) {
  .Call(
    dbarts:::C_dbarts_bartcore_runWithCallback,
    sampler$getPointer(),
    numBurnIn,
    numSamples,
    list(),
    function(i) i == stopAt,
    environment()
  )
}
treesOfDraws <- function(trees, draws) {
  kept <- trees[trees$sample %in% draws, ]
  kept$sample <- match(kept$sample, draws)
  rownames(kept) <- NULL
  kept
}

# a full store, and a run stopped before its kept sweep 9: the nine slots it
# wrote are dropped, run 1's draws 10 to 20 are served
sampler <- dbarts(x, y, n.trees = 4L, control = control)
invisible(sampler$run(10L, numDraws))
before <- sampler$predict(x.test)
trees.before <- sampler$getTrees()
expect_null(runStoppedAt(sampler, 0L, numDraws, 9L))
expect_identical(sampler$predict(x.test), before[, 10L:numDraws])
trees.after <- sampler$getTrees()
rownames(trees.after) <- NULL
expect_identical(trees.after, treesOfDraws(trees.before, 10L:numDraws))

# a stop in burn-in writes nothing and loses nothing
expect_null(runStoppedAt(sampler, 10L, numDraws, 3L))
expect_identical(sampler$predict(x.test), before[, 10L:numDraws])

# after the stop, a stored state and a reload serve the live sampler's draws
sampler$storeState()
stateFile <- tempfile(fileext = ".rds")
saveRDS(sampler, stateFile)
reloaded <- readRDS(stateFile)
unlink(stateFile)
expect_identical(reloaded$predict(x.test), sampler$predict(x.test))

# a full run afterwards writes from where the stopped run found the cursor
invisible(sampler$run(0L, numDraws))
expect_identical(dim(sampler$predict(x.test)), c(5L, numDraws))

# a store not full keeps every draw a short stopped run leaves unwritten
half <- dbarts(x, y, n.trees = 4L, control = control)
invisible(half$run(10L, numDraws / 2L))
half.before <- half$predict(x.test)
expect_null(runStoppedAt(half, 0L, numDraws, 3L))
expect_identical(half$predict(x.test), half.before)
