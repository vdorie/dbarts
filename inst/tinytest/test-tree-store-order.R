# The saved-tree store is a circular buffer and Sampler::run bases each chain's
# writes at the sampler's own cursor, which CARRIES across run() calls. Every
# reader indexes from that cursor: output draw i is slot
# (cursor + capacity - filled + i) mod capacity over the
# filled = min(recorded draws, capacity) most recent draws, oldest first. So a
# hand-driven sampler's replayed draws line up with the run channels that
# recorded them, a store short of capacity reports the draws it holds rather
# than padding with slots nothing wrote, and one holding no draws refuses.
# Readers that walked raw slot order returned a second recorded run's draws
# rotated by the first run's count.

set.seed(22L)
n <- 60L
p <- 6L
x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
y <- 2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.5)

keepSampler <- function(n.samples, n.chains = 1L, n.thin = 1L, n.trees = 5L) {
  control <- dbarts::dbartsControl(
    keepTrees = TRUE,
    n.samples = n.samples,
    n.chains = n.chains,
    n.thin = n.thin,
    n.threads = 1L,
    n.trees = n.trees,
    updateState = FALSE,
    verbose = FALSE
  )
  set.seed(7L)
  dbarts::dbarts(x, y, control = control)
}

# no offset and a gaussian response, so the replayed surface IS the recorded
# train channel, draw for draw
agrees <- function(sampler, expected) {
  max(abs(sampler$predict(x) - expected)) < 1e-12
}

# rotation: a recorded burn-in run of nb draws left the store's cursor at nb,
# and every draw of the sampling run came back rotated by it
for (nb in 1:4) {
  sampler <- keepSampler(6L)
  invisible(sampler$run(0L, nb))
  samples <- sampler$run(0L, 6L)
  expect_true(agrees(sampler, samples$train))
}
rm(nb)

# the single-run case that was always correct stays correct
sampler.one <- keepSampler(6L)
samples.one <- sampler.one$run(2L, 6L)
expect_true(agrees(sampler.one, samples.one$train))

# a run longer than the store keeps its last capacity draws, and the next run
# overwrites them all
sampler.wrap <- keepSampler(4L)
samples.wrap1 <- sampler.wrap$run(0L, 7L)
expect_true(agrees(sampler.wrap, samples.wrap1$train[, 4:7]))
samples.wrap2 <- sampler.wrap$run(0L, 4L)
expect_true(agrees(sampler.wrap, samples.wrap2$train))

# runs that do not divide the capacity: the store holds the tail of one run
# and the head of the next, chronologically
sampler.mix <- keepSampler(6L)
samples.mix1 <- sampler.mix$run(2L, 4L)
expect_true(agrees(sampler.mix, samples.mix1$train))
samples.mix2 <- sampler.mix$run(2L, 4L)
expect_true(agrees(
  sampler.mix,
  cbind(samples.mix1$train[, 3:4], samples.mix2$train)
))
samples.mix3 <- sampler.mix$run(2L, 4L)
expect_true(agrees(
  sampler.mix,
  cbind(samples.mix2$train[, 3:4], samples.mix3$train)
))

# one cursor fans out to every chain, so the rotation was per-sampler
sampler.chains <- keepSampler(6L, n.chains = 2L)
invisible(sampler.chains$run(0L, 2L))
samples.chains <- sampler.chains$run(0L, 6L)
expect_equal(dim(sampler.chains$predict(x)), c(n, 6L, 2L))
expect_true(agrees(sampler.chains, samples.chains$train))

# thinning changes the sweeps per draw, not the draws recorded
sampler.thin <- keepSampler(4L, n.thin = 3L)
samples.thin1 <- sampler.thin$run(0L, 4L)
expect_true(agrees(sampler.thin, samples.thin1$train))
samples.thin2 <- sampler.thin$run(0L, 2L)
expect_true(agrees(
  sampler.thin,
  cbind(samples.thin1$train[, 3:4], samples.thin2$train)
))

# partial fill: k real draws, in order, and no phantom slots behind them
sampler.partial <- keepSampler(6L)
samples.partial <- sampler.partial$run(0L, 3L)
expect_equal(dim(sampler.partial$predict(x)), c(n, 3L))
expect_true(agrees(sampler.partial, samples.partial$train))

# the count rides the state, so a restore replays the same draws in the same
# order - and so does a copy
sampler.state <- keepSampler(6L)
invisible(sampler.state$run(0L, 2L))
samples.state <- sampler.state$run(0L, 6L)
predicted.state <- sampler.state$predict(x)
sampler.state$storeState()
expect_equal(attr(sampler.state$state, "recordedDraws"), 6L)

sampler.restored <- keepSampler(6L)
sampler.restored$setState(sampler.state$state)
expect_equal(sampler.restored$predict(x), predicted.state)
expect_true(agrees(sampler.restored, samples.state$train))
expect_equal(sampler.state$copy()$predict(x), predicted.state)

# a partly filled store restores as partly filled, rather than being promoted
sampler.partialState <- keepSampler(6L)
samples.partialState <- sampler.partialState$run(0L, 3L)
sampler.partialState$storeState()
expect_equal(attr(sampler.partialState$state, "recordedDraws"), 3L)
sampler.partialRestored <- keepSampler(6L)
sampler.partialRestored$setState(sampler.partialState$state)
expect_equal(dim(sampler.partialRestored$predict(x)), c(n, 3L))
expect_true(agrees(sampler.partialRestored, samples.partialState$train))

# a store nothing has been recorded into refuses at every reader: its slots
# hold zero-leaf trees, which replay as a legitimate constant surface
sampler.empty <- keepSampler(4L)
expect_error(sampler.empty$predict(x), pattern = "no recorded draws")
expect_error(sampler.empty$getTrees(), pattern = "no recorded draws")
expect_error(
  sampler.empty$printTrees(treeNums = 1L),
  pattern = "no recorded draws"
)
# ... and so does one whose store was emptied by a resize
sampler.resized <- keepSampler(4L)
invisible(sampler.resized$run(0L, 4L))
control.resized <- sampler.resized$control
control.resized@n.samples <- 6L
sampler.resized$setControl(control.resized)
expect_error(sampler.resized$predict(x), pattern = "no recorded draws")

# the count cannot be inferred from a state that lacks it, so the encoding
# floor refuses one written before it existed
expect_equal(attr(sampler.state$state, "formatVersion"), 3L)
state.old <- sampler.state$state
attr(state.old, "formatVersion") <- 2L
sampler.floor <- keepSampler(6L)
expect_error(
  sampler.floor$setState(state.old),
  pattern = "encoding version 2.*oldest this dbarts \\(3\\)"
)
state.stripped <- sampler.state$state
attr(state.stripped, "recordedDraws") <- NULL
expect_error(
  sampler.floor$setState(state.stripped),
  pattern = "malformed recorded draw count"
)

# the per-forest replay reads the same store on the same terms
sampler.forests <- local({
  control <- dbarts::dbartsControl(
    keepTrees = TRUE,
    n.samples = 6L,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    updateState = FALSE,
    verbose = FALSE
  )
  set.seed(7L)
  z <- rbinom(n, 1L, 0.5)
  dbarts::dbarts(
    x,
    y,
    forests = list(dbarts::forest(), dbarts::forest(basis = ~ factor(z))),
    control = control
  )
})
invisible(sampler.forests$run(0L, 2L))
samples.forests <- sampler.forests$run(0L, 6L)
expect_true(
  max(abs(sampler.forests$predictForests(x) - samples.forests$forestFits)) <
    1e-12
)

# getTrees and printTrees address the same oldest-first draws predict does.
# With one tree every predicted column is one affine image of that draw's leaf
# values, so the leaves returned for draw d identify the draw predict reports
# in column d - the cross-check that catches a reader left on slot order.
sampler.trees <- keepSampler(6L, n.trees = 1L)
invisible(sampler.trees$run(20L, 2L))
invisible(sampler.trees$run(0L, 6L))
predicted.trees <- sampler.trees$predict(x)
leavesOf <- function(sampler, draw) {
  nodes <- sampler$getTrees(treeNums = 1L, sampleNums = draw)
  sort(nodes$value[nodes$var == -1L])
}
leaves.first <- leavesOf(sampler.trees, 1L)
values.first <- sort(unique(round(predicted.trees[, 1L], 12)))
expect_true(length(leaves.first) > 1L)
expect_equal(length(values.first), length(leaves.first))
scale.trees <- diff(range(values.first)) / diff(range(leaves.first))
shift.trees <- values.first[1L] - scale.trees * leaves.first[1L]
drawsMatch <- logical(6L)
for (draw in seq_len(6L)) {
  leaves.draw <- leavesOf(sampler.trees, draw)
  values.draw <- sort(unique(round(predicted.trees[, draw], 12)))
  drawsMatch[draw] <- length(values.draw) == length(leaves.draw) &&
    max(abs(values.draw - (shift.trees + scale.trees * leaves.draw))) < 1e-10
}
expect_true(all(drawsMatch))
rm(draw)

# printTrees names the draw getTrees returns, not the slot holding it
printed <- capture.output(
  sampler.trees$printTrees(treeNums = 1L, sampleNums = 1L)
)
printed.leaves <- sort(as.numeric(
  sub(".*pred: ", "", printed[grepl("pred: ", printed, fixed = TRUE)])
))
expect_true(max(abs(printed.leaves - leaves.first)) < 1e-5)

# a warm start selects donor draws on the same axis: 'samples' 1 is the oldest
# recorded draw, and the trees it installs are that draw's own
donor.trees <- keepSampler(6L, n.trees = 1L)
invisible(donor.trees$run(20L, 2L))
invisible(donor.trees$run(0L, 4L))
donor.oldest <- donor.trees$getTrees(treeNums = 1L, sampleNums = 1L)
dest.trees <- keepSampler(6L, n.trees = 1L)
dest.trees$installTrees(donor.trees, samples = 1L)
installed <- dest.trees$getTrees(treeNums = 1L, current = TRUE)
expect_equal(installed$var, donor.oldest$var)
expect_equal(installed$value, donor.oldest$value)
