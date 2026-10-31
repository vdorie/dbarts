# What each forest's record and the model's other records hold reaches the run:
# a forest's own tree prior and blocks, the variance forest's power and columns,
# a forest's label in a refusal, each model record across $setModel, and each
# cell's tree count in a cross-validation.

forest <- dbartsForests$forest
varianceForest <- dbartsForests$varianceForest
blocks <- dbartsForests$blocks

set.seed(41)
n <- 80L
frame <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
frame$z <- rbinom(n, 1L, 0.5)
frame$zn <- rnorm(n)
frame$y <- with(frame, 2 * x1 + x2 * zn + rnorm(n, sd = 0.2))
x <- as.matrix(frame[c("x1", "x2", "x3")])
y <- frame$y
run <- function(...) {
  dbartsControl(
    n.samples = 20L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 9L,
    keepTrees = TRUE,
    ...
  )
}
secondForestTrees <- function(...) {
  sampler <- dbarts(
    y ~ x1 + x2 + x3,
    frame,
    forests = list(forest(), forest(basis = ~zn, n.trees = 20L, ...)),
    control = run()
  )
  invisible(sampler$run(40L, 20L))
  sampler$getTrees(forest = 2L)
}
# the deepest a tree's depth-first var column reaches, -1 marking a leaf
depth <- function(trees) {
  deepest <- function(tree) {
    at <- 0L
    most <- 0L
    stack <- integer()
    for (var in tree$var) {
      most <- max(most, at)
      if (var > 0L) {
        stack <- c(stack, at + 1L)
        at <- at + 1L
      } else if (length(stack) > 0L) {
        at <- stack[length(stack)]
        stack <- stack[-length(stack)]
      }
    }
    most
  }
  max(vapply(split(trees, paste(trees$tree, trees$sample)), deepest, 0L))
}

# a forest's own base: near 0 keeps every tree a stump
trees <- secondForestTrees(base = 1e-8, power = 1)
expect_true(mean(trees$var != -1L) < 0.05)
# its own power: a huge one stops at depth 1
expect_true(depth(secondForestTrees(base = 0.95, power = 40)) <= 1L)
# its own blocks: each tree splits on one group
trees <- secondForestTrees(blocks = blocks(groups = list("x1", c("x2", "x3"))))
mixed <- vapply(
  split(trees$var, paste(trees$tree, trees$sample)),
  function(var) (1L %in% var) && any(var %in% 2:3),
  NA
)
expect_false(any(mixed))

# the variance forest's own power: one tree at a huge one has at most two leaves
frame$level <- round(frame$x3 * 3) / 3
leavesPerDraw <- dbarts(
  y ~ x1 + x2 + level,
  frame,
  variance = varianceForest(n.trees = 1L, base = 0.95, power = 40),
  control = run()
)
invisible(leavesPerDraw$run(40L, 20L))
for (draw in seq_len(20L)) {
  invisible(leavesPerDraw$run(1L, 1L))
  expect_true(length(unique(leavesPerDraw$getVariance())) <= 2L, info = draw)
}
# its columns: the surface is a function of the one column it may split on
onLevel <- dbarts(
  y ~ x1 + x2 + level,
  frame,
  variance = varianceForest(n.trees = 5L, vars = ~level),
  control = run()
)
invisible(onLevel$run(40L, 20L))
levelsSeen <- tapply(onLevel$getVariance(), frame$level, function(values) {
  length(unique(values))
})
expect_true(all(levelsSeen == 1L))

# a refusal naming a forest gives its label
frame$const <- 1
refused <- tryCatch(
  dbarts(
    y ~ x1 + x2,
    frame,
    forests = list(forest(), treatment = forest(basis = ~const)),
    control = run()
  ),
  error = conditionMessage
)
expect_true(grepl("treatment", refused))

# a model made by hand keeps the sampler's categories, shape and periods
keeps <- function(sampler, name) {
  model <- sampler$model
  held <- slot(model, name)
  sampler$setModel(new(
    "dbartsModel",
    model@tree.prior,
    model@leaf.prior,
    model@leaf.hyperprior,
    model@resid.prior,
    leaf.scale = model@leaf.scale
  ))
  expect_identical(slot(sampler$model, name), held, info = name)
}
ordered <- cut(
  y,
  quantile(y, c(0, 0.3, 0.7, 1)),
  include.lowest = TRUE,
  ordered_result = TRUE
)
keeps(dbarts(x, ordered, family = "ordinal", control = run()), "n.categories")
keeps(dbarts(x, rpois(n, 2), family = "nbinom", control = run()), "shape")
keeps(
  dbarts(
    x,
    cbind(ceiling(rexp(n) * 3), rbinom(n, 1L, 0.7)),
    family = "hazard",
    control = run()
  ),
  "hazard.periods"
)
# and a variance forest's moves are installed as a mean forest's are
moved <- function(probs) {
  dbarts(
    y ~ x1 + x2 + x3,
    frame,
    variance = varianceForest(proposal.probs = probs),
    forests = list(forest()),
    control = run()
  )
}
sampler <- moved(c(birth_death = 0.3, swap = 0.2, change = 0.5, birth = 0.5))
other <- moved(c(birth_death = 0.6, swap = 0.1, change = 0.3, birth = 0.5))
given <- sampler$model
given@variance@moves <- other$model@variance@moves
sampler$setModel(given)
expect_identical(sampler$model@variance@moves, other$model@variance@moves)

# a sampler's data taken by $setData keeps the sampler's status: none on a plain
# one, given a survival fit's
aftFit <- dbarts(
  x,
  cbind(rexp(n), rbinom(n, 1L, 0.7)),
  family = "aft",
  control = run()
)
plain <- dbarts(x, y, control = run())
plain$setData(aftFit$data)
expect_null(plain$data@status)

# each cell of a cross-validation has the grid's tree count: a loss counting the
# distinct fits of a draw sees at most a stump's two at one tree, whichever order
# the grid is in, and many at forty
set.seed(11)
xSmall <- matrix(runif(180L), 60L, 3L)
ySmall <- xSmall[, 1L] + rnorm(60L, sd = 0.2)
distinctFits <- function(y.test, y.test.hat, weights) {
  max(apply(y.test.hat, 2L, function(column) length(unique(column))))
}
for (grid in list(c(1L, 40L), c(40L, 1L))) {
  losses <- xbart(
    xSmall,
    ySmall,
    n.trees = grid,
    k = 2,
    power = 30,
    base = 0.95,
    n.samples = 3L,
    method = "k-fold",
    n.test = 3L,
    n.reps = 2L,
    n.burn = c(3L, 1L),
    n.threads = 1L,
    seed = 7L,
    loss = distinctFits,
    n.cuts = 10L
  )
  expect_true(all(losses[, "1"] <= 2), info = paste(grid, collapse = " "))
  expect_true(all(losses[, "40"] > 2), info = paste(grid, collapse = " "))
}

## The rule a fit runs is the one its verbose summary reports.
summaryOf <- function(...) {
  capture.output(dbarts(x, y, ..., control = run(verbose = TRUE)))
}
expect_true(any(grepl(
  "use quantiles for rule cut points: true",
  summaryOf(useQuantiles = TRUE),
  fixed = TRUE
)))
expect_true(any(grepl(
  "use quantiles for rule cut points: false",
  summaryOf(),
  fixed = TRUE
)))

## A fit's margins take the forests' labels, a forest held by amplitude = fixed()
## is named, and a fit's printed count is its record's.
fixed <- dbartsPriors$fixed
labeled <- bart(
  y ~ forest(x1 + x2) + forest(x3, basis = zn, amplitude = fixed()),
  frame,
  n.samples = 20L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 9L,
  keepTrees = TRUE,
  verbose = FALSE
)
labels <- vapply(labeled$fit$model@forests, methods::slot, "", "label")
expect_false(identical(labels, c("forest1", "forest2")))
expect_identical(names(labeled$leaf.prior), labels)
expect_identical(names(labeled$fixed$amplitude), labels[2L])
printed <- capture.output(print(bart(
  y ~ x1 + x2 + x3,
  frame,
  n.trees = 30L,
  n.samples = 20L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 9L,
  keepTrees = TRUE,
  verbose = FALSE
)))
expect_true("n.trees: 30" %in% printed)
