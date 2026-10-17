# A warm start from a donor of another tree count. Into a forest of more
# trees the donor's go in all, in its order, and the rest start as single
# leaves at 0, or at 1 in a variance forest: the donor as it stands, TRUE.
# Into a forest of fewer the donor is declined, FALSE, and forced takes the
# trees whose fitted values vary most over this sampler's rows, each row
# weighing by its case weight, in the donor's order and as installed; the
# mean and the variance forest choose apart, the variance forest by the log
# of its factors, and each chain among the trees of the donor sample it
# takes.

visibly <- function(value) list(value = value, visible = TRUE)
answerOf <- function(expr) tryCatch(withVisible(expr), error = conditionMessage)
# the state the engine holds now, the sampler's state field left alone
stateOf <- function(sampler) {
  .Call(dbarts:::C_dbarts_bartcore_storeState, sampler$getPointer())
}
cutPointsOf <- function(sampler) attr(stateOf(sampler), "cutPoints")
treeFields <- c("tree.vars", "tree.values", "tree.sizes", "tree.flags")
# A tree in pre-order: a column, counted from one, for a split and -1 for a
# leaf, and a threshold or a leaf value. A row goes left at or below the
# threshold.
tree <- function(vars, values) {
  list(vars = as.integer(vars), values = as.double(values))
}
stump <- function(value = 0) tree(-1L, value)
split1 <- function(column, point, low, high) {
  tree(c(column, -1L, -1L), c(point, low, high))
}
# a state whose live trees of a chain's mean forest, or of its variance
# forest, are the hand-built ones given, as many as they are
withTrees <- function(state, trees, chain = 1L, variance = FALSE) {
  vars <- unlist(lapply(trees, `[[`, "vars"))
  built <- list(
    vars,
    writeBin(unlist(lapply(trees, `[[`, "values")), raw()),
    lengths(lapply(trees, `[[`, "vars")),
    as.raw(ifelse(vars > 0L, 2L, 0L))
  )
  if (variance) {
    state[[chain]][paste0(
      "variance.",
      c("vars", "values", "sizes", "flags")
    )] <-
      built
  } else {
    state[[chain]]$forests[[1L]][treeFields] <- built
  }
  state
}
# the live trees of a chain's mean or variance forest, as `tree` writes them
treesOf <- function(state, chain = 1L, variance = FALSE) {
  held <- if (variance) {
    unname(state[[chain]][paste0("variance.", c("vars", "values", "sizes"))])
  } else {
    unname(state[[chain]]$forests[[1L]][treeFields[1L:3L]])
  }
  values <- readBin(held[[2L]], "double", length(held[[1L]]))
  ends <- cumsum(held[[3L]])
  lapply(seq_along(ends), function(t) {
    nodes <- (ends[t] - held[[3L]][t] + 1L):ends[t]
    tree(held[[1L]][nodes], values[nodes])
  })
}
# the trees of a kept draw, from getTrees
drawOf <- function(sampler, sample) {
  held <- sampler$getTrees(chainNums = 1L, sampleNums = sample)
  unname(lapply(split(held, held$tree), function(t) tree(t$var, t$value)))
}
# a tree's value at each row of `data`
fitsOf <- function(t, data) {
  fits <- numeric(nrow(data))
  at <- 0L
  walk <- function(rows) {
    at <<- at + 1L
    node <- at
    if (t$vars[node] < 0L) {
      fits[rows] <<- t$values[node]
      return(invisible())
    }
    left <- data[rows, t$vars[node]] <= t$values[node]
    walk(rows[left])
    walk(rows[!left])
  }
  walk(seq_len(nrow(data)))
  fits
}
# the variance over the rows that weigh more than 0, each by its weight
rowVariance <- function(values, weights) {
  keep <- weights > 0
  values <- values[keep]
  weights <- weights[keep]
  center <- sum(weights * values) / sum(weights)
  sum(weights * (values - center)^2) / sum(weights)
}
# The trees of `trees` a forest of `slots` keeps, by hand: the most varying
# over the rows of `data`, a tie to the earlier, in the donor's order. A
# variance tree varies by the log of its factors.
keptByHand <- function(
  trees,
  slots,
  data = x,
  weights = rep(1, nrow(data)),
  scale = identity
) {
  varies <- vapply(
    trees,
    function(t) rowVariance(scale(fitsOf(t, data)), weights),
    0
  )
  sort(order(-varies)[seq_len(slots)])
}
lineOf <- function(sampler, donor, samples = NULL) {
  dbarts:::warmStartLine(
    dbarts:::installDonorTrees(sampler, donor, samples, TRUE)
  )
}

set.seed(47L, sample.kind = "Rejection")
n <- 300L
columns <- paste0("x", 1L:4L)
x <- matrix(runif(4L * n), n, 4L, dimnames = list(NULL, columns))
truth <- function(x) 3 * sin(3 * x[, 1L]) + 2 * x[, 2L]^2 + x[, 3L]
y <- truth(x) + rnorm(n, 0, 0.3)
xOther <- matrix(runif(4L * n), n, 4L, dimnames = list(NULL, columns))
yOther <- truth(xOther) + rnorm(n, 0, 0.3)
controlWith <- function(...) {
  settings <- list(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 4L,
    n.samples = 3L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 7L
  )
  do.call(dbarts::dbartsControl, utils::modifyList(settings, list(...)))
}
make <- function(..., data = x, response = y, sweeps = 0L, control = list()) {
  sampler <- dbarts::dbarts(
    data,
    response,
    control = do.call(controlWith, control),
    ...
  )
  if (sweeps > 0L) {
    invisible(sampler$run(sweeps, 3L))
  }
  sampler
}
trees <- function(count) list(n.trees = count)
grid <- cutPointsOf(make())
at <- function(column, k) grid[[column]][k]
# a donor of the hand trees given, as a state
handDonor <- function(given, ...) {
  withTrees(stateOf(make(..., control = trees(length(given)))), given)
}
# forced into `slots` trees: the trees it then holds
forcedInto <- function(donor, slots, ...) {
  recipient <- make(..., control = trees(slots))
  expect_identical(answerOf(recipient$installTrees(donor)), visibly(FALSE))
  expect_null(recipient$installTrees(donor, forceUpdate = TRUE))
  treesOf(stateOf(recipient))
}

# ---- more trees: all of the donor's, then single leaves at 0 --------------------

donor50 <- make(sweeps = 200L, control = trees(50L))
given50 <- stateOf(donor50)
into200 <- make(control = trees(200L))
expect_identical(answerOf(into200$installTrees(donor50)), visibly(TRUE))
held <- treesOf(stateOf(into200))
expect_identical(held[1L:50L], treesOf(given50))
expect_true(all(vapply(held[51L:200L], identical, NA, stump(0))))
expect_equal(into200$predict(x), donor50$predict(x), tolerance = 1e-12)

# ---- fewer trees: declined, untouched ----------------------------------------

donor200 <- make(sweeps = 50L, control = trees(200L))
into50 <- make(control = trees(50L))
before <- stateOf(into50)
expect_identical(answerOf(into50$installTrees(donor200)), visibly(FALSE))
expect_identical(stateOf(into50), before)
expect_null(into50$installTrees(donor200, forceUpdate = TRUE))
expect_identical(
  treesOf(stateOf(into50)),
  treesOf(stateOf(donor200))[keptByHand(treesOf(stateOf(donor200)), 50L)]
)

# ---- the trees that vary most over the rows ----------------------------------

# six trees into three, the leaves holding unequal shares of the rows: by
# the rows trees 6, 1 and 5 vary most, by the leaves alone 1, 4 and 6, and
# they sit in the donor's order, not the order of their variation
six <- list(
  split1(1L, at(1L, 10L), 2.0, 0.0),
  split1(1L, at(1L, 50L), 0.5, -0.5),
  split1(1L, at(1L, 50L), 0.4, -0.4),
  split1(1L, at(1L, 90L), 0.0, 1.5),
  split1(1L, at(1L, 30L), 0.6, -0.6),
  split1(1L, at(1L, 50L), 0.7, -0.7)
)
expect_identical(keptByHand(six, 3L), c(1L, 5L, 6L))
byLeaves <- order(-vapply(six, function(t) var(t$values[-1L]), 0))
expect_false(identical(sort(byLeaves[1L:3L]), c(1L, 5L, 6L)))
expect_identical(forcedInto(handDonor(six), 3L), six[c(1L, 5L, 6L)])

# by variance, not by size: two leaves at 5.0 and 5.1 vary least
offset <- list(
  split1(1L, at(1L, 50L), 5.0, 5.1),
  split1(1L, at(1L, 50L), 0.1, -0.1),
  split1(1L, at(1L, 50L), 0.12, -0.12),
  split1(1L, at(1L, 50L), 0.15, -0.15)
)
expect_identical(forcedInto(handDonor(offset), 3L), offset[2L:4L])

# seven into three: the round of trees 4 to 6 holds the two that vary most
# and the last round, tree 7 and two single leaves, the third
seven <- lapply(
  c(0.1, 0.2, 0.15, 0.9, 0.8, 0.05, 0.5),
  function(a) split1(1L, at(1L, 50L), a, -a)
)
expect_identical(forcedInto(handDonor(seven), 3L), seven[c(4L, 5L, 7L)])

# a tie goes to the earlier tree
tied <- list(
  split1(1L, at(1L, 50L), 0.3, -0.3),
  split1(1L, at(1L, 50L), -0.3, 0.3)
)
expect_identical(forcedInto(handDonor(tied), 1L), tied[1L])

# ---- case weights rank, the working weights do not ------------------------------

# rows below the 20th point of x3 weigh 0 and rows below the 50th of x2
# weigh 10: tree 1 then varies most, tree 2 does with the rows of weight 0
# out and the others at 1, and tree 3 does at unit weights
weights <- ifelse(
  x[, 3L] <= at(3L, 20L),
  0,
  ifelse(x[, 2L] <= at(2L, 50L), 10, 1)
)
weighed <- list(
  split1(1L, at(1L, 50L), 0.1, -0.1),
  split1(2L, at(2L, 50L), 0.12, -0.12),
  split1(3L, at(3L, 20L), 3.0, 0.0)
)
expect_identical(keptByHand(weighed, 1L, weights = weights), 1L)
expect_identical(keptByHand(weighed, 1L, weights = as.numeric(weights > 0)), 2L)
expect_identical(keptByHand(weighed, 1L), 3L)
expect_warning(make(weights = weights), "weights' of 0", fixed = TRUE)
expect_identical(
  suppressWarnings(forcedInto(handDonor(weighed), 1L, weights = weights)),
  weighed[1L]
)

# a logistic sampler that has run, whose working weights are its latent
# draws: the choice is the one unit weights give, on every seed. Two trees
# close at unit weights, each splitting a fifth of the rows from the rest
yBinary <- as.integer(y > median(y))
close <- list(
  split1(1L, at(1L, 20L), 0.300, -0.300),
  split1(2L, at(2L, 20L), 0.300, -0.300)
)
byUnit <- keptByHand(close, 1L)
for (seed in 1L:6L) {
  logistic <- make(
    response = yBinary,
    family = "logistic",
    control = list(n.trees = 1L, seed = seed)
  )
  invisible(logistic$run(10L, 1L))
  expect_null(logistic$installTrees(handDonor(close), forceUpdate = TRUE))
  expect_identical(treesOf(stateOf(logistic)), close[byUnit], info = seed)
}

# ---- ranked as installed --------------------------------------------------------

# under a limit of one variable the first tree is cut back to its root,
# whose two leaves sit near 0: it varies least as installed, most as stored
crossed <- list(
  tree(
    c(1L, 2L, -1L, -1L, 2L, -1L, -1L),
    c(at(1L, 50L), at(2L, 50L), 1, -1, at(2L, 50L), -1, 1)
  ),
  split1(1L, at(1L, 50L), 0.2, -0.2)
)
expect_identical(keptByHand(crossed, 1L), 1L)
expect_identical(
  forcedInto(
    handDonor(crossed),
    1L,
    interactions = interactions(max.order = 1)
  ),
  crossed[2L]
)

# under blocks a tree is ranked as cut back for the slot it is tried in and
# installed in the slot its rank gives, cut back for that one: slots 1 and
# 2 split on x1 alone, 3 and 4 on the others. Trees 3 and 4, tried in
# slots 3 and 4, and 5 and 6, tried in 1 and 2, are kept and land in the
# other block's slots, where nothing they split on is allowed
blocked <- list(
  split1(1L, at(1L, 50L), 0.01, -0.01),
  split1(1L, at(1L, 50L), 0.01, -0.01),
  split1(2L, at(2L, 50L), 0.5, -0.5),
  split1(2L, at(2L, 50L), 0.5, -0.5),
  split1(1L, at(1L, 50L), 0.4, -0.4),
  split1(1L, at(1L, 50L), 0.4, -0.4)
)
held <- forcedInto(
  handDonor(blocked),
  4L,
  blocks = blocks(
    groups = list("x1", c("x2", "x3", "x4")),
    trees.per.group = c(2L, 2L)
  )
)
splitsOn <- lapply(held, function(t) t$vars[t$vars > 0L])
expect_true(all(unlist(splitsOn[1L:2L]) == 1L))
expect_true(all(unlist(splitsOn[3L:4L]) != 1L))

# across cut points each tree is tried as it is moved: onto other rows the
# trees kept are the most varying of the donor's trees as a recipient of as
# many holds them after the move
moved <- make(data = xOther, response = yOther, control = trees(50L))
expect_null(moved$installTrees(donor50, forceUpdate = TRUE))
movedTrees <- treesOf(stateOf(moved))
expect_false(identical(movedTrees, treesOf(given50)))
expect_identical(
  forcedInto(donor50, 20L, data = xOther, response = yOther),
  movedTrees[keptByHand(movedTrees, 20L, data = xOther)]
)

# ---- the donor sample each chain takes -----------------------------------------

keeper <- make(
  sweeps = 100L,
  control = list(n.trees = 10L, keepTrees = TRUE, seed = 10L)
)
draws <- lapply(1L:3L, function(s) drawOf(keeper, s))
twoChains <- make(control = list(n.trees = 4L, n.chains = 2L))
expect_null(
  twoChains$installTrees(keeper, samples = c(1L, 3L), forceUpdate = TRUE)
)
for (chain in 1L:2L) {
  draw <- draws[[c(1L, 3L)[chain]]]
  expect_identical(
    treesOf(stateOf(twoChains), chain),
    draw[keptByHand(draw, 4L)],
    info = chain
  )
}
expect_false(identical(
  keptByHand(draws[[1L]], 4L),
  keptByHand(draws[[3L]], 4L)
))

# a kept draw of a donor of other count is that draw's trees, into fewer
# trees and into more
fewer <- make(control = trees(6L))
expect_null(fewer$installTrees(keeper, samples = 2L, forceUpdate = TRUE))
expect_identical(
  treesOf(stateOf(fewer)),
  draws[[2L]][keptByHand(draws[[2L]], 6L)]
)
more <- make(control = trees(13L))
expect_identical(
  answerOf(more$installTrees(keeper, samples = 2L)),
  visibly(TRUE)
)
expect_identical(treesOf(stateOf(more))[1L:10L], draws[[2L]])

# ---- the variance forest ------------------------------------------------------

heteroscedastic <- function(meanTrees, varianceTrees, ...) {
  make(
    ...,
    variance = varianceForest(n.trees = varianceTrees),
    control = trees(meanTrees)
  )
}
# ranked by the log of the factors: (0.25, 1) varies more than (1, 2.2)
# there, less on their own scale
factors <- list(
  split1(1L, at(1L, 50L), 0.25, 1),
  split1(1L, at(1L, 50L), 1, 2.2)
)
expect_identical(keptByHand(factors, 1L, scale = log), 1L)
expect_identical(keptByHand(factors, 1L), 2L)
logDonor <- withTrees(
  stateOf(heteroscedastic(4L, 2L)),
  factors,
  variance = TRUE
)
logRecipient <- heteroscedastic(4L, 1L)
expect_null(logRecipient$installTrees(logDonor, forceUpdate = TRUE))
expect_identical(
  treesOf(stateOf(logRecipient), variance = TRUE),
  factors[1L]
)

# each forest its own choice: the mean trees 1 and 3, the variance trees 1
# and 2
meanThree <- lapply(c(0.5, 0.1, 0.4), function(a) {
  split1(1L, at(1L, 50L), a, -a)
})
varianceThree <- lapply(c(2.0, 3.0, 1.1), function(h) {
  split1(1L, at(1L, 50L), 1, h)
})
both <- withTrees(
  withTrees(stateOf(heteroscedastic(3L, 3L)), meanThree),
  varianceThree,
  variance = TRUE
)
twoAndTwo <- heteroscedastic(2L, 2L)
expect_null(twoAndTwo$installTrees(both, forceUpdate = TRUE))
expect_identical(treesOf(stateOf(twoAndTwo)), meanThree[c(1L, 3L)])
expect_identical(
  treesOf(stateOf(twoAndTwo), variance = TRUE),
  varianceThree[c(1L, 2L)]
)

# four variance trees into six: the rest are leaves at 1, so the variance is
# the donor's
fourDonor <- heteroscedastic(4L, 4L, sweeps = 50L)
intoSix <- heteroscedastic(4L, 6L)
expect_identical(answerOf(intoSix$installTrees(fourDonor)), visibly(TRUE))
expect_equal(intoSix$getVariance(), fourDonor$getVariance(), tolerance = 1e-12)
expect_identical(
  treesOf(stateOf(intoSix), variance = TRUE)[5L:6L],
  list(stump(1), stump(1))
)

# ---- what bart says --------------------------------------------------------------

expect_identical(
  lineOf(make(control = trees(200L)), donor50),
  "warm start: 50 of 50 donor trees taken; 150 trees started as single leaves"
)
expect_identical(
  lineOf(make(control = trees(50L)), donor200),
  "warm start: 50 of 200 donor trees taken"
)
sixDonor <- heteroscedastic(4L, 6L, sweeps = 50L)
expect_identical(
  lineOf(heteroscedastic(4L, 4L), sixDonor),
  "warm start: 4 of 4 donor trees taken; variance trees: 4 of 6 taken"
)
expect_identical(
  lineOf(heteroscedastic(4L, 6L), fourDonor),
  paste0(
    "warm start: 4 of 4 donor trees taken; variance trees: 4 of 4 taken, ",
    "2 started as single leaves"
  )
)
