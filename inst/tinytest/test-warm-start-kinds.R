# A warm start from a donor of another tree count, kind of leaf or set of
# parts. Into a forest of more
# trees the donor's go in all, in its order, and the rest start as single
# leaves at 0, or at 1 in a variance forest: the donor as it stands, TRUE.
# Into a forest of fewer the donor is declined, FALSE, and forced takes the
# trees whose fitted values vary most over this sampler's rows, each row
# weighing by its case weight, in the donor's order and as installed; the
# mean and the variance forest choose apart, the variance forest by the log
# of its factors, and each chain among the trees of the donor sample it
# takes. DART or a variance forest on one side only, and constant against
# linear leaves, are taken as the last sections say.

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

# linear and gp leaves start their single leaves with slopes and fits of 0:
# six trees into nine are the donor's fit
leaves <- function(kind, ...) {
  if (kind == "linear") {
    make(..., leaf.prior = linear(columns = 2L))
  } else {
    make(..., leaf.prior = gp(columns = 2L))
  }
}
for (kind in c("linear", "gp")) {
  sixLeaves <- leaves(kind, control = trees(6L), sweeps = 30L)
  nineLeaves <- leaves(kind, control = trees(9L))
  expect_identical(
    answerOf(nineLeaves$installTrees(sixLeaves)),
    visibly(TRUE),
    info = kind
  )
  expect_equal(
    nineLeaves$predict(x),
    sixLeaves$predict(x),
    tolerance = 1e-10,
    info = kind
  )
}

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

# a row the active-row mask leaves out weighs 0, gaussian and probit: tree 1
# then varies most, tree 2 with the rows under x3's 20th point in
yBinary <- as.integer(y > median(y))
masked <- weighed[c(1L, 3L)]
mask <- as.numeric(x[, 3L] > at(3L, 20L))
expect_identical(keptByHand(masked, 1L, weights = mask), 1L)
expect_identical(keptByHand(masked, 1L), 2L)
for (family in c("gaussian", "probit")) {
  activeRows <- make(
    response = if (family == "gaussian") y else yBinary,
    family = family,
    control = trees(1L)
  )
  activeRows$setActiveRows(mask)
  expect_null(activeRows$installTrees(handDonor(masked), forceUpdate = TRUE))
  expect_identical(treesOf(stateOf(activeRows)), masked[1L], info = family)
}

# a Student-t sampler weighs rows by its weights and a logistic one by its
# counts: the rows of weight 0 at a count of 1, tree 1 against tree 2 flips
expect_identical(
  suppressWarnings(
    forcedInto(handDonor(weighed), 1L, family = "student", weights = weights)
  ),
  weighed[1L]
)
counts <- ifelse(weights == 0, 1, weights)
counted <- weighed[1L:2L]
expect_identical(keptByHand(counted, 1L, weights = counts), 1L)
expect_identical(keptByHand(counted, 1L), 2L)
expect_identical(
  forcedInto(
    handDonor(counted),
    1L,
    response = yBinary,
    family = "logistic",
    weights = counts
  ),
  counted[1L]
)

# a logistic sampler that has run, whose working weights are its latent
# draws: the choice is the one unit weights give, on every seed. The
# response is 1 wherever x1 is above its 50th point, so the fit is large
# there and the draws small, under which tree 2 would vary most; at unit
# weights tree 1 does
far <- x[, 1L] > at(1L, 50L)
set.seed(5L)
ySure <- ifelse(far, 1L, rbinom(n, 1L, 0.5))
close <- list(
  split1(1L, at(1L, 80L), 0.30, -0.30),
  split1(2L, at(2L, 20L), 0.32, -0.32)
)
expect_identical(keptByHand(close, 1L), 1L)
expect_identical(keptByHand(close, 1L, weights = ifelse(far, 0.5, 1)), 2L)
for (seed in 1L:2L) {
  logistic <- make(
    response = ySure,
    family = "logistic",
    control = list(n.trees = 1L, seed = seed)
  )
  invisible(logistic$run(200L, 1L))
  draws <- stateOf(logistic)[[1L]]$latents
  expect_true(mean(draws[far]) < 0.7 * mean(draws[!far]), info = seed)
  expect_null(logistic$installTrees(handDonor(close), forceUpdate = TRUE))
  expect_identical(treesOf(stateOf(logistic)), close[1L], info = seed)
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

# without blocks a tree installed as a single leaf is ranked like any
# other: into a forest increasing in x1, tree 2, decreasing, collapses to one
# leaf and ties at 0 with tree 3, a single leaf, and the earlier, tree 2, is
# kept at the mean of its fits over the rows
falling <- list(
  split1(1L, at(1L, 50L), -0.5, 0.5),
  split1(1L, at(1L, 50L), 0.6, 0.2),
  stump(0.2)
)
held <- forcedInto(handDonor(falling), 2L, monotone = c(x1 = "increasing"))
expect_identical(held[[1L]], falling[[1L]])
expect_identical(held[[2L]]$vars, -1L)
expect_equal(held[[2L]]$values, mean(fitsOf(falling[[2L]], x)))

# under blocks each tree is tried in the slots of every block, and each
# block in turn keeps its count of the trees that vary most there among
# those no earlier block kept: slots 1 and 2 split on x1 alone, 3 and 4 on
# the others. Trees 5 and 6 vary most on x1 and go to the first block, trees
# 3 and 4 on x2 to the second, each as stored
blocked <- list(
  split1(1L, at(1L, 50L), 0.01, -0.01),
  split1(1L, at(1L, 50L), 0.01, -0.01),
  split1(2L, at(2L, 50L), 0.5, -0.5),
  split1(2L, at(2L, 50L), 0.5, -0.5),
  split1(1L, at(1L, 50L), 0.4, -0.4),
  split1(1L, at(1L, 50L), 0.4, -0.4)
)
twoBlocks <- function(first, second, groups = list("x1", c("x2", "x3", "x4"))) {
  dbartsForests$blocks(groups = groups, trees.per.group = c(first, second))
}
expect_identical(
  forcedInto(handDonor(blocked), 4L, blocks = twoBlocks(2L, 2L)),
  blocked[c(5L, 6L, 3L, 4L)]
)

# a tree a block would hold as a single leaf is not kept there, and a block
# with too few others starts the rest as single leaves at 0: into 1 and 3 the
# first block keeps tree 1, and the second trees 3 and 4 and then a leaf at
# 0, not tree 2 cut back to a leaf
blocked[c(1L, 2L, 5L, 6L)] <- list(
  split1(1L, at(1L, 50L), 0.6, 0.0),
  split1(1L, at(1L, 50L), 0.01, -0.01),
  split1(1L, at(1L, 50L), 0.02, -0.02),
  split1(1L, at(1L, 50L), 0.02, -0.02)
)
expect_identical(
  forcedInto(handDonor(blocked), 4L, blocks = twoBlocks(1L, 3L)),
  c(blocked[c(1L, 3L, 4L)], list(stump(0)))
)

# trees that fit no block, each splitting on x1 and then on x2: declined,
# and forced the first block keeps the tree that varies most as cut back
# there, tree 2 and not tree 3, which varies most as stored, its leaves on
# x2 meaning 0 over their rows; the second block, where each would be a
# single leaf, starts as single leaves
right <- x[, 1L] > at(1L, 50L)
low <- mean(x[right, 2L] <= at(2L, 50L))
across <- lapply(list(c(0.1, 1), c(0.4, 1), c(0.2, 4)), function(v) {
  tree(
    c(1L, -1L, 2L, -1L, -1L),
    c(at(1L, 50L), v[1L], at(2L, 50L), v[2L] * (1 - low), -v[2L] * low)
  )
})
expect_identical(keptByHand(across, 1L), 3L)
held <- forcedInto(handDonor(across), 3L, blocks = twoBlocks(1L, 2L))
expect_identical(held[[1L]]$vars, c(1L, -1L, -1L))
expect_identical(held[[1L]]$values[2L], 0.4)
expect_identical(held[2L:3L], list(stump(0), stump(0)))
expect_identical(
  lineOf(
    make(blocks = twoBlocks(1L, 2L), control = trees(3L)),
    handDonor(across)
  ),
  paste(
    "warm start: 1 of 3 donor trees taken; 2 trees started as single leaves;",
    "1 cut back to fit"
  )
)

# chain 1 holds trees each block holds as stored and chain 2 the trees
# above, which fit no block: declined, and forced chain 1 goes in as placed
# and chain 2 by rank, its first block the tree that splits on x1 and then
# on x2 cut back, and tree 3, its second trees 1 and 4
placeable <- list(
  split1(1L, at(1L, 50L), 0.5, -0.5),
  stump(0.1),
  split1(2L, at(2L, 50L), 0.4, -0.4),
  split1(3L, at(3L, 50L), 0.2, -0.2)
)
mixed <- tree(
  c(1L, 2L, -1L, -1L, -1L),
  c(at(1L, 50L), at(2L, 50L), 0.3, -0.1, 0.5)
)
unplaceable <- list(
  split1(2L, at(2L, 40L), 0.3, -0.3),
  mixed,
  split1(1L, at(1L, 30L), 0.2, -0.2),
  split1(3L, at(3L, 50L), 0.25, -0.25)
)
twoChains <- list(n.trees = 4L, n.chains = 2L)
chained <- withTrees(
  withTrees(stateOf(make(control = twoChains)), placeable, 1L),
  unplaceable,
  2L
)
recipient <- make(blocks = twoBlocks(2L, 2L), control = twoChains)
expect_identical(answerOf(recipient$installTrees(chained)), visibly(FALSE))
expect_null(recipient$installTrees(chained, forceUpdate = TRUE))
expect_identical(treesOf(stateOf(recipient), 1L), placeable)
held <- treesOf(stateOf(recipient), 2L)
expect_identical(held[2L:4L], unplaceable[c(3L, 1L, 4L)])
expect_identical(held[[1L]]$vars, c(1L, -1L, -1L))
expect_identical(held[[1L]]$values[3L], 0.5)
# where chain 2's trees each split on x1 and then on x2, its first block
# keeps two of them cut back and its second starts as single leaves: the line
# gives the counts that differ between chains as ranges
splitTwice <- lapply(1L:4L, function(i) {
  tree(
    c(1L, 2L, -1L, -1L, -1L),
    c(at(1L, 50L), at(2L, 50L), 0.1 * i, -0.1 * i, 0.2 * i)
  )
})
expect_identical(
  lineOf(
    make(blocks = twoBlocks(2L, 2L), control = twoChains),
    withTrees(chained, splitTwice, 2L)
  ),
  paste(
    "warm start: 2 to 4 of 4 donor trees taken; 0 to 2 trees started as",
    "single leaves; 2 cut back to fit, over 2 chains"
  )
)

# a donor whose single leaves come early in its order: placed in the order
# given the leaf would take the first block's only slot, so the trees that
# split are placed first, the leaf after them, each block's in the donor's
# order: TRUE, every tree as stored
early <- list(
  stump(0.3),
  split1(1L, at(1L, 50L), 0.5, -0.5),
  split1(2L, at(2L, 50L), 0.4, -0.4)
)
intoEarly <- make(blocks = twoBlocks(1L, 2L), control = trees(3L))
expect_identical(
  answerOf(intoEarly$installTrees(handDonor(early))),
  visibly(TRUE)
)
expect_identical(treesOf(stateOf(intoEarly)), early[c(2L, 1L, 3L)])

# a tree a block keeps is not kept again by a later one: into 2 and 2 the
# first block keeps tree 1 and the single leaf tree 2, the earliest of those
# that vary by 0 there, and the second tree 3 and then tree 4, the earliest
# such left
stumped <- list(
  split1(1L, at(1L, 50L), 0.6, 0.0),
  stump(0.3),
  split1(2L, at(2L, 50L), 0.5, -0.5),
  stump(0.7),
  stump(0.9)
)
expect_identical(
  forcedInto(handDonor(stumped), 4L, blocks = twoBlocks(2L, 2L)),
  stumped[1L:4L]
)

# a donor of 100 and 100 trees on x1 and x2 and on x3 and x4: into 25 and 25,
# and into 40 and 10, each block keeps the most varying of the donor's trees
# on its own columns, from both halves of the donor
halves <- list(c("x1", "x2"), c("x3", "x4"))
blockDonor <- make(
  sweeps = 200L,
  blocks = twoBlocks(100L, 100L, halves),
  control = trees(200L)
)
blockTrees <- treesOf(stateOf(blockDonor))
for (counts in list(c(25L, 25L), c(40L, 10L))) {
  kept <- c(
    keptByHand(blockTrees[1L:100L], counts[1L]),
    100L + keptByHand(blockTrees[101L:200L], counts[2L])
  )
  expect_identical(
    forcedInto(
      blockDonor,
      sum(counts),
      blocks = twoBlocks(counts[1L], counts[2L], halves)
    ),
    blockTrees[kept],
    info = paste(counts, collapse = "/")
  )
}

# a donor of the blocks and counts it goes into goes in as it stands, each
# tree in its own slot
same <- make(blocks = twoBlocks(100L, 100L, halves), control = trees(200L))
expect_identical(answerOf(same$installTrees(blockDonor)), visibly(TRUE))
expect_identical(treesOf(stateOf(same)), blockTrees)

# a donor of 20 and 20 into 25 and 25 is taken whole: the first block takes
# the donor's first 20 trees and then, in its order, single leaves of the
# second half, up to 5, the second block the rest, each block then single
# leaves at 0
small <- make(
  sweeps = 100L,
  blocks = twoBlocks(20L, 20L, halves),
  control = trees(40L)
)
smallTrees <- treesOf(stateOf(small))
into25 <- make(blocks = twoBlocks(25L, 25L, halves), control = trees(50L))
expect_identical(answerOf(into25$installTrees(small)), visibly(TRUE))
leaves <- which(lengths(lapply(smallTrees, `[[`, "vars")) == 1L)
first <- c(1L:20L, head(intersect(leaves, 21L:40L), 5L))
second <- setdiff(21L:40L, first)
padded <- function(kept) {
  c(smallTrees[kept], rep(list(stump(0)), 25L - length(kept)))
}
expect_identical(treesOf(stateOf(into25)), c(padded(first), padded(second)))

# 25 and 25 into 24 and 26, the first half holding a single leaf at tree 10:
# in the donor's order the first block fills with trees 1 to 24 and tree 25
# fits nowhere, so the first block takes the first half's trees that split
# and the second the rest, the single leaves among them: TRUE, taken whole
halfTrees <- treesOf(stateOf(make(
  sweeps = 100L,
  blocks = twoBlocks(25L, 25L, halves),
  control = trees(50L)
)))
halfTrees[[10L]] <- stump(0.05)
leaves <- which(lengths(lapply(halfTrees, `[[`, "vars")) == 1L)
expect_false(25L %in% leaves)
into24 <- make(blocks = twoBlocks(24L, 26L, halves), control = trees(50L))
expect_identical(
  answerOf(into24$installTrees(handDonor(halfTrees))),
  visibly(TRUE)
)
splitting <- setdiff(1L:25L, leaves)
first <- sort(c(splitting, head(leaves, 24L - length(splitting))))
second <- sort(c(setdiff(26L:50L, leaves), setdiff(leaves, first)))
expect_identical(
  treesOf(stateOf(into24)),
  c(
    halfTrees[first],
    halfTrees[second],
    rep(list(stump(0)), 26L - length(second))
  )
)

# into three blocks of 60, on x1, on x2 and on x3 and x4, the trees on x3
# and x4 would be single leaves in the first two and are not kept there, so
# the third keeps the 60 of them that vary most
held <- forcedInto(
  blockDonor,
  180L,
  blocks = dbartsForests$blocks(
    groups = list("x1", "x2", c("x3", "x4")),
    trees.per.group = c(60L, 60L, 60L)
  )
)
expect_identical(
  held[121L:180L],
  blockTrees[100L + keptByHand(blockTrees[101L:200L], 60L)]
)

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

# and a tree that varies most as stored, its narrow band at 3 and -3 merging
# where the move leaves its split no point, varies least as moved
finePoints <- coarse <- grid
finePoints[[1L]] <- c(0.25, 0.5, 0.6, 0.65, 0.7, 0.75)
coarse[[1L]] <- c(0.25, 0.5, 0.75)
band <- list(
  tree(c(1L, -1L, 1L, 1L, -1L, -1L, -1L), c(0.6, 0, 0.7, 0.65, 3, -3, 0)),
  split1(1L, 0.5, 0.5, -0.5)
)
expect_identical(keptByHand(band, 1L), 1L)
bandSource <- make(control = trees(2L))
bandSource$setCutPoints(finePoints)
bandDonor <- withTrees(stateOf(bandSource), band)
onCoarse <- function(count) {
  sampler <- make(control = trees(count))
  sampler$setCutPoints(coarse)
  sampler
}
bandMoved <- onCoarse(2L)
expect_null(bandMoved$installTrees(bandDonor, forceUpdate = TRUE))
movedBand <- treesOf(stateOf(bandMoved))
expect_true(length(movedBand[[1L]]$vars) < length(band[[1L]]$vars))
expect_identical(keptByHand(movedBand, 1L), 2L)
bandOne <- onCoarse(1L)
expect_null(bandOne$installTrees(bandDonor, forceUpdate = TRUE))
expect_identical(treesOf(stateOf(bandOne)), movedBand[2L])

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

# a live forest whose count does not divide its kept draws names no draw
# whole: refused, not a slice across two draws
uneven <- withTrees(stateOf(keeper), lapply(1L:4L, function(i) stump(0.01 * i)))
expect_error(
  make(control = trees(4L))$installTrees(
    uneven,
    samples = 2L,
    forceUpdate = TRUE
  ),
  "not shape-compatible",
  fixed = TRUE
)

# a pooled factor's masks go with the trees chosen: the masks of the four
# kept are the donor's for those trees, in their order
set.seed(3L)
pooled <- data.frame(
  x,
  f = factor(sample(paste0("l", 1L:70L), n, replace = TRUE))
)
pooledDonor <- make(
  data = pooled,
  response = y + as.integer(pooled$f) %% 3L,
  sweeps = 60L,
  control = trees(8L)
)
pooledFour <- make(data = pooled, response = y, control = trees(4L))
expect_null(pooledFour$installTrees(pooledDonor, forceUpdate = TRUE))
masksOf <- function(sampler) {
  held <- stateOf(sampler)[[1L]]$forests[[1L]]$tree.masks
  splits <- vapply(treesOf(stateOf(sampler)), function(t) sum(t$vars == 5L), 0L)
  words <- length(held) / sum(splits)
  split(held, rep(seq_along(splits), splits * words))
}
donorTrees <- treesOf(stateOf(pooledDonor))
heldTrees <- treesOf(stateOf(pooledFour))
which4 <- vapply(heldTrees, function(t) match(list(t), donorTrees), 0L)
donorMasks <- masksOf(pooledDonor)
expect_false(anyNA(which4))
expect_true(sum(vapply(heldTrees, function(t) sum(t$vars == 5L), 0L)) > 0L)
expect_identical(
  unname(unlist(masksOf(pooledFour))),
  unname(unlist(donorMasks[as.character(which4)]))
)

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

# the variance forest needs rounds where the mean forest needs none: three
# variance trees into one, the third varying most
varianceOnly <- lapply(c(1.1, 1.2, 3), function(h) {
  split1(1L, at(1L, 50L), 1, h)
})
roundsDonor <- withTrees(
  stateOf(heteroscedastic(4L, 3L)),
  varianceOnly,
  variance = TRUE
)
roundsRecipient <- heteroscedastic(4L, 1L)
expect_null(roundsRecipient$installTrees(roundsDonor, forceUpdate = TRUE))
expect_identical(
  treesOf(stateOf(roundsRecipient), variance = TRUE),
  varianceOnly[3L]
)

# the variance trees rank by case weights: rows at or below x1's 50th point
# weigh 10, under which the tree on x2 varies least
heavy <- ifelse(x[, 1L] <= at(1L, 50L), 10, 1)
weighedFactors <- list(
  split1(1L, at(1L, 50L), 1, 2),
  split1(2L, at(2L, 50L), 1, 1.9)
)
expect_identical(keptByHand(weighedFactors, 1L, scale = log), 1L)
expect_identical(
  keptByHand(weighedFactors, 1L, weights = heavy, scale = log),
  2L
)
weighedDonor <- withTrees(
  stateOf(heteroscedastic(4L, 2L, weights = heavy)),
  weighedFactors,
  variance = TRUE
)
weighedRecipient <- heteroscedastic(4L, 1L, weights = heavy)
expect_null(weighedRecipient$installTrees(weighedDonor, forceUpdate = TRUE))
expect_identical(
  treesOf(stateOf(weighedRecipient), variance = TRUE),
  weighedFactors[2L]
)

# nor do the precisions a variance forest serves: one at 100 and 0.01 on the
# two halves of x1 would choose tree 1, the rows at their case weights choose
# tree 2
pair <- list(
  split1(1L, at(1L, 20L), 0.30, -0.30),
  split1(2L, at(2L, 20L), 0.31, -0.31)
)
precision <- ifelse(x[, 1L] <= at(1L, 50L), 100, 0.01)
expect_identical(keptByHand(pair, 1L), 2L)
expect_identical(keptByHand(pair, 1L, weights = precision), 1L)
pairDonor <- withTrees(
  withTrees(stateOf(heteroscedastic(2L, 1L)), pair),
  list(split1(1L, at(1L, 50L), 0.01, 100)),
  variance = TRUE
)
pairRecipient <- heteroscedastic(1L, 1L)
expect_null(pairRecipient$installTrees(pairDonor, forceUpdate = TRUE))
expect_identical(treesOf(stateOf(pairRecipient)), pair[2L])

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
# single leaves are said before the changes, as in the variance clause
expect_true(grepl(
  paste0(
    "^warm start: 50 of 50 donor trees taken; ",
    "150 trees started as single leaves; [0-9]+ moved"
  ),
  lineOf(make(data = xOther, response = yOther, control = trees(200L)), donor50)
))
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

# ---- parts one side has ----------------------------------------------------------

# DART's weights go only between two DART fits and a variance forest only
# where both have one, neither judged: a part the donor alone has is left
# behind, and one this sampler alone has starts as a new fit's once the
# donor goes in, never on a decline. DART starts with its delay spent, so
# that the next sweep draws its weights from the trees taken
forest <- dbarts::dbartsForests$forest
dartOf <- function(state) {
  unname(state[[1L]][c(
    "dart.probabilities",
    "dart.alpha",
    "dart.updates.skipped"
  )])
}
onThree <- function(...) {
  make(
    ...,
    tree.prior = dart(update.delay = 5L),
    forests = list(forest(vars = 1L:3L))
  )
}
created <- dartOf(stateOf(onThree()))
expect_identical(created[[1L]], c(1, 1, 1, 0) / 3)
plainDonor <- make(sweeps = 50L)
ran <- onThree(sweeps = 30L)
expect_false(identical(dartOf(stateOf(ran))[1L:2L], created[1L:2L]))
expect_null(ran$installTrees(plainDonor, forceUpdate = TRUE))
expect_identical(dartOf(stateOf(ran)), c(created[1L:2L], list(5L)))
invisible(ran$run(0L, 1L))
weights <- dartOf(stateOf(ran))[[1L]]
expect_true(weights[4L] == 0 && length(unique(weights[1L:3L])) == 3L)
# not judged: a clean plain donor into a DART fit is taken as it stands
expect_identical(
  answerOf(make(tree.prior = dart())$installTrees(plainDonor)),
  visibly(TRUE)
)

# a DART donor into a plain fit: taken, and nothing of its DART read, so the
# sampler is one that took the donor without it
dartDonor <- make(tree.prior = dart(update.delay = 0L), sweeps = 50L)
noDart <- stateOf(dartDonor)
noDart[[1L]][c("dart.probabilities", "dart.alpha", "dart.updates.skipped")] <-
  list(NULL)
intoPlain <- make()
twin <- make()
expect_identical(answerOf(intoPlain$installTrees(dartDonor)), visibly(TRUE))
expect_true(twin$installTrees(noDart))
expect_identical(stateOf(intoPlain), stateOf(twin))
expect_null(stateOf(intoPlain)[[1L]]$dart.probabilities)
expect_identical(intoPlain$run(0L, 3L), twin$run(0L, 3L))
# both DART: the donor's weights, concentration and counter
bothDart <- make(tree.prior = dart(update.delay = 5L), sweeps = 10L)
expect_identical(answerOf(bothDart$installTrees(dartDonor)), visibly(TRUE))
expect_identical(dartOf(stateOf(bothDart)), dartOf(stateOf(dartDonor)))

# a variance forest this sampler alone has: every tree a root at the
# starting variance and every row's variance the starting variance, as a new
# fit's, after a run that moved both
fresh <- heteroscedastic(4L, 4L)
startTrees <- treesOf(stateOf(fresh), variance = TRUE)
ranVariance <- heteroscedastic(4L, 4L, sweeps = 30L)
expect_false(identical(ranVariance$getVariance(), fresh$getVariance()))
expect_identical(
  answerOf(ranVariance$installTrees(plainDonor)),
  visibly(TRUE)
)
expect_identical(treesOf(stateOf(ranVariance), variance = TRUE), startTrees)
expect_identical(ranVariance$getVariance(), fresh$getVariance())
expect_identical(treesOf(stateOf(ranVariance)), treesOf(stateOf(plainDonor)))
# one the donor alone has: the mean trees taken, the variance trees nowhere
intoHomoscedastic <- make()
expect_identical(
  answerOf(intoHomoscedastic$installTrees(fourDonor)),
  visibly(TRUE)
)
expect_identical(
  treesOf(stateOf(intoHomoscedastic)),
  treesOf(stateOf(fourDonor))
)
expect_null(stateOf(intoHomoscedastic)[[1L]]$variance.vars)

# a decline leaves a sampler's own parts as they were
bothParts <- make(
  tree.prior = dart(update.delay = 5L),
  variance = varianceForest(n.trees = 4L),
  sweeps = 30L
)
offRows <- make(data = xOther, response = yOther, sweeps = 50L)
before <- stateOf(bothParts)
expect_identical(answerOf(bothParts$installTrees(offRows)), visibly(FALSE))
expect_identical(stateOf(bothParts), before)

# the starting variance is read on the response transform in force: after
# setResponse re-anchors it, and after setData, the variance forest starts
# where a new fit's on the response in force does
yMoved <- 3 * y + 2
movedFresh <- heteroscedastic(4L, 4L, sigest = 0.4, response = yMoved)
for (move in c("setResponse", "setData")) {
  moved <- heteroscedastic(4L, 4L, sigest = 0.4, sweeps = 30L)
  if (move == "setResponse") {
    moved$setResponse(yMoved, updateScale = TRUE)
  } else {
    moved$setData(dbarts::dbartsData(x, yMoved))
  }
  expect_true(moved$installTrees(plainDonor))
  expect_identical(moved$getVariance(), movedFresh$getVariance())
}

# and where the trees are chosen, from a donor of more trees forced
partsOnSix <- function(sweeps = 0L) {
  heteroscedastic(6L, 4L, tree.prior = dart(update.delay = 5L), sweeps = sweeps)
}
createdParts <- partsOnSix()
chosen <- partsOnSix(sweeps = 30L)
twelveDonor <- make(control = trees(12L), sweeps = 40L)
expect_false(identical(chosen$getVariance(), createdParts$getVariance()))
expect_false(chosen$installTrees(twelveDonor))
expect_null(chosen$installTrees(twelveDonor, forceUpdate = TRUE))
expect_identical(
  dartOf(stateOf(chosen)),
  c(dartOf(stateOf(createdParts))[1L:2L], list(5L))
)
expect_identical(chosen$getVariance(), createdParts$getVariance())

# ---- kinds of leaf -----------------------------------------------------------------

linearOn <- function(columns, ...) {
  make(..., leaf.prior = linear(columns = columns))
}
# a constant donor seeds linear leaves exactly: its trees, slopes of 0
intoLinear <- linearOn(2L)
expect_identical(
  answerOf(intoLinear$installTrees(plainDonor)),
  visibly(TRUE)
)
expect_identical(treesOf(stateOf(intoLinear)), treesOf(stateOf(plainDonor)))
expect_true(all(stateOf(intoLinear)[[1L]]$forests[[1L]]$tree.params == 0))
# a linear donor on another number of covariates, standardized on other
# columns: declined; forced, each leaf the intercept it stores, slopes of 0
oneCovariate <- linearOn(3L, sweeps = 40L)
expect_true(any(stateOf(oneCovariate)[[1L]]$forests[[1L]]$tree.params != 0))
twoCovariates <- linearOn(2L:3L)
expect_identical(
  answerOf(twoCovariates$installTrees(oneCovariate)),
  visibly(FALSE)
)
expect_null(twoCovariates$installTrees(oneCovariate, forceUpdate = TRUE))
expect_identical(
  treesOf(stateOf(twoCovariates)),
  treesOf(stateOf(oneCovariate))
)
expect_true(all(stateOf(twoCovariates)[[1L]]$forests[[1L]]$tree.params == 0))

# gp leaves go only to gp leaves: an error in both forms, by name, the
# sampler untouched
onRows <- seq_len(60L)
small <- function(...) make(..., data = x[onRows, ], response = y[onRows])
refusedBoth <- function(sampler, donor, message) {
  before <- stateOf(sampler)
  for (force in c(FALSE, TRUE)) {
    expect_identical(
      answerOf(sampler$installTrees(donor, forceUpdate = force)),
      message
    )
  }
  expect_identical(stateOf(sampler), before)
}
gpSmall <- function(...) small(..., leaf.prior = gp(columns = 2L))
gpDonor <- gpSmall(sweeps = 10L)
linearDonor <- small(leaf.prior = linear(columns = 2L), sweeps = 10L)
refusedBoth(
  small(),
  gpDonor,
  "warm-start donor has gp leaves and this sampler constant ones"
)
refusedBoth(
  small(leaf.prior = linear(columns = 2L)),
  gpDonor,
  "warm-start donor has gp leaves and this sampler linear ones"
)
refusedBoth(
  gpSmall(),
  small(sweeps = 10L),
  "warm-start donor has constant leaves and this sampler gp ones"
)
refusedBoth(
  gpSmall(),
  linearDonor,
  "warm-start donor has linear leaves and this sampler gp ones"
)
refusedBoth(
  make(data = x[1L:50L, ], response = y[1L:50L], leaf.prior = gp(columns = 2L)),
  gpDonor,
  "warm-start donor's gp leaves were fit to 60 rows and this sampler holds 50"
)
# a numeric donor column where this sampler holds an unordered factor
withFactor <- data.frame(x[, 1L:3L], x4 = factor(cut(x[, 4L], 3L)))
refusedBoth(
  make(data = withFactor),
  plainDonor,
  "warm-start donor's column 4 is numeric and this sampler's an unordered factor"
)
# DART weights not one a predictor are the donor's shape, the length 1 and
# a long vector included
for (length in c(1L, 3L, 5L, 100000L)) {
  misfit <- stateOf(dartDonor)
  misfit[[1L]]$dart.probabilities <- rep(1 / length, length)
  refusedBoth(
    make(tree.prior = dart(update.delay = 0L), sweeps = 5L),
    misfit,
    paste(
      "warm-start donor is not shape-compatible with this sampler (number of",
      "forests or predictors differ, the donor's kept draws do not match its",
      "tree count, or its DART weights are not one per predictor)"
    )
  )
}
# a variance leaf that is not finite is refused on the way across grids too,
# where the build on the donor's points meets it first
acrossDonor <- function(leaf) {
  withTrees(
    stateOf(
      heteroscedastic(4L, 4L, data = xOther, response = yOther, sweeps = 20L)
    ),
    list(stump(leaf), stump(1), stump(1), stump(1)),
    variance = TRUE
  )
}
expect_null(
  heteroscedastic(4L, 4L)$installTrees(acrossDonor(2), forceUpdate = TRUE)
)
refusedBoth(
  heteroscedastic(4L, 4L),
  acrossDonor(Inf),
  paste(
    "warm-start donor's variance trees cannot be installed on this sampler's",
    "data (a scale leaf is not positive and finite, or a flat tree failed to",
    "rebuild); the donor's variance surface is incompatible with this",
    "sampler, and nothing was installed"
  )
)

# ---- what bart says of parts and kinds ------------------------------------------

taken <- "warm start: 4 of 4 donor trees taken; "
expect_identical(
  lineOf(make(tree.prior = dart()), plainDonor),
  paste0(taken, "DART weights started fresh")
)
expect_identical(
  lineOf(make(), dartDonor),
  paste0(taken, "donor's DART weights left behind")
)
expect_identical(
  lineOf(heteroscedastic(4L, 4L), plainDonor),
  paste0(taken, "variance forest started fresh")
)
expect_identical(
  lineOf(make(), fourDonor),
  paste0(taken, "donor's variance forest left behind")
)
expect_identical(
  lineOf(linearOn(2L), plainDonor),
  paste0(taken, "slopes started at 0")
)
expect_identical(
  lineOf(linearOn(2L:3L), oneCovariate),
  paste0(taken, "slopes left behind")
)
