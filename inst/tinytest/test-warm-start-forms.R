# installTrees has two forms. Unforced, a donor every tree of which stands in
# the sampler as the donor holds it is installed, TRUE; one that would have
# to be changed is not, FALSE, and the sampler is as it was. With
# forceUpdate = TRUE the donor goes in with the changes and NULL comes back
# invisibly: a split off the cut points moved or merged as a forced setState
# moves it, a leaf no row reaches and a split its forest does not allow
# merged at the mean of the leaves beneath by the rows each held, a side
# dropped, a monotone tree out of order set to zero. bart's warm.start forces
# and, under verbose, says in one line what was taken and what was changed.

visibly <- function(value) list(value = value, visible = TRUE)
forcedValue <- list(value = NULL, visible = FALSE)
answerOf <- function(expr) tryCatch(withVisible(expr), error = conditionMessage)
# the state the engine holds now, the sampler's state field left alone
stateOf <- function(sampler) {
  .Call(dbarts:::C_dbarts_bartcore_storeState, sampler$getPointer())
}
cutPointsOf <- function(sampler) attr(stateOf(sampler), "cutPoints")
treeFields <- c("tree.vars", "tree.values", "tree.sizes", "tree.flags")
# the live mean trees of every chain as stored, and with their leaf blocks
liveOf <- function(state, fields = treeFields) {
  lapply(state, function(chain) chain$forests[[1L]][fields])
}
leavesOf <- function(state) liveOf(state, c(treeFields, "tree.params"))
# A tree in pre-order: a column, counted from one, for a split and -1 for a
# leaf; a threshold or a leaf value; and whether a split sends a missing
# value right. `stump` is a tree of one leaf.
tree <- function(vars, values, right = logical(length(vars))) {
  list(vars = as.integer(vars), values = as.double(values), right = right)
}
stump <- function(value = 0) tree(-1L, value)
# a state with the live trees of a chain's mean forest, or of its variance
# forest, replaced by hand-built ones, as many as the forest holds
withTrees <- function(state, trees, chain = 1L, variance = FALSE) {
  vars <- unlist(lapply(trees, `[[`, "vars"))
  right <- unlist(lapply(trees, `[[`, "right"))
  built <- list(
    vars,
    writeBin(unlist(lapply(trees, `[[`, "values")), raw()),
    lengths(lapply(trees, `[[`, "vars")),
    as.raw(ifelse(vars > 0L, 2L, 0L) + right)
  )
  if (variance) {
    fields <- paste0("variance.", c("vars", "values", "sizes", "flags"))
    stopifnot(length(trees) == length(state[[chain]]$variance.sizes))
    state[[chain]][fields] <- built
  } else {
    stopifnot(
      length(trees) == length(state[[chain]]$forests[[1L]]$tree.sizes)
    )
    state[[chain]]$forests[[1L]][treeFields] <- built
  }
  state
}
# tree `index` of a chain's mean or variance forest, as `tree` writes one
treeOf <- function(state, index = 1L, chain = 1L, variance = FALSE) {
  held <- if (variance) {
    unname(state[[chain]][paste0("variance.", c("vars", "values", "sizes"))])
  } else {
    unname(state[[chain]]$forests[[1L]][treeFields[1L:3L]])
  }
  nodes <- sum(held[[3L]][seq_len(index - 1L)]) + seq_len(held[[3L]][index])
  values <- readBin(held[[2L]], "double", length(held[[1L]]))
  list(vars = held[[1L]][nodes], values = values[nodes])
}
# the line bart prints for a forced warm start, the donor installed
lineOf <- function(sampler, donor, samples = NULL) {
  dbarts:::warmStartLine(
    dbarts:::installDonorTrees(sampler, donor, samples, TRUE)
  )
}
# A donor that is not clean: declined as it stands, the sampler untouched,
# and installed by force, which returns the line and the forced sampler. The
# line says why: it never ends at the trees taken.
expectChanged <- function(build, donor, info, samples = NULL) {
  sampler <- build()
  before <- stateOf(sampler)
  expect_identical(
    answerOf(sampler$installTrees(donor, samples)),
    visibly(FALSE),
    info = info
  )
  expect_identical(stateOf(sampler), before, info = info)
  line <- lineOf(sampler, donor, samples)
  expect_true(grepl("trees taken; ", line, fixed = TRUE), info = info)
  list(sampler = sampler, line = line, state = stateOf(sampler), was = before)
}
# an error in both forms that leaves the sampler as it was
expectRefused <- function(sampler, donor, pattern, info) {
  before <- stateOf(sampler)
  for (force in c(FALSE, TRUE)) {
    expect_error(
      sampler$installTrees(donor, forceUpdate = force),
      pattern,
      fixed = TRUE,
      info = info
    )
  }
  expect_identical(stateOf(sampler), before, info = info)
}

set.seed(31L, sample.kind = "Rejection")
n <- 300L
columns <- paste0("x", 1L:4L)
x <- matrix(runif(4L * n), n, 4L, dimnames = list(NULL, columns))
xOther <- matrix(runif(4L * n), n, 4L, dimnames = list(NULL, columns))
truth <- function(x) 3 * sin(3 * x[, 1L]) + 2 * x[, 2L]^2 + x[, 3L]
y <- truth(x) + rnorm(n, 0, 0.3)
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
# a sampler, run `sweeps` sweeps; built twice the same way, two are one
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
onOther <- function(...) make(..., data = xOther, response = yOther)
grid <- cutPointsOf(make())
fine <- grid[[1L]]
finer <- lapply(grid, function(g) sort(c(g, (g[-1L] + g[-length(g)]) / 2)))
# rows the hand trees below split: at or below the 30th point of x1, between
# it and the 50th, and above
low <- x[, 1L] <= fine[30L]
mid <- !low & x[, 1L] <= fine[50L]
byRows <- function(values, rows, weights = rep(1, n)) {
  sum(values * vapply(rows, function(r) sum(weights[r]), 0)) /
    sum(weights[Reduce(`|`, rows)])
}
twoSplits <- tree(
  c(1L, 1L, -1L, -1L, -1L),
  c(fine[50L], fine[30L], 0.10, -0.20, 0.05)
)
stumps <- function(first, ...) {
  c(list(first, ...), list(stump(), stump(), stump()))[1L:4L]
}

donor <- make(sweeps = 80L)
given <- stateOf(donor)
big <- list(n.trees = 50L)
bigDonor <- make(sweeps = 200L, control = big)

# ---- the two forms: clean is TRUE, the default declines, force is NULL ----

# a clean donor, onto a sampler holding kept draws and a stale state field
clean <- make(control = list(keepTrees = TRUE))
clean$storeState()
invisible(clean$run(0L, 3L))
expect_true(attr(stateOf(clean), "recordedDraws") > 0L)
expect_identical(answerOf(clean$installTrees(donor)), visibly(TRUE))
expect_identical(liveOf(stateOf(clean)), liveOf(given))
expect_identical(clean$state, stateOf(clean))
expect_identical(attr(clean$state, "recordedDraws"), 0L)
expect_identical(lineOf(make(), donor), "warm start: 4 of 4 donor trees taken")
expect_identical(
  answerOf(make()$installTrees(donor, forceUpdate = TRUE)),
  forcedValue
)

# 40 of the donor's rows on the donor's cut points: a leaf with no row
few <- function(control = list()) {
  sampler <- make(data = x[1L:40L, ], response = y[1L:40L], control = control)
  sampler$setCutPoints(grid)
  sampler
}
small <- few(big)
expect_identical(answerOf(small$installTrees(bigDonor)), visibly(FALSE))
expect_identical(
  answerOf(small$installTrees(bigDonor, forceUpdate = TRUE)),
  forcedValue
)
viaState <- few(big)
expect_null(viaState$setState(stateOf(bigDonor), forceUpdate = TRUE))
expect_identical(liveOf(stateOf(small)), liveOf(stateOf(viaState)))
expect_true(
  length(stateOf(small)[[1L]]$forests[[1L]]$tree.vars) <
    length(stateOf(bigDonor)[[1L]]$forests[[1L]]$tree.vars)
)

# a decline leaves the trees, the generators, the kept draws and the state
# field, here one stored before the last run, as a twin that never asked
kept <- c(big, keepTrees = TRUE)
stale <- few(kept)
twin <- few(kept)
for (sampler in list(stale, twin)) {
  sampler$storeState()
  invisible(sampler$run(0L, 3L))
}
field <- stale$state
expect_false(identical(field, stateOf(stale)))
expect_identical(stale$installTrees(bigDonor), FALSE)
expect_identical(stale$state, field)
expect_identical(stateOf(stale), stateOf(twin))
expect_identical(stale$predict(x[1L:40L, ]), twin$predict(x[1L:40L, ]))
expect_identical(stale$run(0L, 5L), twin$run(0L, 5L))

for (bad in list(1, NA, c(TRUE, TRUE), "partial")) {
  expect_error(
    make()$installTrees(donor, forceUpdate = bad),
    "'forceUpdate' must be TRUE or FALSE",
    fixed = TRUE
  )
}

# ---- cut points ------------------------------------------------------------

# a donor fitted to other rows: declined, and forced its splits move onto the
# sampler's cut points, which with their weights are as they were. Sigma and
# k go in where the sampler draws them, and its kept draws are dropped
drawsK <- function(...) make(..., leaf.prior = normal(k = chi(1.5, 2)))
kDonor <- drawsK(sweeps = 80L)
own <- stateOf(onOther())
away <- expectChanged(
  function() {
    drawsK(
      data = xOther,
      response = yOther,
      sweeps = 3L,
      control = list(keepTrees = TRUE)
    )
  },
  kDonor,
  "other rows"
)
expect_identical(attr(away$state, "cutPoints"), attr(own, "cutPoints"))
expect_identical(attr(away$state, "cutMass"), attr(own, "cutMass"))
expect_true(grepl("moved to this fit's cut points", away$line, fixed = TRUE))
expect_identical(
  c(attr(away$was, "recordedDraws"), attr(away$state, "recordedDraws")),
  c(3L, 0L)
)
expect_false(identical(
  away$was[[1L]]$forests[[1L]]$k,
  stateOf(kDonor)[[1L]]$forests[[1L]]$k
))
expect_identical(away$state[[1L]]$sigma, stateOf(kDonor)[[1L]]$sigma)
expect_identical(
  away$state[[1L]]$forests[[1L]]$k,
  stateOf(kDonor)[[1L]]$forests[[1L]]$k
)
expect_false(identical(liveOf(away$state), liveOf(stateOf(kDonor))))

# a grid that holds the donor's points and others: every split stands
onFiner <- make()
onFiner$setCutPoints(finer)
expect_identical(answerOf(onFiner$installTrees(donor)), visibly(TRUE))
expect_identical(liveOf(stateOf(onFiner)), liveOf(given))

# The forced move is the donor installed on its own points and then
# setCutPoints(splits = "value") back to the sampler's: the whole state and
# the next draws of a twin made that way, on 30 rows, where a leaf with no
# row merges before any split moves, with constant and with linear leaves.
# Cut points derived from 30 rows carry weights, which the twin's last call
# leaves behind and the forced start keeps, and the next draws read them: so
# the pair starts with its points set to the donor's and back, which leaves
# them unweighted, and a third sampler with its weights is held to the same
# trees and to its weights kept
x30 <- xOther[1L:30L, ]
y30 <- yOther[1L:30L]
linear2 <- function(...) make(..., leaf.prior = linear(columns = 2L))
for (case in list(
  list(build = make, donor = bigDonor, control = big, info = "constant"),
  list(
    build = linear2,
    donor = linear2(sweeps = 100L, control = list(n.trees = 8L)),
    control = list(n.trees = 8L),
    info = "linear"
  )
)) {
  on30 <- function(unweighted = TRUE) {
    sampler <- case$build(data = x30, response = y30, control = case$control)
    if (unweighted) {
      points <- cutPointsOf(sampler)
      sampler$setCutPoints(cutPointsOf(case$donor))
      sampler$setCutPoints(points)
    }
    sampler
  }
  weighted <- on30(FALSE)
  own <- stateOf(weighted)
  ownPoints <- attr(own, "cutPoints")
  expect_false(is.null(attr(own, "cutMass")), info = case$info)
  forced <- on30()
  expect_null(attr(stateOf(forced), "cutMass"), info = case$info)
  hand <- on30()
  for (sampler in list(weighted, forced)) {
    expect_null(sampler$installTrees(case$donor, forceUpdate = TRUE))
  }
  hand$setCutPoints(cutPointsOf(case$donor))
  expect_null(hand$installTrees(case$donor, forceUpdate = TRUE))
  hand$setCutPoints(ownPoints, splits = "value")
  expect_identical(stateOf(forced), stateOf(hand), info = case$info)
  expect_identical(cutPointsOf(forced), ownPoints, info = case$info)
  moved <- stateOf(weighted)
  expect_identical(leavesOf(moved), leavesOf(stateOf(forced)), info = case$info)
  for (record in c("cutPoints", "cutMass")) {
    expect_identical(attr(moved, record), attr(own, record), info = case$info)
  }
  expect_identical(forced$run(0L, 5L), hand$run(0L, 5L), info = case$info)
}
# and the trees of a forced setState
viaState <- make(data = x30, response = y30, control = big)
expect_null(viaState$setState(stateOf(bigDonor), forceUpdate = TRUE))
forced <- make(data = x30, response = y30, control = big)
forced$installTrees(bigDonor, forceUpdate = TRUE)
expect_identical(liveOf(stateOf(forced)), liveOf(stateOf(viaState)))
expect_true(
  length(stateOf(forced)[[1L]]$forests[[1L]]$tree.vars) <
    length(stateOf(bigDonor)[[1L]]$forests[[1L]]$tree.vars)
)

# A split left no point merges by rows. The hand tree splits x1 at its 50th
# point and, beneath, at its 30th; on the 50th and 75th points alone the
# lower split has no point left and becomes one leaf at the mean of its two
# leaves weighted by the rows each held, and by their case weights
caseWeights <- runif(n, 0.2, 3)
handDonor <- withTrees(given, stumps(twoSplits))
for (weights in list(NULL, caseWeights)) {
  onTwo <- function() {
    sampler <- make(weights = weights)
    sampler$setCutPoints(fine[c(50L, 75L)], 1L)
    sampler
  }
  moved <- expectChanged(onTwo, handDonor, "two points")
  merged <- byRows(
    c(0.10, -0.20),
    list(low, mid),
    if (is.null(weights)) rep(1, n) else weights
  )
  got <- treeOf(moved$state)
  expect_identical(got$vars, c(1L, -1L, -1L))
  expect_identical(got$values[c(1L, 3L)], c(fine[50L], 0.05))
  expect_true(abs(got$values[2L] - merged) < 1e-12)
  expect_true(abs(merged + 0.05) > 0.01 && abs(merged - 0.10) > 0.01)
  expect_identical(
    moved$line,
    "warm start: 4 of 4 donor trees taken; 1 cut back to fit"
  )
}
# the weights tell: the last value is not the one unit weights give
expect_true(abs(byRows(c(0.10, -0.20), list(low, mid)) - merged) > 1e-3)

# ---- splits a forest does not allow ----------------------------------------

# x1 and, beneath it on the right, x2: under a limit of one variable a tree,
# and in a tree confined to x1, the lower split is cut back to one leaf at
# the mean by rows and the split above it stays
x2Point <- grid[[2L]][40L]
high <- !low & !mid
twoColumns <- tree(
  c(1L, -1L, 2L, -1L, -1L),
  c(fine[50L], 0.10, x2Point, -0.20, 0.05)
)
cutBack <- byRows(
  c(-0.20, 0.05),
  list(high & x[, 2L] <= x2Point, high & x[, 2L] > x2Point)
)
additive <- function(...) make(..., interactions = interactions(max.order = 1))
blocked <- function(...) {
  make(
    ...,
    blocks = blocks(
      groups = list("x1", c("x2", "x3", "x4")),
      trees.per.group = c(2L, 2L)
    )
  )
}
for (build in list(additive, blocked)) {
  back <- expectChanged(
    build,
    withTrees(given, stumps(twoColumns)),
    "forbidden"
  )
  got <- treeOf(back$state)
  expect_identical(got$vars, c(1L, -1L, -1L))
  expect_identical(got$values[1L:2L], c(fine[50L], 0.10))
  expect_true(abs(got$values[3L] - cutBack) < 1e-12)
  expect_true(abs(cutBack) > 0.01)
  expect_identical(
    back$line,
    "warm start: 4 of 4 donor trees taken; 1 cut back to fit"
  )
}
# a fitted donor, on its own rows and across grids: no tree is left with two
# variables, and on its own rows the fit's mean over the rows is unchanged
oneVariable <- function(state) {
  forest <- state[[1L]]$forests[[1L]]
  index <- rep(seq_along(forest$tree.sizes), forest$tree.sizes)
  splits <- forest$tree.vars > 0L
  all(
    tapply(forest$tree.vars[splits], index[splits], function(v) {
      length(unique(v))
    }) ==
      1L
  )
}
fitOf <- function(sampler, rows = x) as.vector(sampler$predict(rows))
expect_false(oneVariable(stateOf(bigDonor)))
sameRows <- expectChanged(function() additive(control = big), bigDonor, "limit")
expect_true(oneVariable(sameRows$state))
expect_true(
  abs(mean(fitOf(sameRows$sampler)) - mean(fitOf(bigDonor))) < 1e-10
)
across <- expectChanged(
  function() additive(data = xOther, response = yOther, control = big),
  bigDonor,
  "limit, across grids"
)
expect_true(oneVariable(across$state))
expect_true(grepl("cut back to fit", across$line, fixed = TRUE))

# ---- the variance forest ----------------------------------------------------

# a variance threshold off the grid alone declines the donor; forced, that
# split moves onto the grid and the mean trees stand
heteroscedastic <- function(...) make(..., variance = TRUE)
ownState <- stateOf(heteroscedastic(sweeps = 20L))
between <- mean(grid[[2L]][40L:41L])
offGrid <- withTrees(
  ownState,
  c(
    list(tree(c(2L, -1L, -1L), c(between, 1.2, 0.8))),
    rep(list(stump(1)), length(ownState[[1L]]$variance.sizes) - 1L)
  ),
  variance = TRUE
)
attr(offGrid, "cutPoints")[[2L]] <- sort(c(grid[[2L]], between))
moved <- expectChanged(heteroscedastic, offGrid, "variance threshold")
expect_identical(liveOf(moved$state), liveOf(ownState))
got <- treeOf(moved$state, variance = TRUE)
expect_identical(got$vars, c(2L, -1L, -1L))
expect_true(got$values[1L] %in% grid[[2L]])
expect_true(
  grepl("; variance trees: [0-9]+ of [0-9]+ taken, 1 moved to", moved$line)
)
# a split a restricted variance forest may not use is cut back to one leaf at
# the geometric mean by rows
restricted <- function(...) {
  make(..., variance = varianceForest(vars = 2L, n.trees = 3L))
}
ownState <- stateOf(restricted())
onX1 <- withTrees(
  ownState,
  list(tree(c(1L, -1L, -1L), c(fine[50L], 1.2, 0.8)), stump(1), stump(1)),
  variance = TRUE
)
back <- expectChanged(restricted, onX1, "variance restriction")
got <- treeOf(back$state, variance = TRUE)
expect_identical(got$vars, -1L)
expect_true(
  abs(log(got$values) - byRows(log(c(1.2, 0.8)), list(low | mid, high))) < 1e-12
)
expect_identical(
  back$line,
  paste0(
    "warm start: 4 of 4 donor trees taken; ",
    "variance trees: 3 of 3 taken, 1 cut back to fit"
  )
)

# ---- a monotone forest -------------------------------------------------------

# x1 increasing: a tree whose leaves fall along x1 is out of order, and one
# that also loses a split is still one tree reset. The trees in order are the
# donor's
rising <- function(...) make(..., monotone = c(x1 = "increasing"))
falling <- tree(c(1L, -1L, -1L), c(fine[50L], 0.5, -0.5))
inOrder <- tree(c(1L, -1L, -1L), c(fine[40L], -0.3, 0.3))
# the lower split sits below the one over it: its left leaf holds no row,
# and the leaf it merges into is below the one on its left
stacked <- tree(
  c(1L, -1L, 1L, -1L, -1L),
  c(fine[50L], 0.5, fine[30L], 0.2, -0.4)
)
outOfOrder <- withTrees(given, list(falling, inOrder, stacked, stump(0.1)))
reset <- expectChanged(rising, outOfOrder, "monotone")
expect_identical(treeOf(reset$state, 1L)$values, c(fine[50L], 0, 0))
expect_identical(treeOf(reset$state, 2L), inOrder[c("vars", "values")])
expect_identical(treeOf(reset$state, 3L)$values, c(fine[50L], 0, 0))
expect_identical(treeOf(reset$state, 4L)$values, 0.1)
expect_identical(
  reset$line,
  "warm start: 4 of 4 donor trees taken; 2 reset to 0"
)
expect_identical(
  answerOf(rising()$installTrees(withTrees(given, stumps(inOrder)))),
  visibly(TRUE)
)

# ---- sides for missing values ------------------------------------------------

# A donor that could hold no missing value in x2, onto rows that hold some:
# each rule on x2 draws its side before the trees are judged, the start is
# clean and the generator stands past its coin. A start then declined, here
# for a leaf with no row, puts the generator back
xHoled <- x
xHoled[seq(7L, n, by = 9L), 2L] <- NA
onX2 <- tree(c(2L, -1L, -1L), c(x2Point, 0.1, -0.1))
holed <- make(data = xHoled)
before <- stateOf(holed)
expect_identical(cutPointsOf(holed), grid)
expect_identical(
  answerOf(holed$installTrees(withTrees(given, stumps(onX2)))),
  visibly(TRUE)
)
expect_false(identical(stateOf(holed)[[1L]]$rng.state, before[[1L]]$rng.state))
empty <- tree(
  c(1L, 1L, -1L, -1L, -1L),
  c(fine[30L], fine[50L], 0.10, -0.20, 0.05)
)
holed <- make(data = xHoled)
expect_identical(
  holed$installTrees(withTrees(given, stumps(onX2, empty))),
  FALSE
)
expect_identical(stateOf(holed), before)
expect_identical(holed$run(0L, 3L), make(data = xHoled)$run(0L, 3L))
# a side on a column that has never held a missing value here is dropped
sided <- tree(c(2L, -1L, -1L), c(x2Point, 0.1, -0.1), c(TRUE, FALSE, FALSE))
dropped <- expectChanged(make, withTrees(given, stumps(sided)), "side")
expect_identical(
  dropped$state[[1L]]$forests[[1L]]$tree.flags[1L:3L],
  as.raw(c(2L, 0L, 0L))
)
expect_identical(
  dropped$line,
  paste0(
    "warm start: 4 of 4 donor trees taken; ",
    "1 with a side for missing values dropped"
  )
)

# ---- the draw `samples` names ------------------------------------------------

# The verdict is on the trees each chain takes. A donor that keeps draws
# gives those and never its live trees: fitted under a limit of one variable
# a tree its draws are clean for such a sampler, and fitted free they are not
keeping <- list(n.trees = 10L, keepTrees = TRUE)
limited <- stateOf(additive(sweeps = 60L, control = keeping))
free <- stateOf(make(sweeps = 60L, control = keeping))
savedFields <- sub("^tree", "saved", treeFields)
# the live trees of one state over the kept draws of another
mixed <- function(live, kept) {
  live[[1L]]$forests[[1L]][savedFields] <-
    kept[[1L]]$forests[[1L]][savedFields]
  live
}
recipient <- function(chains = 1L) {
  additive(control = list(n.trees = 10L, n.chains = chains))
}
for (draw in 1L:3L) {
  expect_identical(recipient()$installTrees(mixed(limited, free), draw), FALSE)
  expect_identical(recipient()$installTrees(mixed(free, limited), draw), TRUE)
}
# two chains: the second chain's draw alone declines the donor, and nothing
# of the first chain is touched
both <- limited
both[[2L]] <- free[[1L]]
pair <- recipient(2L)
before <- stateOf(pair)
expect_identical(pair$installTrees(both, samples = c(1L, 4L)), FALSE)
expect_identical(stateOf(pair), before)
expect_identical(pair$installTrees(both, samples = c(1L, 2L)), TRUE)

# ---- leaf values that are not the donor's ------------------------------------

# Gaussian-process fits are values by row: copied where the donor's recorded
# points are the sampler's, and started at zero, the donor declined, where
# they are not, though every split is a point of the grid
rows <- seq_len(60L)
processes <- function(...) {
  make(
    ...,
    data = x[rows, ],
    response = y[rows],
    leaf.prior = gp(columns = 2L)
  )
}
gpDonor <- processes(sweeps = 20L)
gpGiven <- stateOf(gpDonor)
copied <- processes()
expect_identical(answerOf(copied$installTrees(gpDonor)), visibly(TRUE))
fitsOf <- function(state) {
  state[[1L]]$forests[[1L]][c("tree.vars", "tree.params")]
}
expect_identical(fitsOf(stateOf(copied)), fitsOf(gpGiven))
cold <- expectChanged(
  function() {
    sampler <- processes()
    points <- cutPointsOf(sampler)
    sampler$setCutPoints(lapply(points, function(g) {
      sort(c(g, (g[-1L] + g[-length(g)]) / 2))
    }))
    sampler
  },
  gpDonor,
  "gp, other points"
)
expect_true(any(gpGiven[[1L]]$forests[[1L]]$tree.params != 0))
expect_true(all(cold$state[[1L]]$forests[[1L]]$tree.params == 0))
expect_identical(
  cold$state[[1L]]$forests[[1L]]$tree.vars,
  gpGiven[[1L]]$forests[[1L]]$tree.vars
)
expect_identical(
  cold$line,
  "warm start: 4 of 4 donor trees taken; gp fits started at 0"
)

# leaf blocks a sampler of constant leaves does not read: each leaf goes in
# at the value the donor stores for it
linearDonor <- linear2(sweeps = 40L)
for (case in list(
  list(donor = linearDonor, build = make, said = "slopes left behind"),
  list(
    donor = gpDonor,
    build = function() make(data = x[rows, ], response = y[rows]),
    said = "gp fits left behind"
  )
)) {
  unread <- expectChanged(case$build, case$donor, case$said)
  expect_identical(liveOf(unread$state), liveOf(stateOf(case$donor)))
  expect_identical(
    unread$line,
    paste0("warm start: 4 of 4 donor trees taken; ", case$said)
  )
}

# A linear leaf's covariate that holds one value. On one side only the
# coefficients are restated and the start is the donor's function at the
# sampler's rows; at the same value on both sides nothing is restated; at
# another value on each side the slopes arrive as built, which is not the
# donor's function there, and the donor is declined. The donor's trees are
# hand-built so that none splits on the covariate
withX2 <- function(value) {
  frame <- x
  frame[, 2L] <- value
  frame
}
spread <- runif(n, 0, 10)
slopedOn <- function(value) {
  sampler <- linear2(data = withX2(value))
  state <- withTrees(
    stateOf(sampler),
    stumps(tree(c(1L, -1L, -1L), c(fine[50L], 0.1, -0.2)))
  )
  state[[1L]]$forests[[1L]]$tree.params <- c(0.3, -0.4, 0, 0, 0)
  expect_true(sampler$setState(state))
  sampler
}
for (case in list(
  list(donor = 0, here = spread, info = "one value in the donor"),
  list(donor = spread, here = 3, info = "one value here"),
  list(donor = 0, here = 0, info = "one value, the same")
)) {
  sloped <- slopedOn(case$donor)
  recipient <- linear2(data = withX2(case$here))
  expect_identical(recipient$installTrees(sloped), TRUE, info = case$info)
  expect_true(
    max(abs(
      fitOf(recipient, withX2(case$here)) - fitOf(sloped, withX2(case$here))
    )) <
      1e-12,
    info = case$info
  )
}
sloped <- slopedOn(0)
keptSlopes <- expectChanged(
  function() linear2(data = withX2(3)),
  sloped,
  "one value on each side"
)
expect_identical(leavesOf(keptSlopes$state), leavesOf(stateOf(sloped)))
expect_identical(
  keptSlopes$line,
  paste0(
    "warm start: 4 of 4 donor trees taken; ",
    "slopes kept across a covariate of one value"
  )
)

# ---- an ordered factor ---------------------------------------------------------

# Its level count fixes its cut points, so a donor recorded on others there,
# a numeric column or a factor of another level count, is another column: an
# error in both forms. The same levels are clean; and a donor's factor onto a
# numeric column is moved like any other, its splits above every value here
# merged
framed <- function(x4, ...) {
  frame <- data.frame(y = y, x[, 1L:3L], x4 = x4)
  dbarts::dbarts(y ~ ., frame, control = controlWith(...))
}
levels4 <- ordered(cut(x[, 4L], 4L, labels = FALSE))
levels6 <- ordered(cut(x[, 4L], 6L, labels = FALSE))
structure <- "does not support the warm-start donor's predictor structure"
expectRefused(framed(levels4), framed(x[, 4L]), structure, "numeric")
expectRefused(framed(levels4), framed(levels6), structure, "six levels")
factorDonor <- framed(levels4)
points4 <- cutPointsOf(factorDonor)[[4L]]
onFactor <- withTrees(
  stateOf(factorDonor),
  stumps(tree(c(4L, -1L, -1L), c(points4[length(points4)], 0.1, -0.1)))
)
expect_identical(framed(levels4)$installTrees(onFactor), TRUE)
expect_true(points4[length(points4)] > max(x[, 4L]))
merged <- expectChanged(function() framed(x[, 4L]), onFactor, "factor donor")
expect_identical(treeOf(merged$state)$vars, -1L)

# ---- errors stay errors in both forms ------------------------------------------

shape <- "not shape-compatible"
blocksOf <- "malformed leaf parameters"
refusals <- list(
  list(make(), make(data = x[, 1L:3L]), "malformed cut points", "columns"),
  list(
    make(),
    framed(factor(as.integer(levels4))),
    structure,
    "an unordered factor"
  ),
  list(make(), make(control = list(n.trees = 6L)), shape, "tree count"),
  list(make(), make(tree.prior = dart()), "DART", "DART in the donor"),
  list(make(tree.prior = dart()), make(), "DART", "DART here"),
  list(make(), heteroscedastic(), shape, "variance in the donor"),
  list(heteroscedastic(), make(), shape, "variance here"),
  list(linear2(), make(), blocksOf, "constant into linear"),
  list(
    processes(),
    make(data = x[rows, ], response = y[rows]),
    blocksOf,
    "constant into gp"
  ),
  list(
    processes(),
    linear2(data = x[rows, ], response = y[rows]),
    blocksOf,
    "linear into gp"
  ),
  list(
    linear2(data = x[rows, ], response = y[rows]),
    processes(),
    blocksOf,
    "gp into linear"
  ),
  list(
    make(leaf.prior = linear(columns = 2L:3L)),
    linear2(),
    blocksOf,
    "covariates"
  ),
  list(
    make(
      data = x[1L:40L, ],
      response = y[1L:40L],
      leaf.prior = gp(columns = 2L)
    ),
    gpDonor,
    blocksOf,
    "gp on other rows"
  ),
  list(make(), list(1, 2, 3), "must be a dbarts sampler", "not a donor")
)
for (case in refusals) {
  expectRefused(case[[1L]], case[[2L]], case[[3L]], case[[4L]])
}

# ---- what bart says --------------------------------------------------------------

# One line, after the setup lines and before the loop, at every warm start
# under verbose and at no other time; no warning. The donor goes in by
# force: under a limit of one variable a tree it is cut back
fit <- function(...) {
  dbarts::bart(
    x,
    y,
    n.trees = 50L,
    n.samples = 2L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 3L,
    ...
  )
}
said <- function(...) {
  warnings <- 0L
  output <- withCallingHandlers(
    suppressMessages(utils::capture.output(invisible(fit(...)))),
    warning = function(w) warnings <<- warnings + 1L
  )
  expect_identical(warnings, 0L)
  output
}
output <- said(
  warm.start = bigDonor,
  verbose = TRUE,
  interactions = interactions(max.order = 1)
)
line <- grep("^warm start:", output)
expect_identical(length(line), 1L)
expect_true(line < grep("^Running mcmc loop:", output))
expect_true(
  grepl(
    "^warm start: 50 of 50 donor trees taken; [0-9]+ cut back to fit$",
    output[line]
  )
)
expect_identical(output[line], sameRows$line)
expect_identical(length(grep("warm start", said(verbose = TRUE))), 0L)
expect_identical(
  length(said(warm.start = bigDonor, verbose = FALSE)),
  0L
)

# ---- what the line counts --------------------------------------------------------

# Each changed tree once, under the first name that holds. A donor recorded
# on other points of x2: a tree as stored; one moved, which also loses a
# side; one cut back, whose splits are off the grid too; one that only loses
# its side
shifted <- grid[[2L]][-length(grid[[2L]])] + diff(grid[[2L]]) / 2
asStored <- tree(c(1L, -1L, -1L), c(fine[40L], 0.1, -0.1))
movedOne <- function(right = FALSE) {
  tree(c(2L, -1L, -1L), c(shifted[30L], 0.2, -0.2), c(right, FALSE, FALSE))
}
cutOne <- tree(
  c(2L, 2L, -1L, -1L, -1L),
  c(shifted[20L], shifted[60L], 0.3, 0.5, -0.3)
)
elsewhere <- function(trees, state = given) {
  state <- withTrees(state, trees)
  attr(state, "cutPoints")[[2L]] <- shifted
  state
}
sidedOnX1 <- tree(
  c(1L, -1L, -1L),
  c(fine[50L], 0.05, -0.05),
  c(TRUE, FALSE, FALSE)
)
counted <- expectChanged(
  make,
  elsewhere(list(asStored, movedOne(TRUE), cutOne, sidedOnX1)),
  "counts"
)
expect_identical(
  counted$line,
  paste0(
    "warm start: 4 of 4 donor trees taken; ",
    "1 moved to this fit's cut points, 1 cut back to fit, ",
    "1 with a side for missing values dropped"
  )
)
expect_identical(treeOf(counted$state)[c("vars", "values")], asStored[1L:2L])
# two chains, each given the other donor chain's trees: trees are counted a
# chain and changes summed over the chains
two <- function() make(control = list(n.chains = 2L))
twoChains <- elsewhere(
  list(movedOne(), asStored, asStored, asStored),
  stateOf(two())
)
twoChains[[2L]] <- elsewhere(list(movedOne(), movedOne(), cutOne, asStored))[[
  1L
]]
summed <- expectChanged(two, twoChains, "two chains", samples = c(2L, 1L))
expect_identical(
  summed$line,
  paste0(
    "warm start: 4 of 4 donor trees taken; ",
    "3 moved to this fit's cut points, 1 cut back to fit, over 2 chains"
  )
)
expect_identical(liveOf(summed$state)[[2L]]$tree.sizes, rep(3L, 4L))
expect_identical(liveOf(summed$state)[[1L]]$tree.sizes, rep(3L, 4L))

# ---- sparse and mixed stores -------------------------------------------------------

# a forced start across grids re-quantizes a store that keeps no dense matrix
# and puts its cut points back
mixedFrame <- function(rows) {
  data.frame(rows, f = factor(rep_len(c("p", "q", "r"), n)))
}
stores <- list(
  mixed = function(rows = x, response = y) {
    dbarts::dbarts(
      mixedFrame(rows),
      response,
      control = controlWith(n.trees = 50L)
    )
  }
)
if (requireNamespace("Matrix", quietly = TRUE)) {
  sparsely <- function(rows) {
    rows[, 2L] <- ifelse(seq_len(n) %% 5L == 0L, round(rows[, 2L], 1) + 0.1, 0)
    Matrix::Matrix(rows, sparse = TRUE)
  }
  stores$sparse <- function(rows = x, response = y) {
    dbarts::dbarts(
      sparsely(rows),
      response,
      control = controlWith(n.trees = 50L)
    )
  }
}
for (store in names(stores)) {
  storeDonor <- stores[[store]]()
  invisible(storeDonor$run(100L, 1L))
  ownPoints <- cutPointsOf(stores[[store]](xOther, yOther))
  moved <- expectChanged(
    function() stores[[store]](xOther, yOther),
    storeDonor,
    store
  )
  expect_identical(attr(moved$state, "cutPoints"), ownPoints, info = store)
  expect_true(all(is.finite(moved$sampler$run(0L, 3L)$train)), info = store)
}
