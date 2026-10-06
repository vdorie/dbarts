# One declared forest honours forest(vars = ): the forest splits on the named
# columns and no others for the sampler's whole life, and the fit is, seed for
# seed and given the same residual scale estimate, the fit on those columns
# alone. The engine's own gates live in tests/cpp/test_state.cpp
# (testSingleForestColumnRestriction).

forest <- dbarts::dbartsForests$forest
blocks <- dbarts::dbartsForests$blocks
interactions <- dbarts::dbartsForests$interactions
cgm <- dbarts::dbartsPriors$cgm
dart <- dbarts::dbartsPriors$dart
linear <- dbarts::dbartsPriors$linear
gp <- dbarts::dbartsPriors$gp
student <- dbarts::dbartsFamilies$student

set.seed(31L)
n <- 150L
x <- cbind(a = rnorm(n), b = rnorm(n), c = rnorm(n))
# every response leans hardest on c, the column a restricted forest may not
# use, so a forest that ignores its restriction splits on it at once
signal <- x[, "a"] + 2 * x[, "c"]
y <- signal + rnorm(n, 0, 0.5)
y.binary <- as.double(signal + rnorm(n) > 0)
allowed <- c("a", "b")

makeControl <- function(...) {
  settings <- list(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    n.samples = 1L,
    updateState = FALSE,
    verbose = FALSE
  )
  do.call(dbartsControl, modifyList(settings, list(...)))
}
control <- makeControl()
# a chain's generator is seeded when its sampler is created, so two samplers
# created under one seed draw the same stream. sigest is fixed because the
# residual scale estimate is the data object's, from a linear fit on every
# column, and so differs between a design and its sub-matrix
sampler <- function(design, response, ..., settings = control) {
  set.seed(7L)
  dbarts(design, response, ..., control = settings, sigest = 1)
}
restrictTo <- function(vars, ...) list(forest(vars = vars, ...))
# total splits per column over a run
splits <- function(draws) rowSums(draws$varcount)
# the split variables of the live trees, one set per tree
treeVars <- function(trees) {
  internal <- trees$var > 0L
  lapply(
    split(trees$var[internal], trees$tree[internal]),
    function(v) sort(unique(v))
  )
}

# ---- the defect: a single declared forest ignored its 'vars' ----------------

# read three ways after every sweep, beside the same fit without 'vars'
restricted <- sampler(x, y, forests = restrictTo(allowed))
plain <- sampler(x, y)
expect_identical(attr(restricted$model, "forest.columns"), 1:2)
expect_null(attr(plain$model, "forest.columns"))
outside <- c(trees = 0, counts = 0, run = 0)
inside <- outside
plainOutside <- outside
for (sweep in seq_len(300L)) {
  draw <- restricted$run(0L, 1L)
  vars <- restricted$getTrees()$var
  counts <- restricted$getForestVariableCounts()
  outside <- outside +
    c(sum(vars == 3L), counts[3L, 1L], draw$varcount[3L, 1L])
  inside <- inside +
    c(sum(vars %in% 1:2), sum(counts[1:2, 1L]), sum(draw$varcount[1:2, 1L]))
  draw <- plain$run(0L, 1L)
  plainOutside <- plainOutside +
    c(
      sum(plain$getTrees()$var == 3L),
      plain$getForestVariableCounts()[3L, 1L],
      draw$varcount[3L, 1L]
    )
}
expect_identical(unname(outside), c(0, 0, 0))
expect_true(all(inside > 0))
expect_true(all(plainOutside > 0))

# ---- the rule: the fit on the named columns alone, seed for seed, at an ------
# ---- equal residual scale estimate -------------------------------------------

# left to be estimated, the scale comes from a linear fit on every column of
# each design, and the two fits differ
set.seed(7L)
estimatedOnList <- dbarts(
  x,
  y,
  forests = restrictTo(allowed),
  control = control
)
set.seed(7L)
estimatedOnColumns <- dbarts(x[, allowed], y, control = control)
expect_false(estimatedOnList$data@sigma == estimatedOnColumns$data@sigma)
expect_false(identical(
  estimatedOnList$run(0L, 20L)$train,
  estimatedOnColumns$run(0L, 20L)$train
))


expectSubMatrixFit <- function(response, ..., alone = list(...)) {
  onList <- sampler(x, response, forests = restrictTo(allowed), ...)$run(
    0L,
    60L
  )
  onColumns <- do.call(sampler, c(list(x[, allowed], response), alone))$run(
    0L,
    60L
  )
  expect_identical(onList$train, onColumns$train)
  expect_identical(onList$sigma, onColumns$sigma)
  expect_identical(unname(onList$varcount[1:2, ]), unname(onColumns$varcount))
  expect_true(all(onList$varcount[3L, ] == 0L))
  expect_true(sum(onList$varcount) > 0L)
  invisible(list(onList = onList, onColumns = onColumns))
}
expectSubMatrixFit(y)
expectSubMatrixFit(y.binary)
expectSubMatrixFit(y, interactions = interactions(max.order = 1L))
# the caller's ratios hold among the allowed columns
expectSubMatrixFit(
  y,
  tree.prior = cgm(split.probs = c(0.1, 0.3, 0.6)),
  alone = list(tree.prior = cgm(split.probs = c(0.25, 0.75)))
)
# DART's Dirichlet is laid over the allowed columns: the same probabilities
# there, and none anywhere else. No delay, so every sweep draws them
dartFits <- expectSubMatrixFit(y, tree.prior = dart(update.delay = 0))
expect_identical(
  unname(dartFits$onList$varprobs[1:2, ]),
  unname(dartFits$onColumns$varprobs)
)
expect_true(all(dartFits$onList$varprobs[3L, ] == 0))
expect_true(all(dartFits$onList$varprobs[1:2, ] != 0.5))

# allowed columns that run out along a path: two 0/1 columns held to one cut
# each, which the default grid of 100 cuts would never exhaust. Below a split
# on each, no allowed column is left. The column list still confines the
# forest there; zero split probabilities on the same fixture do not, the first
# column available being proposed whatever its probability
x.binary <- cbind(a = rbinom(n, 1L, 0.5), b = rbinom(n, 1L, 0.5), c = x[, "c"])
y.runout <- x.binary[, "a"] + 2 * x.binary[, "c"] + rnorm(n, 0, 0.5)
runOut <- function(...) {
  fit <- sampler(
    x.binary,
    y.runout,
    ...,
    settings = makeControl(n.cuts = c(1L, 1L, 100L))
  )
  expect_equal(fit$data@n.cuts, c(1, 1, 100))
  splits(fit$run(0L, 2000L))
}
listed <- runOut(forests = restrictTo(allowed))
expect_identical(unname(listed[3L]), 0)
expect_true(all(listed[1:2] > 0))
zeroed <- runOut(tree.prior = cgm(split.probs = c(0.5, 0.5, 0)))
expect_true(zeroed[3L] > 0)
expect_true(runOut()[3L] > 0)

# ---- how 'vars' is stated ---------------------------------------------------

# by index, and through dbartsSpec and the direct door
byIndex <- sampler(x, y, forests = restrictTo(1:2))
expect_identical(byIndex$model, restricted$model)
spec <- dbartsSpec(
  dbartsData(x, y),
  control,
  forests = restrictTo(allowed),
  sigest = 1
)
expect_identical(attr(spec$model, "forest.columns"), 1:2)
door <- new("dbartsSampler", spec$control, spec$model, spec$data)
expect_identical(unname(splits(door$run(0L, 100L))[3L]), 0)
# the attribute alone restricts a model built without it, and the bridge
# checks what it is handed
plainSpec <- dbartsSpec(dbartsData(x, y), control, sigest = 1)
byHand <- plainSpec$model
attr(byHand, "forest.columns") <- 1:2
door <- new("dbartsSampler", plainSpec$control, byHand, plainSpec$data)
expect_identical(unname(splits(door$run(0L, 100L))[3L]), 0)
attr(byHand, "forest.columns") <- 4L
expect_error(
  new("dbartsSampler", plainSpec$control, byHand, plainSpec$data),
  pattern = "forest column index out of range"
)
attr(byHand, "forest.columns") <- c(1, 2)
expect_error(
  new("dbartsSampler", plainSpec$control, byHand, plainSpec$data),
  pattern = "forest columns must be resolved integer indices"
)

# a factor is one column, named as the formula names it
frame <- data.frame(x, g = factor(sample(c("u", "v", "w"), n, TRUE)))
frame$y <- y + 2 * (frame$g == "v")
set.seed(7L)
withFactor <- dbarts(
  y ~ a + g + b + c,
  frame,
  forests = restrictTo(c("a", "g")),
  control = control
)
expect_identical(attr(withFactor$model, "forest.columns"), 1:2)
factorSplits <- splits(withFactor$run(0L, 200L))
expect_identical(unname(factorSplits[3:4]), c(0, 0))
expect_true(all(factorSplits[1:2] > 0))

# naming every column restricts nothing and stores nothing
everyColumn <- sampler(x, y, forests = restrictTo(c("a", "b", "c")))
expect_null(attr(everyColumn$model, "forest.columns"))
expect_identical(
  everyColumn$run(0L, 30L)$train,
  sampler(x, y)$run(0L, 30L)$train
)

# refused as on any forest, and nothing is dropped unread
expect_error(
  sampler(x, y, forests = restrictTo("z")),
  pattern = "'vars' name not found in the design's column names"
)
expect_error(
  sampler(x, y, forests = restrictTo(character(0L))),
  pattern = "'vars' is empty"
)
expect_error(
  sampler(x, y, forests = restrictTo(4L)),
  pattern = "'vars' column index out of range"
)
# a missing value, of any type, on a single forest and on a forest of several
for (missing in list(NA, NA_integer_, c(1L, NA), c("a", NA))) {
  expect_error(
    sampler(x, y, forests = restrictTo(missing)),
    pattern = "'vars' contains missing values"
  )
}
expect_error(
  sampler(
    x,
    y,
    forests = list(forest(), forest(basis = rep_len(c(0, 1), n), vars = NA))
  ),
  pattern = "'vars' contains missing values"
)
expect_error(
  sampler(
    x,
    y,
    forests = restrictTo(allowed),
    tree.prior = cgm(split.probs = c(0, 0, 1))
  ),
  pattern = "'split.probs' gives no positive probability to any column"
)

# ---- every family and leaf --------------------------------------------------

kinds <- list(
  probit = list(y.binary, family = "probit"),
  logistic = list(y.binary, family = "logistic"),
  ordinal = list(
    as.integer(cut(y, quantile(y, 0:3 / 3), include.lowest = TRUE)),
    family = "ordinal"
  ),
  nbinom = list(rpois(n, exp(0.5 * signal)), family = "nbinom"),
  aft = list(cbind(exp(y), rbinom(n, 1L, 0.8)), family = "aft"),
  student = list(y, family = student(df = 5)),
  linear = list(y, leaf.prior = linear("a")),
  gp = list(y, leaf.prior = gp("a")),
  monotone = list(y, monotone = c(a = 1))
)
for (kind in names(kinds)) {
  arguments <- c(list(x), kinds[[kind]])
  kindSplits <- splits(
    do.call(sampler, c(arguments, list(forests = restrictTo(allowed))))$run(
      0L,
      150L
    )
  )
  expect_identical(unname(kindSplits[3L]), 0, info = kind)
  expect_true(sum(kindSplits[1:2]) > 0, info = kind)
  expect_true(splits(do.call(sampler, arguments)$run(0L, 150L))[3L] > 0)
}

# every category forest of a multinomial fit
y.category <- cut(
  x[, "c"] + rnorm(n, 0, 0.3),
  c(-Inf, -0.5, 0.5, Inf),
  labels = c("low", "mid", "high")
)
categoryCounts <- function(...) {
  set.seed(7L)
  fit <- dbarts(x, y.category, family = "multinomial", ..., control = control)
  apply(fit$run(0L, 150L)$varcount, 1:2, sum)
}
restrictedCategories <- categoryCounts(forests = restrictTo(allowed))
expect_identical(dim(restrictedCategories), c(3L, 3L))
expect_true(all(restrictedCategories[3L, ] == 0))
expect_true(all(colSums(restrictedCategories[1:2, ]) > 0))
expect_true(all(categoryCounts()[3L, ] > 0))

# a hazard fit's own period column stays allowed: the caller did not supply
# it, and 'vars' restricts the caller's columns
time <- pmin(1L + rpois(n, exp(0.5 * x[, "c"])), 5L)
status <- rbinom(n, 1L, 0.8)
set.seed(7L)
hazard <- dbarts(
  x,
  cbind(time, status),
  family = "hazard",
  forests = restrictTo(allowed),
  control = control
)
expect_identical(colnames(hazard$data@x), c("a", "b", "c", "period"))
expect_identical(attr(hazard$model, "forest.columns"), c(1L, 2L, 4L))
hazardSplits <- splits(hazard$run(0L, 200L))
expect_identical(unname(hazardSplits[3L]), 0)
expect_true(hazardSplits[4L] > 0)
# with several forests each forest's 'vars' is taken as written: period is
# split on only where it is named
hazardFirstForest <- function(vars) {
  set.seed(7L)
  fit <- dbarts(
    x,
    cbind(time, status),
    family = "hazard",
    forests = list(
      forest(vars = vars),
      forest(basis = rep_len(c(0, 1), sum(time)))
    ),
    control = control
  )
  apply(fit$run(0L, 200L)$varcount, 1:2, sum)[, 1L]
}
unnamedPeriod <- hazardFirstForest(allowed)
expect_true(all(unnamedPeriod[3:4] == 0))
expect_true(all(unnamedPeriod[1:2] > 0))
namedPeriod <- hazardFirstForest(c("a", "period"))
expect_true(all(namedPeriod[2:3] == 0))
expect_true(all(namedPeriod[c(1L, 4L)] > 0))

# ---- beside the other constraints -------------------------------------------

# blocks partition the allowed columns: each tree within its group
inBlocks <- sampler(
  x,
  y,
  forests = restrictTo(allowed, blocks = blocks(groups = list("a", "b")))
)
expect_identical(attr(inBlocks$model, "block.of.column"), c(0L, 1L, -1L))
expect_identical(unname(splits(inBlocks$run(0L, 200L))[3L]), 0)
blockVars <- treeVars(inBlocks$getTrees())
expect_true(all(lengths(blockVars) == 1L))
expect_true(all(unlist(blockVars[as.integer(names(blockVars)) <= 10L]) == 1L))
expect_true(all(unlist(blockVars[as.integer(names(blockVars)) > 10L]) == 2L))
expect_true(setequal(unlist(blockVars), 1:2))
expect_error(
  sampler(
    x,
    y,
    forests = restrictTo(
      allowed,
      blocks = blocks(groups = list("a", c("b", "c")))
    )
  ),
  pattern = "names column\\(s\\) not among the forest's available predictors: c"
)

# and on the first forest of two, which was refused for not naming the columns
# it may not split on
z <- rbinom(n, 1L, 0.5)
twoForests <- function(groups) {
  set.seed(7L)
  dbarts(
    x,
    y + z * x[, "b"],
    forests = list(
      forest(vars = allowed, blocks = blocks(groups = groups)),
      forest(basis = z)
    ),
    control = control
  )
}
firstOfTwo <- twoForests(list("a", "b"))
expect_null(attr(firstOfTwo$model, "forest.columns"))
expect_identical(
  attr(firstOfTwo$control, "bartcore.forests")$vars,
  list(1:2, NULL)
)
firstCounts <- apply(firstOfTwo$run(0L, 200L)$varcount, 1:2, sum)
expect_true(firstCounts[3L, 1L] == 0)
expect_true(firstCounts[3L, 2L] > 0)
firstTrees <- firstOfTwo$getTrees(forest = 1L)
firstVars <- treeVars(firstTrees)
expect_true(all(lengths(firstVars) == 1L))
expect_true(setequal(unlist(firstVars), 1:2))
expect_error(
  twoForests(list("a", c("b", "c"))),
  pattern = "names column\\(s\\) not among the forest's available predictors: c"
)

# a variance forest keeps its own columns
withVariance <- sampler(x, y, forests = restrictTo(allowed), variance = "c")
expect_identical(unname(splits(withVariance$run(0L, 150L))[3L]), 0)
withVariance$storeState()
varianceVars <- withVariance$state[[1L]]$variance.vars
expect_true(any(varianceVars == 3L))
expect_true(all(varianceVars[varianceVars > 0L] == 3L))

# ---- for the sampler's whole life -------------------------------------------

staysWithin <- function(s, sweeps = 100L) {
  counts <- splits(s$run(0L, sweeps))
  counts[3L] == 0 && sum(counts[1:2]) > 0
}
live <- sampler(x, y, forests = restrictTo(allowed))
invisible(live$run(0L, 50L))
donor <- sampler(x, y)
invisible(donor$run(0L, 50L))
live$storeState()
donor$storeState()

# a copy and a reload rebuild the restriction from the model, and agree
duplicate <- live$copy()
expect_identical(attr(duplicate$model, "forest.columns"), 1:2)
stored <- tempfile(fileext = ".rds")
saveRDS(live, stored)
reloaded <- readRDS(stored)
unlink(stored)
copyDraws <- duplicate$run(0L, 100L)
reloadDraws <- reloaded$run(0L, 100L)
expect_identical(copyDraws$train, reloadDraws$train)
expect_identical(unname(splits(copyDraws)[3L]), 0)
expect_true(sum(splits(copyDraws)[1:2]) > 0)

# the sampler's own state installs; a state or a donor that splits on an
# excluded column is refused with the restricted variance forest's messages,
# and the sampler is left as it was
expect_true(live$setState(live$state))
stateBefore <- live$state
treesBefore <- live$getTrees()
expect_true(any(donor$getTrees()$var == 3L))
expect_error(
  live$setState(donor$state),
  pattern = "state holds a tree that splits on a variable outside this forest's allowed column set"
)
expect_error(
  live$installTrees(donor),
  pattern = "warm-start donor holds a tree that splits on a variable outside this forest's allowed column set"
)
expect_identical(live$state, stateBefore)
expect_identical(live$getTrees(), treesBefore)
live$installTrees(duplicate)
expect_identical(live$getTrees()$var, duplicate$getTrees()$var)
expect_true(staysWithin(live))

# new predictors and new data
expect_true(live$setPredictor(x[, "a"] + 0.1, 1L))
expect_true(staysWithin(live))
live$setPredictor(x[n:1, ])
expect_true(staysWithin(live))
live$setData(dbartsData(x, 2 * x[, "c"] + rnorm(n, 0, 0.5)))
expect_true(staysWithin(live))

# trees drawn from the prior and grown from the root
priorOutside <- 0L
plainPriorOutside <- 0L
for (draw in seq_len(30L)) {
  live$sampleTreesFromPrior()
  priorOutside <- priorOutside + sum(live$getTrees()$var == 3L)
  donor$sampleTreesFromPrior()
  plainPriorOutside <- plainPriorOutside + sum(donor$getTrees()$var == 3L)
}
expect_identical(priorOutside, 0L)
expect_true(plainPriorOutside > 0L)
live$growFromRoot(2L)
expect_identical(sum(live$getTrees()$var == 3L), 0L)
expect_true(sum(live$getTrees()$var %in% 1:2) > 0L)
donor$growFromRoot(2L)
expect_true(sum(donor$getTrees()$var == 3L) > 0L)
expect_true(staysWithin(live))

# ---- setModel changes parameters, never the restriction ---------------------

# the sampler's own model, edited, is accepted and keeps the restriction
edited <- live$model
edited@tree.prior@power <- 3
live$setModel(edited)
expect_identical(live$model@tree.prior@power, 3)
expect_identical(attr(live$copy()$model, "forest.columns"), 1:2)
expect_true(staysWithin(live))

# split probabilities are a parameter setModel may change, held to what
# creation holds them to: positive on some allowed column
edited@tree.prior@splitProbabilities <- c(0.25, 0.25, 0.5)
live$setModel(edited)
expect_true(staysWithin(live))
modelBefore <- live$model
edited@tree.prior@splitProbabilities <- c(0, 0, 1)
expect_error(
  live$setModel(edited),
  pattern = "'split.probs' gives no positive probability to any column"
)
expect_identical(live$model, modelBefore)
expect_true(all(splits(live$run(0L, 100L))[1:2] > 0))

# a model with another restriction or none is refused by name, in both
# directions, before anything is stored
otherModel <- sampler(x, y, forests = restrictTo("a"))$model
modelRefusal <- "cannot change the columns a forest may split on, its 'vars'"
for (target in list(live, donor)) {
  given <- if (identical(target, live)) {
    list(donor$model, otherModel)
  } else {
    list(live$model, otherModel)
  }
  target$storeState()
  modelBefore <- target$model
  stateBefore <- target$state
  for (model in given) {
    expect_error(target$setModel(model), pattern = modelRefusal)
  }
  expect_identical(target$model, modelBefore)
  target$storeState()
  expect_identical(target$state, stateBefore)
}
expect_true(staysWithin(live))
expect_true(splits(donor$run(0L, 100L))[3L] > 0)
