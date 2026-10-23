# A linear or Gaussian-process leaf's covariate standardization - each
# covariate's center and scale, and a gp leaf's heuristic lengthscale - is
# data its sampler derived. It is recorded on the sampler's data object as
# standardization, written from the engine after creation and after every
# call that changes it, and taken by every creation from the object, so a
# copy and a reload hold it and no state brings one. An install, setState
# and a warm start alike, keeps the sampler's standardization and puts a
# linear leaf's slopes in as stored (dec-B457); kept gp draws made under
# another center, scale or lengthscale are refused, a kept gp draw replaying
# only under its own. xbart's folds derive their own.

set.seed(57L, sample.kind = "Rejection")
n <- 120L
x <- cbind(a = runif(n), b = runif(n))
y <- 2 * x[, 1L] + sin(3 * x[, 2L]) + rnorm(n, 0, 0.2)
x.new <- cbind(a = runif(30L, -0.2, 1.2), b = runif(30L))
# the same rows with the leaf covariate on another scale
stretched <- x
stretched[, 1L] <- 3 * x[, 1L] + 5
controlWith <- function(...) {
  settings <- list(
    n.chains = 2L,
    n.threads = 1L,
    n.trees = 10L,
    n.samples = 4L,
    n.burn = 0L,
    keepTrees = TRUE,
    updateState = FALSE,
    seed = 7L
  )
  do.call(dbarts::dbartsControl, utils::modifyList(settings, list(...)))
}
# a sampler of linear leaves on column a unless another leaf prior is named;
# control settings by name
make <- function(
  ...,
  data = x,
  response = y,
  sweeps = 0L,
  control = list()
) {
  named <- names(match.call(expand.dots = FALSE)$...)
  if (!"leaf.prior" %in% named) {
    sampler <- dbarts::dbarts(
      data,
      response,
      leaf.prior = linear(columns = "a"),
      control = do.call(controlWith, control),
      ...
    )
  } else {
    sampler <- dbarts::dbarts(
      data,
      response,
      control = do.call(controlWith, control),
      ...
    )
  }
  if (sweeps > 0L) {
    invisible(sampler$run(sweeps, 4L))
  }
  sampler
}
stateOf <- function(sampler) {
  .Call(dbarts:::C_dbarts_bartcore_storeState, sampler$getPointer())
}
readerOf <- function(sampler) {
  .Call(
    dbarts:::C_dbarts_bartcore_derivedRecords,
    sampler$getPointer(),
    TRUE
  )$standardization
}
# the standardization a chain's leaves hold, as its stored state records it
inForce <- function(sampler, chain = 1L) {
  forest <- stateOf(sampler)[[chain]]$forests[[1L]]
  held <- list(
    center = forest$leaf.covariate.center,
    scale = forest$leaf.covariate.scale
  )
  if (!is.null(forest$leaf.lengthscales)) {
    held$lengthscale <- forest$leaf.lengthscales
  }
  held
}
# the sampler's record without the covariates' names
recordOf <- function(sampler) {
  lapply(sampler$data@standardization[[1L]], unname)
}
liveOf <- function(state) {
  lapply(state, function(chain) {
    chain$forests[[1L]][c("tree.vars", "tree.values", "tree.params")]
  })
}
reloadOf <- function(sampler) {
  path <- tempfile(fileext = ".rds")
  saveRDS(sampler, path)
  reloaded <- readRDS(path)
  unlink(path)
  reloaded
}
noSplitsOnA <- function(...) {
  make(..., tree.prior = cgm(split.probs = c(0, 1)))
}

# ---- written at creation ---------------------------------------------------

# the record is the engine's, named by the covariates, and the object handed
# in carries none
handed <- dbarts::dbartsData(x, y)
linearFit <- dbarts::dbarts(
  handed,
  leaf.prior = linear(columns = "a"),
  control = controlWith()
)
record <- linearFit$data@standardization
expect_identical(length(record), 1L)
expect_identical(names(record[[1L]]), c("center", "scale"))
expect_identical(names(record[[1L]]$center), "a")
expect_identical(names(record[[1L]]$scale), "a")
expect_equal(unname(record[[1L]]$center), mean(x[, 1L]), tolerance = 1e-14)
expect_equal(unname(record[[1L]]$scale), sd(x[, 1L]), tolerance = 1e-14)
expect_identical(recordOf(linearFit), inForce(linearFit))
expect_identical(recordOf(linearFit), inForce(linearFit, 2L))
expect_identical(recordOf(linearFit), readerOf(linearFit)[[1L]])
expect_null(handed@standardization)
# a gp leaf records its heuristic lengthscale beside them, and one whose
# model supplies it does not
gpFit <- make(leaf.prior = gp(columns = "a"))
expect_identical(
  names(gpFit$data@standardization[[1L]]),
  c(
    "center",
    "scale",
    "lengthscale"
  )
)
expect_identical(recordOf(gpFit), inForce(gpFit))
suppliedFit <- make(leaf.prior = gp(columns = "a", lengthscale = 0.5))
expect_identical(
  names(suppliedFit$data@standardization[[1L]]),
  c(
    "center",
    "scale"
  )
)
# two covariates in the leaf's order, a column without spread NA, columns
# without names an unnamed record
flat <- x
flat[, 2L] <- 3
twoFit <- make(data = flat, leaf.prior = linear(columns = c("b", "a")))
expect_identical(names(twoFit$data@standardization[[1L]]$center), c("b", "a"))
expect_identical(twoFit$data@standardization[[1L]]$scale[["b"]], NA_real_)
expect_identical(twoFit$data@standardization[[1L]]$center[["b"]], 3)
expect_identical(recordOf(twoFit), inForce(twoFit))
unnamedFit <- make(data = unname(x), leaf.prior = linear(columns = 1L))
expect_null(names(unnamedFit$data@standardization[[1L]]$center))
expect_identical(recordOf(unnamedFit), recordOf(linearFit))
# leaves that take no covariate record none
expect_null(make(leaf.prior = normal())$data@standardization)
expect_null(
  dbarts::dbarts(
    x,
    y,
    forests = list(
      dbarts::dbartsForests$forest(),
      dbarts::dbartsForests$forest(basis = factor(rep_len(0:1, n)))
    ),
    control = controlWith()
  )$data@standardization
)

# ---- same frame exact ------------------------------------------------------

# under the standardization in force a copy, a reload and a twin handed the
# stored state hold one state, the record and the kept draws, and draw alike.
# A gp tree's flat leaf values are read off its fits and can move in the last
# place when a tree is rebuilt, as they did before the record, so a gp state
# is compared by its fits
fixtures <- list(
  linear = function(seed = 7L) {
    make(sweeps = 20L, control = list(seed = seed))
  },
  gp = function(seed = 7L) {
    make(
      leaf.prior = gp(columns = "a"),
      sweeps = 20L,
      control = list(seed = seed)
    )
  },
  two = function(seed = 7L) {
    make(
      leaf.prior = linear(columns = c("a", "b")),
      sweeps = 20L,
      control = list(seed = seed)
    )
  }
)
for (name in names(fixtures)) {
  sampler <- fixtures[[name]]()
  sampler$storeState()
  stored <- sampler$state
  copied <- sampler$copy()
  reloaded <- reloadOf(sampler)
  twin <- fixtures[[name]](8L)
  expect_true(twin$setState(stored), info = name)
  for (door in list(copied, reloaded, twin)) {
    expect_identical(door$data@standardization, sampler$data@standardization)
    expect_identical(inForce(door), inForce(sampler), info = name)
    expect_identical(door$predict(x.new), sampler$predict(x.new), info = name)
  }
  for (door in list(copied, reloaded)) {
    if (name == "gp") {
      expect_identical(
        lapply(stateOf(door), function(chain) chain$forests[[1L]]$tree.params),
        lapply(stored, function(chain) chain$forests[[1L]]$tree.params),
        info = name
      )
    } else {
      expect_identical(stateOf(door), stored, info = name)
    }
  }
  draws <- twin$run(0L, 10L)
  expect_identical(copied$run(0L, 10L), draws, info = name)
  expect_identical(reloaded$run(0L, 10L), draws, info = name)
}

# ---- keeps its own ---------------------------------------------------------

# A state stored before setData gave the leaf covariate another scale: forced
# or not, setState keeps the sampler's standardization, the one setData
# derived, and puts the slopes in as stored. Column a is never split, so the
# trees fit both grids and the call unforced installs them
keeper <- noSplitsOnA(sweeps = 30L)
keeper$storeState()
stored <- keeper$state
storedOwn <- inForce(keeper)
keeper$setData(dbarts::dbartsData(stretched, y))
own <- inForce(keeper)
expect_equal(own$center, mean(stretched[, 1L]), tolerance = 1e-14)
expect_equal(own$scale, sd(stretched[, 1L]), tolerance = 1e-14)
expect_true(abs(own$center - storedOwn$center) > 4)
record <- keeper$data@standardization
expect_true(keeper$setState(stored))
expect_identical(inForce(keeper), own)
expect_identical(inForce(keeper, 2L), own)
expect_identical(keeper$data@standardization, record)
expect_identical(liveOf(stateOf(keeper)), liveOf(stored))
expect_null(keeper$setState(stored, forceUpdate = TRUE))
expect_identical(inForce(keeper), own)
expect_identical(liveOf(stateOf(keeper)), liveOf(stored))
# the kept draws predict under it, as in a sampler made on those values
madeThere <- noSplitsOnA(data = stretched, control = list(seed = 8L))
expect_true(madeThere$setState(stored))
expect_identical(inForce(madeThere), own)
expect_identical(madeThere$predict(x.new), keeper$predict(x.new))

# ---- copy holds it ---------------------------------------------------------

# setData derives a standardization over the stretched column and setPredictor
# then hands the old values back, which a derivation would read as another:
# a copy and a reload hold the record's, and predict the source's kept draws
holderOf <- function(...) {
  sampler <- make(...)
  sampler$setData(dbarts::dbartsData(stretched, y))
  sampler$setPredictor(x[, 1L], 1L, forceUpdate = TRUE)
  invisible(sampler$run(20L, 4L))
  sampler
}
for (leaf in c("linear", "gp")) {
  holder <- if (leaf == "linear") {
    holderOf()
  } else {
    holderOf(leaf.prior = gp(columns = "a"))
  }
  held <- holder$data@standardization
  expect_equal(
    unname(held[[1L]]$center),
    mean(stretched[, 1L]),
    tolerance = 1e-14,
    info = leaf
  )
  expect_identical(holder$data@x, x, info = leaf)
  holder$storeState()
  for (door in list(holder$copy(), reloadOf(holder))) {
    expect_identical(door$data@standardization, held, info = leaf)
    expect_identical(inForce(door), recordOf(holder), info = leaf)
    expect_identical(door$predict(x.new), holder$predict(x.new), info = leaf)
  }
}

# ---- carried at creation ---------------------------------------------------

# new(), dbarts() and bart() handed that data object take its record, where a
# fresh object of the same values derives another
holder <- holderOf()
carried <- holder$data
doors <- list(
  new = new("dbartsSampler", holder$control, holder$model, carried),
  dbarts = dbarts::dbarts(
    carried,
    leaf.prior = linear(columns = "a"),
    control = controlWith()
  ),
  bart = dbarts::bart(
    carried,
    leaf.prior = linear(columns = "a"),
    n.trees = 10L,
    n.samples = 2L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    verbose = FALSE
  )$fit
)
for (door in names(doors)) {
  expect_identical(inForce(doors[[door]]), recordOf(holder), info = door)
  expect_identical(
    doors[[door]]$data@standardization,
    carried@standardization,
    info = door
  )
}
derived <- make(data = carried@x)
expect_false(identical(inForce(derived), recordOf(holder)))
# xbart's folds derive their own: the carried record changes nothing, and
# folds whose leaves take no covariate are never handed it
crossValidated <- function(data, leafPrior = quote(linear(columns = "a"))) {
  set.seed(9L)
  eval(bquote(dbarts::xbart(
    data,
    leaf.prior = .(leafPrior),
    n.samples = 3L,
    n.burn = 2L,
    n.reps = 1L,
    n.trees = 5L,
    k = 2,
    n.threads = 1L,
    verbose = FALSE
  )))
}
for (leafPrior in list(quote(linear(columns = "a")), quote(normal()))) {
  expect_identical(
    crossValidated(carried, leafPrior),
    crossValidated(dbarts::dbartsData(carried@x, carried@y), leafPrior)
  )
}

# ---- gp kept refused -------------------------------------------------------

# Kept gp draws made under another center and scale, or another heuristic
# lengthscale, are refused in both forms, the sampler untouched; a state
# that keeps none installs its live fits under the sampler's own
gpDonor <- make(leaf.prior = gp(columns = "a"), sweeps = 20L)
gpDonor$storeState()
gpStored <- gpDonor$state
gpOther <- make(
  leaf.prior = gp(columns = "a"),
  data = stretched,
  sweeps = 5L,
  control = list(seed = 8L)
)
before <- stateOf(gpOther)
refusedGp <- "another covariate standardization or other lengthscales"
for (force in c(FALSE, TRUE)) {
  expect_error(
    gpOther$setState(gpStored, forceUpdate = force),
    refusedGp,
    fixed = TRUE
  )
  expect_identical(stateOf(gpOther), before)
}
# the center and scale the same, the lengthscale another: a record handed at
# creation
longer <- dbarts::dbartsData(x, y)
longer@standardization <- gpDonor$data@standardization
longer@standardization[[1L]]$lengthscale <- 2 *
  longer@standardization[[1L]]$lengthscale
gpLonger <- dbarts::dbarts(
  longer,
  leaf.prior = gp(columns = "a"),
  control = controlWith(seed = 9L)
)
expect_identical(
  inForce(gpLonger)[c("center", "scale")],
  inForce(gpDonor)[c("center", "scale")]
)
expect_identical(
  inForce(gpLonger)$lengthscale,
  2 * inForce(gpDonor)$lengthscale
)
before <- stateOf(gpLonger)
expect_error(gpLonger$setState(gpStored), refusedGp, fixed = TRUE)
expect_identical(stateOf(gpLonger), before)
# the scale alone another, by a record with it doubled
wider <- dbarts::dbartsData(x, y)
wider@standardization <- gpDonor$data@standardization
wider@standardization[[1L]]$scale <- 2 * wider@standardization[[1L]]$scale
gpWider <- dbarts::dbarts(
  wider,
  leaf.prior = gp(columns = "a"),
  control = controlWith(seed = 9L)
)
expect_identical(
  inForce(gpWider)[c("center", "lengthscale")],
  inForce(gpDonor)[c("center", "lengthscale")]
)
expect_identical(inForce(gpWider)$scale, 2 * inForce(gpDonor)$scale)
expect_error(gpWider$setState(gpStored), refusedGp, fixed = TRUE)
# and by data: a covariate centered exactly at 0 and then doubled keeps its
# center and its heuristic lengthscale, read off the standardized values
centered <- x
centered[, 1L] <- c(seq_len(n / 2L), -seq_len(n / 2L)) / 64
doubled <- centered
doubled[, 1L] <- 2 * centered[, 1L]
centeredDonor <- noSplitsOnA(
  leaf.prior = gp(columns = "a"),
  data = centered,
  sweeps = 5L
)
centeredDonor$storeState()
doubledSampler <- noSplitsOnA(
  leaf.prior = gp(columns = "a"),
  data = doubled,
  control = list(seed = 8L)
)
expect_identical(inForce(centeredDonor)$center, 0)
expect_identical(
  inForce(doubledSampler)[c("center", "lengthscale")],
  inForce(centeredDonor)[c("center", "lengthscale")]
)
expect_identical(
  inForce(doubledSampler)$scale,
  2 * inForce(centeredDonor)$scale
)
before <- stateOf(doubledSampler)
expect_error(
  doubledSampler$setState(centeredDonor$state),
  refusedGp,
  fixed = TRUE
)
expect_identical(stateOf(doubledSampler), before)
# live fits alone, the state keeping no draw: installed, the sampler's own
# standardization and lengthscale standing
gpLive <- make(
  leaf.prior = gp(columns = "a"),
  sweeps = 20L,
  control = list(keepTrees = FALSE)
)
gpLive$storeState()
own <- inForce(gpLonger)
expect_true(gpLonger$setState(gpLive$state))
expect_identical(inForce(gpLonger), own)
own <- inForce(gpOther)
expect_null(gpOther$setState(gpLive$state, forceUpdate = TRUE))
expect_identical(inForce(gpOther), own)

# ---- supplied lengthscale --------------------------------------------------

# a lengthscale the model supplies is the model's: a record's is not read
givenLengthscale <- dbarts::dbartsData(x, y)
givenLengthscale@standardization <- list(list(
  center = c(a = 0.5),
  scale = c(a = 0.3),
  lengthscale = c(a = 3)
))
supplied <- dbarts::dbarts(
  givenLengthscale,
  leaf.prior = gp(columns = "a", lengthscale = 0.5),
  control = controlWith()
)
expect_identical(
  inForce(supplied),
  list(
    center = 0.5,
    scale = 0.3,
    lengthscale = 0.5
  )
)
expect_identical(
  supplied$data@standardization,
  list(list(center = c(a = 0.5), scale = c(a = 0.3)))
)

# ---- warm start as stored --------------------------------------------------

# A linear donor under another standardization, live and from a kept draw:
# clean, the recipient's slopes the donor's stored numbers and its
# standardization its own. A twin made on the donor's standardization takes
# the same numbers, which there are the draw's function. A column without
# spread on each side, at another value, is no change either
donor <- noSplitsOnA(sweeps = 30L)
donor$storeState()
newest <- c(4L, 8L)
for (samples in list(newest, c(2L, 5L))) {
  info <- paste(samples, collapse = " ")
  recipient <- noSplitsOnA(
    data = stretched,
    control = list(keepTrees = FALSE, seed = 9L)
  )
  own <- inForce(recipient)
  expect_true(recipient$installTrees(donor, samples = samples), info = info)
  expect_identical(inForce(recipient), own, info = info)
  twin <- noSplitsOnA(control = list(keepTrees = FALSE, seed = 10L))
  expect_true(twin$installTrees(donor, samples = samples), info = info)
  expect_identical(liveOf(stateOf(recipient)), liveOf(stateOf(twin)))
  # one function of the standardized covariate, two of the raw one
  expect_true(
    max(abs(recipient$predict(x.new) - twin$predict(x.new))) > 1e-3,
    info = info
  )
  if (identical(samples, newest)) {
    # warm start, same frame: the twin holds the donor's live trees and, to
    # the rounding of the donor's running sums, its fits
    expect_identical(liveOf(stateOf(twin)), liveOf(donor$state))
    expect_equal(
      twin$getFitsWithoutOffset(),
      donor$getFitsWithoutOffset(),
      tolerance = 1e-12
    )
  }
}
constantAt <- function(value) {
  frame <- x
  frame[, 1L] <- value
  frame
}
flatDonor <- noSplitsOnA(data = constantAt(0), sweeps = 30L)
flatRecipient <- noSplitsOnA(
  data = constantAt(3),
  control = list(keepTrees = FALSE, seed = 9L)
)
expect_identical(inForce(flatRecipient)$scale, NA_real_)
expect_true(flatRecipient$installTrees(flatDonor, samples = newest))
flatDonor$storeState()
expect_identical(liveOf(stateOf(flatRecipient)), liveOf(flatDonor$state))
expect_identical(inForce(flatRecipient)$center, 3)

# ---- matched to the leaves -------------------------------------------------

# Which covariates the leaves take is the model's, so a record is matched to
# them by name, never refused for them: an entry the leaves do not take is
# ignored, and a covariate without one derived
withRecord <- function(standardization) {
  data <- dbarts::dbartsData(x, y)
  data@standardization <- standardization
  data
}
fits <- list(list(center = c(a = 0.5), scale = c(a = 0.3)))
fresh <- function(leafPrior) {
  inForce(eval(bquote(make(leaf.prior = .(leafPrior)))))
}
# leaves taking no covariate, at creation and at setData
constantFit <- dbarts::dbarts(carried, control = controlWith())
expect_null(constantFit$data@standardization)
constantBart <- dbarts::bart(
  carried,
  n.trees = 10L,
  n.samples = 2L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)$fit
expect_null(constantBart$data@standardization)
constantFit$setData(carried)
expect_null(constantFit$data@standardization)
# leaves taking another covariate derive it, and taking another beside the
# recorded one use the record's for it
otherCovariate <- dbarts::dbarts(
  carried,
  leaf.prior = linear(columns = "b"),
  control = controlWith()
)
expect_identical(inForce(otherCovariate), fresh(quote(linear(columns = "b"))))
expect_identical(names(otherCovariate$data@standardization[[1L]]$center), "b")
both <- dbarts::dbarts(
  carried,
  leaf.prior = linear(columns = c("b", "a")),
  control = controlWith()
)
expect_identical(
  lapply(recordOf(both), `[`, 2L),
  recordOf(holder)
)
expect_identical(
  lapply(recordOf(both), `[`, 1L),
  fresh(quote(linear(columns = "b")))
)
# an entry for another column, an unnamed one of another count and a second
# forest's: the first two derive, the third is not read
served <- list(
  "another column" = list(
    record = list(list(center = c(b = 0.5), scale = c(b = 0.3))),
    held = fresh(quote(linear(columns = "a")))
  ),
  "two covariates" = list(
    record = list(list(center = c(0.5, 0.5), scale = c(0.3, 0.3))),
    held = fresh(quote(linear(columns = "a")))
  ),
  "a second forest's" = list(
    record = c(fits, fits),
    held = list(center = 0.5, scale = 0.3)
  ),
  # a logical NA is a scale without spread, as NA_real_ is
  "a logical NA" = list(
    record = list(list(center = c(a = 0.5), scale = c(a = NA))),
    held = inForce(dbarts::dbarts(
      withRecord(list(list(center = c(a = 0.5), scale = c(a = NA_real_)))),
      leaf.prior = linear(columns = "a"),
      control = controlWith()
    ))
  )
)
sampler <- make(sweeps = 5L)
for (case in names(served)) {
  taken <- dbarts::dbarts(
    withRecord(served[[case]]$record),
    leaf.prior = linear(columns = "a"),
    control = controlWith()
  )
  expect_identical(inForce(taken), served[[case]]$held, info = case)
  sampler$setData(withRecord(served[[case]]$record))
}
# one that fits, unnamed or named, is taken
for (record in list(fits, list(lapply(fits[[1L]], unname)))) {
  taken <- dbarts::dbarts(
    withRecord(record),
    leaf.prior = linear(columns = "a"),
    control = controlWith()
  )
  expect_identical(inForce(taken), list(center = 0.5, scale = 0.3))
  expect_identical(taken$data@standardization, fits)
}

# setData reads the columns by name, else by position, the record taking the
# new names: a sampler's data object reordered, stripped of its names or
# renamed is served, linear and gp, and the record is the reader's under the
# names the data then carry
variants <- list(
  reordered = function(data) {
    data@x <- data@x[, c("b", "a")]
    data
  },
  unnamed = function(data) {
    colnames(data@x) <- NULL
    data
  },
  renamed = function(data) {
    colnames(data@x) <- c("c", "d")
    data
  }
)
for (leaf in c("linear", "gp")) {
  for (variant in names(variants)) {
    swapped <- if (leaf == "linear") {
      make(control = list(keepTrees = FALSE))
    } else {
      make(leaf.prior = gp(columns = "a"), control = list(keepTrees = FALSE))
    }
    swapped$setData(variants[[variant]](swapped$data))
    info <- paste(leaf, variant)
    expect_identical(
      names(swapped$data@standardization[[1L]]$center),
      colnames(swapped$data@x)[1L],
      info = info
    )
    expect_identical(recordOf(swapped), readerOf(swapped)[[1L]], info = info)
  }
}

# ---- malformed -------------------------------------------------------------

# a record not of the slot's form is refused at creation, whatever the leaves
# take, and at setData, naming the slot, the sampler untouched
malformed <- list(
  "not a list" = list(
    record = c(a = 0.5),
    message = "must be a list with one entry per forest"
  ),
  "an entry not a list" = list(
    record = list(c(0.5, 0.3)),
    message = "must be a list with one entry per forest"
  ),
  "no entry" = list(
    record = list(),
    message = "must be a list with one entry per forest"
  ),
  "a scale of zero" = list(
    record = list(list(center = c(a = 0.5), scale = c(a = 0))),
    message = "must give a finite 'center', a positive or NA 'scale'"
  ),
  "a center not finite" = list(
    record = list(list(center = c(a = NA), scale = c(a = 0.3))),
    message = "must give a finite 'center', a positive or NA 'scale'"
  ),
  "a center not numeric" = list(
    record = list(list(center = c(a = "0.5"), scale = c(a = 0.3))),
    message = "must give a finite 'center', a positive or NA 'scale'"
  ),
  "lengths that differ" = list(
    record = list(list(center = c(a = 0.5), scale = c(0.3, 0.3))),
    message = "must give a finite 'center', a positive or NA 'scale'"
  ),
  "a lengthscale of zero" = list(
    record = list(list(
      center = c(a = 0.5),
      scale = c(a = 0.3),
      lengthscale = c(a = 0)
    )),
    message = "must give a finite 'center', a positive or NA 'scale'"
  )
)
before <- stateOf(sampler)
for (case in names(malformed)) {
  message <- paste0(
    "the data's 'standardization' record ",
    malformed[[case]]$message
  )
  data <- withRecord(malformed[[case]]$record)
  for (leafPrior in list(quote(linear(columns = "a")), quote(normal()))) {
    expect_error(
      eval(bquote(dbarts::dbarts(
        data,
        leaf.prior = .(leafPrior),
        control = controlWith()
      ))),
      message,
      fixed = TRUE,
      info = case
    )
  }
  expect_error(sampler$setData(data), message, fixed = TRUE, info = case)
  expect_identical(stateOf(sampler), before, info = case)
}

# ---- record-less object ----------------------------------------------------

# an object saved before the slot existed records nothing and derives, as a
# fresh object does
old <- holder$data
attr(old, "standardization") <- NULL
expect_false(methods::.hasSlot(old, "standardization"))
fromOld <- dbarts::dbarts(
  old,
  leaf.prior = linear(columns = "a"),
  control = controlWith()
)
expect_identical(inForce(fromOld), inForce(derived))
