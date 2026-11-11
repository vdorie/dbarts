# A data object built from a formula and a data frame keeps the frame, the
# positions of the rows it kept and what reads the formula again; x, y, weights
# and offset hold a stand-in that stops, naming model.frame(o). Every door that
# takes a data object reads it through one resolution to the block-holding
# object keep.data = FALSE gives, so a fit from the object is the fit from the
# formula.

set.seed(51)
n <- 40L
d <- data.frame(
  y = rnorm(n),
  x1 = rnorm(n),
  x2 = runif(n),
  f = factor(sample(c("a", "b", "c"), n, TRUE)),
  w = runif(n) + 0.5
)
d$x1[3L] <- NA
d$x2[9L] <- NA

ctl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = 2L,
  updateState = FALSE
)
draws <- function(s) s$run(2L, 5L)$train
fitDraws <- function(data, family = "auto") {
  draws(dbarts(data, n.trees = 5L, control = ctl, seed = 7L, family = family))
}

# the slots -------------------------------------------------------------------

o <- dbartsData(y ~ x1 + x2 + f, d, weights = w)
expect_identical(o@frame, d)
expect_identical(o@frame.columns, names(d))
expect_null(o@rows)
expect_true(inherits(o@x, "dbartsStandIn"))
expect_true(inherits(o@weights, "dbartsStandIn"))
expect_null(dbartsData(y ~ x1 + x2 + f, d)@weights)
expect_null(dbartsData(y ~ x1 + x2 + f, d)@offset)

b <- dbartsData(y ~ x1 + x2 + f, d, weights = w, keep.data = FALSE)
expect_null(b@frame)
expect_identical(b@frame.columns, names(d))
expect_true(is.numeric(b@y))
expect_true(is.numeric(b@weights))
expect_identical(nrow(b@x), nrow(d) - 0L)

expect_silent(same <- dbartsData(o, keep.data = FALSE))
expect_identical(same, o)

# the positions of the rows kept ---------------------------------------------

expect_identical(dbartsData(y ~ x1 + x2, d, subset = 2:10)@rows, 2:10)
expect_identical(dbartsData(y ~ x1 + x2, d, subset = -1)@rows, 2:n)
expect_identical(
  dbartsData(y ~ x1 + x2, d, subset = c(10, 2, 5, 7))@rows,
  c(10L, 2L, 5L, 7L)
)
expect_identical(
  dbartsData(y ~ x1 + x2, d, subset = c(2, 2, 5, 7, 8))@rows,
  c(2L, 2L, 5L, 7L, 8L)
)
expect_identical(
  dbartsData(y ~ x1 + x2, d, subset = f != "b")@rows,
  which(d$f != "b")
)
expect_identical(
  dbartsData(y ~ x1 + x2, d, na.action = na.omit)@rows,
  setdiff(seq_len(n), c(3L, 9L))
)
expect_identical(
  dbartsData(y ~ x1 + x2, d, subset = 2:10, na.action = na.omit)@rows,
  setdiff(2:10, c(3L, 9L))
)
dn <- d
rownames(dn) <- paste0("r", seq_len(n))
expect_identical(
  dbartsData(y ~ x1 + x2, dn, subset = c("r4", "r2", "r4"))@rows,
  c(4L, 2L, 4L)
)
sub <- replace(rep(TRUE, n), c(5L, 6L), NA)
oNA <- dbartsData(y ~ x1 + x2, d, subset = sub, na.action = na.exclude)
expect_identical(oNA@rows, setdiff(seq_len(n), c(3L, 5L, 6L, 9L)))
expect_identical(nrow(model.frame(oNA)), n - 4L)

# an na.action that records nothing of what it drops: the rows kept are the
# model frame's, at every door, as lm's are
noNote <- function(object, ...) {
  object[stats::complete.cases(object), , drop = FALSE]
}
lmNoNote <- lm(y ~ x1 + x2, d, na.action = noNote)
oNoNote <- dbartsData(y ~ x1 + x2, d, na.action = noNote)
expect_identical(oNoNote@rows, setdiff(seq_len(n), c(3L, 9L)))
expect_identical(nobs(oNoNote), nobs(lmNoNote))
expect_identical(nrow(model.frame(oNoNote)), n - 2L)
expect_identical(
  unname(model.response(model.frame(oNoNote))),
  unname(model.response(model.frame(lmNoNote)))
)
bNoNote <- dbartsData(y ~ x1 + x2, d, na.action = noNote, keep.data = FALSE)
expect_identical(bNoNote@rows, oNoNote@rows)
expect_identical(dbarts:::resolveDataObject(oNoNote)@x, bNoNote@x)
expect_identical(
  fitDraws(oNoNote),
  draws(dbarts(
    y ~ x1 + x2,
    d,
    na.action = noNote,
    n.trees = 5L,
    control = ctl,
    seed = 7L
  ))
)
bartNoNote <- function(data, ...) {
  bart(
    data,
    ...,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 7L,
    verbose = FALSE,
    keepTrees = FALSE
  )$yhat.train
}
expect_identical(
  bartNoNote(oNoNote),
  bartNoNote(y ~ x1 + x2, d, na.action = noNote)
)
expect_identical(ncol(bartNoNote(oNoNote)), n - 2L)
expect_identical(
  dbartsData(y ~ x1 + x2, as.list(d), na.action = noNote)@rows,
  oNoNote@rows
)
expect_identical(
  dbartsData(y ~ x1 + x2, d, subset = 2:20, na.action = noNote)@rows,
  setdiff(2:20, c(3L, 9L))
)
# a full-length basis is cut by the same rows
basisNoNote <- dbartsSpec(
  oNoNote,
  forests = list(forest(), forest(basis = d$w))
)
expect_identical(NROW(basisNoNote$data@bases[[2L]]), n - 2L)

# the stand-ins ---------------------------------------------------------------

stops <- function(expr) {
  expect_error(expr, "model.frame\\(o\\)")
}
stops(dim(o@x))
stops(nrow(o@x))
stops(colnames(o@x))
stops(length(o@y))
stops(o@x[1L, ])
stops(o@y[[1L]])
stops(o@x$dense)
stops(as.matrix(o@x))
stops(as.vector(o@y))
stops(o@y + 1)
stops(sum(o@y))
stops(sqrt(o@y))
stops(mtfrm(o@y))
expect_null(attr(o@x, "factor.levels"))
expect_null(attr(o@x, "term.labels"))
expect_silent(capture.output(print(o@y)))
expect_silent(capture.output(str(o)))
expect_silent(capture.output(show(o)))

# all.equal() on an object that keeps its frame: TRUE against itself and a
# copy, and a difference in the data reported, as 0.9-34 did -------------------

oEq <- dbartsData(y ~ x1 + x2, d)
expect_true(all.equal(oEq, oEq))
expect_true(all.equal(oEq, dbartsData(y ~ x1 + x2, d)))
dEq <- d
dEq$y[3L] <- dEq$y[3L] + 10
expect_false(isTRUE(all.equal(oEq, dbartsData(y ~ x1 + x2, dEq))))
expect_false(
  isTRUE(all.equal(oEq, dbartsData(y ~ x1 + x2, d, weights = rep(1, n))))
)
rm(oEq, dEq)

# model.frame() and nobs() ----------------------------------------------------

offv <- d$w / 10
fit <- lm(
  y ~ x1 + x2 + f,
  d,
  subset = 2:30,
  weights = w,
  offset = offv,
  na.action = na.omit
)
oOff <- dbartsData(
  y ~ x1 + x2 + f,
  d,
  subset = 2:30,
  weights = w,
  offset = offv,
  na.action = na.omit
)
mf <- model.frame(oOff)
expect_identical(model.response(mf), model.response(model.frame(fit)))
expect_identical(model.weights(mf), model.weights(model.frame(fit)))
expect_identical(model.offset(mf), model.offset(model.frame(fit)))
expect_true(is.factor(mf$f))
expect_identical(nobs(oOff), nrow(mf))
expect_identical(nobs(oOff), nobs(fit))

dz <- d
dz$w[c(4L, 5L)] <- 0
oZ <- suppressWarnings(dbartsData(y ~ x1 + x2, dz, weights = w))
# the rows with a positive weight, as lm counts them (the data keep their NAs,
# which lm drops, so the comparison is on complete data)
expect_identical(nobs(oZ), n - 2L)
dzc <- dz
dzc$x1[3L] <- 0
dzc$x2[9L] <- 0
expect_identical(
  nobs(suppressWarnings(dbartsData(y ~ x1 + x2, dzc, weights = w))),
  nobs(lm(y ~ x1 + x2, dzc, weights = w))
)
# the same count where the block holds the data
expect_identical(
  nobs(suppressWarnings(
    dbartsData(y ~ x1 + x2, dz, weights = w, keep.data = FALSE)
  )),
  n - 2L
)
rm(dzc)
expect_identical(nobs(lm(y ~ x1 + x2, dz, weights = w)), n - 4L)

# nobs() reads the count the build recorded and leaves the random-number
# stream alone, whatever the frame and its response
streamAfter <- function(expr) {
  set.seed(5L)
  force(expr)
  runif(1L)
}
listO <- dbartsData(y ~ x1, list(y = d$y, x1 = d$x1))
expect_false(is.null(listO@frame))
expect_identical(streamAfter(nobs(listO)), streamAfter(NULL))
expect_identical(nobs(listO), n)
drawnO <- dbartsData(rnorm(20L) ~ x1, list(x1 = d$x1[1:20]))
expect_identical(streamAfter(nobs(drawnO)), streamAfter(NULL))
expect_identical(nobs(drawnO), 20L)
expect_identical(streamAfter(nobs(oOff)), streamAfter(NULL))

# the kind of data the column names came from is a slot, and the response's
# name is read from the terms record; no loose attribute carries either
envData <- list2env(as.list(d[c("y", "x1", "x2")]))
envO <- dbartsData(y ~ x1 + x2, data = envData)
expect_identical(envO@frame.source, "environment")
expect_identical(dbartsData(y ~ x1 + x2, d)@frame.source, "frame")
expect_identical(
  dbartsData(y ~ x1 + x2, as.list(d), keep.data = FALSE)@frame.source,
  "frame"
)
expect_null(dbartsData(as.matrix(d[c("x1", "x2")]), d$y)@frame.source)
for (attribute in c("frame.environment", "response.label")) {
  expect_null(attr(envO, attribute, exact = TRUE))
  expect_null(attr(
    dbartsData(log(abs(y) + 1) ~ x1, d, keep.data = FALSE),
    attribute,
    exact = TRUE
  ))
}
labeled <- dbartsData(log(abs(y) + 1) ~ x1, d, keep.data = FALSE)
expect_identical(names(model.frame(labeled))[1L], "log(abs(y) + 1)")
expect_identical(names(model.frame(envO))[1L], "y")

# the record slots hold what they say
typed <- dbartsData(y ~ x1 + x2, d, subset = 2:10, keep.data = FALSE)
expect_error(typed@rows <- "junk", "rows")
expect_error(typed@rows.total <- 3.5, "rows.total")
expect_error(typed@frame.columns <- 3, "frame.columns")
expect_error(typed@build <- 7, "build")
expect_error(typed@basis.written <- list(1), "basis.written")
expect_error(typed@frame.source <- 1, "frame.source")

# a block-holding object decodes to the model's own shape where it is exact
decoded <- function(o) {
  mf <- model.frame(o)
  expect_true(is.factor(mf$y))
  expect_true(is.factor(mf$g))
  expect_identical(levels(mf$y), c("no", "yes"))
  expect_identical(levels(mf$g), c("p", "q", "r"))
}
dd <- data.frame(
  y = factor(sample(c("no", "yes"), n, TRUE)),
  x = rnorm(n),
  g = factor(sample(c("p", "q", "r"), n, TRUE))
)
decoded(dbartsData(y ~ x + g, dd, keep.data = FALSE))
decoded(dbarts(y ~ x + g, dd, n.trees = 5L, control = ctl)$data)
xm <- dd[c("x", "g")]
decodedX <- model.frame(dbartsData(xm, dd$y))
expect_true(is.factor(decodedX$y) && is.factor(decodedX$g))
e <- list2env(as.list(dd))
decoded(dbartsData(y ~ x + g, e))

# a change after the build leaves the object unchanged -------------------------

dc <- d
oc <- dbartsData(y ~ x1 + x2, dc)
before <- draws(dbarts(oc, n.trees = 5L, control = ctl, seed = 7L))
dc$x2[1L] <- 9
expect_identical(oc@frame$x2[1L], d$x2[1L])
expect_identical(
  draws(dbarts(oc, n.trees = 5L, control = ctl, seed = 7L)),
  before
)
ee <- new.env()
ee$y <- d$y
ee$x1 <- d$x1
oe <- dbartsData(y ~ x1, ee)
ee$x1[1L] <- 9
expect_identical(nrow(oe@x), n)
expect_false(oe@x[1L, 1L] == 9)
ee$later <- 1
ee$.hidden <- 2
expect_identical(
  dbartsData(y ~ x1, ee)@frame.columns,
  c(".hidden", "later", "x1", "y")
)
expect_null(dbartsData(y ~ x1, ee)@frame)

# a build that draws from R's generator holds the drawn values ----------------

moved <- function(expr) {
  set.seed(9L)
  before <- .Random.seed
  force(expr)
  !identical(before, .Random.seed)
}
expect_false(moved(dbartsData(y ~ x1, d)))
jittered <- dbartsData(y ~ jitter(x2), d)
expect_null(jittered@frame)
expect_true(is.numeric(jittered@y))
expect_identical(
  model.frame(
    dbarts(jittered, n.trees = 5L, control = ctl)$data
  )[[2L]],
  model.frame(dbarts(jittered, n.trees = 5L, control = ctl)$data)[[2L]]
)
drawnResponse <- dbartsData(rnorm(n) ~ x1, d)
expect_null(drawnResponse@frame)
expect_identical(
  dbarts(drawnResponse, n.trees = 5L, control = ctl)$data@y,
  drawnResponse@y
)
expect_null(dbartsData(y ~ x1, d, weights = runif(n))@frame)
# the formula and the object built from it draw the same values under a seed
set.seed(4L)
fromFormula <- dbarts(y ~ jitter(x2), d, n.trees = 5L, control = ctl)
set.seed(4L)
fromObject <- dbarts(dbartsData(y ~ jitter(x2), d), n.trees = 5L, control = ctl)
expect_identical(fromObject$data@x, fromFormula$data@x)
# 'data' and 'subset' are kept, so drawing them keeps the frame
drawnRows <- dbartsData(y ~ x1, d, subset = sample(n, 20L))
expect_false(is.null(drawnRows@frame))
expect_identical(nobs(drawnRows), 20L)
drawnData <- dbartsData(y ~ x1, d[sample(n), ])
expect_false(is.null(drawnData@frame))

# makeTestModelMatrix ----------------------------------------------------------

newdata <- d[1:6, c("x1", "x2", "f")]
newdata$x1[1L] <- 0.1
expect_identical(
  makeTestModelMatrix(o, newdata),
  makeTestModelMatrix(b, newdata)
)

# the missing-value flags ride the record: a column complete in training that
# the new rows leave missing is refused as it is beside the block
newNA <- newdata
newNA$f[2L] <- NA
messageOf <- function(x) {
  tryCatch(makeTestModelMatrix(x, newNA), error = conditionMessage)
}
expect_true(is.character(messageOf(b)))
expect_identical(messageOf(o), messageOf(b))
newOK <- newdata
newOK$x1[2L] <- NA
expect_identical(makeTestModelMatrix(o, newOK), makeTestModelMatrix(b, newOK))

# a caller's value in the formula is held as it was at the build ---------------

v <- rnorm(n)
ref <- dbartsData(y ~ x1 + v, d, subset = 2:30, keep.data = FALSE)
ov <- dbartsData(y ~ x1 + v, d, subset = 2:30)
v <- rnorm(n)
expect_identical(dbarts:::resolveDataObject(ov)@x, ref@x)
expect_identical(nrow(model.frame(ov)), 29L)

# a package function held under a name only the formula's environment knows
ns <- splines::ns
sp <- dbartsData(y ~ ns(w, df = 3), d)
spRef <- dbartsData(y ~ ns(w, df = 3), d, keep.data = FALSE)
rm(ns)
expect_identical(dbarts:::resolveDataObject(sp)@x, spRef@x)

# a function from no package is held by value, with the names it reads: a
# redefinition or an rm after the build changes nothing
wt0 <- runif(n) + 0.5
sq <- function(v) v^power
power <- 2
oSq <- dbartsData(y ~ x1, d, weights = sq(wt0))
built <- model.weights(model.frame(oSq))
expect_identical(built, wt0^2)
sq <- function(v) v^3
power <- 4
expect_identical(model.weights(model.frame(oSq)), built)
rm(sq, power)
expect_identical(model.weights(model.frame(oSq)), built)
expect_identical(
  fitDraws(oSq),
  fitDraws(dbartsData(y ~ x1, d, weights = wt0^2, keep.data = FALSE))
)
# one defined at the top level, in the response
assign("dbartsHeldSq", function(v) v^2, envir = globalenv())
oTop <- dbartsData(dbartsHeldSq(y) ~ x1, d)
assign("dbartsHeldSq", function(v) v^3, envir = globalenv())
expect_identical(unname(model.response(model.frame(oTop))), d$y^2)
rm("dbartsHeldSq", envir = globalenv())
expect_identical(unname(model.response(model.frame(oTop))), d$y^2)
expect_true(is.numeric(dbarts(oTop, n.trees = 5L, control = ctl)$data@y))
# a function local to a caller's function keeps none of that function's frame
localBuild <- function() {
  bulk <- numeric(5e6)
  lsq <- function(v) log1p(v^2)
  dbartsData(lsq(y) ~ x1, d)
}
oLocal <- localBuild()
expect_true(length(serialize(oLocal, NULL)) < 1e6)
expect_identical(unname(model.response(model.frame(oLocal))), log1p(d$y^2))

# a held function holds the names it looks up as variables: not the field after
# $ or @, nor either side of :: or :::, and a name it assigns after reading it
# still holds the workspace's value
local({
  heldNames <- function(object, fn) {
    held <- get(fn, environment(object@build$formula))
    sort(ls(environment(held), all.names = TRUE))
  }
  df <- data.frame(unrelated = 1:3)
  fit <- list(df = 2)
  fieldFn <- function(v) fit$df + v
  x <- 1
  assignFn <- function(v) {
    x <- x + v
    x
  }
  qualifiedFn <- function(v) v + stats::sd(v) + base:::length(v)
  expect_identical(heldNames(dbartsData(fieldFn(y) ~ x1, d), "fieldFn"), "fit")
  expect_identical(heldNames(dbartsData(assignFn(y) ~ x1, d), "assignFn"), "x")
  expect_identical(
    heldNames(dbartsData(qualifiedFn(y) ~ x1, d), "qualifiedFn"),
    character(0L)
  )
})

# the sibling identity: a fit from the object is the fit from the formula ------

sibling <- function(label, call, family = "auto", fit = TRUE) {
  where <- parent.frame()
  build <- function(fn, ...) {
    out <- call
    out[[1L]] <- as.name(fn)
    extra <- list(...)
    for (name in names(extra)) {
      out[[name]] <- extra[[name]]
    }
    eval(out, where)
  }
  o <- build("dbartsData")
  k <- build("dbartsData", keep.data = FALSE)
  r <- dbarts:::resolveDataObject(o)
  expect_identical(r@x, k@x, info = label)
  expect_identical(r@y, k@y, info = label)
  expect_identical(r@weights, k@weights, info = label)
  expect_identical(r@offset, k@offset, info = label)
  expect_identical(r@na.action, k@na.action, info = label)
  expect_identical(r@rowNames, k@rowNames, info = label)
  expect_identical(r@rows, k@rows, info = label)
  for (name in c("survivalStatus", "survivalTime", "survivalTimeOmitted")) {
    expect_identical(attr(r, name), attr(k, name), info = label)
  }
  if (!fit) {
    return(invisible(NULL))
  }
  fromObject <- fitDraws(o, family = family)
  fromFormula <- {
    out <- call
    out[[1L]] <- as.name("dbarts")
    out$n.trees <- 5L
    out$control <- ctl
    out$seed <- 7L
    out$family <- family
    draws(eval(out, where))
  }
  expect_identical(fromObject, fromFormula, info = label)
}
sibling("plain", quote(dbartsData(y ~ x1 + x2 + f, d)))
sibling("subset", quote(dbartsData(y ~ x1 + x2 + f, d, subset = 2:30)))
sibling(
  "repeat",
  quote(dbartsData(y ~ x1 + x2, d, subset = c(2, 2, 5, 7, 8, 10, 11, 12, 13)))
)
sibling("logical subset", quote(dbartsData(y ~ x1 + x2, d, subset = f != "b")))
sibling("na.omit", quote(dbartsData(y ~ x1 + x2, d, na.action = na.omit)))
sibling(
  "na.omit subset",
  quote(dbartsData(y ~ x1 + x2, d, subset = 2:30, na.action = na.omit))
)
sibling("weights", quote(dbartsData(y ~ x1 + x2, d, weights = w)))
sibling("offset", quote(dbartsData(y ~ x1 + x2, d, offset = offv)))
sibling("offset term", quote(dbartsData(y ~ x1 + x2 + offset(w), d)))
sibling(
  "poly and scale",
  quote(dbartsData(y ~ poly(w, 2) + scale(x2), d, subset = 2:30))
)
sibling("indicators", quote(dbartsData(y ~ x1 + f, d, factors = "indicators")))
sibling("callers vector", quote(dbartsData(y ~ x1 + v, d, subset = 2:30)))
sibling(
  "callers vector na.omit",
  quote(dbartsData(y ~ x1 + v, d, na.action = na.omit))
)
sibling(
  "callers weights",
  quote(dbartsData(y ~ x1, d, weights = d$w + 1, subset = 2:30))
)
dp <- d
dp$yb <- as.numeric(dp$y > 0)
sibling(
  "probit numeric",
  quote(dbartsData(yb ~ x1 + x2, dp)),
  family = "probit"
)
dp$yf <- factor(ifelse(dp$y > 0, "up", "down"))
sibling("probit factor", quote(dbartsData(yf ~ x1 + x2, dp)), family = "probit")
dp$yl <- dp$y > 0
sibling(
  "probit logical",
  quote(dbartsData(yl ~ x1 + x2, dp)),
  family = "probit"
)
dp$yo <- cut(
  dp$y,
  quantile(dp$y, c(0, 1 / 3, 2 / 3, 1)),
  include.lowest = TRUE,
  ordered_result = TRUE
)
sibling(
  "ordinal factor",
  quote(dbartsData(yo ~ x1 + x2, dp)),
  family = "ordinal"
)
dp$ym <- factor(sample(c("p", "q", "r"), n, TRUE))
sibling("multinomial factor", quote(dbartsData(ym ~ x1 + x2, dp)), fit = FALSE)
ds <- d
ds$time <- rexp(n)
ds$status <- rbinom(n, 1L, 0.7)
if (requireNamespace("survival", quietly = TRUE)) {
  sibling(
    "aft",
    quote(dbartsData(survival::Surv(time, status) ~ x1 + x2, ds)),
    family = "aft"
  )
  sibling(
    "aft na.omit",
    quote(dbartsData(
      survival::Surv(time, status) ~ x1 + x2,
      ds,
      na.action = na.omit
    )),
    family = "aft"
  )
  sibling(
    "hazard",
    quote(dbartsData(survival::Surv(ceiling(time * 3), status) ~ x1 + x2, ds)),
    family = "hazard"
  )
}

# the other doors take the object and fit as their formula sibling ------------

ob <- dbartsData(y ~ x1 + x2, d, subset = 2:30, weights = w)
obv <- dbartsData(y ~ x1 + x2, d, subset = 2:30, weights = d$w)
expect_identical(
  bart(
    ob,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 7L,
    verbose = FALSE,
    keepTrees = FALSE
  )$yhat.train,
  bart(
    y ~ x1 + x2,
    d,
    subset = 2:30,
    weights = w,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 7L,
    verbose = FALSE,
    keepTrees = FALSE
  )$yhat.train
)

sampler <- dbarts(y ~ x1 + x2, d, subset = 2:30, n.trees = 5L, control = ctl)
expect_silent(sampler$setData(ob))
expect_null(sampler$data@frame)
expect_true(is.numeric(sampler$data@y))
expect_identical(sampler$data@frame.columns, names(d))
spec <- dbartsSpec(ob)
expect_null(spec$data@frame)
expect_true(is.numeric(spec$data@x[, 1L]))
ob100 <- dbartsData(y ~ x1 + x2, d, subset = 2:30, n.cuts = 100L)
ob100@sigma <- spec$data@sigma
made <- new("dbartsSampler", control = ctl, model = spec$model, data = ob100)
expect_null(made$data@frame)
sFit <- dbarts(ob, n.trees = 5L, control = ctl)
expect_null(sFit$data@frame)
expect_identical(nobs(sFit$data), 29L)
xbArgs <- list(
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 2L,
  n.reps = 2L,
  n.test = 2L,
  n.threads = 1L,
  seed = 7L
)
expect_identical(
  do.call(xbart, c(list(obv), xbArgs)),
  do.call(xbart, c(list(y ~ x1 + x2, d, subset = 2:30, weights = d$w), xbArgs))
)

# a build's warnings are given once, at the build -------------------------------

caught <- character()
collect <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      caught <<- c(caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
}
dz$w <- pmin(dz$w, 1)
collect(oz <- dbartsData(y ~ x1 + x2, dz, weights = w))
expect_identical(length(caught), 1L)
caught <- character()
collect({
  dbarts(oz, n.trees = 5L, control = ctl)
  dbartsSpec(oz)
  sampler$setData(oz)
  makeTestModelMatrix(oz, d[1:3, c("x1", "x2")])
  nobs(oz)
  model.frame(oz)
})
expect_identical(caught, character())

# saved size ---------------------------------------------------------------------

nBig <- 4000L
big <- data.frame(y = rnorm(nBig), matrix(rnorm(nBig * 8L), nBig))
sizeOf <- function(x) length(serialize(x, NULL))
kept <- dbartsData(y ~ ., big)
block <- dbartsData(y ~ ., big, keep.data = FALSE)
expect_true(sizeOf(kept) < sizeOf(block) + sizeOf(big))
expect_true(sizeOf(kept@frame) < sizeOf(block@x) + sizeOf(block@y) + 1000)
bs <- dbarts(kept, n.trees = 5L, control = ctl)
expect_null(bs$data@frame)
expect_true(
  abs(
    sizeOf(bs$data) - sizeOf(dbarts(block, n.trees = 5L, control = ctl)$data)
  ) <
    1000
)

# rbart_vi takes the object and fits as its formula sibling ----------------------

grp <- factor(rep(1:4, length.out = n))
rbartArgs <- list(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 12L,
  n.burn = 2L,
  n.thin = 1L,
  verbose = FALSE,
  keepTrees = FALSE,
  seed = 3L
)
rbartDraws <- function(first, ...) {
  withCallingHandlers(
    do.call(rbart_vi, c(list(first, ...), rbartArgs))$yhat.train,
    warning = function(w) invokeRestart("muffleWarning")
  )
}
dr <- d[c("y", "x2")]
dr$x2[is.na(dr$x2)] <- 0.5
expect_identical(
  rbartDraws(dbartsData(y ~ x2, dr), group.by = grp),
  rbartDraws(y ~ x2, dr, group.by = grp)
)

# the na.action record rides model.frame() --------------------------------------

expect_identical(
  attr(model.frame(oOff), "na.action"),
  attr(model.frame(fit), "na.action")
)

# the widened slots still hold their types --------------------------------------

bw <- dbartsData(y ~ x1 + x2, d, weights = w, offset = offv, keep.data = FALSE)
expect_true(validObject(bw))
tooBad <- bw
tooBad@y <- letters[seq_len(n)]
expect_error(validObject(tooBad), "'y' must be numeric")
tooBad <- bw
tooBad@weights <- letters[seq_len(n)]
expect_error(validObject(tooBad), "'weights' must be null or numeric")
tooBad <- bw
tooBad@offset <- letters[seq_len(n)]
expect_error(validObject(tooBad), "'offset' must be null or numeric")
tooBad <- bw
tooBad@x <- letters[seq_len(n)]
expect_error(validObject(tooBad), "'x' must be a matrix")
tooBad <- o
tooBad@y <- letters[seq_len(n)]
expect_error(validObject(tooBad), "'y' must be a stand-in")

# the sampler's data are read-only: a write stops and names $setData ----------

fieldSampler <- dbarts(
  y ~ x1 + x2,
  d,
  subset = 2:30,
  n.trees = 5L,
  control = ctl
)
heldBefore <- fieldSampler$data
readOnly <- "the sampler's data are read-only here; change them with s\\$setData"
expect_error(fieldSampler$data <- ob100, readOnly)
expect_error(fieldSampler$data <- 3, readOnly)
expect_error(
  fieldSampler$data <- dbartsData(y ~ x1 + x2, d, keep.data = FALSE),
  readOnly
)
expect_identical(fieldSampler$data, heldBefore)
# the sampler goes on as it was, and copies
expect_true(is.numeric(fieldSampler$run(0L, 2L)$train))
expect_identical(fieldSampler$copy()$data@y, heldBefore@y)

# a stand-in does not answer for values it does not hold ------------------------

stops(is.na(o@y))
stops(anyNA(o@y))
stops(is.finite(o@y))
stops(is.infinite(o@y))
stops(is.nan(o@y))
stops(mean(o@y))
stops(c(o@y))
stops(rep(o@y, 2L))
stops(unique(o@y))
stops(as.numeric(o@y))
stops(as.double(o@weights))
stops(cbind(o@y, 1))
stops(cbind(1, o@y))
stops(cbind(o@x, o@y))
stops(rbind(o@y, 1))

# a sparse column: model.frame() keeps it as stored, never dense --------------

nsp <- 12L
spf <- factor(rep(c("a", "b", "c"), length.out = nsp))
dsp <- data.frame(y = d$y[seq_len(nsp)], w = d$w[seq_len(nsp)])
dsp$s <- Matrix::sparseVector(c(1, 2), c(3L, 7L), length = nsp)
dsp$sf <- sparseFactor(spf)
spM <- Matrix::sparseMatrix(
  i = c(1L, 4L, 6L, 2L, 9L),
  j = c(1L, 1L, 1L, 2L, 2L),
  x = c(0.5, -1, 2, 3, 1.5),
  dims = c(nsp, 2L),
  dimnames = list(NULL, c("a", "b"))
)
dsp$M <- spM
spFormula <- y ~ poly(w, 2) + s + sf + M
osp <- dbartsData(spFormula, dsp)
expect_false(is.null(osp@frame))
expect_identical(nobs(osp), nsp)
mfsp <- model.frame(osp)
expect_identical(names(mfsp), c("y", "poly(w, 2)", "s", "sf", "M"))
expect_true(methods::is(mfsp$s, "sparseVector"))
expect_true(methods::is(mfsp$sf, "sparseFactor"))
expect_true(methods::is(mfsp$M, "dgCMatrix"))
expect_identical(mfsp$s, dsp$s)
expect_identical(mfsp$sf, dsp$sf)
expect_identical(mfsp$M, spM)
expect_true(inherits(mfsp[["poly(w, 2)"]], "poly"))
expect_false(is.null(attr(mfsp[["poly(w, 2)"]], "coefs")))
expect_identical(unname(model.response(mfsp)), dsp$y)
expect_identical(attr(attr(mfsp, "terms"), "response"), 1L)
# over the rows kept
ospSub <- dbartsData(spFormula, dsp, subset = c(2L, 5L, 5L, 9L))
mfSub <- model.frame(ospSub)
expect_identical(nrow(mfSub), 4L)
expect_identical(as.vector(mfSub$s), as.vector(dsp$s)[c(2L, 5L, 5L, 9L)])
expect_identical(
  as.character(as.factor(mfSub$sf)),
  as.character(spf[c(2L, 5L, 5L, 9L)])
)
expect_identical(mfSub$sf@levels, levels(spf))
expect_identical(as.matrix(mfSub$M), as.matrix(spM)[c(2L, 5L, 5L, 9L), ])

# a sampler's data decode a sparse column sparse too
spSampler <- dbarts(y ~ w + s + sf + M, dsp, n.trees = 5L, control = ctl)
mfSampler <- model.frame(spSampler$data)
expect_identical(names(mfSampler), c("y", "w", "s", "sf", "M"))
expect_true(methods::is(mfSampler$s, "sparseVector"))
expect_identical(as.vector(mfSampler$s), as.vector(dsp$s))
expect_true(methods::is(mfSampler$sf, "sparseFactor"))
expect_identical(as.factor(mfSampler$sf), spf)
expect_true(methods::is(mfSampler$M, "dgCMatrix"))
expect_identical(as.matrix(mfSampler$M), as.matrix(spM))
expect_identical(unname(model.response(mfSampler)), dsp$y)
# the terms record lists each sparse column, so a fitter that cannot read one
# stops at it, as it does on the kept-frame route, and fits no dense part alone
spTerms <- terms(mfSampler)
expect_identical(attr(spTerms, "term.labels"), c("w", "s", "sf", "M"))
expect_identical(attr(spTerms, "response"), 1L)
expect_identical(attr(spTerms, "intercept"), 1L)
expect_identical(
  length(attr(spTerms, "predvars")),
  length(attr(spTerms, "variables"))
)
expect_error(lm(mfSampler), "not a matrix")
expect_error(lm(mfsp), "not a matrix")
# a record with a data-dependent basis keeps its predvars beside the sparse term
spPoly <- dbarts(y ~ poly(w, 2) + s, dsp, n.trees = 5L, control = ctl)
mfPoly <- model.frame(spPoly$data)
expect_identical(
  attr(terms(mfPoly), "term.labels"),
  c("poly(w, 2)", "s")
)
expect_true(any(grepl("coefs", deparse(attr(terms(mfPoly), "predvars")))))
expect_identical(unname(model.response(mfPoly)), dsp$y)
rm(spTerms, spPoly, mfPoly)

# a block-holding object's frame is lm's: the response named as written, a
# matrix term one column under its own name, the formula's order and a terms
# record ------------------------------------------------------------------------

d0 <- d[c("y", "w", "x2")]
d0$x2[is.na(d0$x2)] <- 0.5
d0$M <- matrix(rnorm(2L * n), n, 2L, dimnames = list(NULL, c("a", "b")))
blockFormula <- log(abs(y) + 1) ~ x2 + poly(w, 2) + scale(x2) + M
bd <- dbartsData(blockFormula, d0, keep.data = FALSE)
mfb <- model.frame(bd)
mfRef <- model.frame(blockFormula, d0)
expect_identical(names(mfb), names(mfRef))
expect_identical(
  names(mfb),
  c("log(abs(y) + 1)", "x2", "poly(w, 2)", "scale(x2)", "M")
)
expect_true(is.matrix(mfb[["poly(w, 2)"]]))
expect_identical(colnames(mfb[["poly(w, 2)"]]), c("1", "2"))
expect_true(is.matrix(mfb[["scale(x2)"]]))
expect_identical(dim(mfb[["scale(x2)"]]), c(n, 1L))
# a one-column matrix is named as lm's frame names it: nothing for scale() or
# an unnamed matrix, and the name the user gave otherwise
expect_null(colnames(mfb[["scale(x2)"]]))
expect_identical(colnames(mfRef[["scale(x2)"]]), colnames(mfb[["scale(x2)"]]))
d0$U <- matrix(d0$x2, n, 1L)
d0$Q <- cbind(q = d0$x2)
oneFormula <- y ~ U + Q + M
mfOne <- model.frame(dbartsData(oneFormula, d0, keep.data = FALSE))
mfOneRef <- model.frame(oneFormula, d0)
for (name in c("U", "Q", "M")) {
  expect_identical(
    colnames(mfOne[[name]]),
    colnames(mfOneRef[[name]]),
    info = name
  )
}
expect_identical(colnames(mfOne$Q), "q")
rm(mfOne, mfOneRef, oneFormula)
d0$U <- NULL
d0$Q <- NULL
expect_identical(colnames(mfb$M), c("a", "b"))
for (name in names(mfRef)[-1L]) {
  expect_equal(
    unclass(as.matrix(mfb[[name]])),
    unclass(as.matrix(mfRef[[name]])),
    check.attributes = FALSE,
    info = name
  )
}
expect_equal(mfb[[1L]], mfRef[[1L]])
expect_equal(model.response(mfb), model.response(mfRef))
bdTerms <- attr(mfb, "terms")
expect_identical(attr(bdTerms, "response"), 1L)
expect_identical(
  attr(bdTerms, "term.labels"),
  attr(attr(mfRef, "terms"), "term.labels")
)
expect_true(any(grepl("coefs", deparse(attr(bdTerms, "predvars")))))
expect_true(any(grepl("center", deparse(attr(bdTerms, "predvars")))))
expect_equal(
  unname(coef(lm(mfb))),
  unname(coef(lm(blockFormula, d0)))
)
fittedFrame <- model.frame(
  dbarts(blockFormula, d0, n.trees = 5L, control = ctl)$data
)
expect_identical(names(fittedFrame), names(mfRef))
expect_equal(model.response(fittedFrame), model.response(mfRef))
# an offset() term is counted once, as lm's frame counts it, with and without
# offset =; the block holds the sum
offsetFormula <- y ~ x2 + offset(w)
for (extra in list(NULL, rep(100, n))) {
  lmOffset <- model.offset(model.frame(lm(offsetFormula, d0, offset = extra)))
  offsetSampler <- dbarts(
    offsetFormula,
    d0,
    offset = extra,
    n.trees = 5L,
    control = ctl
  )
  expect_equal(
    unname(model.offset(model.frame(offsetSampler$data))),
    unname(lmOffset),
    info = paste("offset =", !is.null(extra))
  )
  expect_equal(
    unname(model.offset(model.frame(
      dbartsData(offsetFormula, d0, offset = extra, keep.data = FALSE)
    ))),
    unname(lmOffset)
  )
}
# where no formula was, or a hazard fit added its period, the frame has no
# terms record
expect_null(attr(
  model.frame(dbartsData(as.matrix(d0[c("w", "x2")]), d0$y)),
  "terms"
))
hazardFrame <- model.frame(
  dbarts(
    survival::Surv(ceiling(time * 3), status) ~ x1 + x2,
    ds,
    family = "hazard",
    n.trees = 5L,
    control = ctl
  )$data
)
expect_identical(names(hazardFrame)[4L], "period")
expect_null(attr(hazardFrame, "terms"))
