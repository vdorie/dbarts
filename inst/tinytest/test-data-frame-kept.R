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
expect_identical(nobs(oZ), n)
expect_identical(nobs(lm(y ~ x1 + x2, dz, weights = w)), n - 4L)

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

# a field write resolves the object, as $setData does ---------------------------

fieldSampler <- dbarts(
  y ~ x1 + x2,
  d,
  subset = 2:30,
  n.trees = 5L,
  control = ctl
)
fieldSampler$data <- ob100
expect_null(fieldSampler$data@frame)
expect_true(is.numeric(fieldSampler$data@y))

# a stand-in does not answer for values it does not hold ------------------------

stops(is.na(o@y))
stops(anyNA(o@y))
stops(is.finite(o@y))
stops(is.infinite(o@y))
stops(is.nan(o@y))

# a sparse column: model.frame() answers as a sampler's data do -----------------

dsp <- d[c("y", "x2", "w")]
dsp$s <- Matrix::sparseVector(c(1, 2), c(3L, 7L), length = n)
osp <- dbartsData(y ~ x2 + s, dsp)
expect_identical(nobs(osp), n)
mfsp <- model.frame(osp)
expect_identical(nrow(mfsp), n)
expect_true("s" %in% names(mfsp))
expect_identical(
  mfsp$s,
  model.frame(dbartsData(y ~ x2 + s, dsp, keep.data = FALSE))$s
)

# a block-holding object's frame names the response as written and joins a
# matrix term back into one column --------------------------------------------

d0 <- d[c("y", "w", "x2")]
d0$x2[is.na(d0$x2)] <- 0.5
bd <- dbartsData(log(abs(y) + 1) ~ poly(w, 2) + x2, d0, keep.data = FALSE)
mfb <- model.frame(bd)
mfRef <- model.frame(log(abs(y) + 1) ~ poly(w, 2) + x2, d0)
expect_identical(names(mfb), names(mfRef))
expect_true(is.matrix(mfb[["poly(w, 2)"]]))
expect_identical(colnames(mfb[["poly(w, 2)"]]), c("1", "2"))
expect_equal(
  unclass(as.matrix(mfb[["poly(w, 2)"]])),
  unclass(as.matrix(mfRef[["poly(w, 2)"]])),
  check.attributes = FALSE
)
expect_equal(mfb[[1L]], mfRef[[1L]])
expect_identical(
  names(model.frame(
    dbarts(
      log(abs(y) + 1) ~ poly(w, 2) + x2,
      d0,
      n.trees = 5L,
      control = ctl
    )$data
  )),
  names(mfRef)
)
