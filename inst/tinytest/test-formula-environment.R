# A fit stores its formula's terms and its 'offset' expression with every name
# that is no data column frozen at its fit-time value and with no environment
# of its own, so it carries neither its caller's frame nor the global
# environment, and predicts from newdata alone, before and after a reload. A
# formula calling a function from no package is refused.

fitArgs <- list(
  n.trees = 10L,
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 5L,
  keepTrees = TRUE
)

set.seed(31L)
n <- 50L
d <- data.frame(x = runif(n), o = rnorm(n))
d$y <- sin(3 * d$x) + d$o + rnorm(n, 0, 0.1)
nd <- data.frame(x = c(0.2, 0.5, 0.8), o = c(1, 0, -1))

serializedSize <- function(object) length(serialize(object, NULL))

# ---- a fit made inside a function carries none of that function's frame

topLevel <- do.call(bart, c(list(y ~ x + offset(o), d), fitArgs))
insideFunction <- function(data) {
  big <- numeric(5e6) # 40 MB that must not ride along
  K <- 1
  fit <- do.call(bart, c(list(y ~ x + offset(o * K), data), fitArgs))
  invisible(big)
  fit
}
wrapper <- function(data, ...) insideFunction(data, ...)
nested <- wrapper(d)
expect_true(abs(serializedSize(nested) - serializedSize(topLevel)) < 1e5)
expect_identical(environment(attr(nested$fit$data@x, "terms")), baseenv())

# the offset argument's expression is stored the same way
argumentInside <- function(data) {
  big <- numeric(5e6)
  fit <- bart(
    y ~ x,
    data,
    offset = o / 2,
    n.trees = 10L,
    n.samples = 10L,
    n.burn = 10L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    seed = 5L,
    keepTrees = TRUE
  )
  invisible(big)
  fit
}
fit.argument <- argumentInside(d)
expect_true(serializedSize(fit.argument) < serializedSize(topLevel) + 1e5)
expect_identical(
  environment(attr(fit.argument$fit$data, "offset.argument")),
  baseenv()
)
expect_equal(
  unname(apply(
    predict(fit.argument, nd[c(2L, 2L), ] + c(0, 0, 0, 2)),
    1L,
    diff
  )),
  rep(1, 10L)
)

# ---- an offset expression mixing a data column and a local is evaluated as
# model.frame does - the column from data, the local from the calling scope -
# at fit, for the default test offset, and at predict

mixedOffset <- function(data, test) {
  k <- 2
  bart(
    y ~ x,
    data,
    test = test,
    offset = o * k,
    n.trees = 10L,
    n.samples = 10L,
    n.burn = 10L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    seed = 5L,
    keepTrees = TRUE
  )
}
fit.mixed <- mixedOffset(d, nd)
expect_identical(fit.mixed$fit$data@offset, d$o * 2)
expect_identical(fit.mixed$fit$data@offset.test, nd$o * 2)
expect_equal(unname(predict(fit.mixed, nd)), unname(fit.mixed$yhat.test))
fit.literal <- bart(
  y ~ x,
  d,
  offset = o * 2,
  n.trees = 10L,
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 5L,
  keepTrees = TRUE
)
expect_identical(fit.mixed$yhat.train, fit.literal$yhat.train)
# a name found in neither the data nor any calling scope is refused, not
# dropped
expect_error(
  bart(y ~ x, d, offset = o * not.anywhere, verbose = FALSE),
  pattern = "'offset' cannot be evaluated: object 'not.anywhere' not found"
)
rm(mixedOffset, fit.mixed, fit.literal)

# ---- a local constant in a term or offset predicts the same after a reload

localFit <- function(data) {
  K <- 2
  deg <- 2L
  do.call(bart, c(list(y ~ poly(x, deg) + offset(o * K), data), fitArgs))
}
fit.local <- localFit(d)
before <- predict(fit.local, nd)
expect_equal(
  unname(apply(predict(fit.local, nd[c(2L, 2L), ] + c(0, 0, 0, 1)), 1L, diff)),
  rep(2, 10L)
)
file.fit <- tempfile(fileext = ".rds")
fit.local$fit$storeState()
saveRDS(fit.local, file.fit)
expect_equal(predict(readRDS(file.fit), nd), before)
# the stored terms hold no environment, so the reloaded fit predicts the same
# with every local it was fit beside gone
expect_identical(
  environment(attr(fit.local$fit$data@x, "terms")),
  baseenv()
)
rm(localFit, fit.local)
expect_equal(predict(readRDS(file.fit), nd), before)
invisible(file.remove(file.fit))

# ---- package bases keep working, unqualified or qualified

bs <- splines::bs
ns <- splines::ns
for (rhs in c("ns(x, 3)", "bs(x, 3)", "scale(x)", "log(x + 1)", "I(x^2)")) {
  fit.basis <- do.call(
    bart,
    c(list(as.formula(paste("y ~", rhs, "+ offset(o)")), d), fitArgs)
  )
  expect_identical(
    environment(attr(fit.basis$fit$data@x, "terms")),
    baseenv()
  )
  expect_equal(
    unname(predict(fit.basis, d[1:3, ])),
    unname(fit.basis$yhat.train[, 1:3])
  )
}

# ---- a function from no package is refused at fit, naming it

myTransform <- function(v) v^2
expect_error(
  do.call(bart, c(list(y ~ myTransform(x), d), fitArgs)),
  pattern = "the formula calls 'myTransform', a function from no package"
)
expect_error(
  dbarts::xbart(y ~ myTransform(x), d, n.reps = 1L, verbose = FALSE),
  pattern = "the formula calls 'myTransform'"
)
# and in the offset expression, on either interface
expect_error(
  do.call(bart, c(list(y ~ x, d, offset = quote(myTransform(o))), fitArgs)),
  pattern = "the 'offset' expression calls 'myTransform', a function from no package"
)
expect_error(
  bart(
    d["x"],
    d$y,
    offset = myTransform(d$o),
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    verbose = FALSE
  ),
  pattern = "the 'offset' expression calls 'myTransform'"
)

# ---- xbart and bartBT inside a function carry no frame either

xbartInside <- function(data) {
  big <- numeric(5e6)
  K <- 1
  result <- dbarts::xbart(
    y ~ x + offset(o * K),
    data,
    n.reps = 1L,
    n.trees = 10L,
    n.samples = 10L,
    n.burn = c(10L, 5L),
    n.threads = 1L,
    verbose = FALSE
  )
  invisible(big)
  result
}
expect_true(serializedSize(xbartInside(d)) < 1e6)
bartBTInside <- function(data) {
  big <- numeric(5e6)
  fit <- dbarts::bartBT(
    data["x"],
    data$y,
    ntree = 10L,
    ndpost = 10L,
    nskip = 10L,
    verbose = FALSE,
    keeptrees = TRUE
  )
  invisible(big)
  fit
}
expect_true(serializedSize(bartBTInside(d)) < serializedSize(topLevel) + 1e6)

rm(
  fitArgs,
  n,
  d,
  nd,
  serializedSize,
  topLevel,
  insideFunction,
  wrapper,
  nested,
  argumentInside,
  fit.argument,
  before,
  file.fit,
  bs,
  ns,
  rhs,
  fit.basis,
  myTransform,
  xbartInside,
  bartBTInside
)
