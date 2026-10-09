# sigest beside a fixed residual scale: accepted with a message where it
# equals the square root of the fixed variance, refused where it differs,
# the same at every door that reaches the check

library(dbarts)
set.seed(7)
n <- 40L
x <- matrix(runif(n * 2L), n)
y <- x[, 1L] + rnorm(n, 0, 0.3)
ctl <- dbartsControl(
  n.trees = 3L,
  n.samples = 3L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  updateState = FALSE
)
msgKey <- "tombstone.sigest.fixed"
onceState <- dbarts:::onceWarnState
resetMsg <- function() onceState[[msgKey]] <- NULL
agreeMsg <- "has no effect under a fixed residual scale .*error from dbarts 1\\.1-0"
differ <- "no effect under a fixed residual scale: sigma = fixed"
gauss <- dbartsFamilies$gaussian
fixedPrior <- dbartsPriors$fixed
fx <- gauss(sigma = fixedPrior(4))
bartQuick <- function(...) {
  bart(
    x,
    y,
    n.trees = 3L,
    n.samples = 3L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    ...
  )
}

# rbart_vi and bartBT take no family and so cannot reach a fixed residual
# scale: their prior is always chisq

# the loops' exact call
resetMsg()
expect_message(
  s <- dbarts(x, y, resid.prior = fixedPrior(1), sigma = 1, control = ctl),
  agreeMsg
)
expect_true(is(s$model@resid.prior, "dbartsFixedPrior"))
# the message is shown once per session
expect_silent(dbarts(x, y, family = fx, sigest = 2, control = ctl))

# agreeing: a message, no error, at each door
resetMsg()
expect_message(dbarts(x, y, family = fx, sigest = 2, control = ctl), agreeMsg)
resetMsg()
expect_message(
  suppressWarnings(dbarts(x, y, family = fx, sigma = 2, control = ctl)),
  agreeMsg
)
resetMsg()
expect_message(bartQuick(family = fx, sigest = 2), agreeMsg)
resetMsg()
expect_message(
  dbartsSpec(dbartsData(x, y), family = fx, sigest = 2, control = ctl),
  agreeMsg
)
resetMsg()
expect_message(
  xbart(
    x,
    y,
    family = fx,
    sigest = 2,
    n.trees = 3L,
    n.samples = 3L,
    n.burn = 0L,
    n.reps = 1L,
    n.test = 2L,
    n.threads = 1L,
    verbose = FALSE
  ),
  agreeMsg
)
resetMsg()
expect_message(
  suppressWarnings(xbart(
    x,
    y,
    family = fx,
    sigma = 2,
    n.trees = 3L,
    n.samples = 3L,
    n.burn = 0L,
    n.reps = 1L,
    n.test = 2L,
    n.threads = 1L,
    verbose = FALSE
  )),
  agreeMsg
)

# a few ulps off the square root still agrees; a different scale does not
resetMsg()
expect_message(
  dbarts(
    x,
    y,
    family = gauss(sigma = fixedPrior(0.09)),
    sigest = 0.3,
    control = ctl
  ),
  agreeMsg
)
expect_error(
  dbarts(x, y, family = fx, sigest = 2 + 1e-6, control = ctl),
  differ
)

# differing: refused at each door
expect_error(dbarts(x, y, family = fx, sigest = 1.5, control = ctl), differ)
expect_error(
  suppressWarnings(dbarts(x, y, family = fx, sigma = 1.5, control = ctl)),
  "'sigma' has no effect"
)
expect_error(bartQuick(family = fx, sigest = 1.5), differ)
expect_error(
  dbartsSpec(dbartsData(x, y), family = fx, sigest = 1.5, control = ctl),
  differ
)
expect_error(
  xbart(
    x,
    y,
    family = fx,
    sigest = 1.5,
    n.trees = 3L,
    n.samples = 3L,
    n.burn = 0L,
    n.reps = 1L,
    n.test = 2L,
    n.threads = 1L,
    verbose = FALSE
  ),
  differ
)
expect_error(
  suppressWarnings(xbart(
    x,
    y,
    family = fx,
    sigma = 1.5,
    n.trees = 3L,
    n.samples = 3L,
    n.burn = 0L,
    n.reps = 1L,
    n.test = 2L,
    n.threads = 1L,
    verbose = FALSE
  )),
  "'sigma' has no effect"
)

# registered for the 1.1-0 error
reg <- dbarts:::dbartsTombstones
expect_true(
  "sigest beside a fixed residual scale" %in% vapply(reg, `[[`, "", "name")
)
