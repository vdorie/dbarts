# A gp sampler that holds saved draws refuses updateScale = TRUE where its
# family has a response scale to re-derive (gaussian, Student-t, aft and nbinom),
# and takes it where it has none (probit, logistic and ordinal), in every setter
# that carries the argument: setResponse, setOffset, setWeights, setActiveRows
# and setData (dec-B481 reading dec-B330). A saved gp draw replays only under the
# transform it was drawn with; a scale-free family has no transform, so the call
# derives nothing and moves no saved draw.

gp <- dbartsPriors$gp
set.seed(3L)
n <- 80L
x <- cbind(a = runif(n), b = rnorm(n))
x.new <- cbind(a = runif(30L), b = rnorm(30L))
y <- x[, "b"] + sin(4 * x[, "a"]) + rnorm(n, 0, 0.3)
y.binary <- as.double(y > median(y))
y.ordinal <- as.integer(cut(
  y,
  quantile(y, c(0, 1 / 3, 2 / 3, 1)),
  include.lowest = TRUE
))
y.count <- as.double(rpois(n, exp(0.3 * x[, "b"] + 1)))
control.gp <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 8L,
  n.samples = 3L,
  n.burn = 5L,
  keepTrees = TRUE,
  updateState = FALSE,
  verbose = FALSE
)
samplerOf <- function(...) {
  sampler <- dbarts(x, ..., control = control.gp, leaf.prior = gp(columns = 2L))
  invisible(sampler$run())
  sampler
}
# the call each setter makes, on the replacement the family takes: a mask of
# 0/1 for the families whose weights are one, counts for the logistic, precision
# multipliers for the others
mask <- rep_len(c(1, 1, 0), n)
cases <- list(
  gaussian = list(
    make = function() samplerOf(y),
    response = y + 1,
    weights = rep_len(c(1, 1, 2), n)
  ),
  student = list(
    make = function() samplerOf(y, family = dbartsFamilies$student(df = 4)),
    response = y + 1,
    weights = rep_len(c(1, 1, 2), n)
  ),
  aft = list(
    make = function() samplerOf(cbind(exp(y), 1), family = "aft"),
    response = exp(y + 1),
    weights = rep_len(c(1, 1, 2), n)
  ),
  nbinom = list(
    make = function() samplerOf(y.count, family = "nbinom"),
    response = y.count + 1,
    weights = mask
  ),
  probit = list(
    make = function() samplerOf(y.binary, family = "probit"),
    response = 1 - y.binary,
    weights = mask
  ),
  logistic = list(
    make = function() samplerOf(y.binary, family = "logistic"),
    response = 1 - y.binary,
    weights = rep_len(c(1, 1, 2), n)
  ),
  ordinal = list(
    make = function() samplerOf(y.ordinal, family = "ordinal"),
    response = rev(y.ordinal),
    weights = mask
  )
)
refused <- c("gaussian", "student", "aft", "nbinom")
taken <- c("probit", "logistic", "ordinal")
refusal <- function(caller) {
  paste0(
    "^",
    caller,
    ": 'updateScale' cannot be TRUE for a sampler with gp leaves that holds ",
    "saved draws: re-derive the scale during burn-in, before draws are saved ",
    "\\(burn-in sweeps, or keepTrees off until the scale is fixed\\), or make ",
    "a new sampler$"
  )
}
calls <- list(
  setResponse = function(s, case) {
    s$setResponse(case$response, updateScale = TRUE)
  },
  setOffset = function(s, case) s$setOffset(rep(-5, n), updateScale = TRUE),
  setWeights = function(s, case) {
    s$setWeights(case$weights, updateScale = TRUE)
  },
  setActiveRows = function(s, case) s$setActiveRows(mask, updateScale = TRUE),
  setData = function(s, case) s$setData(s$data, updateScale = TRUE)
)
# what a call that is taken has installed
installed <- list(
  setResponse = function(s, case) {
    expect_equal(as.double(s$data@y), as.double(case$response))
  },
  setOffset = function(s, case) expect_equal(s$data@offset, rep(-5, n)),
  setWeights = function(s, case) {
    if (identical(case$weights, mask)) {
      expect_null(s$data@weights)
    } else {
      expect_equal(s$data@weights, case$weights)
    }
  },
  setActiveRows = function(s, case) NULL,
  setData = function(s, case) NULL
)
# everything a refused call could have touched
held <- function(s) {
  s$storeState()
  list(
    state = s$state,
    data = s$data,
    model = s$model,
    kept = s$predict(x.new)
  )
}

for (kind in c(refused, taken)) {
  case <- cases[[kind]]
  for (name in names(calls)) {
    info <- paste(kind, name)
    # aft fixes its censoring structure at creation, which setData refuses on
    # its own ground, before the gp rule
    if (kind == "aft" && name == "setData") {
      next
    }
    sampler <- case$make()
    before <- held(sampler)
    if (kind %in% refused) {
      expect_error(
        calls[[name]](sampler, case),
        if (name == "setData") {
          paste0(
            "$setData: 'newData' cannot change the covariate standardization ",
            "or the response range of a sampler with gp leaves that holds ",
            "saved draws"
          )
        } else {
          # a mask of 0/1 weights is installed as the active rows, which is the
          # call that refuses
          refusal(paste0(
            "\\$",
            if (name == "setWeights" && identical(case$weights, mask)) {
              "setActiveRows"
            } else {
              name
            }
          ))
        },
        fixed = name == "setData",
        info = info
      )
      expect_identical(held(sampler), before, info = info)
    } else {
      # taken: no error, no warning, and the saved draws predict what they did
      expect_silent(calls[[name]](sampler, case), info = info)
      installed[[name]](sampler, case)
      expect_identical(sampler$predict(x.new), before$kept, info = info)
    }
    # and the sampler runs on
    expect_true(all(is.finite(sampler$run(0L, 2L)$train)), info = info)
  }
}

# updateScale left FALSE is taken by every family, as it was
for (kind in c(refused, taken)) {
  case <- cases[[kind]]
  sampler <- case$make()
  expect_silent(sampler$setOffset(rep(1, n), updateScale = FALSE), info = kind)
  # an aft sampler takes no case weights at all
  if (kind != "aft") {
    expect_silent(
      sampler$setWeights(case$weights, updateScale = FALSE),
      info = kind
    )
  }
  expect_silent(
    sampler$setResponse(case$response, updateScale = FALSE),
    info = kind
  )
}

# with nothing saved, a gaussian sampler takes the call too: the refusal is of
# the saved draws, not of the family
bare <- dbarts(
  x,
  y,
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 8L,
    n.samples = 3L,
    n.burn = 5L,
    keepTrees = FALSE,
    updateState = FALSE,
    verbose = FALSE
  ),
  leaf.prior = gp(columns = 2L)
)
invisible(bare$run())
expect_silent(bare$setResponse(y + 1, updateScale = TRUE))
expect_silent(bare$setOffset(rep(-5, n), updateScale = TRUE))

# a scale-free sampler still refuses a re-derived standardization while it
# holds saved draws, which is the standardization's rule and not the scale's
sampler <- cases$probit$make()
expect_error(
  sampler$setWeights(mask, updateStandardization = TRUE),
  "'updateStandardization' cannot be TRUE for a sampler with gp leaves"
)
expect_error(
  sampler$setData(dbartsData(x, y.binary)),
  "cannot change the covariate standardization or the response range"
)
