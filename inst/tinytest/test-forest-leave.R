# A forest sharing the fitting function's leaf prior as a part of one total
# leaves the group when $setForest states its sd (dec-B525): it takes that size,
# its k pinned at 1, and the forests still sharing keep their parts, set once at
# creation and recorded on each forest's record. No creation reaches a group
# with a forest gone, so the oracles are a sampler created from the live
# sampler's records, creation stating the parts in force as sizes, and
# creation itself where every forest has left. Three sharers at unequal tree
# counts beside a factor forest, on two chains.

forest <- dbartsForests$forest
normal <- dbartsPriors$normal
chi <- dbartsPriors$chi

set.seed(20261103)
n <- 120L
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
z <- factor(rbinom(n, 1L, 0.5))
y <- sin(3 * x1) + x2 + x3^2 + (z == "1") * 0.5 + rnorm(n, sd = 0.3)
yb <- as.double(y > stats::median(y))
d <- data.frame(y, yb, x1, x2, x3, z)

ctl <- function(...) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 5L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 7L,
    verbose = FALSE,
    ...
  )
}
# three sharers, 38, 37 and 37 trees, and a factor forest; `sds` states a
# forest's sd by position
build <- function(response = "y", sds = list(), ...) {
  term <- function(index) {
    column <- as.name(paste0("x", index))
    trees <- if (index == 1L) 38L else 37L
    sd <- sds[[as.character(index)]]
    if (is.null(sd)) {
      bquote(forest(.(column), n.trees = .(trees)))
    } else {
      bquote(forest(.(column), n.trees = .(trees), sd = .(sd)))
    }
  }
  formula <- bquote(
    .(as.name(response)) ~
      .(term(1L)) +
      .(term(2L)) +
      .(term(3L)) +
      forest(x1, basis = z, n.trees = 6L)
  )
  dbarts(stats::as.formula(formula), d, control = ctl(), ...)
}
# a model's records with what a formula holds of its frame left out
modelOf <- function(sampler) {
  model <- sampler$model
  for (index in seq_along(model@forests)) {
    model@forests[[index]]@basis.term <- NULL
  }
  model
}
recordedTwin <- function(sampler) {
  new("dbartsSampler", sampler$control, sampler$model, sampler$data)
}
draws <- function(sampler) {
  sampler$run(0L, 5L)[c("sigma", "train", "k")]
}
shares <- function(sampler) {
  vapply(sampler$model@forests, methods::slot, 0, "share")
}
sizes <- function(sampler) {
  vapply(1:4, function(index) unname(sampler$getLeafPrior(index)$sd), 0)
}
reloaded <- function(sampler) {
  sampler$storeState()
  file <- tempfile(fileext = ".rds")
  saveRDS(sampler, file)
  loaded <- readRDS(file)
  unlink(file)
  loaded
}
refusal <- function(expr) {
  tryCatch(
    {
      expr
      "taken"
    },
    error = conditionMessage
  )
}
third <- 1 / 3
# the sd each forest leaves at
leaveSd <- c(0.3, 0.4, 0.5)

# --- the records: a part on every sharer from creation, gone on a leave ---
live <- build()
expect_identical(shares(live), c(third, third, third, NA))
expect_identical(dbarts:::samplerSharers(live), 1:3)
expect_identical(
  shares(dbarts(y ~ forest(x1) + forest(x2, basis = z), d, control = ctl())),
  c(1, NA)
)
expect_identical(
  shares(dbarts(y ~ forest(x1) + forest(x2, sd = 1), d, control = ctl())),
  c(1, NA)
)
before <- sizes(live)
expect_null(live$setForest(1L, sd = 0.4))
expect_identical(shares(live), c(NA, third, third, NA))
expect_identical(live$model@forests[[1L]]@sd, 0.4)
expect_identical(dbarts:::samplerSharers(live), 2:3)
after <- sizes(live)
expect_identical(after[-1L], before[-1L])
expect_identical(after[[1L]], 0.4)
expect_identical(unname(live$getK()[, 1L]), c(1, 2, 2, 1))
expect_identical(live$getLeafPrior(1L)$sd.stated, TRUE)
expect_identical(live$getLeafPrior(2L)$prior.sd.of, "forest, shared")
expect_identical(live$getForest(2L)$sd, before[[2L]])
# a value equal to the one in force is no change: the sharer stays
expect_null(live$setForest(2L, sd = before[[2L]]))
expect_identical(shares(live), c(NA, third, third, NA))
expect_error(
  methods::validObject(local({
    record <- live$model@forests[[1L]]
    record@share <- third
    record
  })),
  "where no 'sd' is stated"
)
# records whose shares name no one group are refused
fresh <- build()
model <- fresh$model
model@forests[[2L]]@share <- NA_real_
expect_error(
  new("dbartsSampler", fresh$control, model, fresh$data),
  "give a share of the prior total to some forests that take it and not to others",
  fixed = TRUE
)
model@forests[[2L]]@share <- 2
expect_error(
  new("dbartsSampler", fresh$control, model, fresh$data),
  "'share' must be"
)

# --- the recorded twin: a leave before any run, bitwise ---
for (response in c("y", "yb")) {
  for (leaving in list(1L, c(3L, 1L), c(2L, 1L, 3L))) {
    info <- paste(response, "leaving", paste(leaving, collapse = " "))
    live <- build(response)
    for (index in leaving) {
      live$setForest(index, sd = leaveSd[[index]])
    }
    # the records of all left create, the leaf prior gone with the group
    twin <- tryCatch(recordedTwin(live), error = conditionMessage)
    expect_true(
      is(twin, "dbartsSampler"),
      info = paste(info, if (is.character(twin)) twin)
    )
    if (is.character(twin)) {
      next
    }
    expect_identical(modelOf(twin), modelOf(live), info = info)
    expect_identical(twin$getK(), live$getK(), info = info)
    liveDraws <- draws(live)
    expect_identical(dim(liveDraws$train), c(n, 5L, 2L), info = info)
    expect_identical(liveDraws, draws(twin), info = paste(info, "draws"))
  }
}

# --- the stated-sizes twin: the parts creation gave, stated ---
for (leaving in list(1L, c(1L, 3L))) {
  info <- paste("stated sizes, leaving", paste(leaving, collapse = " "))
  live <- build()
  parts <- sizes(live)
  for (index in leaving) {
    live$setForest(index, sd = leaveSd[[index]])
  }
  stated <- as.list(parts[1:3])
  names(stated) <- as.character(1:3)
  stated[as.character(leaving)] <- as.list(leaveSd[leaving])
  twin <- build(sds = stated)
  spread <- function(sampler) {
    vapply(
      1:3,
      function(index) {
        prior <- sampler$getLeafPrior(index)
        unname(prior$k.scale / sampler$getK(index)[1L])
      },
      0
    )
  }
  expect_equal(spread(live), spread(twin), tolerance = 1e-12, info = info)
  expect_false(isTRUE(all.equal(spread(live), spread(build()))), info = info)
  liveDraws <- draws(live)
  twinDraws <- draws(twin)
  expect_equal(liveDraws$train, twinDraws$train, tolerance = 1e-8, info = info)
  expect_equal(liveDraws$sigma, twinDraws$sigma, tolerance = 1e-8, info = info)
}

# --- creation where none is left ---
for (response in c("y", "yb")) {
  info <- paste("none left,", response)
  live <- build(response)
  for (index in 1:3) {
    live$setForest(index, sd = leaveSd[[index]])
  }
  twin <- build(response, sds = list(`1` = 0.3, `2` = 0.4, `3` = 0.5))
  expect_identical(modelOf(live), modelOf(twin), info = info)
  expect_identical(draws(live), draws(twin), info = info)
}

# --- forest 1 leaving under a drawn k: the k column stays, the group's ---
live <- build("yb")
live$setForest(1L, sd = 0.5)
run <- live$run()
expect_false(is.null(run$k))
expect_identical(as.numeric(run$k[5L, ]), unname(live$getK(2L)))
expect_identical(unname(live$getK(1L)), c(1, 1))
expect_identical(unname(live$getK(3L)), unname(live$getK(2L)))
# one left of the group, forest 2: its k still recorded
live$setForest(3L, sd = 0.4)
again <- live$run(0L, 4L)
expect_identical(dim(again$k), c(4L, 2L))
expect_identical(as.numeric(again$k[4L, ]), unname(live$getK(2L)))
# none left: no k is drawn, and none recorded
live$setForest(2L, sd = 0.3)
expect_null(live$run(0L, 4L)$k)

# --- after a run: a copy, a reload and the state go on as the live sampler ---
for (when in c("before", "after")) {
  info <- paste("state", when, "the leave")
  live <- build("yb")
  invisible(live$run(0L, 5L))
  if (when == "after") {
    live$setForest(2L, sd = 0.6)
  }
  live$storeState()
  saved <- live$state
  if (when == "before") {
    live$setForest(2L, sd = 0.6)
  }
  twin <- recordedTwin(live)
  twin$setState(saved)
  live$setState(saved)
  copied <- live$copy()
  loaded <- reloaded(live)
  next_ <- draws(twin)
  expect_identical(draws(copied), next_, info = paste(info, "copy"))
  expect_identical(draws(loaded), next_, info = paste(info, "reload"))
  expect_identical(draws(live), next_, info = paste(info, "live"))
  expect_identical(
    modelOf(loaded),
    modelOf(live),
    info = paste(info, "records")
  )
}
# a state from after a leave on a sampler where the leaver still shares: the
# group takes one k
left <- build("yb")
invisible(left$run(0L, 5L))
left$setForest(1L, sd = 0.5)
invisible(left$run(0L, 5L))
left$storeState()
sharing <- build("yb")
sharing$setState(left$state)
ks <- sharing$getK()
expect_identical(ks[1L, ], ks[2L, ])
expect_identical(ks[3L, ], ks[2L, ])
expect_identical(ks[2L, ], left$getK()[2L, ])

# --- a total stated after a leave reaches those left at their parts ---
live <- build()
part <- sizes(live)[[2L]]
live$setForest(1L, sd = 0.4)
live$setLeafPrior(normal(k = 4))
expect_equal(sizes(live)[2:3], rep(part / 2, 2L), tolerance = 1e-14)
expect_identical(sizes(live)[[1L]], 0.4)
expect_identical(shares(live), c(NA, third, third, NA))
twin <- recordedTwin(live)
expect_identical(draws(live), draws(twin))
expect_null(live$setModel(twin$model))
reference <- build(leaf.prior = normal(k = 1))$model
model <- live$model
model@leaf.prior <- reference@leaf.prior
model@leaf.hyperprior <- reference@leaf.hyperprior
model@prior.scale <- reference@prior.scale
live$setModel(model)
expect_equal(sizes(live)[2:3], rep(part * 2, 2L), tolerance = 1e-14)
# under a named total each left has its part of it, the parts no longer adding
live$setLeafPrior(normal(sd = 1.2))
expect_equal(sizes(live)[2:3], rep(1.2 * sqrt(third), 2L), tolerance = 1e-14)
# one left of three, at a held 2: its part, not the whole
one <- build()
one$setForest(1L, sd = 0.4)
one$setForest(3L, sd = 0.3)
oneSize <- sizes(one)[[2L]]
expect_equal(oneSize, part, tolerance = 1e-14)
one$setLeafPrior(normal(k = 3))
expect_equal(sizes(one)[[2L]], part * 2 / 3, tolerance = 1e-14)
one$setLeafPrior(normal())
expect_equal(sizes(one)[[2L]], part, tolerance = 1e-14)
expect_identical(shares(one), c(NA, third, NA, NA))
expect_identical(one$getLeafPrior(2L)$prior.sd.of, "forest, shared")
expect_identical(draws(one), draws(recordedTwin(one)))

# --- the last leave: no forest takes the leaf prior ---
live <- build("yb", leaf.prior = normal(k = chi(1.5, 2)))
for (index in 1:3) {
  live$setForest(index, sd = leaveSd[[index]])
}
expect_identical(
  modelOf(live),
  modelOf(build("yb", sds = list(`1` = 0.3, `2` = 0.4, `3` = 0.5)))
)
loaded <- tryCatch(reloaded(live), error = conditionMessage)
expect_true(
  is(loaded, "dbartsSampler"),
  info = if (is.character(loaded)) loaded
)
liveDraws <- tryCatch(draws(live), error = conditionMessage)
expect_true(is.list(liveDraws), info = if (is.character(liveDraws)) liveDraws)
loadedDraws <- if (is.character(loaded)) {
  loaded
} else {
  tryCatch(draws(loaded), error = conditionMessage)
}
expect_identical(liveDraws, loadedDraws)
expect_identical(refusal(live$setLeafPrior(normal())), "taken")
expect_error(
  live$setLeafPrior(normal(k = 3)),
  "no forest of this sampler takes the fitting function's leaf prior",
  fixed = TRUE
)
expect_error(live$setLeafPrior(normal(k = 3)), "$setForest(", fixed = TRUE)

# --- $setModel: a share is set once ---
live <- build()
saved <- live$model
live$setForest(1L, sd = 0.4)
expect_error(live$setModel(saved), "a new sampler", fixed = TRUE)
other <- build(sds = list(`1` = 0.4))
expect_identical(shares(other), c(NA, 0.5, 0.5, NA))
expect_error(
  live$setModel(other$model),
  "$setModel cannot change forest 2's share of the prior total: a share is set once",
  fixed = TRUE
)
expect_identical(shares(live), c(NA, third, third, NA))
expect_null(live$setModel(live$model))

# --- the one candidate taking the leaf prior records a share of 1 ---
lone <- dbarts(
  y ~ forest(x1, sd = 2) + forest(x2, basis = z),
  d,
  control = ctl()
)
expect_identical(shares(lone), c(NA_real_, NA_real_))
lone$setLeafPrior(normal())
expect_identical(shares(lone), c(1, NA))
expect_identical(draws(lone), draws(recordedTwin(lone)))

# --- print ---
live <- build()
live$setForest(1L, sd = 0.4)
live$setForest(3L, sd = 0.3)
lines <- dbarts:::forestPriorLines(
  live$getLeafPrior(),
  c("x1", "x2", "x3", "z"),
  live$data@bases,
  unit = 1
)
expect_true(
  any(grepl(
    "forest 2 holds its part of one prior total, k 2",
    lines,
    fixed = TRUE
  )),
  info = paste(lines, collapse = "\n")
)
verbose <- utils::capture.output(
  new(
    "dbartsSampler",
    local({
      control <- live$control
      control@verbose <- TRUE
      control
    }),
    live$model,
    live$data
  )
)
expect_true(
  any(grepl(
    "forest sharing one prior total: 2, 1/3 of its variance",
    verbose,
    fixed = TRUE
  )),
  info = paste(verbose, collapse = "\n")
)
twoLeft <- build()
twoLeft$setForest(1L, sd = 0.4)
verbose <- utils::capture.output(
  new(
    "dbartsSampler",
    local({
      control <- twoLeft$control
      control@verbose <- TRUE
      control
    }),
    twoLeft$model,
    twoLeft$data
  )
)
expect_true(
  any(grepl(
    "forests sharing one prior total: 2 3, each 1/3 of its variance",
    verbose,
    fixed = TRUE
  )),
  info = paste(verbose, collapse = "\n")
)

# --- the multinomial model and the variance forest have no group to leave ---
counts <- cbind(a = rpois(n, 2), b = rpois(n, 2), c = rpois(n, 2))
multinomial <- dbarts(
  cbind(x1, x2),
  counts,
  family = "multinomial",
  control = dbartsControl(n.chains = 1L, n.threads = 1L, verbose = FALSE)
)
expect_error(
  multinomial$setForest(sd = 0.5),
  "a named leaf-prior 'sd'",
  fixed = TRUE
)
