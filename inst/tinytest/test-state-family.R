# What a state of one response family meets in a sampler of another. A state
# names no family and none is asked: setState, copy and a reload install a
# state whose blocks fit the sampler, and what the sampler then holds of
# another family's state is not promised. A latent block of precisions holding
# a value that is not positive and finite is refused, which is every pair that
# would leave fits that are not finite or a sweep that does not return. A warm
# start takes any family's trees. The draws that follow a refusal are asked
# only of a sampler that refused and still stores the state it had.

student <- dbartsFamilies$student
nbinom <- dbartsFamilies$nbinom
forest <- dbartsForests$forest

set.seed(20261007L)
n <- 60L
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, c("x1", "x2", "x3")))
f <- 2.5 * x[, 1L] - 1.2
z <- rep(c(0, 1), length.out = n)
status <- rep(c(1, 1, 0), length.out = n)
continuous <- f + rnorm(n, sd = 0.5)
binary <- as.double(rbinom(n, 1L, plogis(2 * f)))
ordered <- cut(f + rnorm(n, sd = 0.5), c(-Inf, -0.4, 0.5, Inf), labels = FALSE)
counts <- as.double(rnbinom(n, size = 3, mu = exp(0.6 * f + 1)))
category <- factor(sample(c("a", "b", "c"), n, TRUE))
survival <- cbind(time = exp(f + rnorm(n, sd = 0.5)), status = status)
periods <- cbind(time = as.double(pmin(3L, 1L + rpois(n, 1))), status = status)
twoForests <- list(forest(), forest(basis = ~z))

# the eight families first, each under its own name, then other forms of them
kinds <- list(
  gaussian = list(continuous, family = "gaussian"),
  student = list(continuous, family = student(5)),
  probit = list(binary, family = "probit"),
  logistic = list(binary, family = "logistic"),
  ordinal = list(ordered, family = "ordinal"),
  nbinom = list(counts, family = "nbinom"),
  aft = list(survival, family = "aft"),
  multinomial = list(category, family = "multinomial"),
  studentDrawn = list(continuous, family = student()),
  nbinomFixed = list(counts, family = nbinom(shape = 3)),
  hazard = list(periods, family = "hazard"),
  hazardLogistic = list(periods, family = "hazard.logistic"),
  probitTwo = list(binary, family = "probit", forests = twoForests),
  logisticTwo = list(binary, family = "logistic", forests = twoForests),
  monotone = list(continuous, family = "gaussian", monotone = c(x1 = 1)),
  probitMonotone = list(binary, family = "probit", monotone = c(x1 = 1))
)
make <- function(kind, seed = 72L, n.chains = 1L, active = NULL) {
  control <- dbartsControl(
    n.chains = n.chains,
    n.threads = 1L,
    n.trees = 5L,
    n.samples = 4L,
    updateState = FALSE,
    seed = seed,
    verbose = FALSE
  )
  sampler <- do.call(dbarts, c(list(x), kinds[[kind]], list(control = control)))
  if (!is.null(active)) {
    sampler$setActiveRows(active)
  }
  sampler$run(3L, 2L)
  sampler
}
stored <- function(sampler) {
  sampler$storeState()
  sampler$state
}
# the state a sampler stores, the generator in it, and its data, as bytes
bytes <- function(sampler) serialize(list(stored(sampler), sampler$data), NULL)
finite <- function(sampler) {
  all(is.finite(unlist(sampler$run(0L, 3L)[c("train", "sigma")])))
}
# The message an install is refused with, or NA where it went through. After a
# refusal the sampler is as it was: its state field, the state it stores, its
# data, and its next three draws, which are those of a twin offered nothing;
# the draws are asked only of a sampler whose state is the one before.
refusal <- function(kind, state, info, ...) {
  sampler <- make(kind, ...)
  twin <- make(kind, ...)
  before <- bytes(sampler)
  field <- sampler$state
  message <- tryCatch(
    {
      sampler$setState(state)
      NA_character_
    },
    error = conditionMessage
  )
  if (!is.na(message)) {
    expect_identical(sampler$state, field, info = info)
    untouched <- identical(bytes(sampler), before)
    expect_true(untouched, info = info)
    if (untouched) {
      expect_identical(sampler$run(0L, 3L), twin$run(0L, 3L), info = info)
    }
  }
  message
}
plain <- "state is not consistent with this sampler"
states <- lapply(setNames(nm = names(kinds)), function(kind) {
  stored(make(kind, 71L))
})

# --- a state names no family ---
severalChains <- stored(make("logistic", 71L, n.chains = 2L))
for (state in c(states, list(severalChains))) {
  expect_null(attr(state, "family"))
}
# and an attribute of that name is not read, whatever it holds: the state
# installs as the one without it does
for (value in list("student", "probit", 3)) {
  labelled <- states$student
  attr(labelled, "family") <- value
  with <- make("student")
  without <- make("student")
  expect_identical(with$setState(labelled), without$setState(states$student))
  expect_identical(with$run(0L, 3L), without$run(0L, 3L))
}

# --- one family, another setting: the state installs ---
# the live trees of the first chain and its latent block
chain <- function(state) {
  blocks <- c("tree.vars", "tree.values", "tree.sizes", "tree.flags")
  c(state[[1L]]$forests[[1L]][c(blocks, "tree.params")], state[[1L]]["latents"])
}
for (pair in list(
  c("student", "studentDrawn"),
  c("studentDrawn", "student"),
  c("nbinom", "nbinomFixed"),
  c("nbinomFixed", "nbinom"),
  c("monotone", "gaussian"),
  c("probitMonotone", "probit")
)) {
  info <- paste(pair, collapse = " into ")
  state <- states[[pair[1L]]]
  sampler <- make(pair[2L])
  expect_true(sampler$setState(state), info = info)
  expect_identical(chain(stored(sampler)), chain(state), info = info)
  expect_true(finite(sampler), info = info)
}

# --- latent responses offered as precisions ---
# a block of precisions is taken from any family that holds one
logistic <- make("logistic")
logistic$setState(states$student)
expect_true(finite(logistic))
# Every pair here whose state has real-valued latents where the sampler holds
# Student-t scales or Polya-Gamma variates and fits the sampler otherwise:
# installed, such a block leaves fits that are not finite or a sweep that does
# not return. Each is refused, and for those values alone: the same state with
# its latents made positive installs.
for (pair in list(
  c("probit", "student"),
  c("ordinal", "student"),
  c("aft", "student"),
  c("probitMonotone", "student"),
  c("probit", "studentDrawn"),
  c("ordinal", "studentDrawn"),
  c("aft", "studentDrawn"),
  c("probitMonotone", "studentDrawn"),
  c("probit", "logistic"),
  c("ordinal", "logistic"),
  c("aft", "logistic"),
  c("probitMonotone", "logistic"),
  c("aft", "nbinom"),
  c("hazard", "hazardLogistic"),
  c("probitTwo", "logisticTwo")
)) {
  info <- paste(pair, collapse = " into ")
  state <- states[[pair[1L]]]
  expect_true(any(state[[1L]]$latents <= 0), info = info)
  expect_identical(refusal(pair[2L], state, info), plain, info = info)
  state[[1L]]$latents <- abs(state[[1L]]$latents) + 1
  expect_identical(refusal(pair[2L], state, info), NA_character_, info = info)
}

# --- a latent block edited by hand ---
withLatent <- function(state, value, row = 7L, chain = 1L) {
  state[[chain]]$latents[row] <- value
  state
}
for (value in c(0, -1, NaN, Inf)) {
  for (kind in c("student", "studentDrawn", "logistic", "nbinom")) {
    info <- paste(kind, value)
    message <- refusal(kind, withLatent(states[[kind]], value), info)
    expect_identical(message, plain, info = info)
  }
  # a latent response is any number
  state <- withLatent(states$probit, value)
  probit <- make("probit")
  expect_true(probit$setState(state), info = paste("probit", value))
  expect_identical(stored(probit)[[1L]]$latents, state[[1L]]$latents)
}
# every chain is asked, the last included
expect_identical(
  refusal(
    "logistic",
    withLatent(severalChains, -1, chain = 2L),
    "last chain",
    n.chains = 2L
  ),
  plain
)
# Every row is asked, a row the mask has out included: a variate that is not
# positive there would leave the sampler's next sweep without return. The
# install is refused before the engine sweeps, and the draws that follow are
# asked only of a sampler that refused, so nothing here sweeps such a state.
active <- rep(c(1, 0, 1), length.out = n)
masked <- stored(make("logistic", 71L, active = active))
expect_true(make("logistic", active = active)$setState(masked))
for (value in c(0, -1)) {
  state <- withLatent(masked, value, row = 2L)
  expect_identical(
    refusal("logistic", state, paste("masked", value), active = active),
    plain
  )
}

# --- copy and a reload install the field by the same rule ---
original <- make("student")
twin <- make("student")
original$state <- states$probit
expect_error(original$copy(), plain, fixed = TRUE)
reloaded <- unserialize(serialize(original, NULL))
expect_error(reloaded$run(0L, 1L), plain, fixed = TRUE)
expect_error(reloaded$run(0L, 1L), plain, fixed = TRUE)
reloaded$state <- states$student
expect_true(finite(reloaded))
expect_identical(original$run(0L, 3L), twin$run(0L, 3L))

# --- chains spliced from two samplers of one family ---
spliced <- stored(make("logistic", 81L))
spliced[[2L]] <- stored(make("logistic", 82L))[[1L]]
both <- make("logistic", n.chains = 2L)
expect_true(both$setState(spliced))
expect_true(finite(both))

# --- a warm start reads no latents ---
warm <- make("student")
warm$installTrees(states$probit)
taken <- stored(warm)[[1L]]$forests[[1L]]
expect_identical(taken$tree.vars, states$probit[[1L]]$forests[[1L]]$tree.vars)
expect_true(all(stored(warm)[[1L]]$latents > 0))
expect_true(finite(warm))
