# A state records the response family of the sampler that stored it. setState,
# copy and a reload refuse a state of another family, naming both, before
# anything of the sampler is touched. A state with no record is judged by its
# blocks, and a latent block of precisions holding a value that is not positive
# and finite is refused whatever the record says. A warm start takes any
# family's trees. No sampler here is swept after an install unless the state
# was of its own family: the draws that follow a refusal are asked only once
# the sampler has refused.

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
  probitMonotone = list(binary, family = "probit", monotone = c(x1 = 1)),
  aftVariance = list(
    survival,
    family = "aft",
    variance = dbartsForests$varianceForest(n.trees = 3L)
  )
)
eight <- names(kinds)[1:8]
familyOf <- c(
  setNames(eight, eight),
  studentDrawn = "student",
  nbinomFixed = "nbinom",
  hazard = "probit",
  hazardLogistic = "logistic",
  probitTwo = "probit",
  logisticTwo = "logistic",
  monotone = "gaussian",
  probitMonotone = "probit",
  aftVariance = "aft"
)

make <- function(kind, seed = 72L, n.chains = 1L) {
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
# data, and its next three draws, which are those of a twin offered nothing.
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
    expect_identical(bytes(sampler), before, info = info)
    expect_identical(sampler$run(0L, 3L), twin$run(0L, 3L), info = info)
  }
  message
}
plain <- "state is not consistent with this sampler"
named <- function(state, own) {
  sprintf(
    "%s: its family is \"%s\" and the sampler's is \"%s\"",
    plain,
    state,
    own
  )
}
edited <- function(state, name, value) {
  attr(state, name) <- value
  state
}
unrecorded <- function(state) edited(state, "family", NULL)

states <- lapply(setNames(nm = names(kinds)), function(kind) {
  stored(make(kind, 71L))
})

# --- what a state carries ---
expect_identical(lapply(states, attr, "family"), as.list(familyOf))

# --- another family's state is refused by name, the sampler untouched ---
pairs <- c(
  unlist(
    lapply(eight, function(own) {
      lapply(setdiff(eight, own), c, own)
    }),
    recursive = FALSE
  ),
  list(c("hazard", "hazardLogistic"), c("hazardLogistic", "hazard")),
  list(c("probitTwo", "logisticTwo"), c("logisticTwo", "probitTwo"))
)
expect_identical(length(pairs), 60L)
for (pair in pairs) {
  info <- paste(pair, collapse = " into ")
  expect_identical(
    refusal(pair[2L], states[[pair[1L]]], info),
    named(familyOf[[pair[1L]]], familyOf[[pair[2L]]]),
    info = info
  )
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
# a plain state into a monotone sampler is judged by its leaf values
for (kind in c("gaussian", "probit")) {
  target <- if (kind == "gaussian") "monotone" else "probitMonotone"
  message <- refusal(target, states[[kind]], paste(kind, "into", target))
  expect_false(grepl("family", message, fixed = TRUE), info = kind)
}

# --- a state from before the record: judged by its blocks ---
for (kind in names(kinds)) {
  with <- make(kind)
  without <- make(kind)
  expect_identical(
    without$setState(unrecorded(states[[kind]])),
    with$setState(states[[kind]]),
    info = kind
  )
  expect_identical(without$run(0L, 3L), with$run(0L, 3L), info = kind)
}
logistic <- make("logistic")
logistic$setState(unrecorded(states$student))
expect_true(finite(logistic))
# latent responses offered as precisions: the pairs that broke the sampler
for (pair in list(
  c("probit", "student"),
  c("ordinal", "student"),
  c("aft", "student"),
  c("probit", "logistic"),
  c("ordinal", "logistic"),
  c("aft", "logistic"),
  c("aft", "nbinom")
)) {
  info <- paste(pair, collapse = " into ")
  expect_true(any(states[[pair[1L]]][[1L]]$latents <= 0), info = info)
  message <- refusal(pair[2L], unrecorded(states[[pair[1L]]]), info)
  expect_identical(message, plain, info = info)
}

# --- malformed records ---
for (value in list(3, c("student", "student"), NA_character_, as.raw(1L))) {
  state <- edited(states$student, "family", value)
  expect_identical(
    refusal("student", state, "malformed"),
    "malformed family in bartcore state"
  )
}
# a name no family has is refused by that name, to a bounded width
for (value in c("", "Student", strrep("x", 500L))) {
  state <- edited(states$student, "family", value)
  expect_identical(
    refusal("student", state, "unknown"),
    named(substr(value, 1L, 32L), "student")
  )
}

# --- a latent block edited by hand ---
withLatent <- function(state, value) {
  state[[1L]]$latents[7L] <- value
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

# --- copy and a reload install the field by the same rule ---
original <- make("student")
twin <- make("student")
original$state <- states$probit
refused <- named("probit", "student")
expect_error(original$copy(), refused, fixed = TRUE)
reloaded <- unserialize(serialize(original, NULL))
expect_error(reloaded$run(0L, 1L), refused, fixed = TRUE)
expect_error(reloaded$run(0L, 1L), refused, fixed = TRUE)
reloaded$state <- states$student
expect_true(finite(reloaded))
expect_identical(original$run(0L, 3L), twin$run(0L, 3L))

# --- chains spliced from two samplers of one family ---
spliced <- stored(make("logistic", 81L))
spliced[[2L]] <- stored(make("logistic", 82L))[[1L]]
both <- make("logistic", n.chains = 2L)
expect_true(both$setState(spliced))
expect_true(finite(both))
expect_identical(
  refusal("probit", spliced, "spliced", n.chains = 2L),
  named("logistic", "probit")
)

# --- a warm start reads no record and no latents ---
for (state in list(states$probit, unrecorded(states$probit))) {
  warm <- make("student")
  warm$installTrees(state)
  expect_true(all(stored(warm)[[1L]]$latents > 0))
  expect_true(finite(warm))
}
