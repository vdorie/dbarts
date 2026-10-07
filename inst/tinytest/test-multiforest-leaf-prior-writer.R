# $setLeafPrior on the two multi-forest samplers: a multinomial sampler takes
# normal(k = ) on every category forest, and a sampler whose forests carry
# amplitudes takes forests = list(forest(sd = ), ...) in creation's channels.
# The main oracle is the twin: creation consumes R's RNG and a write does not,
# so a sampler created with P that writes P' at once is bitwise the sampler
# created with P' under the same seed, on every chain and every channel.

forest <- dbartsForests$forest

set.seed(71)
n <- 150L
p <- 3L
x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
z <- rbinom(n, 1L, 0.5)
g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
y <- 4 * sin(pi * x[, 1L]) + z * (1 + 2 * x[, 2L]) + rnorm(n, sd = 0.3)
yBinary <- as.double(y > median(y))
labels <- factor(sample(0:2, n, replace = TRUE))

writerControl <- function(...) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.trees = 15L,
    n.samples = 4L,
    updateState = FALSE,
    seed = 71L,
    ...
  )
}
numWarnings <- 0L
counted <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    numWarnings <<- numWarnings + 1L
    invokeRestart("muffleWarning")
  })
}
amplitudeSampler <- function(forests, response = y, ...) {
  set.seed(5L)
  counted(dbarts(x, response, forests = forests, control = writerControl(...)))
}
multinomialSampler <- function(...) {
  set.seed(5L)
  counted(dbarts(
    x,
    labels,
    family = "multinomial",
    control = writerControl(),
    ...
  ))
}
forestInfo <- function(sampler) {
  attr(sampler$control, "bartcore.forests", exact = TRUE)
}
expectTwins <- function(a, b, info) {
  expect_identical(a$getLeafPrior(), b$getLeafPrior(), info = info)
  expect_identical(forestInfo(a)$params, forestInfo(b)$params, info = info)
  expect_identical(forestInfo(a)$anchor, forestInfo(b)$anchor, info = info)
  expect_identical(counted(a$run(3L, 4L)), counted(b$run(3L, 4L)), info = info)
}

twoDefault <- list(forest(), forest(basis = ~ factor(z)))
threeDefault <- list(forest(), forest(basis = ~ factor(z)), forest(basis = ~g))

# --- twin identity ---
a <- multinomialSampler()
a$setLeafPrior(normal(k = 3))
b <- multinomialSampler(leaf.prior = normal(k = 3))
expect_identical(a$model, b$model)
expectTwins(a, b, "multinomial k 2 -> 3")
a <- multinomialSampler()
a$setLeafPrior(normal(k = Inf))
expectTwins(
  a,
  multinomialSampler(leaf.prior = normal(k = Inf)),
  "multinomial k -> Inf"
)

# the fixed-variance twin's sd is one whose reported anchor moves under a
# reassociation of the map's leaf-scale expression, found against this data's
# own anchor s, so the twin's bitwise anchor sees a reassociation
reference <- amplitudeSampler(twoDefault)
s <- forestInfo(reference)$anchor
unit <- reference$getLeafPrior(2L)$response.scale * sqrt(50)
reassociationMoves <- function(sd) {
  !identical(
    sd * s / 0.674 / sqrt(50) * unit,
    sd * (s / 0.674) / sqrt(50) * unit
  )
}
fixedSd <- Find(reassociationMoves, seq(0.41, 0.79, by = 0.01))
expect_false(is.null(fixedSd))

twinCases <- list(
  list(
    "gaussian, fixed-variance sd",
    y,
    list(forest(), forest(sd = fixedSd)),
    list(forest(), forest(basis = ~ factor(z), sd = fixedSd))
  ),
  list(
    "gaussian, scale-mixture median",
    y,
    list(forest(sd = 1.7)),
    list(forest(sd = 1.7), forest(basis = ~ factor(z)))
  ),
  list(
    "probit, both at once",
    yBinary,
    list(forest(sd = 1.7), forest(sd = 0.6)),
    list(forest(sd = 1.7), forest(basis = ~ factor(z), sd = 0.6))
  ),
  list(
    "K = 3, every forest",
    y,
    list(forest(sd = 1.7), forest(sd = 0.6), forest(sd = 0.45)),
    list(
      forest(sd = 1.7),
      forest(basis = ~ factor(z), sd = 0.6),
      forest(basis = ~g, sd = 0.45)
    )
  )
)
for (case in twinCases) {
  default <- if (length(case[[4L]]) == 3L) threeDefault else twoDefault
  a <- amplitudeSampler(default, case[[2L]])
  a$setLeafPrior(forests = case[[3L]])
  expectTwins(a, amplitudeSampler(case[[4L]], case[[2L]]), case[[1L]])
}

# --- mid-run half-Cauchy median, which rides no state: both install A's ---
a <- amplitudeSampler(twoDefault)
counted(a$run(5L, 2L))
a$setLeafPrior(forests = list(forest(sd = 0.35)))
a$storeState()
b <- amplitudeSampler(list(forest(sd = 0.35), forest(basis = ~ factor(z))))
a$setState(a$state)
b$setState(a$state)
expect_identical(counted(a$run(3L, 4L)), counted(b$run(3L, 4L)))

# --- round trip: writing the reader back is inert, as is normal() ---
for (build in list(
  function() multinomialSampler(leaf.prior = normal(k = 3)),
  function() amplitudeSampler(twoDefault),
  function() amplitudeSampler(threeDefault, yBinary)
)) {
  a <- build()
  b <- build()
  expect_identical(counted(a$run(3L, 2L)), counted(b$run(3L, 2L)))
  if (dbarts:::samplerCarriesCounts(a)) {
    a$setLeafPrior(a$getLeafPrior(1L)$leaf.prior)
  } else {
    a$setLeafPrior(forests = lapply(a$getLeafPrior(), `[[`, "leaf.prior"))
    a$setLeafPrior(normal())
    a$setLeafPrior(normal(k = 2))
  }
  expectTwins(a, b, "round trip")
}

# --- discrimination: a changed write moves draws and the reader follows ---
a <- multinomialSampler()
b <- multinomialSampler()
a$setLeafPrior(normal(k = 4))
expect_identical(a$getK(), matrix(4, 3L, 2L))
expect_false(identical(counted(a$run(3L, 2L)), counted(b$run(3L, 2L))))
a <- amplitudeSampler(twoDefault)
b <- amplitudeSampler(twoDefault)
a$setLeafPrior(forests = list(forest(), forest(sd = 0.6)))
expect_equal(
  a$getLeafPrior(2L)$k.scale / b$getLeafPrior(2L)$k.scale,
  0.6,
  tolerance = 1e-14
)
expect_identical(a$getLeafPrior(2L)$leaf.prior, forest(sd = 0.6))
expect_identical(a$getLeafPrior(1L), b$getLeafPrior(1L))
expect_false(identical(counted(a$run(3L, 2L)), counted(b$run(3L, 2L))))
a <- amplitudeSampler(twoDefault)
b <- amplitudeSampler(twoDefault)
a$setLeafPrior(forests = list(forest(sd = 0.35)))
expect_identical(a$getLeafPrior(1L)$amplitude.prior.scale, 0.35)
expect_identical(a$getLeafPrior(1L)$k.scale, b$getLeafPrior(1L)$k.scale)
expect_false(identical(counted(a$run(3L, 2L)), counted(b$run(3L, 2L))))

# --- re-creation: the anchor s is carried, so a response swap at
# updateScale = FALSE survives a copy, a save and load, and a basis install ---
live <- amplitudeSampler(twoDefault)
counted(live$run(3L, 2L))
live$setResponse(y + 3 * x[, 3L], updateScale = FALSE)
live$storeState()
anchorOf <- function(sampler) {
  vapply(sampler$getLeafPrior(), `[[`, 0, "k.scale")
}
factorOf <- function(sampler) sampler$getLeafPrior(2L)$leaf.scale.factor
copied <- live$copy()
file <- tempfile(fileext = ".rds")
saveRDS(live, file)
loaded <- readRDS(file)
unlink(file)
rebased <- live$copy()
rebased$setForestBasis(2L, cbind(1 - z, z))
for (other in list(copied, loaded, rebased)) {
  expect_identical(anchorOf(other), anchorOf(live))
  expect_identical(factorOf(other), factorOf(live))
  expect_false(is.na(factorOf(other)))
}
spread <- list(forest(sd = 1.3), forest(sd = 0.9))
live$setLeafPrior(forests = spread)
copied$setLeafPrior(forests = spread)
expect_identical(copied$getLeafPrior(), live$getLeafPrior())

# a write followed by a stateless copy keeps the k, the factor and the median
a <- multinomialSampler()
a$setLeafPrior(normal(k = 3))
expect_null(a$state)
expect_identical(a$copy()$getK(), a$getK())
a <- amplitudeSampler(twoDefault)
a$setLeafPrior(forests = list(forest(sd = 0.35), forest(sd = 0.6)))
expect_null(a$state)
copied <- a$copy()
expect_identical(copied$getLeafPrior(1L)$leaf.prior, forest(sd = 0.35))
expect_identical(copied$getLeafPrior(2L)$leaf.prior, forest(sd = 0.6))
expect_identical(copied$getLeafPrior(), a$getLeafPrior())

# a pre-write state leaves the write standing on both forests: a state holds
# no leaf scale and no fixed amplitude variance, so restoring one undoes
# neither, and the reader and the next draws are a twin's that wrote and then
# restored its own state
a <- amplitudeSampler(twoDefault)
twin <- amplitudeSampler(twoDefault)
a$storeState()
before <- a$state
a$setLeafPrior(forests = list(forest(sd = 0.35), forest(sd = 0.6)))
a$setState(before)
twin$setLeafPrior(forests = list(forest(sd = 0.35), forest(sd = 0.6)))
twin$storeState()
twin$setState(twin$state)
expect_identical(a$getLeafPrior(2L)$leaf.prior, forest(sd = 0.6))
expect_identical(a$getLeafPrior(1L)$leaf.prior, forest(sd = 0.35))
expect_identical(a$getLeafPrior(), twin$getLeafPrior())
expect_identical(a$run(0L, 3L), twin$run(0L, 3L))

# a basis install on the scale-mixture forest keeps its median and channel
a <- amplitudeSampler(twoDefault)
a$setLeafPrior(forests = list(forest(sd = 0.35)))
a$setForestBasis(1L, matrix(2, n, 1L))
expect_identical(a$getLeafPrior(1L)$leaf.prior, forest(sd = 0.35))
expect_identical(a$getLeafPrior(1L)$prior.sd.of, "amplitude scale")
expect_identical(a$getLeafPrior(2L)$prior.sd.of, "forest total")

# --- refusals ---
m <- multinomialSampler()
expect_error(m$setLeafPrior(normal(sd = 1)), "named 'sd' has nowhere to land")
expect_error(m$setLeafPrior(normal(k = chi())), "'k' hyperprior is not suppo")
expect_error(m$setLeafPrior(normal(k = "chi(1.5)")), "'k' hyperprior")
expect_error(
  m$setLeafPrior(forests = list(forest())),
  "its forests are its categories"
)
a <- amplitudeSampler(twoDefault)
expect_error(a$setLeafPrior(normal(k = 3)), "forests = list\\(forest\\(sd")
expect_error(a$setLeafPrior(normal(sd = 1)), "multi-forest calibration map")
expect_error(
  a$setLeafPrior(forests = list(forest(n.trees = 5L))),
  "'n.trees' is fixed at creation"
)
expect_error(
  a$setLeafPrior(forests = list(forest(), forest(basis = ~ factor(z)))),
  "change it with \\$setForestBasis"
)
expect_error(
  a$setLeafPrior(forests = list(forest(), forest(), forest())),
  "names 3 forests; this sampler has 2"
)
expect_error(
  a$setLeafPrior(forests = list(first = forest(sd = 1))),
  "created unnamed"
)
expect_error(
  a$setLeafPrior(forests = list(dbartsPriors$normal())),
  "forest\\(\\) specif"
)
for (bad in list(0, -1, Inf)) {
  expect_error(
    a$setLeafPrior(forests = list(forest(sd = bad))),
    "forest 'sd' must be positive and finite"
  )
}
for (bad in list(NaN, NA)) {
  expect_error(
    a$setLeafPrior(forests = list(forest(sd = bad))),
    "forest 'sd' must not be NA"
  )
}
expect_error(
  a$setLeafPrior(forests = list(forest(sd = c(1, 2)))),
  "forest 'sd' must be a single number, not a vector of length 2"
)
expect_error(
  a$setLeafPrior(normal(), forests = list(forest())),
  "not both"
)
labelled <- amplitudeSampler(list(
  mu = forest(),
  tau = forest(basis = ~ factor(z))
))
labelled$setLeafPrior(forests = list(mu = forest(), tau = forest(sd = 0.6)))
expect_error(
  labelled$setLeafPrior(forests = list(tau = forest(sd = 0.6))),
  "created as 'mu'"
)
# every refusal came before any write
expect_identical(a$getLeafPrior(), amplitudeSampler(twoDefault)$getLeafPrior())
single <- counted(dbarts(x, y, control = writerControl()))
expect_error(
  single$setLeafPrior(forests = list(forest(sd = 1))),
  "only a sampler whose forests carry amplitudes"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(sd = Inf), forest(basis = ~ factor(z))),
    control = writerControl()
  ),
  "forest 'sd' must be positive and finite"
)

expect_identical(numWarnings, 0L)
