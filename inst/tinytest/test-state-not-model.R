# A saved state holds the chain, not the model. It carries trees, leaf
# values, the quantities the sampler draws, the generator and the units its
# numbers are stored in; it carries no prior parameter and no value the
# sampler holds fixed. An install - setState, copy, a reload, a warm start -
# never changes the sampler's model, and a state stored in other response
# units is converted into the sampler's own. The oracles: a store, write,
# restore matches a twin that wrote and restored its own state, reader and
# draws; a state carrying a block for a value the recipient holds fixed
# installs as one without it; and a converted state replays the function it
# held, to rounding.

# the vocabularies, for argument lists built ahead of the call
for (name in c("normal", "chi", "invchi", "linear", "gp", "dart", "fixed")) {
  assign(name, dbartsPriors[[name]])
}
student <- dbartsFamilies$student
nbinom <- dbartsFamilies$nbinom
gaussian <- dbartsFamilies$gaussian
forest <- dbartsForests$forest

set.seed(61)
n <- 80L
p <- 3L
x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
z <- rbinom(n, 1L, 0.5)
y <- 4 * sin(pi * x[, 1L]) + 2 * x[, 2L] + rnorm(n, sd = 0.3)
counts <- rpois(n, exp(0.5 + x[, 1L]))
# a response whose range is symmetric about zero, exactly, so a stretch of it
# moves the range and leaves the shift, the midpoint, at zero to the last bit
centred <- y - (min(y) + max(y)) / 2
centred[which.min(centred)] <- -max(centred)

stateControl <- function(...) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.trees = 12L,
    n.samples = 4L,
    updateState = FALSE,
    seed = 17L,
    ...
  )
}
make <- function(..., response = y, control = stateControl()) {
  set.seed(5L)
  sampler <- dbarts(x, response, control = control, ...)
  invisible(sampler$run(10L, 0L))
  sampler
}
stored <- function(sampler) {
  sampler$storeState()
  sampler$state
}
sweeps <- function(sampler) sampler$run(0L, 3L)
# the units a sampler's numbers are in, the sampler's transform
units <- function(sampler) stored(sampler)[[1L]]$fit.scale

# --- store, write, restore: the write stands, and the reader and the next
# draws are those of a twin that wrote and restored its own state ---
storeWriteRestore <- function(label, build, write) {
  a <- build()
  twin <- build()
  before <- stored(a)
  write(a)
  a$setState(before)
  write(twin)
  twin$setState(stored(twin))
  expect_identical(a$getLeafPrior(), twin$getLeafPrior(), info = label)
  expect_identical(a$getK(), twin$getK(), info = label)
  expect_identical(a$getSigmas(), twin$getSigmas(), info = label)
  expect_identical(sweeps(a), sweeps(twin), info = label)
}
storeWriteRestore("k", function() make(), function(s) {
  s$setLeafPrior(normal(k = 5))
})
storeWriteRestore(
  "named sd",
  function() make(leaf.prior = normal(sd = 1)),
  function(s) s$setLeafPrior(normal(sd = 0.3))
)
storeWriteRestore("sd law", function() make(), function(s) {
  s$setLeafPrior(normal(sd = invchi(3, 0.5)))
})
storeWriteRestore(
  "fixed sigma",
  function() {
    make(family = gaussian(sigma = fixed(1)))
  },
  function(s) s$setSigma(2)
)
storeWriteRestore(
  "forest spread",
  function() make(forests = list(forest(), forest(basis = ~ factor(z)))),
  function(s) {
    s$setLeafPrior(forests = list(forest(sd = 0.35), forest(sd = 0.6)))
  }
)

# a sigma written to a sampler that holds it fixed is recorded on the model,
# so a copy keeps it with no state stored at all
fixedSigma <- make(family = gaussian(sigma = fixed(1)))
fixedSigma$setSigma(2)
expect_equal(fixedSigma$model@resid.prior@value, 4)
expect_equal(unname(fixedSigma$copy()$getSigmas()), c(2, 2))

# --- a state carrying a block for a value the recipient holds fixed, as one
# written before did, installs as one without it ---
withBlocks <- function(label, build, edit) {
  a <- build()
  state <- stored(a)
  plain <- build()
  plain$setState(state)
  edited <- build()
  edited$setState(edit(state))
  expect_identical(edited$getLeafPrior(), plain$getLeafPrior(), info = label)
  expect_identical(sweeps(edited), sweeps(plain), info = label)
}
editChains <- function(state, f) {
  for (chain in seq_along(state)) {
    state[[chain]] <- f(state[[chain]])
  }
  state
}
withBlocks("k and leaf scale", function() make(), function(state) {
  editChains(state, function(chain) {
    chain$forests[[1L]]$k <- 4
    chain$forests[[1L]]$leaf.scale <- 0.01
    chain
  })
})
withBlocks(
  "sigma",
  function() {
    make(family = gaussian(sigma = fixed(1)))
  },
  function(state) {
    editChains(state, function(chain) {
      chain$sigma <- 3
      chain
    })
  }
)
withBlocks("df", function() make(family = student(df = 3)), function(state) {
  editChains(state, function(chain) {
    chain$resid.df <- 30
    chain
  })
})
withBlocks(
  "shape",
  function() {
    make(response = counts, family = nbinom(shape = 5))
  },
  function(state) {
    editChains(state, function(chain) {
      chain$shape <- 20
      chain
    })
  }
)
withBlocks(
  "concentration",
  function() {
    make(
      tree.prior = dart(alpha = 0.5, update.alpha = FALSE, update.delay = 0L)
    )
  },
  function(state) {
    editChains(state, function(chain) {
      chain$dart.alpha <- 5
      chain
    })
  }
)
# the glue is [K, q_1..q_K, amplitudes, K prior variances]: forest 2 holds a
# fixed variance in one build and fixed amplitudes in the other
twoFixed <- function(...) {
  make(forests = list(forest(), forest(basis = ~ factor(z), ...)))
}
withBlocks(
  "fixed amplitude variance",
  function() {
    twoFixed(amplitude.prior.variance = 2)
  },
  function(state) {
    editChains(state, function(chain) {
      chain$glue[[length(chain$glue)]] <- 50
      chain
    })
  }
)
withBlocks(
  "fixed amplitudes",
  function() {
    twoFixed(amplitude = dbartsPriors$fixed())
  },
  function(state) {
    editChains(state, function(chain) {
      chain$glue[5:6] <- c(3, -2)
      chain
    })
  }
)
# none of those blocks is written where the value is fixed
expect_null(stored(make(family = gaussian(sigma = fixed(1))))[[1L]]$sigma)
expect_null(stored(make(family = student(df = 3)))[[1L]]$resid.df)
expect_null(
  stored(make(response = counts, family = nbinom(shape = 5)))[[1L]]$shape
)
expect_null(
  stored(make(tree.prior = dart(alpha = 0.5, update.alpha = FALSE)))[[
    1L
  ]]$dart.alpha
)
expect_null(stored(make())[[1L]]$forests[[1L]]$k)

# --- drawn values still install ---
drawn <- make(leaf.prior = normal(k = chi(1.5, 2)), family = student())
state <- editChains(stored(drawn), function(chain) {
  chain$sigma <- 0.7
  chain$forests[[1L]]$k <- 3.5
  chain$resid.df <- 10
  chain
})
recipient <- make(leaf.prior = normal(k = chi(1.5, 2)), family = student())
recipient$setState(state)
expect_equal(unname(recipient$getSigmas()), c(0.7, 0.7))
expect_identical(recipient$getK(), c(3.5, 3.5))
expect_identical(stored(recipient)[[2L]]$resid.df, 10)
drawnShape <- make(response = counts, family = nbinom())
shapeState <- editChains(stored(drawnShape), function(chain) {
  chain$shape <- 30
  chain
})
drawnShape$setState(shapeState)
expect_identical(drawnShape$getShape(), c(30, 30))
drawnAlpha <- make(tree.prior = dart())
alphaState <- editChains(stored(drawnAlpha), function(chain) {
  chain$dart.alpha <- 0.25
  chain
})
drawnAlpha$setState(alphaState)
expect_identical(stored(drawnAlpha)[[1L]]$dart.alpha, 0.25)

# a gaussian state leaves a probit sampler's sigma, pinned at 1, where it is
yb <- as.double(y > median(y))
probit <- make(response = yb, family = "binomial")
probit$setState(stored(make()))
expect_identical(unname(probit$getSigmas()), c(1, 1))

# --- a state from one model under another: the recipient's reader stands ---
underAnother <- function(label, recipientArgs, donorArgs) {
  recipient <- do.call(make, recipientArgs)
  priorBefore <- recipient$getLeafPrior()
  recipient$setState(stored(do.call(make, donorArgs)))
  expect_identical(recipient$getLeafPrior(), priorBefore, info = label)
}
underAnother(
  "k",
  list(leaf.prior = normal(k = 2)),
  list(leaf.prior = normal(k = chi(1.5, 2)))
)
underAnother(
  "sd",
  list(leaf.prior = normal(sd = 0.5)),
  list(leaf.prior = normal(sd = 2))
)
underAnother(
  "df",
  list(family = student(df = 10)),
  list(family = student())
)
withDf <- make(family = student(df = 10))
withDf$setState(stored(make(family = student())))
expect_true(all(withDf$run(0L, 2L)$resid.df == 10))
withShape <- make(response = counts, family = nbinom(shape = 10))
withShape$setState(stored(make(response = counts, family = nbinom())))
expect_identical(withShape$getShape(), c(10, 10))

# --- the anchor, and conversion on install ---

# the reader is the recipient's across an install from a sampler on another
# response, in other units, for every naming and leaf model
rescaled <- 3 * y + 10
stretched <- 3 * centred # another range, the same shift
acrossUnits <- function(label, args, response = rescaled, base = y) {
  recipient <- do.call(make, c(args, list(response = base)))
  priorBefore <- recipient$getLeafPrior()
  unitsBefore <- units(recipient)
  donor <- do.call(make, c(args, list(response = response)))
  expect_false(identical(units(donor), unitsBefore), info = label)
  # a converted state is not the stored one
  expect_false(recipient$setState(stored(donor)), info = label)
  expect_identical(recipient$getLeafPrior(), priorBefore, info = label)
  expect_identical(units(recipient), unitsBefore, info = label)
}
acrossUnits("k-named", list())
acrossUnits("sd-named", list(leaf.prior = normal(sd = 1)))
acrossUnits("drawn k", list(leaf.prior = normal(k = chi(1.5, 2))))
acrossUnits("sd law", list(leaf.prior = normal(sd = invchi(3, 1))))
acrossUnits("linear", list(leaf.prior = linear("x2")))
acrossUnits("monotone", list(monotone = c(1L, 0L, 0L)))
acrossUnits("gp", list(leaf.prior = gp("x2")), stretched, centred)
acrossUnits("variance forest", list(variance = TRUE))
acrossUnits("count", list(family = nbinom()), 2L * counts, counts)

# a constant response's transform is units too, the window of width 1
# centred on its value, recorded as (c, c): a
# sampler on one keeps them across an install from a normal response, a state
# stored in them replays in a normal recipient, and its own state is untouched
makeConstant <- function(...) {
  sampler <- NULL
  expect_warning(
    sampler <- make(response = rep(2, n), ...),
    "indistinguishable"
  )
  sampler
}
constantRecipient <- makeConstant()
constantPrior <- constantRecipient$getLeafPrior()
constantRecipient$setState(stored(make()))
expect_identical(constantRecipient$getLeafPrior(), constantPrior)
expect_identical(units(constantRecipient), c(2, 2))
keepConstant <- stateControl(keepTrees = TRUE)
constantDonor <- makeConstant(control = keepConstant)
invisible(constantDonor$run(0L, 3L))
normalRecipient <- make(control = keepConstant)
normalRecipient$setState(stored(constantDonor))
expect_equal(
  normalRecipient$predict(x),
  constantDonor$predict(x),
  tolerance = 1e-12
)
constantSelf <- makeConstant()
constantState <- stored(constantSelf)
constantSelf$setState(constantState)
constantTwin <- makeConstant()
constantTwin$setState(constantState)
expect_identical(sweeps(constantSelf), sweeps(constantTwin))
# a stored (c, c) read back is the window centred on c: a sampler made on a
# constant and given a varying response under the pinned transform still
# records (2, 2), and a reload, which re-creates on the varying response and
# installs the stored pair, reports it, shifts by 2 and continues as the
# sampler it was saved from (a pair read as c to c + 1 put the shift at 2.5)
pinnedConstant <- makeConstant()
pinnedConstant$setResponse(y, updateScale = FALSE)
invisible(pinnedConstant$run(5L, 0L))
expect_identical(units(pinnedConstant), c(2, 2))
expect_identical(pinnedConstant$getLeafPrior()$response.shift, 2)
pinnedFits <- pinnedConstant$getFitsWithoutOffset()
pinnedFile <- tempfile(fileext = ".rds")
saveRDS(pinnedConstant, pinnedFile)
pinnedReloaded <- readRDS(pinnedFile)
unlink(pinnedFile)
expect_identical(units(pinnedReloaded), c(2, 2))
expect_identical(pinnedReloaded$getLeafPrior()$response.shift, 2)
expect_equal(
  pinnedReloaded$getFitsWithoutOffset(),
  pinnedFits,
  tolerance = 1e-12
)
# the same streams, to the last digits a state install keeps
expect_equal(
  sweeps(pinnedReloaded),
  sweeps(pinnedConstant),
  tolerance = 1e-12
)

# a record that is not a transform is refused when a re-creation reads it: a
# non-finite or decreasing pair anywhere, and an equal one on the count family
withRecord <- function(sampler, record) {
  model <- sampler$model
  attr(model, "response.range") <- record
  sampler$model <- model
  sampler
}
for (record in list(c(NaN, NaN), c(2, 1))) {
  expect_error(
    withRecord(make(), record)$copy(),
    "response.range record must be two finite numbers"
  )
}
expect_error(
  withRecord(make(response = counts, family = nbinom()), c(1, 1))$copy(),
  "the second above the first"
)

# prior-only draws are at the recipient's anchor after it, not the donor's
priorRecipient <- make()
priorRecipient$setState(stored(make(response = rescaled)))
lp <- priorRecipient$getLeafPrior()
set.seed(9)
priorDraws <- vapply(
  seq_len(300L),
  function(draw) {
    priorRecipient$sampleTreesFromPrior(updateState = FALSE)
    priorRecipient$sampleLeafParametersFromPrior(updateState = FALSE)
    priorRecipient$predict(x[1L, , drop = FALSE])[[1L]]
  },
  numeric(1L)
)
spread <- lp$k.scale / priorRecipient$getK()[[1L]]
expect_true(abs(mean(priorDraws) - lp$prior.mean) < 4 * spread / sqrt(300))
expect_true(abs(sd(priorDraws) / spread - 1) < 0.15)

# a converted state replays the function its donor held, kept draws and live
# fit both, to rounding; a gp leaf's replay solves its kernel system again, so
# its rounding is the solve's
replays <- function(
  label,
  args,
  response = rescaled,
  base = y,
  tolerance = 1e-12
) {
  control <- stateControl(keepTrees = TRUE)
  donor <- do.call(make, c(args, list(response = response, control = control)))
  invisible(donor$run(0L, 4L))
  recipient <- do.call(make, c(args, list(response = base, control = control)))
  recipient$setState(stored(donor))
  expect_false(identical(units(donor), units(recipient)), info = label)
  expect_equal(
    recipient$predict(x),
    donor$predict(x),
    tolerance = tolerance,
    info = label
  )
  expect_equal(
    recipient$getFitsWithoutOffset(),
    donor$getFitsWithoutOffset(),
    tolerance = tolerance,
    info = label
  )
  recipient
}
invisible(replays("constant", list()))
invisible(replays("monotone", list(monotone = c(1L, 0L, 0L))))
invisible(replays("linear", list(leaf.prior = linear("x2"))))
invisible(replays(
  "gp, equal shift",
  list(leaf.prior = gp("x2")),
  stretched,
  centred,
  tolerance = 1e-11
))
invisible(replays("count", list(family = nbinom()), 2L * counts, counts))
varianceDonor <- make(response = rescaled, variance = TRUE)
varianceRecipient <- make(variance = TRUE)
varianceRecipient$setState(stored(varianceDonor))
expect_equal(
  varianceRecipient$getVariance(),
  varianceDonor$getVariance(),
  tolerance = 1e-12
)

# the two refusals: a gp leaf and forests with amplitudes cannot carry another
# shift, and the sampler is left as it was
gpRecipient <- make(leaf.prior = gp("x2"))
gpBefore <- stored(gpRecipient)
expect_error(
  gpRecipient$setState(stored(make(
    leaf.prior = gp("x2"),
    response = rescaled
  ))),
  "response shift cannot be converted"
)
expect_identical(stored(gpRecipient), gpBefore)
twoRecipient <- make(forests = list(forest(), forest(basis = ~ factor(z))))
twoBefore <- stored(twoRecipient)
expect_error(
  twoRecipient$setState(stored(make(
    forests = list(forest(), forest(basis = ~ factor(z))),
    response = rescaled
  ))),
  "response shift cannot be converted"
)
expect_identical(stored(twoRecipient), twoBefore)

# a supplied gp lengthscale is the sampler's: saved draws made under another
# are refused, and without saved draws the state installs and the sampler
# keeps its own
keep <- stateControl(keepTrees = TRUE)
shortScale <- make(leaf.prior = gp("x2", lengthscale = 0.5), control = keep)
invisible(shortScale$run(0L, 4L))
longScale <- make(leaf.prior = gp("x2", lengthscale = 2), control = keep)
longBefore <- stored(longScale)
expect_error(
  longScale$setState(stored(shortScale)),
  "other lengthscales than"
)
expect_identical(stored(longScale), longBefore)
liveOnly <- make(leaf.prior = gp("x2", lengthscale = 2))
liveOnly$setState(stored(make(leaf.prior = gp("x2", lengthscale = 0.5))))
expect_identical(stored(liveOnly)[[1L]]$forests[[1L]]$leaf.lengthscales, 2)

# chains from two samplers on different ranges, combined into one state, are
# installed in one set of units and run as one posterior: weak data, so a
# chain left in its own units would run under its own prior and show it
set.seed(4L)
nWeak <- 25L
xWeak <- matrix(runif(nWeak * 2L), nWeak, 2L)
yWeak <- 2 * xWeak[, 1L] + rnorm(nWeak, sd = 3)
weak <- function() {
  set.seed(5L)
  dbarts(
    xWeak,
    yWeak,
    control = dbartsControl(
      n.chains = 2L,
      n.threads = 1L,
      n.trees = 15L,
      n.samples = 1500L,
      n.burn = 500L,
      updateState = FALSE,
      seed = 71L
    )
  )
}
own <- weak()
invisible(own$run(2L, 2L))
other <- weak()
other$setResponse(3 * yWeak + 10, updateScale = TRUE)
invisible(other$run(2L, 2L))
other$setResponse(yWeak, updateScale = FALSE)
invisible(other$run(20L, 2L))
combined <- stored(own)
combined[[2L]] <- stored(other)[[2L]]
expect_false(identical(combined[[1L]]$fit.scale, combined[[2L]]$fit.scale))
pooled <- weak()
pooled$setState(combined)
expect_identical(pooled$getLeafPrior(), own$getLeafPrior())
reference <- weak()
reference$setState(stored(own))
pooledRun <- pooled$run(500L, 1500L)
referenceRun <- reference$run(500L, 1500L)
rowSpread <- function(run, chain) sd(rowMeans(run$train[,, chain]))
expect_true(
  abs(mean(pooledRun$train[,, 2L]) - mean(referenceRun$train[,, 2L])) < 0.05
)
expect_true(
  abs(rowSpread(pooledRun, 2L) / rowSpread(referenceRun, 2L) - 1) < 0.3
)

# a re-creation after a swap without the scale update: with its state, a
# reload is bitwise the live sampler's own restore; without one, a copy runs
# under the saver's prior, anchored where its model records
frozen <- make(leaf.prior = normal(sd = 1))
frozen$setResponse(rescaled, updateScale = FALSE)
invisible(frozen$run(0L, 2L))
frozenState <- stored(frozen)
reloaded <- unserialize(serialize(frozen, NULL))
expect_identical(sweeps(reloaded), sweeps(frozen$copy()))
expect_identical(reloaded$getLeafPrior(), frozen$getLeafPrior())
stateless <- make()
stateless$setResponse(rescaled, updateScale = FALSE)
expect_null(stateless$state)
copied <- stateless$copy()
expect_identical(copied$getLeafPrior(), stateless$getLeafPrior())
expect_identical(units(copied), units(stateless))

# a copy and a reload whose stored state predates a re-anchor: the state is
# converted into the re-anchored units, and its function is kept
anchored <- make()
anchored$storeState()
heldFits <- anchored$getFitsWithoutOffset()
heldShift <- anchored$getLeafPrior()$response.shift
anchored$setResponse(rescaled, updateScale = TRUE)
expect_true(anchored$getLeafPrior()$response.shift != heldShift)
reanchoredCopy <- anchored$copy()
expect_identical(reanchoredCopy$getLeafPrior(), anchored$getLeafPrior())
expect_equal(reanchoredCopy$getFitsWithoutOffset(), heldFits, tolerance = 1e-12)
reanchoredReload <- unserialize(serialize(anchored, NULL))
expect_identical(reanchoredReload$getLeafPrior(), anchored$getLeafPrior())
expect_equal(
  reanchoredReload$getFitsWithoutOffset(),
  heldFits,
  tolerance = 1e-12
)

# the rollback: a re-anchor is a model change a restore does not undo, so a
# rejected re-anchoring proposal is rolled back by re-anchoring to the old
# response and then restoring, which leaves reader and draws as a twin's that
# never proposed
proposer <- make()
untouched <- make()
saved <- stored(proposer)
proposer$setResponse(rescaled, updateScale = TRUE)
invisible(proposer$run(0L, 2L))
expect_false(proposer$setState(saved))
expect_false(identical(proposer$getLeafPrior(), untouched$getLeafPrior()))
proposer$setResponse(y, updateScale = TRUE)
expect_true(proposer$setState(saved))
expect_true(untouched$setState(saved))
expect_identical(proposer$getLeafPrior(), untouched$getLeafPrior())
# the two re-anchors carry the sigma prior out and back, which need not
# return to its last bit
expect_equal(sweeps(proposer), sweeps(untouched), tolerance = 1e-8)

# a warm start from a donor on another range seeds the function the donor
# held, in the recipient's own units
warmDonor <- make(response = rescaled)
warmed <- make()
warmBefore <- warmed$getLeafPrior()
warmed$installTrees(warmDonor)
expect_identical(warmed$getLeafPrior(), warmBefore)
expect_equal(
  warmed$getFitsWithoutOffset(),
  warmDonor$getFitsWithoutOffset(),
  tolerance = 1e-12
)
