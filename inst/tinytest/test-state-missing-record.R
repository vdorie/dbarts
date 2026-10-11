# A stored state records which predictor columns could hold a missing value
# when it was stored. setState, copy() and a reload that put it into a
# sampler where a further column can draw the side of every rule on that
# column, in the current trees and the kept draws, as the column's first
# missing value would have, before the state is judged. The draw itself -
# its coins, their order and the generator they come from - is held in
# tests/cpp (testStateMissingRecord).

set.seed(20261010L)
n <- 150L
complete <- data.frame(
  x1 = runif(n),
  x2 = runif(n),
  f = factor(sample(letters[1:4], n, replace = TRUE))
)
y <- 2 *
  (complete$x1 > 0.5) +
  complete$x2 +
  as.integer(complete$f) / 2 +
  rnorm(n, sd = 0.1)
controlOf <- function(n.chains = 2L, n.trees = 15L, n.samples = 20L) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = 1L,
    n.trees = n.trees,
    n.burn = 30L,
    n.samples = n.samples,
    keepTrees = TRUE,
    updateState = FALSE,
    verbose = FALSE
  )
}
make <- function(data = complete, seed = 11L, control = controlOf(), ...) {
  sampler <- dbarts(y ~ x1 + x2 + f, data, control = control, seed = seed, ...)
  invisible(sampler$run())
  sampler
}
seen <- function(sampler) dbarts:::dataMissingSeen(sampler$data)
stored <- function(sampler) {
  sampler$storeState()
  sampler$state
}
recorded <- function(state) attr(state, "missing.columns")
generators <- function(sampler) lapply(stored(sampler), `[[`, "rng.state")
# the rules, kept and current, with the side each sends a missing value to
rules <- function(sampler) {
  lapply(c(kept = FALSE, current = TRUE), function(current) {
    as.list(sampler$getTrees(current = current))
  })
}
# the sides of the rules on a column, of the kept draws or the current trees
sides <- function(sampler, current = FALSE, column = 1L) {
  trees <- sampler$getTrees(current = current)
  trees$missing[trees$var == column]
}
bothSides <- function(side) all(c("L", "R") %in% side)
visibly <- function(value) list(value = value, visible = TRUE)
gone <- c(3L, 40L, 77L)
holed <- complete
holed$x1[gone] <- NA
holedF <- complete
holedF$f[gone] <- NA
asked <- complete[1:2, ]
asked$x1[1L] <- NA
inconsistent <- "state is not consistent with this sampler"

# ---- the record -----------------------------------------------------------

sampler <- make()
before <- stored(sampler)
expect_identical(recorded(before), c(FALSE, FALSE, FALSE))
paths <- list(
  "forced column" = function(s) {
    s$setPredictor(holed$x1, "x1", forceUpdate = TRUE)
  },
  "unforced column" = function(s) s$setPredictor(holed$x1, "x1"),
  "whole frame" = function(s) s$setPredictor(holed, forceUpdate = FALSE),
  "row by row" = function(s) {
    s$setPredictor(holed$x1, "x1", forceUpdate = "partial")
  },
  "jointly" = function(s) {
    updatePredictorPerObservationJointly(list(s), holed$x1, "x1")
  },
  "setData" = function(s) s$setData(dbartsData(y ~ x1 + x2 + f, holed))
)
for (path in names(paths)) {
  sampler <- make()
  paths[[path]](sampler)
  expect_identical(
    recorded(stored(sampler)),
    c(TRUE, FALSE, FALSE),
    info = path
  )
  # the flag and not the content: filled, the column can still hold one
  sampler$setPredictor(complete$x1, "x1", forceUpdate = TRUE)
  expect_identical(
    recorded(stored(sampler)),
    c(TRUE, FALSE, FALSE),
    info = path
  )
}
path <- tempfile(fileext = ".rds")
saveRDS(sampler$state, path)
expect_identical(recorded(readRDS(path)), c(TRUE, FALSE, FALSE))
unlink(path)
# a record of another length is another sampler's; of another type, or
# holding NA, it is malformed: an error with and without force
malformed <- "malformed missing-value columns in bartcore state"
for (bad in list(c(TRUE, FALSE), c(1, 0, 0), c(TRUE, NA, FALSE), "x1")) {
  state <- before
  attr(state, "missing.columns") <- bad
  refusal <- if (is.logical(bad) && !anyNA(bad)) inconsistent else malformed
  for (force in c(FALSE, TRUE)) {
    expect_error(sampler$setState(state, forceUpdate = force), refusal)
  }
}

# ---- a state from before the first missing value, installed after it ------

# The update is unforced and taken. Installing the earlier state then is
# clean, and leaves what the update left: every rule and its side, kept and
# current, the generators, and predict's answer for a missing value. A
# second update draws nothing.
more <- holed$x1
more[c(5L, 90L)] <- NA
reviewCase <- function(sampler, info, answers = TRUE) {
  earlier <- stored(sampler)
  expect_true(sampler$setPredictor(holed$x1, "x1"), info = info)
  between <- list(rules(sampler), generators(sampler))
  expect_true(
    all(vapply(
      between[[1L]],
      function(trees) {
        bothSides(trees$missing[trees$var == 1L])
      },
      NA
    )),
    info = info
  )
  if (answers) {
    answer <- sampler$predict(asked)
  }
  expect_identical(withVisible(sampler$setState(earlier)), visibly(TRUE), info)
  expect_identical(list(rules(sampler), generators(sampler)), between, info)
  if (answers) {
    expect_identical(sampler$predict(asked), answer, info = info)
  }
  expect_true(sampler$setPredictor(more, "x1"), info = info)
  expect_identical(
    list(rules(sampler)$kept$missing, generators(sampler)),
    list(between[[1L]]$kept$missing, between[[2L]]),
    info = info
  )
}
reviewCase(make(), "one forest")
z <- as.double(seq_len(n) %% 2L)
reviewCase(
  make(forests = list(forest(), forest(basis = z))),
  "two forests",
  answers = FALSE
)
# a ring of kept draws that has wrapped, with rules on x1 on both sides of
# its seam: draws 1 to 13 lie in slots 7 to 19
wrapped <- make()
invisible(wrapped$run(0L, 7L))
expect_identical(attr(stored(wrapped), "currentSampleNum"), 7L)
trees <- wrapped$getTrees()
expect_true(all(c(TRUE, FALSE) %in% (trees$sample[trees$var == 1L] <= 13L)))
reviewCase(wrapped, "a wrapped ring")

# 20 kept draws into a store of 5: every one draws, and the newest 5 go in
# as drawn
source <- make()
earlier <- stored(source)
expect_true(source$setPredictor(holed$x1, "x1"))
small <- make(control = controlOf(n.samples = 5L))
expect_true(small$setPredictor(holed$x1, "x1"))
expect_true(small$setState(earlier))
trees <- source$getTrees()
newest <- as.list(trees[trees$sample > 15L, ])
newest$sample <- newest$sample - 15L
expect_identical(as.list(small$getTrees()), newest)
expect_identical(rules(small)$current, rules(source)$current)
expect_identical(generators(small), generators(source))

# a copy and a reload install the cached state, stored before the value,
# into an engine re-created with the data object's record: the source's
# sides, and its next draws to the rounding a re-creation leaves
expect_identical(recorded(source$state), c(TRUE, FALSE, FALSE))
source$state <- earlier
expect_identical(seen(source), c(TRUE, FALSE, FALSE))
duplicate <- source$copy()
path <- tempfile(fileext = ".rds")
saveRDS(source, path)
reloaded <- readRDS(path)
unlink(path)
held <- rules(source)
draws <- source$run(0L, 5L)
for (route in list(duplicate, reloaded)) {
  expect_identical(
    lapply(rules(route), `[[`, "missing"),
    lapply(held, `[[`, "missing")
  )
  expect_equal(route$run(0L, 5L), draws, tolerance = 1e-12)
}

# ---- nothing to draw ------------------------------------------------------

# a state stored after the value keeps its sides and its generators
sampler <- make()
expect_true(sampler$setPredictor(holed$x1, "x1"))
after <- stored(sampler)
invisible(sampler$run(0L, 10L))
later <- stored(sampler)
held <- list(rules(sampler), generators(sampler))
invisible(sampler$run(0L, 3L))
expect_true(sampler$setState(later))
expect_identical(list(rules(sampler), generators(sampler)), held)
# a state with no record, as one stored before the record was, installs as
# it is: every side left, and its own generators
unrecorded <- before
attr(unrecorded, "missing.columns") <- NULL
expect_true(sampler$setState(unrecorded))
expect_true(all(sides(sampler) == "L") && all(sides(sampler, TRUE) == "L"))
expect_identical(generators(sampler), lapply(before, `[[`, "rng.state"))
# a state that holds sides on x1 does not fit a sampler whose x1 never could
never <- make()
expect_identical(withVisible(never$setState(after)), visibly(FALSE))
expect_null(seen(never))
expect_identical(recorded(stored(never)), c(FALSE, FALSE, FALSE))

# ---- a factor column ------------------------------------------------------

# One tree of a complete sampler's state: f at the root, a rule on f on
# each side. A missing value reaches one of the two, by the root's side, and
# the other holds none. A level mask is 64 bits in the machine's byte order.
maskBytes <- function(bits) {
  words <- c(as.integer(bits), 0L)
  writeBin(if (.Platform$endian == "big") rev(words) else words, raw())
}
nested <- stored(make())
forest <- nested[[1L]]$forests[[1L]]
forest$tree.vars <- rep(c(3L, 3L, -1L, -1L, 3L, -1L, -1L), 15L)
forest$tree.values <- rep(
  c(
    maskBytes(3L),
    maskBytes(4L),
    writeBin(c(-0.01, 0.01), raw()),
    maskBytes(1L),
    writeBin(c(0.02, -0.02), raw())
  ),
  15L
)
forest$tree.sizes <- rep(7L, 15L)
forest$tree.flags <- as.raw(rep(c(4L, 4L, 0L, 0L, 4L, 0L, 0L), 15L))
nested[[1L]]$forests[[1L]] <- forest
created <- make(holedF)
status <- tryCatch(created$setState(nested), error = conditionMessage)
expect_identical(status, TRUE)
if (isTRUE(status)) {
  trees <- created$getTrees(chainNums = 1L, current = TRUE)
  side <- matrix(trees$missing, 7L)
  expect_true(bothSides(side[1L, ]))
  # the rule on the side the root sends no missing value down
  expect_true(all(side[cbind(ifelse(side[1L, ] == "R", 2L, 5L), 1:15)] == "L"))
  expect_true(bothSides(sides(created, column = 3L)))
}

# Kept rules that send no level right and a missing value right, which a
# sampler whose f has never held one can hold only in kept draws taken from
# another's state. They keep their side through f's first missing value and
# through an install that draws, where a side drawn left would leave a rule
# that sends nothing right, and the sampler could not install its own state.
aloneRight <- function(sampler) {
  lapply(stored(sampler), function(chain) {
    forest <- chain$forests[[1L]]
    noLevel <- colSums(matrix(forest$saved.values != as.raw(0L), 8L)) == 0L
    # 5 tags a rule on a factor that sends a missing value right
    which(forest$saved.vars == 3L & forest$saved.flags == as.raw(5L) & noLevel)
  })
}
dataF <- dbartsData(y ~ x1 + x2 + f, holedF)
found <- FALSE
for (seed in 1:40) {
  donor <- stored(make(holedF, seed = seed))
  holder <- make()
  taken <- tryCatch(
    is.null(holder$setState(donor, forceUpdate = TRUE)),
    error = function(e) FALSE
  )
  if (taken && length(unlist(aloneRight(holder))) > 0L) {
    found <- TRUE
    break
  }
}
expect_true(found)
alone <- aloneRight(holder)
held <- holder$state
expect_identical(recorded(held), c(FALSE, FALSE, FALSE))
holder$setData(dataF)
expect_identical(seen(holder), c(FALSE, FALSE, TRUE))
expect_identical(aloneRight(holder), alone)
expect_true(bothSides(sides(holder, column = 3L)))
reinstalls <- function(sampler) {
  tryCatch(sampler$setState(stored(sampler)), error = conditionMessage)
}
expect_identical(reinstalls(holder), TRUE)
# the state from before the value, installed after it
status <- tryCatch(holder$setState(held), error = conditionMessage)
expect_true(is.logical(status))
expect_identical(aloneRight(holder), alone)
expect_true(bothSides(sides(holder, column = 3L)))
expect_identical(reinstalls(holder), TRUE)

# ---- a drawn side that leaves the state not clean, or not a state ---------

# A declined state leaves the sampler the twin that made no call: its state
# field, its trees, the state its engine stores and its next 20 draws.
# Forced, it goes in.
expectDeclined <- function(make, state, info) {
  sampler <- make()
  twin <- make()
  field <- sampler$state
  expect_identical(withVisible(sampler$setState(state)), visibly(FALSE), info)
  expect_identical(sampler$state, field, info = info)
  expect_identical(rules(sampler)$current, rules(twin)$current, info = info)
  expect_identical(stored(sampler), stored(twin), info = info)
  expect_identical(sampler$run(0L, 20L), twin$run(0L, 20L), info = info)
  forced <- make()
  expect_identical(
    withVisible(forced$setState(state, forceUpdate = TRUE)),
    list(value = NULL, visible = FALSE),
    info = info
  )
  forced
}
single <- controlOf(1L, 1L)
single@keepTrees <- FALSE
# One tree, cut on x1 below every value but the first row's, which is alone
# in the left leaf, in a state stored while x1 could hold no missing value.
# The sampler is then given one in that row, forced. The state's side for it
# is the coin that update took: sent right the left leaf is empty, the update
# merged it and the state is declined; sent left it is clean.
lone <- complete
lone$x1[1L] <- -5
brought <- replace(lone$x1, 1L, NA)
loneLeft <- function(seed, bring = TRUE) {
  sampler <- dbarts(y ~ x1 + x2 + f, lone, control = single, seed = seed)
  state <- stored(sampler)
  forest <- state[[1L]]$forests[[1L]]
  forest$tree.vars <- c(1L, -1L, -1L)
  cut <- attr(state, "cutPoints")[[1L]][1L]
  forest$tree.values <- writeBin(c(cut, -0.1, 0.1), raw())
  forest$tree.sizes <- 3L
  forest$tree.flags <- as.raw(c(2L, 0L, 0L))
  state[[1L]]$forests[[1L]] <- forest
  stopifnot(isTRUE(sampler$setState(state)))
  if (bring) {
    sampler$setPredictor(brought, "x1", forceUpdate = TRUE)
  }
  sampler
}
met <- c(declined = FALSE, clean = FALSE)
for (seed in 1:40) {
  state <- loneLeft(seed, FALSE)$state
  merged <- nrow(loneLeft(seed)$getTrees(current = TRUE)) == 1L
  arm <- if (merged) "declined" else "clean"
  if (met[[arm]]) {
    next
  }
  met[[arm]] <- TRUE
  if (merged) {
    forced <- expectDeclined(
      function() loneLeft(seed),
      state,
      "an emptied leaf"
    )
    expect_identical(nrow(forced$getTrees(current = TRUE)), 1L)
  } else {
    sampler <- loneLeft(seed)
    expect_identical(withVisible(sampler$setState(state)), visibly(TRUE))
    expect_identical(sides(sampler, TRUE), "L")
  }
  if (all(met)) {
    break
  }
}
expect_identical(met, c(declined = TRUE, clean = TRUE))

# A state edited to contradict its record: f recorded as never missing, and
# beneath a rule on f one that sends a missing value alone to the right. The
# draw leaves that rule its side. Where the root's coin sends missing values
# its way the state is one the sampler could have reached; where it sends
# them the other way the rule splits nothing, and the state is refused, with
# and without force, the sampler left as it was.
contradicting <- function(seed) {
  sampler <- dbarts(y ~ x1 + x2 + f, holedF, control = single, seed = seed)
  state <- stored(sampler)
  forest <- state[[1L]]$forests[[1L]]
  forest$tree.vars <- c(3L, -1L, 3L, -1L, -1L)
  forest$tree.values <- c(
    maskBytes(3L),
    writeBin(-0.1, raw()),
    maskBytes(0L),
    writeBin(c(0.1, 0.2), raw())
  )
  forest$tree.sizes <- 5L
  forest$tree.flags <- as.raw(c(4L, 0L, 5L, 0L, 0L))
  state[[1L]]$forests[[1L]] <- forest
  attr(state, "missing.columns") <- c(FALSE, FALSE, FALSE)
  list(sampler = sampler, state = state)
}
met <- c(refused = FALSE, clean = FALSE)
for (seed in 1:40) {
  case <- contradicting(seed)
  twin <- contradicting(seed)$sampler
  status <- tryCatch(
    case$sampler$setState(case$state),
    error = conditionMessage
  )
  arm <- if (isTRUE(status)) "clean" else "refused"
  if (met[[arm]]) {
    next
  }
  met[[arm]] <- TRUE
  if (arm == "clean") {
    expect_identical(
      case$sampler$getTrees(current = TRUE)$missing[c(1L, 3L)],
      c("R", "R")
    )
  } else {
    expect_identical(status, inconsistent)
    expect_error(
      case$sampler$setState(case$state, forceUpdate = TRUE),
      inconsistent
    )
    expect_identical(stored(case$sampler), stored(twin))
    expect_identical(case$sampler$run(0L, 20L), twin$run(0L, 20L))
  }
  if (all(met)) {
    break
  }
}
expect_identical(met, c(refused = TRUE, clean = TRUE))

# ---- a warm start ---------------------------------------------------------

# A donor fitted where x1 could hold no missing value, its trees installed
# where it can: each receiving chain draws the side of every rule on x1 from
# its own generator, before its first sweep.
receiving <- function(seed, data = holed) {
  dbarts(y ~ x1 + x2 + f, data, control = controlOf(), seed = seed)
}
donor <- make()
recipient <- receiving(3L)
recipient$installTrees(donor)
expect_true(bothSides(sides(recipient, TRUE)))
# fair coins: 40 receiving seeds, a binomial band of 1e-3 around one half
right <- total <- 0L
for (seed in 1:40) {
  recipient <- receiving(seed)
  recipient$installTrees(donor)
  right <- right + sum(sides(recipient, TRUE) == "R")
  total <- total + length(sides(recipient, TRUE))
}
expect_true(total >= 200L)
expect_true(abs(right - total / 2) < qnorm(1 - 5e-4) * sqrt(total) / 2)
# two chains given one draw of the donor hold its rules and their own sides
recipient <- receiving(3L)
recipient$installTrees(donor, samples = c(1L, 1L))
trees <- recipient$getTrees(current = TRUE)
byChain <- split(trees[c("tree", "var", "value", "missing")], trees$chain)
expect_identical(lapply(byChain[[1L]][1:3], c), lapply(byChain[[2L]][1:3], c))
expect_false(identical(byChain[[1L]]$missing, byChain[[2L]]$missing))
# nothing to draw: a donor state with no record leaves every side left, a
# donor whose x1 could hold a missing value leaves the sides it learned, and
# neither moves a generator
learned <- make(holed)
learnedTrees <- learned$getTrees(chainNums = 1L, sampleNums = 1L)
unrecorded <- stored(donor)
attr(unrecorded, "missing.columns") <- NULL
for (case in list(
  list(donor = unrecorded, sides = "L"),
  list(donor = learned, sides = learnedTrees$missing[learnedTrees$var == 1L])
)) {
  recipient <- receiving(3L)
  before <- generators(recipient)
  recipient$installTrees(case$donor, samples = c(1L, 1L))
  expect_identical(generators(recipient), before)
  drawn <- split(
    sides(recipient, TRUE),
    rep(1:2, each = length(sides(recipient, TRUE)) / 2L)
  )
  # each chain holds rules on x1
  expect_identical(length(drawn), 2L)
  for (chain in drawn) {
    expect_true(all(chain == case$sides))
  }
}
# a donor state whose record is not one is refused, where a state with no
# record is taken
for (bad in list(c(1, 0, 0), c(TRUE, NA, FALSE), "x1")) {
  misrecorded <- stored(donor)
  attr(misrecorded, "missing.columns") <- bad
  expect_error(
    receiving(3L)$installTrees(misrecorded),
    "malformed missing-value columns in warm-start donor",
    fixed = TRUE
  )
}

# The draw comes before the install merges what the rows leave empty. A
# donor of one tree, cut on x1 above every value the receiving sampler holds:
# the one row above it there is missing here. Sent right, the missing row
# keeps the right leaf and the split stands; sent left, the split is merged.
# On the donor's cut points and on others, the receiving seed searched for
# each side.
aboveAll <- function(top, seed, bring) {
  data <- complete
  data$x1[1L] <- top
  sampler <- dbarts(y ~ x1 + x2 + f, data, control = single, seed = seed)
  if (bring) {
    sampler$setPredictor(replace(data$x1, 1L, NA), "x1", forceUpdate = TRUE)
  }
  sampler
}
handDonor <- stored(aboveAll(5, 1L, FALSE))
forest <- handDonor[[1L]]$forests[[1L]]
forest$tree.vars <- c(1L, -1L, -1L)
cuts <- attr(handDonor, "cutPoints")[[1L]]
forest$tree.values <- writeBin(c(cuts[length(cuts)], -0.1, 0.1), raw())
forest$tree.sizes <- 3L
forest$tree.flags <- as.raw(c(2L, 0L, 0L))
handDonor[[1L]]$forests[[1L]] <- forest
for (top in c(5, 6)) {
  met <- c(stands = FALSE, merged = FALSE)
  for (seed in 1:40) {
    recipient <- aboveAll(top, seed, TRUE)
    recipient$installTrees(handDonor)
    trees <- recipient$getTrees(current = TRUE)
    arm <- if (nrow(trees) == 3L) "stands" else "merged"
    if (met[[arm]]) {
      next
    }
    met[[arm]] <- TRUE
    if (arm == "stands") {
      expect_identical(trees$missing[1L], "R", info = top)
      expect_identical(trees$n, c(n, n - 1L, 1L), info = top)
    } else {
      expect_identical(trees$n, n, info = top)
    }
    if (all(met)) {
      break
    }
  }
  expect_identical(met, c(stands = TRUE, merged = TRUE), info = top)
}
