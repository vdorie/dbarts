# n.perturb.cuts, the most cut positions a perturb proposal moves a split
# either way. The default, 1, is the window the engine always had; a whole
# number up to the cut cap is taken, Inf is the window over the node's whole
# range, and a finite value above the cap is refused, naming both.

set.seed(31L)
n <- 80L
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, c("a", "b", "c")))
y <- x[, 1L] - x[, 3L] + rnorm(n, 0, 0.2)
perturbing <- c(birth_death = 0.5, change = 0.2, perturb = 0.3)
widthControl <- function(...) {
  dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    updateState = FALSE,
    proposal.probs = perturbing,
    ...
  )
}
drawsAt <- function(control) {
  set.seed(7L)
  dbarts::dbarts(x, y, control = control)$run(20L, 20L)$train
}

# ---- the default changes no draw ---------------------------------------------

expect_identical(dbarts::dbartsControl()@n.perturb.cuts, 1)
atDefault <- drawsAt(widthControl())
expect_identical(drawsAt(widthControl(n.perturb.cuts = 1L)), atDefault)
expect_identical(drawsAt(widthControl(n.perturb.cuts = 1)), atDefault)

# ---- refusals, by name -------------------------------------------------------

for (bad in list(0, -1, 1.5, NA, NA_real_, c(1, 2))) {
  expect_error(
    dbarts::dbartsControl(n.perturb.cuts = bad),
    "'n.perturb.cuts' must be a single positive whole number or Inf",
    fixed = TRUE,
    info = deparse(bad)
  )
}
# the bridge holds the same rule for a slot written past the constructor
backstop <- dbarts::dbarts(x, y, control = widthControl())
backstop$control@n.perturb.cuts <- 0
expect_error(backstop$copy(), "'n.perturb.cuts' must be a positive whole")
rm(backstop)
# a whole integer written into the slot installs as its value
integerWidth <- widthControl()
integerWidth@n.perturb.cuts <- 3L
expect_silent(validObject(integerWidth))

# ---- the cut cap and Inf are taken, a finite value above the cap is refused -

atCap <- drawsAt(widthControl(n.perturb.cuts = 65533))
expect_false(identical(atCap, atDefault))
# Inf is the explicit window over the node's whole range, the cap's
expect_identical(drawsAt(widthControl(n.perturb.cuts = Inf)), atCap)
expect_identical(
  dbarts::dbartsControl(n.perturb.cuts = Inf)@n.perturb.cuts,
  Inf
)
expect_identical(
  dbarts::dbartsControl(n.perturb.cuts = 65533L)@n.perturb.cuts,
  65533
)
# a value over the cap is refused where it is written, the message naming the
# cap and both ways to ask for the widest window
overCap <- paste0(
  "is over the cap of 65533, the most cuts a predictor's grid can carry; ",
  "give 65533 or Inf for a window over the node's whole range"
)
for (wide in list(65534, 65534L, 70000, .Machine$integer.max + 1, 1e12)) {
  expect_error(
    dbarts::dbartsControl(n.perturb.cuts = wide),
    paste0(
      "'n.perturb.cuts' of ",
      format(wide, scientific = FALSE),
      " ",
      overCap
    ),
    fixed = TRUE,
    info = format(wide)
  )
}
# by a slot written past the constructor, at validation, at creation, and at
# setControl, whose sampler is left as it was
written <- widthControl()
written@n.perturb.cuts <- 70000
expect_error(validObject(written), overCap, fixed = TRUE)
expect_error(dbarts::dbarts(x, y, control = written), overCap, fixed = TRUE)
live <- dbarts::dbarts(x, y, control = widthControl(n.perturb.cuts = 3L))
expect_error(live$setControl(written), overCap, fixed = TRUE)
expect_identical(live$control@n.perturb.cuts, 3)
for (taken in list(65533, Inf)) {
  taking <- live$control
  taking@n.perturb.cuts <- taken
  live$setControl(taking)
  expect_identical(live$control@n.perturb.cuts, taken)
}
rm(written, live, taking, taken)
# the other refusals are as they were, a value over the cap among them only
# when finite and whole
expect_error(
  dbarts::dbartsControl(n.perturb.cuts = 70000.5),
  "'n.perturb.cuts' must be a single positive whole number or Inf",
  fixed = TRUE
)
# a value under the cap that exceeds a column's own cut count is taken, the
# setting being a maximum: the window stops at the node's range, so it is the
# whole range, as the cap is
smallGrid <- function(width) {
  set.seed(7L)
  dbarts::dbarts(
    x,
    y,
    control = widthControl(n.cuts = 5L, n.perturb.cuts = width)
  )$run(20L, 20L)$train
}
expect_identical(smallGrid(60000), smallGrid(5L))
expect_identical(smallGrid(5L), smallGrid(65533))
expect_identical(smallGrid(5L), smallGrid(Inf))
# the bridge still clamps a slot that reaches it unchecked, as a backstop
copiedAt <- function(width) {
  written <- dbarts::dbarts(x, y, control = widthControl())
  written$control@n.perturb.cuts <- width
  set.seed(7L)
  written$copy()$run(20L, 20L)$train
}
expect_identical(copiedAt(70000), copiedAt(65533))
expect_identical(copiedAt(Inf), copiedAt(65533))
expect_false(identical(copiedAt(70000), copiedAt(1)))

# ---- width 3 reaches every sampler kind --------------------------------------

atThree <- drawsAt(widthControl(n.perturb.cuts = 3L))
expect_false(identical(atThree, atDefault))
expect_identical(drawsAt(integerWidth), atThree)

z <- rbinom(n, 1L, 0.5)
bcfForests <- list(
  dbarts::dbartsForests$forest(),
  dbarts::dbartsForests$forest(basis = ~ factor(z))
)
bcfDraws <- function(width) {
  set.seed(7L)
  dbarts::dbarts(
    x,
    y + z,
    forests = bcfForests,
    control = widthControl(n.perturb.cuts = width)
  )$run(20L, 20L)$train
}
bcfAtThree <- bcfDraws(3L)
expect_true(all(is.finite(bcfAtThree)))
expect_false(identical(bcfAtThree, bcfDraws(1L)))

counts <- matrix(0L, n, 3L)
counts[cbind(seq_len(n), 1L + (x[, 1L] > 0.5) + (x[, 2L] > 0.5))] <- 1L
multinomialDraws <- function(width) {
  set.seed(7L)
  dbarts::dbarts(
    dbarts::dbartsData(x, counts = counts),
    family = "multinomial",
    control = widthControl(n.perturb.cuts = width)
  )$run(20L, 20L)$train
}
multinomialAtThree <- multinomialDraws(3L)
expect_true(all(is.finite(multinomialAtThree)))
expect_false(identical(multinomialAtThree, multinomialDraws(1L)))

varianceDraws <- function(width) {
  set.seed(7L)
  dbarts::dbarts(
    x,
    y,
    variance = varianceForest(n.trees = 5L),
    control = widthControl(n.perturb.cuts = width)
  )$run(20L, 20L)$train
}
varianceAtThree <- varianceDraws(3L)
expect_true(all(is.finite(varianceAtThree)))
expect_false(identical(varianceAtThree, varianceDraws(1L)))

# ---- kept by copy and reload -------------------------------------------------

# a copy and a reload re-create the engine from the stored control and state:
# both carry width 3, and the same reload with its control's width written
# back to 1 before re-creation draws differently, so the width is what the
# re-created engine reads
keptSampler <- dbarts::dbarts(
  x,
  y,
  control = widthControl(n.perturb.cuts = 3L, n.samples = 10L)
)
set.seed(9L)
invisible(keptSampler$run(20L, 10L))
keptSampler$storeState()
duplicate <- keptSampler$copy()
expect_identical(duplicate$control@n.perturb.cuts, 3)
reloadFile <- tempfile(fileext = ".rds")
saveRDS(keptSampler, reloadFile)
reloaded <- readRDS(reloadFile)
reloadedNarrow <- readRDS(reloadFile)
unlink(reloadFile)
expect_identical(reloaded$control@n.perturb.cuts, 3)
reloadedNarrow$control@n.perturb.cuts <- 1
set.seed(5L)
continued <- duplicate$run(0L, 10L)$train
set.seed(5L)
expect_identical(reloaded$run(0L, 10L)$train, continued)
set.seed(5L)
expect_false(identical(reloadedNarrow$run(0L, 10L)$train, continued))

# ---- taken by setControl between runs ----------------------------------------

widened <- dbarts::dbarts(x, y, control = widthControl())
narrow <- dbarts::dbarts(x, y, control = widthControl())
set.seed(3L)
invisible(widened$run(20L, 1L))
set.seed(3L)
invisible(narrow$run(20L, 1L))
widerControl <- widened$control
widerControl@n.perturb.cuts <- 3
widened$setControl(widerControl)
expect_identical(widened$control@n.perturb.cuts, 3)
set.seed(4L)
widenedDraws <- widened$run(0L, 20L)$train
set.seed(4L)
expect_false(identical(widenedDraws, narrow$run(0L, 20L)$train))

# a refused install rolls the control back: a multinomial sampler's priors are
# fixed at creation, so the install that carries the width refuses
fixedPrior <- dbarts::dbarts(
  dbarts::dbartsData(x, counts = counts),
  family = "multinomial",
  control = widthControl()
)
refusedControl <- fixedPrior$control
refusedControl@n.perturb.cuts <- 3
expect_error(fixedPrior$setControl(refusedControl), "\\$setModel")
expect_identical(fixedPrior$control@n.perturb.cuts, 1)
