# The empty-leaf veto counts MEMBERS, whatever weight they carry.
# A row of weight zero is in the design and not in the likelihood: it occupies
# its leaf, so a leaf all of whose rows carry weight zero is legal, takes no
# likelihood term and has its value drawn from the prior. The oracle is
# structural and needs no tree walk: the live trees report, per node, the rows
# that reach it (getTrees), and routing ONLY the positive-weight rows through
# them (getTrees(newdata = )) the rows the likelihood can see. The rule says
# every leaf holds a row, and that some leaves hold none the likelihood sees.
# The zero-weight region is a half-space on one predictor so that whole leaves
# can fall inside it, which is what makes the check bite.

set.seed(20260812L)
n <- 400L
x <- matrix(runif(n * 3L), n, 3L)
w <- as.double(x[, 1L] <= 0.5) # the x1 > 0.5 half-space carries no weight
kept <- w > 0
signal <- 2 * x[, 1L] - x[, 2L]
y <- signal + rnorm(n, 0, 0.5)

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE,
  seed = 7L
)

leavesOf <- function(sampler) {
  members <- sampler$getTrees(current = TRUE)
  weighted <- sampler$getTrees(
    current = TRUE,
    newdata = x[kept, , drop = FALSE]
  )
  isLeaf <- members$var == -1L
  data.frame(members = members$n[isLeaf], weighted = weighted$n[isLeaf])
}

fitAndReport <- function(weights, ...) {
  sampler <- dbarts::dbarts(x, y, weights = weights, control = control, ...)
  draws <- sampler$run(100L, 50L)
  list(leaves = leavesOf(sampler), fits = draws$train)
}

# gaussian: every live leaf holds a row, and some hold only zero-weight rows
gaussian <- fitAndReport(w)
expect_true(nrow(gaussian$leaves) > 50L)
expect_true(all(gaussian$leaves$members > 0L))
expect_true(any(gaussian$leaves$weighted == 0L))

# a zero-weight row is still partitioned and still reported a fit, which
# tracks the signal through the leaves it shares with weighted rows
expect_true(all(is.finite(gaussian$fits)))
expect_true(cor(rowMeans(gaussian$fits)[!kept], signal[!kept]) > 0.5)

# Student-t: the composed weight w * lambda is zero wherever w is, so the same
# rule governs the family that carries the other shipped weight channel
student <- fitAndReport(w, family = dbarts:::student(df = 4))
expect_true(all(student$leaves$members > 0L))
expect_true(any(student$leaves$weighted == 0L))
expect_true(all(is.finite(student$fits)))

# Weights do not ride the tree state, so a vector installed BETWEEN samples on
# an already-grown forest leaves leaves that hold only zero-weight rows. That
# is a legal state under every vector, and the forest keeps moving in it.
structureOf <- function(sampler) {
  live <- sampler$getTrees(current = TRUE)
  paste0(live$tree, ":", live$n, ":", live$var, collapse = "|")
}

# every row zeroed: every leaf scores nothing, so the forest sits at its prior -
# which is a distribution over structures, not a fixed structure
grown <- dbarts::dbarts(x, y, control = control)
invisible(grown$run(100L, 1L))
beforeInstall <- structureOf(grown)
grown$setWeights(rep(0, n))
zeroedDraws <- grown$run(50L, 10L)
expect_false(identical(structureOf(grown), beforeInstall))
expect_true(all(is.finite(zeroedDraws$train)))
expect_true(all(is.finite(zeroedDraws$sigma) & zeroedDraws$sigma > 0))
expect_true(all(leavesOf(grown)$members > 0L))

# partially zeroed: the forest goes on moving, every leaf holds a row, and
# leaves inside the zero-weight half-space stay legal
partial <- dbarts::dbarts(x, y, control = control)
invisible(partial$run(100L, 1L))
partial$setWeights(w)
beforeInstall <- structureOf(partial)
invisible(partial$run(200L, 1L))
expect_false(identical(structureOf(partial), beforeInstall))
partialLeaves <- leavesOf(partial)
expect_true(all(partialLeaves$members > 0L))
expect_true(any(partialLeaves$weighted == 0L))

# restoring positive weights needs no repair: the state was legal throughout
partial$setWeights(rep(1, n))
restored <- partial$run(50L, 10L)
expect_true(all(is.finite(restored$train)))

rm(
  n,
  x,
  w,
  kept,
  signal,
  y,
  control,
  leavesOf,
  fitAndReport,
  gaussian,
  student,
  structureOf,
  grown,
  beforeInstall,
  zeroedDraws,
  partial,
  partialLeaves,
  restored
)
