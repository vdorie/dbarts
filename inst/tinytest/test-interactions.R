# Interaction constraints: the
# interactions() surface and its fit-time validation, a constrained fit obeying
# max.order / forbid / groups, the per-forest BCF split, and warm-start refusal.
# The exactness and containment gates live in tests/cpp/test_interaction.cpp and
# test_state.cpp.

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

set.seed(11L)

# ---- the whole-tree walk over the flat getTrees representation ----------------

# maximum distinct predictors on any root-to-leaf path of one flat (pre-order,
# var = -1 marks a leaf) tree
pathMaxOrder <- function(vars) {
  i <- 0L
  rec <- function(ancestors) {
    i <<- i + 1L
    v <- vars[i]
    if (is.na(v) || v < 0) {
      return(length(ancestors))
    }
    a <- unique(c(ancestors, v))
    l <- rec(a)
    r <- rec(a)
    max(l, r)
  }
  rec(integer(0))
}
# does any path carry all of the given 1-based variable indices together?
pathHasAll <- function(vars, want) {
  i <- 0L
  bad <- FALSE
  rec <- function(ancestors) {
    i <<- i + 1L
    v <- vars[i]
    if (is.na(v) || v < 0) {
      if (all(want %in% ancestors)) {
        bad <<- TRUE
      }
      return(invisible())
    }
    a <- c(ancestors, v)
    rec(a)
    rec(a)
  }
  rec(integer(0))
  bad
}
# group a getTrees data.frame into one tree per (chain, sample, tree), using
# whichever of those columns the format carries (bart2 extract has all three,
# the BCF current-tree query only 'tree')
source(
  system.file("common", "splitTrees.R", package = "dbarts"),
  local = TRUE
)
worstOrder <- function(trees) {
  max(vapply(splitTrees(trees), function(k) pathMaxOrder(k$var), integer(1)))
}
anyCoOccur <- function(trees, want) {
  any(vapply(
    splitTrees(trees),
    function(k) pathHasAll(k$var, want),
    logical(1)
  ))
}

# an interactive signal so an unconstrained fit reaches for several variables
n <- 200L
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- ifelse(x1 > 0.5 & x2 > 0.5, 2, -1) + 0.7 * x3 + rnorm(n, 0, 0.3)
df <- data.frame(y, x1, x2, x3)

fitArgs <- list(
  n.trees = 25L,
  n.samples = 15L,
  n.burn = 40L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
doFitInteractions <- function(interactions = NULL) {
  do.call(
    bart2,
    c(list(y ~ x1 + x2 + x3, df, interactions = interactions), fitArgs)
  )
}

# ---- unconstrained baseline actually uses multiple variables ------------------

expect_true(worstOrder(extract(doFitInteractions(), type = "trees")) >= 2L)

# ---- max.order caps the distinct predictors on every path ---------------------

fit1 <- doFitInteractions(interactions(max.order = 1))
expect_equal(worstOrder(extract(fit1, type = "trees")), 1L)

fit2 <- doFitInteractions(interactions(max.order = 2))
expect_true(worstOrder(extract(fit2, type = "trees")) <= 2L)

# ---- forbid bars a named pair from ever co-occurring --------------------------

fitF <- doFitInteractions(interactions(forbid = list(c("x1", "x2"))))
expect_false(anyCoOccur(extract(fitF, type = "trees"), c(1L, 2L)))

# ---- groups: named columns co-occur only with their group-mates ---------------

fitG <- doFitInteractions(interactions(groups = list(c("x1", "x3"), "x2")))
treesG <- extract(fitG, type = "trees")
expect_false(anyCoOccur(treesG, c(1L, 2L))) # different groups
expect_false(anyCoOccur(treesG, c(2L, 3L))) # different groups

# ---- fit-time validation (safe over fast) -------------------------------------

expect_error(interactions(), "at least one of")
expect_error(doFitInteractions(interactions(max.order = 0)), "max.order")
expect_error(doFitInteractions(interactions(max.order = -1L)), "max.order")
expect_error(
  doFitInteractions(interactions(forbid = list(c("x1", "nope")))),
  "unrecognized variable name 'nope'"
)
expect_error(
  doFitInteractions(interactions(groups = list(character(0)))),
  "at least one column"
)
expect_error(
  doFitInteractions(interactions(forbid = list("x1"))),
  "two or more columns"
)

# ---- warm-start refusal: an unconstrained donor's order-2 trees cannot seed a
#      max.order = 1 fit ---------------------------------------------------------

donorArgs <- modifyList(fitArgs, list(keepSampler = TRUE, n.burn = 80L))
donor <- do.call(
  bart2,
  c(list(y ~ x1 + x2 + x3, df, interactions = NULL), donorArgs)
)
expect_error(
  do.call(
    bart2,
    c(
      list(
        y ~ x1 + x2 + x3,
        df,
        interactions = interactions(max.order = 1),
        warm.start = donor
      ),
      fitArgs
    )
  ),
  "interaction constraint"
)

# ---- per-forest BCF: tau capped additive, mu free -----------------------------

p <- 4L
xb <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
z <- rbinom(n, 1L, 0.5)
tau <- ifelse(xb[, 1] > 0.5 & xb[, 2] > 0.5, 3, -1) # wants an x1*x2 interaction
yb <- 2 * sin(pi * xb[, 1]) + xb[, 3] + z * tau + rnorm(n, sd = 0.2)
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 30L,
  updateState = FALSE
)
bcSampler <- dbarts(xb, yb, control = control)
bc <- dbarts:::bartcoreBCFSampler(
  bcSampler,
  z,
  n.trees.treatment = 30L,
  mu.interactions = interactions(max.order = 3),
  tau.interactions = interactions(max.order = 1)
)
invisible(bartcoreRun(bc, 150L, 0L))
tauTrees <- bartcoreGetTrees(
  bc,
  chainNums = 1L,
  treeNums = 1:30,
  current = TRUE,
  forest = 1L
)
muTrees <- bartcoreGetTrees(
  bc,
  chainNums = 1L,
  treeNums = 1:30,
  current = TRUE,
  forest = 0L
)
expect_equal(worstOrder(tauTrees), 1L) # tau forest honors max.order = 1
expect_true(worstOrder(muTrees) >= 2L) # mu forest is unrestricted and uses more

# ---- warm-start refusal: a BCF donor whose treatment forest splits on a non-
#      moderator cannot seed a target whose moderators forbid that column ---
# The tau signal above forces splits on both x1 and x2, so an unrestricted donor
# tau splits on x2; a target restricting tau to x1 must refuse the transplant (a
# moderator forest would otherwise silently score an out-of-mask split).
donorSampler <- dbarts(xb, yb, control = control)
donorBC <- dbarts:::bartcoreBCFSampler(donorSampler, z, n.trees.treatment = 30L)
invisible(bartcoreRun(donorBC, 150L, 0L))
donorState <- bartcoreStoreState(donorBC)

targetSampler <- dbarts(xb, yb, control = control)
targetBC <- dbarts:::bartcoreBCFSampler(
  targetSampler,
  z,
  n.trees.treatment = 30L,
  moderators = 1L
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_installForests,
    targetBC$ptr,
    donorState,
    NULL
  ),
  "column restriction"
)
