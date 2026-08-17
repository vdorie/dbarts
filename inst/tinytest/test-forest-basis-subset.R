# Formula-path forests = basis declarations, evaluated against the fit's own
# (pre-subset) data and then aligned to 'subset' (R/model.R
# resolveFormulaBasisSubset/alignForestBasisToSubset). With no 'subset', or
# one that keeps every row, a basis's row count must match the data's, as
# always. With 'subset' present the basis must instead cover every row of the
# FULL data and is restricted to the same rows the model frame keeps; a basis
# already at the subset's row count - matching only by coincidence, not
# because it names the pre-subset data - is refused by name. The x/y matrix
# interface, dbartsData(bases = ) called directly, and the forest() formula
# term (':') route are all out of scope and checked here only as regressions.

n <- 40L
set.seed(11)
a <- runif(n)
z <- rbinom(n, 1L, 0.5)
y <- a + z * (1 + a) + rnorm(n, sd = 0.2)
d <- data.frame(y = y, a = a, z = z)
zBasis <- unname(cbind(1 - z, z))
idx <- 1:20

seededControl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.samples = 4L,
    n.burn = 2L,
    updateState = FALSE,
    seed = 7L,
    ...
  )
}
ambiguousPattern <-
  "forest 2's 'basis' has 20 rows, matching 'subset' \\(20\\) but not the full data \\(40 rows\\)"

## --- (i) no subset, basis at the data's row count: unchanged --------------
noSubset <- dbarts(
  y ~ a,
  d,
  forests = list(forest(), forest(basis = zBasis)),
  control = seededControl()
)
expect_equal(dim(noSubset$data@bases[[2L]]), c(n, 2L))

## --- (ii) subset present, basis at the FULL data's row count: now works,
## matching a fit of the manually pre-subset data and basis -----------------
subsetFull <- dbarts(
  y ~ a,
  d,
  subset = idx,
  forests = list(forest(), forest(basis = zBasis)),
  control = seededControl()
)
expect_identical(subsetFull$data@bases[[2L]], zBasis[idx, ])
manual <- dbarts(
  y ~ a,
  d[idx, ],
  forests = list(forest(), forest(basis = zBasis[idx, ])),
  control = seededControl()
)
expect_identical(subsetFull$run(0L, 4L)$train, manual$run(0L, 4L)$train)

## --- (iii) subset present, basis at the SUBSET's row count instead of the
## full data's: refused by name, naming the forest and both counts. This is
## also the mutation-target cell - reverting the ambiguous-count refusal to
## an accept makes exactly this expect_error stop seeing a condition. -------
expect_error(
  dbarts(
    y ~ a,
    d,
    subset = idx,
    forests = list(forest(), forest(basis = zBasis[idx, ])),
    control = seededControl()
  ),
  ambiguousPattern
)

## --- (iv) the same three shapes through a basis FORMULA, where distinguishable
formulaNoSubset <- dbarts(
  y ~ a,
  d,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = seededControl()
)
expect_equal(dim(formulaNoSubset$data@bases[[2L]]), c(n, 2L))

# subset present, the formula evaluates against the full data by
# construction: now works, matching the manual pre-subset fit
formulaSubset <- dbarts(
  y ~ a,
  d,
  subset = idx,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = seededControl()
)
expect_identical(formulaSubset$data@bases[[2L]], zBasis[idx, ])
manualFormula <- dbarts(
  y ~ a,
  d[idx, ],
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = seededControl()
)
expect_identical(
  formulaSubset$run(0L, 4L)$train,
  manualFormula$run(0L, 4L)$train
)

# subset present, a formula whose OWN evaluation slices down to exactly the
# subset's row count instead of the full data's: refused by name, the same
# ambiguous shape a raw matrix hits above
expect_error(
  dbarts(
    y ~ a,
    d,
    subset = idx,
    forests = list(forest(), forest(basis = ~ zBasis[idx, 2L])),
    control = seededControl()
  ),
  ambiguousPattern
)

## --- regression: the forest() formula term (':') route already implements
## the post-subset rule and is unaffected by this slice ---------------------
termFit <- dbarts(
  y ~ a + z:forest(a),
  d,
  subset = idx,
  control = seededControl()
)
expect_equal(dim(termFit$data@bases[[2L]]), c(length(idx), 1L))

## --- regression: the x/y matrix interface keeps its own, already-correct,
## contract - a basis must be full-data length, subset automatically -------
xyFull <- dbarts(
  cbind(a),
  y,
  subset = idx,
  forests = list(forest(), forest(basis = zBasis)),
  control = seededControl()
)
expect_identical(xyFull$data@bases[[2L]], zBasis[idx, ])
expect_error(
  dbarts(
    cbind(a),
    y,
    subset = idx,
    forests = list(forest(), forest(basis = zBasis[idx, ])),
    control = seededControl()
  ),
  "length of 'basis'"
)

## --- regression: dbartsData(bases = ) called directly keeps its own,
## unmodified (count-only) contract - a subset-length basis still passes ---
directData <- dbartsData(
  y ~ a,
  d,
  subset = idx,
  bases = list(NULL, zBasis[idx, ])
)
expect_identical(directData@bases[[2L]], zBasis[idx, ])
