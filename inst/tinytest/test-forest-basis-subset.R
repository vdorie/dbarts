# Formula-path forests = basis declarations, evaluated against the fit's own
# (pre-subset) data and then aligned to 'subset' (R/model.R
# resolveFormulaBasisSubset/alignForestBasisToSubset, and R/data.R
# validateForestBases's 'subsetRows' branch - the same rule dbartsData(bases =
# ) now applies directly on its own formula path). With no 'subset', or one
# that keeps every row, a basis's row count must match the data's, as always.
# With 'subset' present the basis must instead cover every row of the FULL
# data and is restricted to the same rows the model frame keeps; a basis
# already at the subset's row count - matching only by coincidence, not
# because it names the pre-subset data - is refused by name, at both
# dbarts(forests = ) and a direct dbartsData(bases = ) call. At an EQUAL row
# count (subset selects/reorders exactly as many rows as the full data has),
# the two are not ambiguous: a full-data match takes priority over a
# subset-count match, so both read the basis as full-data and reorder it by
# 'subset'. The x/y matrix interface (already implementing the full-data
# rule) and the forest() formula term (':') route (evaluated post-subset by
# construction, and unaffected by this rule) are checked here only as
# regressions.

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
## the post-subset rule and is unaffected by this change --------------------
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
  "'basis' must have the same length"
)

## --- dbartsData(bases = ) called directly now applies the SAME rule as
## dbarts(forests = ): a basis at the FULL data's row count is accepted and
## aligned to 'subset' ---------------------------------------------------
directFull <- dbartsData(y ~ a, d, subset = idx, bases = list(NULL, zBasis))
expect_identical(directFull@bases[[2L]], zBasis[idx, ])

## a basis at the SUBSET's row count instead is refused by name, naming the
## forest and both counts - the same message dbarts(forests = ) gives, since
## dbartsData() is now the single place that raises it. Mutation-target cell
## for this entry: reverting validateForestBases's 'subsetRows' branch back
## to a bare count check makes exactly this expect_error stop seeing a
## condition. ---------------------------------------------------------------
expect_error(
  dbartsData(y ~ a, d, subset = idx, bases = list(NULL, zBasis[idx, ])),
  ambiguousPattern
)

## --- equal row count (subset selects/reorders exactly as many rows as the
## full data has): dbartsData(bases = ) and dbarts(forests = ) resolve it
## IDENTICALLY rather than ambiguously - a full-data-length match takes
## priority over a subset-count match, so both reorder the basis by 'subset'
recycled <- rep(seq_len(20L), 2L)
equalCountDirect <- dbartsData(
  y ~ a,
  d,
  subset = recycled,
  bases = list(NULL, zBasis)
)
equalCountForests <- dbarts(
  y ~ a,
  d,
  subset = recycled,
  forests = list(forest(), forest(basis = zBasis)),
  control = seededControl()
)
expect_identical(equalCountDirect@bases[[2L]], zBasis[recycled, ])
expect_identical(equalCountForests$data@bases[[2L]], zBasis[recycled, ])
