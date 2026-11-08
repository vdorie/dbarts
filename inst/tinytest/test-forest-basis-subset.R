# A forest's basis covers every row of the data the fit was given and is cut
# to the rows 'subset' and the na.action keep. With no 'subset', or one that
# keeps every row, a basis's row count must match the data's, as always. With
# 'subset' present the basis must cover every row of the FULL data and is
# restricted to the same rows the model frame keeps; a basis already at the
# subset's row count - matching only by coincidence, not because it names the
# pre-subset data - is refused by name. A basis written as code says so in
# its own text, at either door; one handed over as a value is refused naming
# the forest and both counts, at the data and at a data object. At an EQUAL row count (subset selects/reorders exactly as many rows
# as the full data has), the two are not ambiguous: a full-data match takes
# priority over a subset-count match, so the basis is read as full-data and
# reordered by 'subset'.

n <- 40L
set.seed(11)
a <- runif(n)
z <- rbinom(n, 1L, 0.5)
y <- a + z * (1 + a) + rnorm(n, sd = 0.2)
d <- data.frame(y = y, a = a, z = z)
zBasis <- unname(cbind(1 - z, z))
zFactor <- factor(z)
idx <- 1:20

seededControlForestBasisSubset <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 4L,
    n.burn = 2L,
    updateState = FALSE,
    seed = 7L,
    ...
  )
}
ambiguousPattern <-
  "forest 2's 'basis' has 20 rows, matching 'subset' \\(20\\) but not the full data \\(40 rows\\)"
# a basis written as code that does not cover the data
shortCode <- function(code, rows) {
  paste0(
    "'basis' (",
    code,
    ") must have the same length as the data: it has ",
    rows,
    " rows and the data 40; a basis covers every row of the data and is cut ",
    "by 'subset' and the na.action with it"
  )
}
# the columns of a matrix written as code are named as lm() names them
named <- function(basis, names) {
  colnames(basis) <- names
  basis
}
valueForests <- function(value) {
  list(
    dbartsForests$forest(),
    do.call(dbartsForests$forest, list(basis = value))
  )
}

## --- (i) no subset, basis at the data's row count: unchanged --------------
noSubset <- dbarts(
  y ~ a,
  d,
  forests = list(forest(), forest(basis = zFactor)),
  control = seededControlForestBasisSubset(),
  n.trees = 5L
)
expect_equal(dim(noSubset$data@bases[[2L]]), c(n, 2L))

## --- (ii) subset present, basis at the FULL data's row count: now works,
## matching a fit of the manually pre-subset data and basis -----------------
subsetFull <- dbarts(
  y ~ a,
  d,
  subset = idx,
  forests = list(forest(), forest(basis = zFactor)),
  control = seededControlForestBasisSubset(),
  n.trees = 5L
)
expect_identical(
  subsetFull$data@bases[[2L]],
  named(zBasis[idx, ], c("zFactor0", "zFactor1"))
)
# the same basis handed over as a value is cut the same way
subsetValue <- dbarts(
  y ~ a,
  d,
  subset = idx,
  forests = valueForests(zFactor),
  control = seededControlForestBasisSubset(),
  n.trees = 5L
)
expect_identical(subsetValue$data@bases[[2L]], zBasis[idx, ])
expect_identical(
  subsetValue$run(0L, 4L)$train,
  dbarts(
    y ~ a,
    d,
    subset = idx,
    forests = list(forest(), forest(basis = zFactor)),
    control = seededControlForestBasisSubset(),
    n.trees = 5L
  )$run(0L, 4L)$train
)
manual <- dbarts(
  y ~ a,
  d[idx, ],
  forests = list(forest(), forest(basis = zFactor[idx])),
  control = seededControlForestBasisSubset(),
  n.trees = 5L
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
    forests = list(forest(), forest(basis = zFactor[idx])),
    control = seededControlForestBasisSubset(),
    n.trees = 5L
  ),
  shortCode("zFactor[idx]", 20L),
  fixed = TRUE
)
expect_error(
  dbarts(
    y ~ a,
    d,
    subset = idx,
    forests = valueForests(zFactor[idx]),
    control = seededControlForestBasisSubset(),
    n.trees = 5L
  ),
  ambiguousPattern
)

## --- (iv) the same three shapes through a basis FORMULA, where distinguishable
formulaNoSubset <- dbarts(
  y ~ a,
  d,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = seededControlForestBasisSubset(),
  n.trees = 5L
)
expect_equal(dim(formulaNoSubset$data@bases[[2L]]), c(n, 2L))

# subset present, the formula evaluates against the full data by
# construction: now works, matching the manual pre-subset fit
formulaSubset <- dbarts(
  y ~ a,
  d,
  subset = idx,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = seededControlForestBasisSubset(),
  n.trees = 5L
)
expect_identical(
  formulaSubset$data@bases[[2L]],
  named(zBasis[idx, ], c("factor(z)0", "factor(z)1"))
)
manualFormula <- dbarts(
  y ~ a,
  d[idx, ],
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = seededControlForestBasisSubset(),
  n.trees = 5L
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
    control = seededControlForestBasisSubset(),
    n.trees = 5L
  ),
  shortCode("zBasis[idx, 2L]", 20L),
  fixed = TRUE
)

## --- a forest() term of the formula is cut the same way ---------------------
termFit <- dbarts(
  y ~ a + forest(a, basis = ~z),
  d,
  subset = idx,
  control = seededControlForestBasisSubset(),
  n.trees = 5L
)
expect_identical(
  termFit$data@bases[[2L]],
  matrix(as.double(z[idx]), dimnames = list(NULL, "z"))
)
expect_error(
  dbarts(
    y ~ a + forest(a, basis = zFactor[idx]),
    d,
    subset = idx,
    control = seededControlForestBasisSubset(),
    n.trees = 5L
  ),
  shortCode("zFactor[idx]", 20L),
  fixed = TRUE
)

## --- regression: the x/y matrix interface keeps its own, already-correct,
## contract - a basis must be full-data length, subset automatically -------
xyFull <- dbarts(
  cbind(a),
  y,
  subset = idx,
  forests = list(forest(), forest(basis = zFactor)),
  control = seededControlForestBasisSubset(),
  n.trees = 5L
)
expect_identical(
  xyFull$data@bases[[2L]],
  named(zBasis[idx, ], c("zFactor0", "zFactor1"))
)
expect_error(
  dbarts(
    cbind(a),
    y,
    subset = idx,
    forests = list(forest(), forest(basis = zFactor[idx])),
    control = seededControlForestBasisSubset(),
    n.trees = 5L
  ),
  shortCode("zFactor[idx]", 20L),
  fixed = TRUE
)
expect_error(
  dbarts(
    cbind(a),
    y,
    subset = idx,
    forests = valueForests(zFactor[idx]),
    control = seededControlForestBasisSubset(),
    n.trees = 5L
  ),
  ambiguousPattern
)

## --- a data object applies the SAME rule as dbarts(forests = ): a basis at
## the FULL data's row count is accepted and aligned to 'subset' -------------
subsetObject <- function(value, rows = idx) {
  dbarts(
    dbartsData(y ~ a, d, subset = rows),
    forests = valueForests(value),
    control = seededControlForestBasisSubset(),
    n.trees = 5L
  )
}
directFull <- subsetObject(zFactor)
expect_identical(unname(directFull$data@bases[[2L]]), zBasis[idx, ])

## a basis at the SUBSET's row count instead is refused by name, naming the
## forest and both counts - the same message dbarts(forests = ) gives
expect_error(subsetObject(zFactor[idx]), ambiguousPattern)

## --- equal row count (subset selects/reorders exactly as many rows as the
## full data has): a data object and dbarts(forests = ) resolve it
## IDENTICALLY rather than ambiguously - a full-data-length match takes
## priority over a subset-count match, so both reorder the basis by 'subset'
recycled <- rep(seq_len(20L), 2L)
equalCountDirect <- subsetObject(zFactor, recycled)$data
equalCountForests <- dbarts(
  y ~ a,
  d,
  subset = recycled,
  forests = list(forest(), forest(basis = zFactor)),
  control = seededControlForestBasisSubset(),
  n.trees = 5L
)
expect_identical(unname(equalCountDirect@bases[[2L]]), zBasis[recycled, ])
expect_identical(
  equalCountForests$data@bases[[2L]],
  named(zBasis[recycled, ], c("zFactor0", "zFactor1"))
)

## --- a basis already shortened by the na.action, with no 'subset': the
## refusal names the rows the na.action kept, not a 'subset' never given ----
dMissing <- d
dMissing$y[3L] <- NA
expect_error(
  dbarts(
    y ~ a,
    dMissing,
    forests = list(forest(), forest(basis = zFactor[-3L])),
    control = seededControlForestBasisSubset(),
    n.trees = 5L
  ),
  shortCode("zFactor[-3L]", 39L),
  fixed = TRUE
)
expect_error(
  dbarts(
    y ~ a,
    dMissing,
    forests = valueForests(zFactor[-3L]),
    control = seededControlForestBasisSubset(),
    n.trees = 5L
  ),
  paste0(
    "forest 2's 'basis' has 39 rows, matching the rows the na.action kept ",
    "\\(39\\) but not the full data \\(40 rows\\)"
  )
)
