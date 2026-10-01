# Arguments the sampler's methods refuse at the R and C boundary, and what a
# refusal leaves behind: the sampler a refused call reached must be the one it
# was before the call, mirrors included.

set.seed(1)
n <- 20L
x <- matrix(rnorm(n * 2L), n)
y <- rnorm(n)
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  keepTrees = TRUE,
  updateState = FALSE,
  verbose = FALSE
)
sampler <- dbarts(x, y, control = control)

# --- run counts: a negative count wrapped to a huge size and reported draws
# the run never recorded; 0 + 0 returned NULL. Both refused as in 0.9-x, at
# the method and at the C entry beneath it
expect_error(
  sampler$run(-2L, 5L),
  "number of burn-in steps must be greater than or equal to 0"
)
expect_error(
  sampler$run(3L, -1L),
  "number of samples must be greater than or equal to 0"
)
expect_error(
  sampler$run(0L, 0L),
  "either number of burn-in or samples must be positive"
)
expect_error(sampler$run(2.5, 1L), "whole number")
expect_error(sampler$run(c(1L, 2L), 1L), "'numBurnIn' must be a single")
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_run,
    sampler$getPointer(),
    -1L,
    2L,
    NULL,
    NULL,
    TRUE
  ),
  "number of burn-in steps must be greater than or equal to 0"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_run,
    sampler$getPointer(),
    0L,
    0L,
    NULL,
    NULL,
    TRUE
  ),
  "either number of burn-in or samples must be positive"
)
# a refused run records nothing: the store still holds no draws
expect_error(sampler$predict(x), "holds no recorded draws")
r <- sampler$run(2L, 5L)
expect_true(all(is.finite(r$sigma)) && all(r$sigma > 1e-100))

# --- missing indices are refused by name, not by a bare missing condition
expect_error(
  sampler$getTrees(chainNums = NA_integer_),
  "'chainNums' contains missing values"
)
expect_error(
  sampler$getTrees(sampleNums = NA_integer_),
  "'sampleNums' contains missing values"
)
expect_error(
  sampler$getTrees(treeNums = NA_integer_),
  "'treeNums' contains missing values"
)
expect_error(
  sampler$setPredictor(x[, 1L], NA_integer_),
  "'column' contains missing values"
)

# --- non-finite swaps: one infinite value left sigma and every fit NaN for
# good, even after the value was put back; refused as creation refuses one,
# and the sampler is unchanged
expect_error(
  sampler$setResponse(replace(y, 1L, Inf)),
  "response contains non-finite values"
)
expect_error(
  sampler$setOffset(replace(numeric(n), 1L, -Inf)),
  "'offset' contains non-finite values"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setResponse,
    sampler$getPointer(),
    replace(y, 1L, Inf),
    FALSE,
    NULL
  ),
  "response contains non-finite values"
)
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setOffset,
    sampler$getPointer(),
    replace(numeric(n), 2L, NaN),
    FALSE
  ),
  "offset contains non-finite values"
)
expect_identical(sampler$data@y, y)
expect_null(sampler$data@offset)
r <- sampler$run(0L, 2L)
expect_true(all(is.finite(r$sigma)) && all(is.finite(r$train)))

# the data conduits state the same rule: dbartsData refuses an infinite
# offset as it refuses an infinite response, and the whole-data swap's entry
# refuses either on a data object that skipped that check
expect_error(
  dbartsData(x, y, offset = replace(numeric(n), 3L, Inf)),
  "'offset' contains non-finite values"
)
badData <- sampler$data
badData@y <- replace(y, 2L, Inf)
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setData, sampler$getPointer(), badData),
  "\\$setData: response contains non-finite values"
)
badData <- sampler$data
badData@offset <- replace(numeric(n), 2L, -Inf)
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setData, sampler$getPointer(), badData),
  "\\$setData: offset contains non-finite values"
)
expect_error(sampler$setSigma(Inf), "'sigma' must be finite and positive")
rm(badData)

# --- updateScale is a single TRUE or FALSE: 1 rescaled the engine while the R
# side skipped its own rescaling steps
expect_error(
  sampler$setOffset(numeric(n), updateScale = 1),
  "'updateScale' must be TRUE or FALSE"
)
expect_error(
  sampler$setResponse(y, updateScale = NA),
  "'updateScale' must be TRUE or FALSE"
)

# --- a refused offset swap leaves the mirror a save and load re-creates from
set.seed(3)
nb <- 60L
xb <- matrix(runif(nb * 2L), nb)
zb <- rbinom(nb, 1L, 0.5)
yb <- xb[, 1L] + zb + rnorm(nb, sd = 0.2)
bcf <- dbarts(
  xb,
  yb,
  forests = list(forest(), forest(basis = ~ factor(zb), n.trees = 10L)),
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    updateState = FALSE,
    verbose = FALSE
  )
)
bcf$setOffset(rep(0.5, nb), updateScale = FALSE)
expect_error(
  bcf$setOffset(rep(5, nb), updateScale = NA),
  "'updateScale' must be TRUE or FALSE"
)
expect_identical(bcf$data@offset, rep(0.5, nb))

# --- allocation failures raise rather than abort: a tree store far past what
# the process can hold, at creation and through setControl. Only at home and
# not under AddressSanitizer: ASan aborts on the deliberate oversize request
# instead of letting it throw, and valgrind flags it as a silly argument, so
# CRAN's memory-check flavors would report the request itself.
if (at_home() && !nzchar(Sys.getenv("ASAN_OPTIONS"))) {
  expect_error(
    dbarts(
      x,
      y,
      control = dbartsControl(
        n.trees = 200000L,
        n.samples = .Machine$integer.max,
        keepTrees = TRUE,
        n.chains = 1L,
        n.threads = 1L,
        verbose = FALSE
      )
    ),
    "sampler creation failed"
  )
  wide <- dbarts(
    x,
    y,
    control = dbartsControl(
      n.trees = 200000L,
      n.samples = 4L,
      keepTrees = TRUE,
      n.chains = 2L,
      n.threads = 1L,
      updateState = FALSE,
      verbose = FALSE
    )
  )
  wideControl <- wide$control
  wideControl@n.samples <- .Machine$integer.max
  expect_error(
    wide$setControl(wideControl),
    "saved-tree storage for 2147483647 samples cannot be allocated"
  )
  expect_identical(wide$control@n.samples, 4L)
  # the new store is built aside, so the old one keeps its capacity and the
  # draws it recorded
  invisible(wide$run(0L, 4L))
  widePredict <- wide$predict(x[1:3, , drop = FALSE])
  expect_error(wide$setControl(wideControl), "cannot be allocated")
  expect_identical(wide$predict(x[1:3, , drop = FALSE]), widePredict)
  expect_identical(dim(widePredict), c(3L, 4L, 2L))
  rm(wide, wideControl, widePredict)
}

# --- a state with another store capacity. A live $setState keeps the
# sampler's own capacity, so the state is refused and the draws stay; the
# re-creation path takes the state's capacity only once the state is
# accepted, so a refused one leaves the store and its draws alone too
storeControl <- function(n.samples, n.trees = 5L) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = n.trees,
    n.samples = n.samples,
    keepTrees = TRUE,
    updateState = FALSE,
    verbose = FALSE,
    seed = 9L
  )
}
held <- dbarts(x, y, control = storeControl(5L))
invisible(held$run(2L, 5L))
heldPredict <- held$predict(x)
other <- dbarts(x, y, control = storeControl(3L))
invisible(other$run(2L, 3L))
other$storeState()
expect_error(held$setState(other$state), "not consistent")
expect_identical(held$predict(x), heldPredict)
expect_identical(dim(held$predict(x)), c(n, 5L))
wrongTrees <- dbarts(x, y, control = storeControl(3L, n.trees = 6L))
invisible(wrongTrees$run(2L, 3L))
wrongTrees$storeState()
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setState,
    held$getPointer(),
    wrongTrees$state,
    x,
    TRUE
  ),
  "not consistent"
)
expect_identical(held$predict(x), heldPredict)
rm(held, heldPredict, other, wrongTrees, storeControl)

# --- a refused warm start touches no chain. A donor whose trees no longer
# route onto the predictors it now holds fails to rebuild part way through
# the install; every chain installed before it is put back
set.seed(1)
nw <- 12L
xw <- data.frame(
  a = rnorm(nw),
  b = factor(sample(letters[1:3], nw, TRUE)),
  c = rnorm(nw)
)
yw <- xw$a + rnorm(nw)
warmControl <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 3L,
  keepTrees = TRUE,
  updateState = FALSE,
  verbose = FALSE,
  seed = 3L
)
donor <- dbarts(yw ~ ., xw, control = warmControl)
invisible(donor$run(2L, 2L))
invisible(donor$setPredictor(xw$c * 2, 3L))
donor$setData(donor$data)
recipientControl <- warmControl
recipientControl@keepTrees <- FALSE
recipientControl@seed <- 7L
recipient <- dbarts(yw ~ ., xw, control = recipientControl)
twin <- dbarts(yw ~ ., xw, control = recipientControl)
invisible(recipient$run(5L, 1L))
invisible(twin$run(5L, 1L))
predictBefore <- recipient$predict(xw)
treesBefore <- recipient$getTrees()
expect_error(
  recipient$installTrees(donor, samples = c(1L, 3L)),
  "cannot be rebuilt on this sampler's data"
)
expect_identical(recipient$predict(xw), predictBefore)
expect_identical(recipient$getTrees(), treesBefore)
# the rebuild is judged on scratch trees before anything is replaced, so the
# continuation is the untouched twin's, bitwise
expect_identical(recipient$run(5L, 5L), twin$run(5L, 5L))

# --- a heteroscedastic donor whose variance forest holds another tree count
# is refused before any mean forest is replaced, and says so
set.seed(0)
nh <- 80L
xh <- cbind(x1 = runif(nh), x2 = runif(nh))
yh <- 2 * xh[, 1L] + rnorm(nh, 0, 0.2 + xh[, 2L])
hetero <- function(keepTrees, seed, numVarianceTrees, numChains) {
  dbarts(
    xh,
    yh,
    variance = dbartsForests$varianceForest(n.trees = numVarianceTrees),
    control = dbartsControl(
      n.trees = 6L,
      n.chains = numChains,
      n.threads = 1L,
      seed = seed,
      keepTrees = keepTrees,
      updateState = TRUE,
      verbose = FALSE
    )
  )
}
heteroDonor <- hetero(TRUE, 1L, 4L, 1L)
invisible(heteroDonor$run(50L, 3L))
heteroRecipient <- hetero(FALSE, 2L, 5L, 2L)
heteroTwin <- hetero(FALSE, 2L, 5L, 2L)
invisible(heteroRecipient$run(30L, 1L))
invisible(heteroTwin$run(30L, 1L))
fitsBefore <- heteroRecipient$getForestFits()
varianceBefore <- heteroRecipient$getVariance()
expect_error(
  heteroRecipient$installTrees(heteroDonor),
  "variance forest has a different number of trees"
)
expect_identical(heteroRecipient$getForestFits(), fitsBefore)
expect_identical(heteroRecipient$getVariance(), varianceBefore)
expect_identical(heteroRecipient$run(5L, 5L), heteroTwin$run(5L, 5L))
