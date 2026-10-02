source(system.file("common", "pdData.R", package = "dbarts"), local = TRUE)
source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

# test that pdbart gives same results when run with different x.train argument types
x <- testData$x
y <- testData$y

set.seed(0L)
pdb1 <- dbarts::pdbart(
  x,
  y,
  xind = c(1, 2),
  pl = FALSE,
  levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2)),
  ntree = 5L,
  ndpost = 10L,
  nskip = 5L,
  verbose = FALSE
)

bartFit <- dbarts::bart(
  x,
  y,
  ntree = 5L,
  ndpost = 10L,
  nskip = 5L,
  verbose = FALSE
)
set.seed(0L)
warnings.pdb2 <- captureWarnings(
  pdb2 <- dbarts::pdbart(
    bartFit,
    xind = c(1, 2),
    pl = FALSE,
    levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2))
  )
)
expect_true(any(vapply(
  warnings.pdb2,
  inherits,
  logical(1L),
  "dbartsFallbackWarning"
)))

set.seed(0)
bartFit <- dbarts::bart(
  x,
  y,
  ntree = 5L,
  ndpost = 10L,
  nskip = 5L,
  verbose = FALSE,
  keeptrees = TRUE
)
pdb3 <- dbarts::pdbart(
  bartFit,
  xind = c(1, 2),
  pl = FALSE,
  levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2))
)

# bartBT's tree-move mixture, which the BayesTree-spelled fits above run under
control <- dbarts::dbartsControl(
  n.trees = 5L,
  n.samples = 10L,
  n.burn = 5L,
  verbose = FALSE,
  n.chains = 1L,
  proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
)
set.seed(0L)
sampler <- dbarts::dbarts(x, y, control = control)
invisible(sampler$run(0L, 5L))
pdb4 <- suppressWarnings(dbarts::pdbart(
  sampler,
  xind = c(1, 2),
  pl = FALSE,
  levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2))
))


control@keepTrees <- TRUE
set.seed(0L)
sampler <- dbarts::dbarts(x, y, control = control)
invisible(sampler$run())
pdb5 <- dbarts::pdbart(
  sampler,
  xind = c(1, 2),
  pl = FALSE,
  levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2))
)


expect_equal(pdb1$fd, pdb2$fd)
expect_equal(pdb1$fd, pdb3$fd)
expect_equal(pdb1$fd, pdb4$fd)
expect_equal(pdb1$fd, pdb5$fd)

# the plot method renders each requested predictor: one device page per
# entry of xind, each carrying that predictor's own axis label
psFile <- tempfile(fileext = ".ps")
postscript(psFile, onefile = TRUE)
expect_silent(plot(pdb1))
dev.off()
psLines <- readLines(psFile, warn = FALSE)
expect_equal(length(grep("^%%Page:", psLines)), 2L)
labelHits <- vapply(
  pdb1$xlbs,
  function(lab) length(grep(paste0("(", lab, ")"), psLines, fixed = TRUE)),
  integer(1L)
)
expect_equivalent(labelHits, c(1L, 1L))

# and xind selects WHICH predictor, not just how many: the second alone
postscript(psFile, onefile = TRUE)
expect_silent(plot(pdb1, xind = 2L))
dev.off()
psLines <- readLines(psFile, warn = FALSE)
expect_equal(length(grep("^%%Page:", psLines)), 1L)
expect_equal(
  length(grep(paste0("(", pdb1$xlbs[2L], ")"), psLines, fixed = TRUE)),
  1L
)
unlink(psFile)
rm(psFile, psLines, labelHits)

# sampleronly is set internally (pdbart always needs the sampler, not just
# its result); an explicit user value must error rather than be silently
# clobbered
expect_error(
  dbarts::pdbart(
    x,
    y,
    xind = c(1, 2),
    pl = FALSE,
    levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2)),
    ntree = 5L,
    ndpost = 10L,
    nskip = 5L,
    sampleronly = TRUE,
    verbose = FALSE
  ),
  pattern = "set internally"
)

rm(pdb5, sampler, pdb4, control, pdb3, bartFit, pdb2, pdb1, y, x, warnings.pdb2)


# test that pd2bart gives same results when run with different x.train argument types
x <- testData$x
y <- testData$y

set.seed(0L)
pdb1 <- dbarts::pd2bart(
  x,
  y,
  xind = c(2, 3),
  pl = FALSE,
  levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
  ntree = 5L,
  ndpost = 10L,
  nskip = 5L,
  verbose = FALSE
)

bartFit <- dbarts::bart(
  x,
  y,
  ntree = 5L,
  ndpost = 10L,
  nskip = 5L,
  verbose = FALSE
)
set.seed(0L)
pdb2 <- suppressWarnings(dbarts::pd2bart(
  bartFit,
  xind = c(2, 3),
  pl = FALSE,
  levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
))

set.seed(0L)
bartFit <- dbarts::bart(
  x,
  y,
  ntree = 5L,
  ndpost = 10L,
  nskip = 5L,
  verbose = FALSE,
  keeptrees = TRUE
)
pdb3 <- dbarts::pd2bart(
  bartFit,
  xind = c(2, 3),
  pl = FALSE,
  levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
)

# bartBT's tree-move mixture, which the BayesTree-spelled fits above run under
control <- dbarts::dbartsControl(
  n.trees = 5L,
  n.samples = 10L,
  n.burn = 5L,
  verbose = FALSE,
  n.chains = 1L,
  proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
)
set.seed(0L)
sampler <- dbarts::dbarts(x, y, control = control)
invisible(sampler$run(0, 5))
pdb4 <- suppressWarnings(dbarts::pd2bart(
  sampler,
  xind = c(2, 3),
  pl = FALSE,
  levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
))

control@keepTrees <- TRUE
set.seed(0L)
sampler <- dbarts::dbarts(x, y, control = control)
invisible(sampler$run())
pdb5 <- dbarts::pd2bart(
  sampler,
  xind = c(2, 3),
  pl = FALSE,
  levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
)

expect_equal(pdb1$fd, pdb2$fd)
expect_equal(pdb1$fd, pdb3$fd)
expect_equal(pdb1$fd, pdb4$fd)
expect_equal(pdb1$fd, pdb5$fd)

# the contour plot renders into a null device
pdf(NULL)
expect_silent(plot(pdb1))
dev.off()

# same guard as pdbart: sampleronly is set internally
expect_error(
  dbarts::pd2bart(
    x,
    y,
    xind = c(2, 3),
    pl = FALSE,
    levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
    ntree = 5L,
    ndpost = 10L,
    nskip = 5L,
    sampleronly = TRUE,
    verbose = FALSE
  ),
  pattern = "set internally"
)

rm(pdb5, sampler, pdb4, control, pdb3, bartFit, pdb2, pdb1, y, x)

# the four own-class fits (bartMultinomial/Ordinal/Negbin/Hurdle) reach
# neither the "bart" nor the dbartsSampler branch, and used to
# fall through to the generic "'x.train' must be a matrix, ..." message,
# naming neither the fit nor why it fails; refused by name instead
negbinFit <- dbarts::bart(
  testData$x,
  rpois(nrow(testData$x), 3L),
  family = "nbinom",
  n.trees = 5L,
  n.samples = 10L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  dbarts::pdbart(negbinFit, xind = c(1, 2), pl = FALSE),
  "pdbart does not support a bartNegbin fit",
  fixed = TRUE
)
rm(negbinFit)


# A factor predictor is evaluated at every level by default, levels are given
# and reported by name, and each value equals a direct prediction at the level.
set.seed(7)
factorFrame <- data.frame(
  y = rnorm(150L),
  x1 = runif(150L),
  f = factor(sample(LETTERS[1:15], 150L, TRUE))
)
factorFrame$y <- factorFrame$y + as.integer(factorFrame$f)
factorFit <- suppressMessages(bart(
  y ~ .,
  factorFrame,
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 2L,
  n.trees = 10L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
pd <- pdbart(factorFit, xind = "f", pl = FALSE)
expect_identical(pd$levs[[1L]], LETTERS[1:15])
atLevel <- factorFrame
atLevel$f <- factor("C", levels = levels(factorFrame$f))
expect_equal(pd$fd[[1L]][, 3L], rowMeans(predict(factorFit, atLevel)))
pd <- pdbart(factorFit, xind = "f", levs = list(c("D", "A")), pl = FALSE)
expect_identical(pd$levs[[1L]], c("D", "A"))
atLevel$f <- factor("A", levels = levels(factorFrame$f))
expect_equal(pd$fd[[1L]][, 2L], rowMeans(predict(factorFit, atLevel)))
expect_error(
  pdbart(factorFit, xind = "f", levs = list(c("A", "ZZ")), pl = FALSE),
  "'ZZ'"
)
expect_error(
  pdbart(factorFit, xind = "f", levs = list(1:2), pl = FALSE),
  "must name its levels"
)
pd2 <- pd2bart(factorFit, xind = c("x1", "f"), pl = FALSE)
expect_identical(pd2$levs[[2L]], LETTERS[1:15])
expect_equal(ncol(pd2$fd), length(pd2$levs[[1L]]) * 15L)
grDevices::pdf(NULL)
plot(pd)
plot(pd2)
grDevices::dev.off()
rm(factorFrame, factorFit, pd, atLevel, pd2)

# With exactly two predictors each grid point is a whole row: the result keeps
# one column per grid point, equal to a direct prediction there, for named
# predictors, a fit object, and either keepTrees setting.
set.seed(8)
xTwo <- matrix(runif(100L), 50L, 2L, dimnames = list(NULL, c("x1", "x2")))
yTwo <- 2 * xTwo[, 1L] + rnorm(50L, sd = 0.1)
for (keepTrees in c(TRUE, FALSE)) {
  pd2 <- pd2bart(
    xTwo,
    yTwo,
    xind = c(1L, 2L),
    pl = FALSE,
    keeptrees = keepTrees,
    ndpost = 10L,
    nskip = 10L,
    ntree = 10L,
    verbose = FALSE
  )
  expect_equal(dim(pd2$fd), c(10L, 121L))
}
grDevices::pdf(NULL)
plot(pd2)
grDevices::dev.off()
twoFit <- suppressMessages(bart(
  xTwo,
  yTwo,
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 2L,
  n.trees = 10L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
pd2 <- pd2bart(twoFit, xind = c(2L, 1L), pl = FALSE)
grid <- as.matrix(expand.grid(pd2$levs[[1L]], pd2$levs[[2L]]))[, c(2L, 1L)]
colnames(grid) <- c("x1", "x2")
expect_equal(pd2$fd, predict(twoFit, grid), check.attributes = FALSE)
rm(xTwo, yTwo, keepTrees, pd2, twoFit, grid)

rm(testData)
