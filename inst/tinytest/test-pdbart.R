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

control <- dbarts::dbartsControl(
  n.trees = 5L,
  n.samples = 10L,
  n.burn = 5L,
  verbose = FALSE,
  n.chains = 1L
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

control <- dbarts::dbartsControl(
  n.trees = 5L,
  n.samples = 10L,
  n.burn = 5L,
  verbose = FALSE,
  n.chains = 1L
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

# rbart_vi and the four own-class fits (bartMultinomial/Ordinal/Negbin/
# Hurdle) reach neither the "bart" nor the dbartsSampler branch, and used to
# fall through to the generic "'x.train' must be a matrix, ..." message,
# naming neither the fit nor why it fails; refused by name instead
negbinFit <- dbarts::bart2(
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

rbartFit <- dbarts::rbart_vi(
  testData$y ~ testData$x,
  group.by = rep(1:2, length.out = nrow(testData$x)),
  n.trees = 5L,
  n.samples = 10L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_error(
  dbarts::pdbart(rbartFit, xind = c(1, 2), pl = FALSE),
  "pdbart does not support a rbart fit",
  fixed = TRUE
)
expect_error(
  dbarts::pd2bart(rbartFit, xind = c(1, 2), pl = FALSE),
  "pd2bart does not support a rbart fit",
  fixed = TRUE
)
rm(rbartFit)

rm(testData)
