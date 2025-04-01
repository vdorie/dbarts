source(system.file("common", "pdData.R", package = "dbarts"), local = TRUE)

# test that pdbart gives same results when run with different x.train argument types
x <- testData$x
y <- testData$y

set.seed(0L)
pdb1 <- dbarts::pdbart(
  x, y, xind = c(1, 2), pl = FALSE,
  levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2)),
  ntree = 5L, ndpost = 10L, nskip = 5L, verbose = FALSE
)

bartFit <- dbarts::bart(
  x, y,
  ntree = 5L, ndpost = 10L, nskip = 5L,
  verbose = FALSE
)
set.seed(0L)
pdb2 <- suppressWarnings(dbarts::pdbart(
  bartFit,
  xind = c(1, 2), pl = FALSE,
  levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2))
))

set.seed(0)
bartFit <- dbarts::bart(
  x, y,
  ntree = 5L, ndpost = 10L, nskip = 5L,
  verbose = FALSE, keeptrees = TRUE
)
pdb3 <- dbarts::pdbart(
  bartFit,
  xind = c(1, 2), pl = FALSE,
  levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2))
)

control <- dbarts::dbartsControl(
  n.trees = 5L, n.samples = 10L, n.burn = 5L,
  verbose = FALSE, n.chains = 1L
)
set.seed(0L)
sampler <- dbarts::dbarts(x, y, control = control)
invisible(sampler$run(0L, 5L))
pdb4 <- suppressWarnings(dbarts::pdbart(
  sampler, xind = c(1, 2), pl = FALSE,
 levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2))
))


control@keepTrees <- TRUE
set.seed(0L)
sampler <- dbarts::dbarts(x, y, control = control)
invisible(sampler$run())
pdb5 <- dbarts::pdbart(
  sampler, xind = c(1, 2), pl = FALSE,
  levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2))
)


expect_equal(pdb1$fd, pdb2$fd)
expect_equal(pdb1$fd, pdb3$fd)
expect_equal(pdb1$fd, pdb4$fd)
expect_equal(pdb1$fd, pdb5$fd)

rm(pdb5, sampler, pdb4, control, pdb3, bartFit, pdb2, pdb1, y, x)


# test that pd2bart gives same results when run with different x.train argument types
x <- testData$x
y <- testData$y

set.seed(0L)
pdb1 <- dbarts::pd2bart(
  x, y, xind = c(2, 3), pl = FALSE,
  levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
  ntree = 5L, ndpost = 10L, nskip = 5L, verbose = FALSE
)

bartFit <- dbarts::bart(
  x, y,
  ntree = 5L, ndpost = 10L, nskip = 5L, verbose = FALSE
)
set.seed(0L)
pdb2 <- suppressWarnings(dbarts::pd2bart(
  bartFit, xind = c(2, 3), pl = FALSE,
  levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
)

set.seed(0L)
bartFit <- dbarts::bart(
  x, y,
  ntree = 5L, ndpost = 10L, nskip = 5L,
  verbose = FALSE, keeptrees = TRUE
)
pdb3 <- dbarts::pd2bart(
  bartFit, xind = c(2, 3), pl = FALSE,
  levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
)

control <- dbarts::dbartsControl(
  n.trees = 5L, n.samples = 10L, n.burn = 5L,
  verbose = FALSE, n.chains = 1L
)
set.seed(0L)
sampler <- dbarts::dbarts(x, y, control = control)
invisible(sampler$run(0, 5))
pdb4 <- suppressWarnings(dbarts::pd2bart(
  sampler, xind = c(2, 3), pl = FALSE,
  levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
))

control@keepTrees <- TRUE
set.seed(0L)
sampler <- dbarts::dbarts(x, y, control = control)
invisible(sampler$run())
pdb5 <- dbarts::pd2bart(
  sampler, xind = c(2, 3), pl = FALSE,
  levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
)

expect_equal(pdb1$fd, pdb2$fd)
expect_equal(pdb1$fd, pdb3$fd)
expect_equal(pdb1$fd, pdb4$fd)
expect_equal(pdb1$fd, pdb5$fd)

rm(pdb5, sampler, pdb4, control, pdb3, bartFit, pdb2, pdb1, y, x)

rm(testData)

