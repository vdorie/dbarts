# $setData reads each predictor column by name and each factor row by its
# label against the level table the sampler was created with, as setPredictor
# and predict read one (dec-B454). Fewer levels, another order and declared
# levels no row holds go in, the record keeping the sampler's table, so a copy
# and a reload continue as the original does; a row at a label the sampler
# never had is refused by name, and so is a factor column of bare numbers.
# Every sampler is seeded before creation: its chain seeds are drawn then.

n.fl <- 200L
set.seed(1L)
a.fl <- runif(n.fl)
effect.fl <- c(a = 0, b = 2, c = 4, d = 0, e = 2, f = 4, g = 1, h = 3)
control.fl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 10L,
  n.burn = 0L,
  keepTrees = TRUE,
  updateState = TRUE,
  verbose = FALSE
)
# rows at 'used' declaring 'levels', the response on the first set's scale
frameOf <- function(levels, used, seed, ordered = FALSE, missing = FALSE) {
  set.seed(seed)
  labels <- sample(used, n.fl, TRUE)
  f <- factor(labels, levels = levels, ordered = ordered)
  if (missing) {
    f[sample(n.fl, 15L)] <- NA
  }
  y <- a.fl + effect.fl[tolower(labels)] + rnorm(n.fl, 0, 0.3)
  data.frame(y = 1 + 5 * (y - min(y)) / diff(range(y)), a = a.fl, f = f)
}
samplerOf <- function(ordered = FALSE, missing = FALSE, frame = NULL) {
  if (is.null(frame)) {
    frame <- frameOf(letters[1:6], letters[1:6], 1L, ordered, missing)
  }
  set.seed(4L)
  sampler <- dbarts::dbarts(y ~ a + f, frame, control = control.fl)
  set.seed(5L)
  invisible(sampler$run(30L, 10L))
  sampler
}
# the same rows, their factor re-coded by label onto the sampler's six levels
byLabel <- function(frame, levels = letters[1:6]) {
  frame$f <- factor(
    as.character(frame$f),
    levels = levels,
    ordered = is.ordered(frame$f)
  )
  frame
}
reloadOf <- function(sampler) {
  path <- tempfile(fileext = ".rds")
  saveRDS(sampler, path)
  reloaded <- readRDS(path)
  unlink(path)
  reloaded
}
fitsOf <- function(sampler) {
  set.seed(77L)
  sampler$run(0L, 5L)$train
}
predictAt <- function(sampler, labels, ordered = FALSE) {
  rowMeans(sampler$predict(data.frame(
    a = 0.5,
    f = factor(labels, levels = letters[1:6], ordered = ordered)
  )))
}
untouched <- function(sampler) list(sampler$data, sampler$state)

# same table: the record is the object handed over, bit for bit
sampler <- samplerOf()
handed <- dbarts::dbartsData(
  y ~ a + f,
  frameOf(letters[1:6], letters[1:6], 2L),
  test = frameOf(letters[1:6], letters[1:3], 3L)[1:10, ]
)
sampler$setData(handed)
expect_identical(sampler$data@x, handed@x)
expect_identical(sampler$data@x.test, handed@x.test)

# fewer levels, the last two dropped, on each kind of factor: the record keeps
# the sampler's six, a copy and a reload continue as the original, and the
# original's stored state goes in everywhere
for (kind in list(c(FALSE, FALSE), c(TRUE, FALSE), c(TRUE, TRUE))) {
  info <- paste("ordered", kind[1L], "missing", kind[2L])
  sampler <- samplerOf(kind[1L], kind[2L])
  sampler$setData(dbarts::dbartsData(
    y ~ a + f,
    frameOf(letters[1:4], letters[1:4], 2L, kind[1L], kind[2L])
  ))
  expect_identical(
    attr(sampler$data@x, "factor.levels")[[2L]],
    letters[1:6],
    info = info
  )
  sampler$storeState()
  own <- sampler$state
  copied <- sampler$copy()
  reloaded <- reloadOf(sampler)
  expect_identical(copied$getTrees(), sampler$getTrees(), info = info)
  original <- fitsOf(sampler)
  expect_equal(fitsOf(copied), original, info = info)
  expect_equal(fitsOf(reloaded), original, info = info)
  expect_true(sampler$setState(own), info = info)
  expect_true(copied$setState(own), info = info)
  expect_true(reloadOf(sampler)$setState(own), info = info)
  expect_true(reloaded$setState(own), info = info)
}
# and the sampler and its copy take all six levels back
copied <- sampler$copy()
again <- dbarts::dbartsData(
  y ~ a + f,
  frameOf(letters[1:6], letters[1:6], 3L, TRUE, TRUE)
)
sampler$setData(again)
copied$setData(again)
expect_equal(fitsOf(copied), fitsOf(sampler))

# fewer levels, b dropped: each row at its own label, as a twin handed the
# rows re-coded by hand, training and test, and the kept draws predict at c
# what they did
for (ordered in c(FALSE, TRUE)) {
  dropped <- frameOf(letters[c(1L, 3:6)], letters[c(1L, 3:6)], 2L, ordered)
  test <- frameOf(letters[c(3L, 5L)], letters[c(3L, 5L)], 3L, ordered)[1:10, ]
  sampler <- samplerOf(ordered)
  twin <- samplerOf(ordered)
  before <- predictAt(sampler, c("c", "e"), ordered)
  sampler$setData(dbarts::dbartsData(y ~ a + f, dropped, test = test))
  twin$setData(dbarts::dbartsData(
    y ~ a + f,
    byLabel(dropped),
    test = byLabel(test)
  ))
  expect_identical(sampler$data@x, twin$data@x, info = ordered)
  expect_equal(predictAt(sampler, c("c", "e"), ordered), before, info = ordered)
  set.seed(9L)
  result <- sampler$run(0L, 2L)
  set.seed(9L)
  expect_identical(result, twin$run(0L, 2L), info = ordered)
}

# a character column rebuilt from a draw that missed b: its levels are those
# present, sorted, and its c rows are c
sampler <- samplerOf()
drawn <- frameOf(letters[1:6], letters[c(1L, 3:6)], 2L)
drawn$f <- as.character(drawn$f)
sampler$setData(dbarts::dbartsData(y ~ a + f, drawn))
expect_identical(
  as.character(sampler$data@x$dense[[2L]]),
  drawn$f
)
expect_identical(levels(sampler$data@x$dense[[2L]]), letters[1:6])

# the levels reordered, unordered and ordered: read by label, the sampler's
# order standing, and a copy continues as the original
for (ordered in c(FALSE, TRUE)) {
  reordered <- frameOf(rev(letters[1:6]), letters[1:6], 2L, ordered)
  sampler <- samplerOf(ordered)
  twin <- samplerOf(ordered)
  sampler$setData(dbarts::dbartsData(y ~ a + f, reordered))
  twin$setData(dbarts::dbartsData(y ~ a + f, byLabel(reordered)))
  expect_identical(sampler$data@x, twin$data@x, info = ordered)
  copied <- sampler$copy()
  expect_equal(fitsOf(copied), fitsOf(sampler), info = ordered)
}

# more levels declared, no row at them: taken, the record keeping six
sampler <- samplerOf()
sampler$setData(dbarts::dbartsData(
  y ~ a + f,
  frameOf(letters[1:8], letters[1:6], 2L)
))
expect_identical(attr(sampler$data@x, "factor.levels")[[2L]], letters[1:6])
copied <- sampler$copy()
expect_equal(fitsOf(copied), fitsOf(sampler))

# a row at a label the sampler never had is refused by name, a level replaced
# at the same count, renamed levels and a test row alike, and the sampler is
# left as it was
sampler <- samplerOf()
sampler$storeState()
held <- untouched(sampler)
expect_error(
  sampler$setData(dbarts::dbartsData(
    y ~ a + f,
    frameOf(letters[c(1L, 3:7)], letters[c(1L, 3:7)], 2L)
  )),
  pattern = "column 'f' has label 'g' not among its training levels"
)
expect_identical(untouched(sampler), held)
expect_error(
  sampler$setData(dbarts::dbartsData(
    y ~ a + f,
    frameOf(letters[1:8], letters[1:8], 2L)
  )),
  pattern = "has labels 'g' and 'h' not among"
)
expect_error(
  sampler$setData(dbarts::dbartsData(
    y ~ a + f,
    frameOf(LETTERS[1:6], LETTERS[1:6], 2L)
  )),
  pattern = "column 'f' has labels '[A-F]'"
)
expect_error(
  sampler$setData(dbarts::dbartsData(
    y ~ a + f,
    frameOf(letters[1:7], letters[1:6], 2L),
    test = frameOf(letters[1:7], "g", 3L)[1:5, ]
  )),
  pattern = "column 'f' has label 'g' not among its training levels"
)
expect_identical(untouched(sampler), held)

# columns are matched by name: y ~ g + f on a sampler of y ~ f + g, labels the
# two share, is read as y ~ f + g. Names that are not the sampler's are read by
# position, as a whole-frame setPredictor reads them
set.seed(1L)
yes.fl <- data.frame(
  f = factor(sample(c("no", "yes"), n.fl, TRUE)),
  g = factor(sample(c("no", "yes"), n.fl, TRUE))
)
yes.fl$y <- 3 * (yes.fl$f == "yes") + rnorm(n.fl, 0, 0.3)
pairOf <- function() {
  set.seed(4L)
  sampler <- dbarts::dbarts(y ~ f + g, yes.fl, control = control.fl)
  set.seed(5L)
  invisible(sampler$run(30L, 10L))
  sampler
}
sampler <- pairOf()
twin <- pairOf()
sampler$setData(dbarts::dbartsData(y ~ g + f, yes.fl, test = yes.fl[1:9, ]))
twin$setData(dbarts::dbartsData(y ~ f + g, yes.fl, test = yes.fl[1:9, ]))
expect_identical(sampler$data@x, twin$data@x)
expect_identical(sampler$data@x.test, twin$data@x.test)
expect_identical(fitsOf(sampler), fitsOf(twin))
renamed <- yes.fl
names(renamed)[2L] <- "h"
sampler$setData(dbarts::dbartsData(y ~ f + h, renamed))
expect_identical(sampler$data@x$dense, twin$data@x$dense)

# an NA level is matched to the NA level, and a missing value stays missing
levels.na <- c("a", "b", "c", NA)
codes.na <- rep_len(c(1L, 2L, 3L, 4L, NA, 4L, 2L), n.fl)
frame.na <- data.frame(
  y = a.fl + (codes.na %in% 4L) * 3 + rnorm(n.fl, 0, 0.3),
  a = a.fl,
  f = structure(codes.na, levels = levels.na, class = "factor")
)
sampler <- samplerOf(frame = frame.na)
reversed <- frame.na
reversed$f <- structure(
  5L - codes.na,
  levels = rev(levels.na),
  class = "factor"
)
sampler$setData(dbarts::dbartsData(y ~ a + f, reversed))
expect_identical(sampler$data@x$dense[[2L]], frame.na$f)

# two factor columns with other tables, handed in another order: each column
# by its name and each row by its label, training and test, as a twin handed
# the rows re-coded by hand, on each kind of factor
twoOf <- function(levels.f, levels.g, seed, ordered) {
  frame <- frameOf(levels.f, levels.f, seed, ordered)
  set.seed(seed + 10L)
  frame$g <- factor(sample(levels.g, n.fl, TRUE), levels.g, ordered = ordered)
  frame
}
for (ordered in c(FALSE, TRUE)) {
  start <- twoOf(letters[1:6], c("u", "v", "w"), 1L, ordered)
  pairs <- lapply(1:2, function(i) {
    set.seed(4L)
    sampler <- dbarts::dbarts(y ~ a + f + g, start, control = control.fl)
    set.seed(5L)
    invisible(sampler$run(30L, 10L))
    sampler
  })
  handed <- twoOf(letters[c(1L, 3:6)], c("w", "u"), 2L, ordered)
  test <- twoOf(letters[c(6L, 4L, 3L)], c("w", "u"), 3L, ordered)[1:10, ]
  hand <- function(frame) {
    frame <- byLabel(frame)
    frame$g <- factor(
      as.character(frame$g),
      c("u", "v", "w"),
      ordered = ordered
    )
    frame
  }
  pairs[[1L]]$setData(dbarts::dbartsData(y ~ g + f + a, handed, test = test))
  pairs[[2L]]$setData(dbarts::dbartsData(
    y ~ a + f + g,
    hand(handed),
    test = hand(test)
  ))
  expect_identical(pairs[[1L]]$data@x, pairs[[2L]]$data@x, info = ordered)
  expect_identical(
    pairs[[1L]]$data@x.test,
    pairs[[2L]]$data@x.test,
    info = ordered
  )
  set.seed(9L)
  result <- pairs[[1L]]$run(0L, 3L)
  set.seed(9L)
  expect_identical(result, pairs[[2L]]$run(0L, 3L), info = ordered)
}

# a numeric and a factor column handed in another order: taken by name
sampler <- samplerOf()
twin <- samplerOf()
reordered <- frameOf(rev(letters[1:6]), letters[1:6], 3L)
sampler$setData(dbarts::dbartsData(y ~ f + a, reordered))
twin$setData(dbarts::dbartsData(y ~ a + f, byLabel(reordered)))
expect_identical(sampler$data@x, twin$data@x)
expect_identical(sampler$data@varTypes, twin$data@varTypes)

# a factor column handed over as bare codes with no level table is refused in
# the words of dec-B360, and the sampler is left as it was
sampler <- samplerOf()
sampler$storeState()
held <- untouched(sampler)
frame <- frameOf(letters[1:6], letters[1:6], 2L)
codes <- cbind(a = frame$a, f = as.double(as.integer(frame$f)) - 1)
attr(codes, "varTypes") <- c(0L, 1L)
expect_error(
  sampler$setData(dbarts::dbartsData(codes, frame$y)),
  pattern = "the data hold factor 'f' as numbers with no level table"
)
expect_identical(untouched(sampler), held)
# a matrix of codes carrying its own table is read by label, and a number that
# is no code of its table is left to the bridge's refusal
attr(codes, "factor.levels") <- list(NULL, rev(letters[1:6]))
sampler$setData(dbarts::dbartsData(codes, frame$y))
expect_identical(
  sampler$data@x[, 2L],
  as.double(6L - as.integer(frame$f))
)
# a missing value there stays missing
missing <- codes
missing[c(3L, 9L), 2L] <- NA
sampler$setData(dbarts::dbartsData(missing, frame$y))
expect_identical(which(is.na(sampler$data@x[, 2L])), c(3L, 9L))
expect_identical(
  sampler$data@x[-c(3L, 9L), 2L],
  as.double(6L - as.integer(frame$f))[-c(3L, 9L)]
)
codes[1L, 2L] <- 7
expect_error(
  sampler$setData(dbarts::dbartsData(codes, frame$y)),
  pattern = "categorical predictor values must be existing category codes"
)

# a sparse source is the bridge's, which refuses it in its own words
sampler <- samplerOf()
sparse <- data.frame(a = a.fl)
sparse$f <- dbarts::sparseFactor(factor(
  sample(letters[1:6], n.fl, TRUE),
  levels = rev(letters[1:6])
))
expect_error(
  sampler$setData(dbarts::dbartsData(sparse, frame$y)),
  pattern = "requires a dense predictor matrix"
)
# and so is a sparse test set beside training rows whose levels are re-coded
dropped <- frameOf(letters[c(1L, 3:6)], letters[c(1L, 3:6)], 2L)
test <- data.frame(a = a.fl[1:10])
test$f <- dbarts::sparseFactor(factor(
  as.character(dropped$f[1:10]),
  levels = levels(dropped$f)
))
expect_error(
  sampler$setData(dbarts::dbartsData(dropped[c("a", "f")], dropped$y, test)),
  pattern = "requires a dense test matrix"
)

# a sampler made from a matrix of codes records no table, and takes codes as
# before
codes <- cbind(a = a.fl, f = as.double(rep_len(0:5, n.fl)))
attr(codes, "varTypes") <- c(0L, 1L)
set.seed(4L)
sampler <- dbarts::dbarts(codes, frame$y, control = control.fl)
codes[, 2L] <- rep_len(c(2, 0, 1), n.fl)
sampler$setData(dbarts::dbartsData(codes, frame$y))
expect_identical(sampler$data@x[, 2L], codes[, 2L])
