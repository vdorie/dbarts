# What the two binary studies' harnesses share: link helpers, the data root and
# its downloads, seeds, atomic and single-precision files, the real causal
# designs (IHDP, ACIC 2016, Twins and the three real causal datasets), the
# ACIC 2016 surface mapping, the four-arm design, and the draw summaries.
# Moved verbatim from benchmarks/R/binary-interval-study.R, which sources it,
# as does benchmarks/R/binary-k-mixing-need.R; a change here changes both
# studies' records. The sourcing script defines harnessDir() and sources this
# file before benchmarks/R/binary-datasets.R.

levelsPointwise <- c(0.5, 0.8, 0.9, 0.95)
interiorRange <- c(0.01, 0.99) # a probability inside it is interior

linkFunction <- function(link) {
  switch(link, probit = pnorm, logistic = plogis, stop("no link ", link))
}
linkQuantile <- function(link) {
  switch(link, probit = qnorm, logistic = qlogis, stop("no link ", link))
}
links <- c("probit", "logistic")

dataRoot <- function() {
  dir <- Sys.getenv("DBARTS_BIS_DATA")
  if (!nzchar(dir)) {
    stop("DBARTS_BIS_DATA must name the data root")
  }
  normalizePath(dir, mustWork = TRUE)
}

if (
  !nzchar(Sys.getenv("DBARTS_BENCH_DATA")) &&
    nzchar(Sys.getenv("DBARTS_BIS_DATA"))
) {
  Sys.setenv(DBARTS_BENCH_DATA = file.path(dataRoot(), "uci"))
}

ensureDir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

# Written once, atomically, so a killed job never leaves a half file that a
# reader takes for a result.
saveAtomic <- function(object, path) {
  ensureDir(dirname(path))
  temp <- paste0(path, ".part", Sys.getpid())
  saveRDS(object, temp)
  file.rename(temp, path)
  invisible(path)
}

# Effect draws are kept as single precision, a flat little-endian file
# beside a small rds holding the dimensions; a path ending in .gz is
# compressed.
openDraws <- function(path, mode) {
  if (grepl("[.]gz$", path)) {
    gzfile(path, mode, compression = 6L)
  } else {
    file(path, mode)
  }
}

writeSingle <- function(x, path) {
  ensureDir(dirname(path))
  con <- openDraws(path, "wb")
  on.exit(close(con))
  writeBin(as.numeric(x), con, size = 4L, endian = "little")
  invisible(dim(x))
}

readSingle <- function(path, dims) {
  con <- openDraws(path, "rb")
  on.exit(close(con))
  x <- readBin(con, "numeric", n = prod(dims), size = 4L, endian = "little")
  array(x, dims)
}

## ------------------------------------------------------------------ seeds

# A stable integer seed from the names of what it seeds, so that every seed in
# the study is a stated function of its cell, replication and role and does
# not depend on the order jobs run in.
seedOf <- function(...) {
  key <- paste(..., sep = "|")
  codes <- utf8ToInt(key)
  h <- 0
  for (code in codes) {
    h <- (h * 131 + code) %% 2147483629
  }
  as.integer(h + 1)
}

withSeed <- function(seed, expr) {
  old <- if (exists(".Random.seed", globalenv())) {
    get(".Random.seed", globalenv())
  } else {
    NULL
  }
  on.exit(
    if (is.null(old)) {
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", old, globalenv())
    }
  )
  set.seed(
    seed,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  expr
}


## ------------------------------------------------------------ downloads

# Each source lands in a new directory of its own under the data root, is
# checked against its recorded sha256, and is read from there by R running
# elsewhere. An archive is extracted into a sibling directory.
bisSources <- list(
  twins = list(
    dir = "twins",
    files = list(
      "twin_pairs_X_3years_samesex.csv" = "6ec6263e678f67205d584ebda083173c1ec44fb46ea7cfaafe3c87a31f612b86",
      "twin_pairs_T_3years_samesex.csv" = "fbe6eb1cad27794c61cbbd8b11b72a26e96549e59d64fb88b07b30aaf2bae07d",
      "twin_pairs_Y_3years_samesex.csv" = "d03428d23ed308673f9151e3934aa7712df27c1fed244ed44ecc80a98a60836b",
      "covar_desc.txt" = "51d126a3a16a7f46bfe42536bd0ae1e2e31dc840b0e5000b595b76b7aee672c0"
    ),
    url = "https://raw.githubusercontent.com/AMLab-Amsterdam/CEVAE/master/datasets/TWINS/"
  ),
  acic2019low = list(
    dir = "acic2019/lowD-test",
    files = list(
      "TestDatasets_lowD_Dec28.zip" = "281c7723026ab81f5abbbe9ee4c9bc17e2c4cc95baccf8cb5c0e4d3f866d39d6"
    ),
    url = "https://www.dropbox.com/s/qaj6fjbzorzmwpp/",
    query = "?dl=1",
    extract = "acic2019/lowD-x"
  ),
  acic2019high = list(
    dir = "acic2019/highD-test",
    files = list(
      "TestDatasets_highD_Dec28.zip" = "9734a8489dbb780c7d373bd853b1862a0a8c233d12826a45b44b07530e69d45c"
    ),
    url = "https://www.dropbox.com/s/7st5ttdihk6dzfz/",
    query = "?dl=1",
    extract = "acic2019/highD-x"
  )
)

sha256Of <- function(path) {
  if (exists("sha256sum", where = asNamespace("tools"))) {
    return(unname(get("sha256sum", asNamespace("tools"))(path)))
  }
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", out[1L])
}

verifiedFile <- function(source, name) {
  path <- file.path(dataRoot(), bisSources[[source]]$dir, name)
  if (!file.exists(path)) {
    stop("missing ", path, "; run `fetch`")
  }
  want <- bisSources[[source]]$files[[name]]
  got <- sha256Of(path)
  if (!identical(got, want)) {
    stop("sha256 mismatch for ", path, ": ", got)
  }
  path
}

readCausal <- function(name) {
  readRDS(file.path(dataRoot(), "causal", paste0(name, ".rds")))
}

# B2's three datasets: binary outcome y, treatment z, baseline covariates x.
causalDatasets <- list(
  lindner = function() {
    d <- readCausal("lindner")
    x <- d[, c(
      "stent",
      "height",
      "female",
      "diabetic",
      "acutemi",
      "ejecfrac",
      "ves1proc"
    )]
    list(x = x, z = as.integer(d$abcix), y = as.integer(!d$sixMonthSurvive))
  },
  actg175 = function() {
    d <- readCausal("ACTG175")
    d <- d[d$arms %in% c(0L, 1L), ]
    # baseline only: the follow-up counts, censoring time and off-treatment
    # flag are outcomes; zprior is constant among these arms
    x <- d[, c(
      "age",
      "wtkg",
      "hemo",
      "homo",
      "drugs",
      "karnof",
      "oprior",
      "z30",
      "preanti",
      "race",
      "gender",
      "str2",
      "strat",
      "symptom",
      "cd40",
      "cd80"
    )]
    x$strat <- factor(x$strat)
    list(x = x, z = as.integer(d$arms == 1L), y = as.integer(d$cens))
  },
  lalonde = function() {
    d <- readCausal("lalonde")
    x <- d[, c(
      "age",
      "educ",
      "black",
      "hispan",
      "married",
      "nodegree",
      "re74",
      "re75"
    )]
    list(x = x, z = as.integer(d$treat), y = as.integer(d$re78 > 0))
  }
)

# npci's IHDP construction: every treated child of a nonwhite mother
# dropped, 25 covariates, the continuous ones standardized.
ihdpData <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) {
      return(cache)
    }
    path <- Sys.getenv("DBARTS_IHDP")
    if (!nzchar(path)) {
      path <- file.path(dataRoot(), "ihdp", "ihdp.RData")
    }
    e <- new.env()
    load(path, envir = e)
    d <- subset(e$ihdp, treat != 1 | momwhite != 0)
    keepNames <- c(
      "bw",
      "b.head",
      "preterm",
      "birth.o",
      "nnhealth",
      "momage",
      "sex",
      "twin",
      "b.marr",
      "mom.lths",
      "mom.hs",
      "mom.scoll",
      "cig",
      "first",
      "booze",
      "drugs",
      "work.dur",
      "prenatal",
      "ark",
      "ein",
      "har",
      "mia",
      "pen",
      "tex",
      "was"
    )
    x <- d[, keepNames]
    for (name in keepNames) {
      if (is.numeric(x[[name]]) && length(unique(x[[name]])) > 2L) {
        x[[name]] <- (x[[name]] - mean(x[[name]])) / sd(x[[name]])
      }
    }
    rownames(x) <- NULL
    cache <<- list(
      x = x,
      z = as.integer(d$treat),
      noComparison = d$momwhite == 0
    )
    cache
  }
})

acic2016Input <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- aciccomp2016::input_2016
    }
    cache
  }
})

# ACIC 2016's generated surfaces for one setting and iteration, on all 4,802
# rows. dgp_2016 seeds R's generator itself, so it runs inside withSeed.
acic2016Surfaces <- function(setting, iteration) {
  withSeed(1L, {
    d <- aciccomp2016::dgp_2016(acic2016Input(), setting, iteration)
    list(
      z = as.integer(d$z),
      mu0 = as.numeric(d$mu.0),
      mu1 = as.numeric(d$mu.1)
    )
  })
}

# Twins (Louizos et al. 2017): same-sex pairs both under 2,000 g. The lighter
# twin is _0 and the heavier _1; the heavier is the treated one. Identifiers
# and the twin-specific birth orders are not covariates.
twinsData <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) {
      return(cache)
    }
    read <- function(name) {
      utils::read.csv(verifiedFile(
        "twins",
        sprintf("twin_pairs_%s_3years_samesex.csv", name)
      ))
    }
    X <- read("X")
    Tw <- read("T")
    Y <- read("Y")
    keep <- Tw$dbirwt_0 < 2000 & Tw$dbirwt_1 < 2000
    drop <- c(
      "X",
      "Unnamed..0",
      "infant_id_0",
      "infant_id_1",
      "bord_0",
      "bord_1"
    )
    x <- X[keep, setdiff(names(X), drop)]
    rownames(x) <- NULL
    cache <<- list(
      x = x,
      y0 = as.integer(Y$mort_0[keep]),
      y1 = as.integer(Y$mort_1[keep])
    )
    cache
  }
})

# Latent surfaces are standardized on the controls (arm 0) by the control
# surface's mean and standard deviation, and mapped to p = invLink(a eta + c), F
# the fitted link, c set so the controls' mean probability is the base rate.
mapSurfaces <- function(etas, control, a, rate, link) {
  invLink <- linkFunction(link)
  center <- mean(etas[[1L]][control])
  scale <- sd(etas[[1L]][control])
  std <- lapply(etas, function(e) (e - center) / scale)
  c0 <- uniroot(
    function(c) mean(invLink(a * std[[1L]][control] + c)) - rate,
    c(-40, 40),
    tol = 1e-10
  )$root
  lapply(std, function(e) invLink(a * e + c0))
}

subsampleRows <- function(n, size, seed) {
  if (size >= n) {
    return(seq_len(n))
  }
  withSeed(seed, sort(sample.int(n, size)))
}

# The high-dimensional cells' 192 permuted copies: copy j is column
# ((j - 1) mod 58) + 1, its rows permuted independently.
permutedCopies <- function(x, count, seed) {
  withSeed(seed, {
    copies <- lapply(seq_len(count), function(j) {
      col <- x[[((j - 1L) %% ncol(x)) + 1L]]
      col[sample.int(length(col))]
    })
    names(copies) <- sprintf("perm_%03d", seq_len(count))
    as.data.frame(copies, stringsAsFactors = FALSE)
  })
}

# The k study's split rule, one seed: 80/20, or 4,000 and 1,000 out of a
# dataset of more than 5,000 rows, redrawn until training holds both classes.
kStudySplit <- function(n, y, seed, large) {
  nTrain <- if (large) subsampleTrain else floor(0.8 * n)
  for (attempt in seq_len(100L)) {
    set.seed(seed + 1000000L * (attempt - 1L))
    if (large) {
      drawn <- sample.int(n, nTrain + subsampleTest)
      train <- drawn[seq_len(nTrain)]
      test <- drawn[nTrain + seq_len(subsampleTest)]
    } else {
      train <- sample.int(n, nTrain)
      test <- setdiff(seq_len(n), train)
    }
    if (length(unique(y[train])) == 2L) break
  }
  list(train = train, test = test)
}

# A6: 2,000 units, four arms from a multinomial logistic model on the
# standardized numeric covariates, six covariates an arm with N(0, 0.35^2)
# coefficients, arm 3 pushed down by 4 on the logit where x_27 is above its
# median (weak overlap there). Surfaces: arm 0 row 31's control surface, arm 1
# row 31's treated surface, arm 2 arm 1 plus the effect of an independent
# draw of row 31's generator (iteration + 50), arm 3 arm 0.
a6Base <- function(rep) {
  rows <- subsampleRows(4802L, 2000L, seedOf("A6-sub", rep))
  x <- acic2016Input()[rows, , drop = FALSE]
  rownames(x) <- NULL
  s <- acic2016Surfaces(31L, rep)
  s2 <- acic2016Surfaces(31L, rep + 50L)
  numericCols <- names(x)[vapply(x, is.numeric, logical(1L))]
  xs <- scale(as.matrix(x[, numericCols]))
  xs[!is.finite(xs)] <- 0
  arm <- withSeed(seedOf("A6-arm", rep), {
    eta <- matrix(0, nrow(xs), 4L)
    for (j in 2:4) {
      cols <- sample.int(ncol(xs), 6L)
      eta[, j] <- drop(xs[, cols] %*% rnorm(6L, 0, 0.35))
    }
    eta[, 4L] <- eta[, 4L] - 4 * (x$x_27 > median(x$x_27))
    prob <- exp(eta) / rowSums(exp(eta))
    u <- runif(nrow(prob))
    rowSums(u > t(apply(prob, 1L, cumsum))) + 0L
  })
  arm <- pmin(arm, 3L)
  list(
    x = x,
    arm = arm,
    mu0 = s$mu0[rows],
    mu1 = s$mu1[rows],
    delta2 = (s2$mu1 - s2$mu0)[rows]
  )
}

a6Probabilities <- function(base, cell, link) {
  par <- cell$par
  s0 <- base$mu0
  s1 <- base$mu1
  s2 <- base$mu1 + base$delta2
  surfaces <- switch(
    as.character(par$arrangement),
    "0" = list(s0, s0, s0, s0),
    "1" = list(s0, s1, s2, s0),
    "2" = list(s0, s1, s1, s0)
  )
  p <- mapSurfaces(surfaces, base$arm == 0L, par$a, 0.3, link)
  do.call(cbind, p)
}

# The six pairwise contrasts, in this order everywhere.
contrastPairs <- rbind(
  c(1L, 0L),
  c(2L, 0L),
  c(3L, 0L),
  c(2L, 1L),
  c(3L, 1L),
  c(3L, 2L)
)

probsPointwise <- sort(c((1 - levelsPointwise) / 2, (1 + levelsPointwise) / 2))

# Pointwise quantiles (type 7), posterior means and the rank-based
# simultaneous band of Besag et al. (1995) at a common tail level, in blocks
# of columns so that a large draw matrix is never sorted whole. Q is draws by
# quantities; `transform` maps it to the reported scale (the band's ranks are
# taken before it, so a probability's band is built on the latent scale);
# `groupSize` adds a band over each run of that many consecutive columns.
summarizeDraws <- function(
  Q,
  transform = identity,
  bandLevels = c(0.9, 0.95),
  groupSize = NULL,
  blockCols = 1000L
) {
  S <- nrow(Q)
  m <- ncol(Q)
  tail <- min(S, ceiling(0.06 * S) + 2L)
  if (!is.null(groupSize)) {
    blockCols <- max(groupSize, (blockCols %/% groupSize) * groupSize)
  }
  means <- numeric(m)
  q <- matrix(0, length(probsPointwise), m)
  low <- matrix(0, tail, m)
  high <- matrix(0, tail, m)
  extreme <- integer(S)
  unitBands <- if (is.null(groupSize)) {
    NULL
  } else {
    lapply(bandLevels, function(l) list(lo = numeric(m), hi = numeric(m)))
  }
  h <- (S - 1) * probsPointwise + 1
  hLow <- floor(h)
  hHigh <- pmin(hLow + 1L, S)
  hFrac <- h - hLow
  for (start in seq.int(1L, m, by = blockCols)) {
    cols <- seq.int(start, min(m, start + blockCols - 1L))
    b <- length(cols)
    B <- Q[, cols, drop = FALSE]
    o <- apply(B, 2L, order)
    if (!is.matrix(o)) {
      o <- matrix(o, S, b)
    }
    offset <- rep((seq_len(b) - 1L) * S, each = S)
    sorted <- transform(matrix(B[o + offset], S, b))
    means[cols] <- colMeans(transform(B))
    q[, cols] <- sorted[hLow, , drop = FALSE] +
      hFrac * (sorted[hHigh, , drop = FALSE] - sorted[hLow, , drop = FALSE])
    low[, cols] <- sorted[seq_len(tail), , drop = FALSE]
    high[, cols] <- sorted[seq.int(S - tail + 1L, S), , drop = FALSE]
    ranks <- matrix(0L, S, b)
    ranks[o + offset] <- rep.int(seq_len(S), b)
    ext <- pmax(ranks, S + 1L - ranks)
    extreme <- pmax(extreme, ext[cbind(seq_len(S), max.col(ext, "first"))])
    if (!is.null(groupSize)) {
      units <- b %/% groupSize
      ext3 <- array(ext, c(S, groupSize, units))
      eg <- ext3[, 1L, ]
      for (k in seq_len(groupSize)[-1L]) {
        eg <- pmax(eg, ext3[, k, ])
      }
      eg <- matrix(eg, S, units)
      for (l in seq_along(bandLevels)) {
        pos <- ceiling(bandLevels[l] * S)
        kStar <- apply(eg, 2L, function(v) sort.int(v, partial = pos)[pos])
        kCol <- rep(kStar, each = groupSize)
        unitBands[[l]]$lo[cols] <- sorted[cbind(S + 1L - kCol, seq_len(b))]
        unitBands[[l]]$hi[cols] <- sorted[cbind(kCol, seq_len(b))]
      }
    }
  }
  bands <- lapply(bandLevels, function(level) {
    kStar <- sort.int(extreme)[ceiling(level * S)]
    if (S + 1L - kStar > tail) {
      stop("band level outside the stored tails")
    }
    list(
      level = level,
      lo = low[S + 1L - kStar, ],
      hi = high[kStar - (S - tail), ],
      collapsed = kStar == S,
      kStar = kStar
    )
  })
  names(bands) <- as.character(bandLevels)
  list(mean = means, q = q, bands = bands, unitBands = unitBands, S = S)
}

interval <- function(summary, level) {
  lo <- match((1 - level) / 2, probsPointwise)
  hi <- match((1 + level) / 2, probsPointwise)
  if (is.na(lo) || is.na(hi)) {
    lo <- which.min(abs(probsPointwise - (1 - level) / 2))
    hi <- which.min(abs(probsPointwise - (1 + level) / 2))
  }
  list(lo = summary$q[lo, ], hi = summary$q[hi, ])
}

intervalScore <- function(lo, hi, truth, level) {
  alpha <- 1 - level
  (hi - lo) +
    (2 / alpha) * pmax(lo - truth, 0) +
    (2 / alpha) * pmax(truth - hi, 0)
}

safeMean <- function(x) if (length(x) == 0L) NA_real_ else mean(x)

# Whether each column's 90 percent interval (5th to 95th percentile) lies
# inside the interior range: a probability that single precision or the
# posterior itself does not pin at 0 or 1.
interiorColumns <- function(P) {
  q <- apply(P, 2L, stats::quantile, probs = c(0.05, 0.95), names = FALSE)
  q[1L, ] >= interiorRange[1L] & q[2L, ] <= interiorRange[2L]
}
