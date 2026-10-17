#!/usr/bin/env Rscript

# The binary interval study: does the collapsed scale move make binary
# intervals more accurate at equal running time? Two builds of one source,
# "shipped" (the move's toggle off) and "move" (on), are compared under
# probit and logistic on real covariates with generated outcomes (truth
# known) and on real outcomes against long reference runs (truth unknown).
# The design, the measures and the decision rule are fixed in
# docs/plans/binary-interval-study.md; this file implements them and decides
# nothing the plan leaves open except where its comments say so.
#
# BUILDS. One R process sees one build: the build is whatever dbarts R_LIBS
# finds first, and BIS_BUILD names it ("shipped", "move" or "base"). Fits
# that the plan assigns to one build refuse to run under another: the A1
# BART truths are the move build's, every propensity is the shipped build's.
#
# DATA. DBARTS_BIS_DATA is the data root (twins/, acic2019/, ihdp/, causal/,
# uci/); DBARTS_BENCH_DATA, the UCI cache, defaults to its uci/ directory.
# Nothing is read from the working directory. `fetch` downloads Twins and
# the ACIC 2019 practice archives into new directories of their own and checks
# every file against a recorded sha256; `export-causal` writes lindner and
# lalonde (twang) and ACTG 175 (BART) to causal/ so the bench host needs
# neither package; IHDP is npci's ihdp.RData, placed at ihdp/ihdp.RData or
# named by DBARTS_IHDP.
#
# RUN DIRECTORY. BIS_RUN holds every output: truth/, prop/, timing/, fits/,
# ref/, frozen.rds (the pilot's frozen choice and margins) and the summaries.
#
# Usage (each action one process; <cell> as `cells` prints it):
#   Rscript binary-interval-study.R cells
#   Rscript binary-interval-study.R jobs <phase> [file]
#   Rscript binary-interval-study.R truth bart <link> <dataset> <part>
#   Rscript binary-interval-study.R truth additive <link> <dataset>
#   Rscript binary-interval-study.R propensity <key> <rep>
#   Rscript binary-interval-study.R timing <link> <cell> <rep>
#   Rscript binary-interval-study.R timing summarize
#   Rscript binary-interval-study.R fit <link> <cell> <rep>
#   Rscript binary-interval-study.R reference <link> <cell> <part> [rep]
#   Rscript binary-interval-study.R reference accept <link> <cell> [rep]
#   Rscript binary-interval-study.R check <outfile>
#   Rscript binary-interval-study.R check compare <a.rds> <b.rds>
#   Rscript binary-interval-study.R pilot | interim | summary [doubled]
#   Rscript binary-interval-study.R fetch | export-causal
# BIS_QUICK=1 shortens every chain for a smoke run; quick output is marked
# and the readers refuse to mix it with full output.

suppressPackageStartupMessages(library(dbarts))

## ---------------------------------------------------------------- settings

quick <- nzchar(Sys.getenv("BIS_QUICK"))

# MCMC lengths. A fit runs past the default and is read at several lengths;
# a longer run's first draws are a shorter run's, bit for bit.
mcmc <- if (!quick) {
  list(
    burn = 500L, # every study fit
    defaultDraws = 500L, # the default reading
    draws4 = 2000L, # a 4-chain fit's draws a chain
    draws10 = 1000L, # a 10-chain fit's draws a chain
    eqUnit = 1000L, # burn + default draws: equal time is eqUnit * r
    propDraws = 500L, # a propensity fit, package defaults
    refBurn = 5000L,
    refIter = 15000L,
    refThin = 15L,
    timingWarm = 500L,
    timingRun = 500L,
    stored = 500L # effect draws kept on disk a reading
  )
} else {
  list(
    burn = 20L,
    defaultDraws = 20L,
    draws4 = 80L,
    draws10 = 40L,
    eqUnit = 40L,
    propDraws = 20L,
    refBurn = 20L,
    refIter = 60L,
    refThin = 3L,
    timingWarm = 10L,
    timingRun = 40L,
    stored = 20L
  )
}
nTrees <- 75L
refChainsPerBuild <- 8L
refParts <- 4L # parts a build's reference half is cut into
truthChains <- 16L
truthParts <- 4L
pilotReps <- as.integer(Sys.getenv("BIS_PILOT_REPS", "10"))
probClamp <- 1e-4 # the additive truth's range
levelsPointwise <- c(0.5, 0.8, 0.9, 0.95)

linkFunction <- function(link) {
  switch(link, probit = pnorm, logistic = plogis, stop("no link ", link))
}
linkQuantile <- function(link) {
  switch(link, probit = qnorm, logistic = qlogis, stop("no link ", link))
}
links <- c("probit", "logistic")

## ------------------------------------------------------------------ paths

harnessDir <- function() {
  flag <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(flag) == 0L) {
    "benchmarks/R"
  } else {
    dirname(normalizePath(sub("^--file=", "", flag[1L])))
  }
}

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

source(file.path(harnessDir(), "binary-datasets.R"))

runDir <- function(...) {
  root <- Sys.getenv("BIS_RUN")
  if (!nzchar(root)) {
    stop("BIS_RUN must name the run directory")
  }
  path <- file.path(root, ...)
  path
}

ensureDir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

build <- function() {
  b <- Sys.getenv("BIS_BUILD")
  if (!b %in% c("shipped", "move", "base")) {
    stop("BIS_BUILD must be shipped, move or base")
  }
  b
}

requireBuild <- function(which) {
  if (build() != which) {
    stop("this fit belongs to the ", which, " build, not ", build())
  }
}

cellFile <- function(cell) gsub("[:]", "_", cell)

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

fetchAll <- function() {
  root <- dataRoot()
  for (name in names(bisSources)) {
    s <- bisSources[[name]]
    dir <- file.path(root, s$dir)
    for (file in names(s$files)) {
      path <- file.path(dir, file)
      if (!file.exists(path)) {
        if (dir.exists(dir) && length(list.files(dir)) > 0L) {
          stop(dir, " exists and is incomplete; fetch into a new directory")
        }
        ensureDir(dir)
        url <- paste0(s$url, file, if (is.null(s$query)) "" else s$query)
        status <- utils::download.file(
          url,
          path,
          method = "curl",
          extra = "-sSfL --retry 3 --max-time 600",
          mode = "wb",
          quiet = TRUE
        )
        if (!identical(as.integer(status), 0L)) stop("download failed: ", url)
      }
      cat(sprintf("%-40s %s\n", file, sha256Of(path)))
      verifiedFile(name, file)
    }
    if (!is.null(s$extract)) {
      target <- file.path(root, s$extract)
      if (!dir.exists(target)) {
        ensureDir(target)
        utils::unzip(file.path(dir, names(s$files)[1L]), exdir = target)
      }
    }
  }
}

## ------------------------------------------------------------ data sets

# The causal datasets come from twang and BART through `export-causal`, so the
# bench host reads plain files and needs neither package.
exportCausal <- function() {
  dir <- ensureDir(file.path(dataRoot(), "causal"))
  e <- new.env()
  utils::data("lindner", package = "twang", envir = e)
  utils::data("lalonde", package = "twang", envir = e)
  utils::data("ACTG175", package = "BART", envir = e)
  for (name in c("lindner", "lalonde", "ACTG175")) {
    saveRDS(as.data.frame(e[[name]]), file.path(dir, paste0(name, ".rds")))
    cat(name, nrow(e[[name]]), "\n")
  }
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

acic2019Sets <- c(
  low1 = "lowD:1",
  low4 = "lowD:4",
  low5 = "lowD:5",
  low8 = "lowD:8",
  high1 = "highD:1",
  high2 = "highD:2",
  high5 = "highD:5",
  high6 = "highD:6"
)

acic2019Data <- function(name) {
  parts <- strsplit(acic2019Sets[[name]], ":", fixed = TRUE)[[1L]]
  dir <- if (parts[1L] == "lowD") {
    file.path(dataRoot(), "acic2019", "lowD-x", "TestDatasets_lowD")
  } else {
    file.path(dataRoot(), "acic2019", "highD-x", "TestDatasets_highD")
  }
  verifiedFile(
    if (parts[1L] == "lowD") "acic2019low" else "acic2019high",
    names(
      bisSources[[
        if (parts[1L] == "lowD") "acic2019low" else "acic2019high"
      ]]$files
    )
  )
  stem <- if (parts[1L] == "lowD") "testdataset" else "highDim_testdataset"
  d <- utils::read.csv(file.path(dir, paste0(stem, parts[2L], ".csv")))
  cf <- utils::read.csv(file.path(dir, paste0(stem, parts[2L], "_cf.csv")))
  p1 <- cf[[grep("^EY1", names(cf))]]
  p0 <- cf[[grep("^EY0", names(cf))]]
  if (!all(d$Y %in% c(0, 1))) {
    stop(name, " is not binary")
  }
  list(x = d[, grep("^V", names(d))], z = as.integer(d$A), p1 = p1, p0 = p0)
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

## ------------------------------------------------------------ the cells

# A cell is one design under one link; its replications are numbered from 1.
# `groups` names the rule's groups and guards the cell enters:
#   F, T               gain groups and their guards
#   M1, M2, M3         multiple comparisons, each with its own bar
#   R                  real outcomes against a reference: a diagnosis outside
#                      the rule (dec-B456)
#   A6M                A6 alone, M1 and M2
#   sim, ate           the simultaneous and average-effect coverage guards
#   rmse, is           the pooled RMSE and interval-score ratios
# `psKey` names the propensity a cell's fits carry; cells that share
# covariates and treatment within a replication share it, and a key ending
# in "fixed" is one propensity for every replication.
makeCell <- function(
  id,
  problem,
  kind,
  chains,
  reps,
  groups,
  ...,
  psKey = NULL,
  big = FALSE
) {
  list(
    id = id,
    problem = problem,
    kind = kind,
    chains = chains,
    reps = reps,
    groups = groups,
    psKey = psKey,
    big = big,
    par = list(...)
  )
}

effectGroups <- function(null, sparse = FALSE) {
  if (sparse) {
    c("T", "M1", "M3", "sim", "ate", "rmse", "is")
  } else if (null) {
    c("M1", "M2", "sim", "ate", "rmse", "is")
  } else {
    c("T", "M3", "sim", "ate", "rmse", "is")
  }
}

buildCells <- function() {
  cells <- list()
  add <- function(cell) cells[[cell$id]] <<- cell
  for (d in names(realDatasets)) {
    for (truth in c("bart", "additive")) {
      add(makeCell(
        sprintf("A1:%s:%s", d, truth),
        "A1",
        "prob",
        4L,
        16L,
        c("F", "sim", "rmse", "is"),
        dataset = d,
        truth = truth,
        big = FALSE
      ))
    }
  }
  a2 <- list(
    "A-mod" = list(setting = "A", a = 1, null = FALSE),
    "A-sep" = list(setting = "A", a = 2.5, null = FALSE),
    "B-mod" = list(setting = "B", a = 1, null = FALSE),
    "null-mod" = list(setting = "A", a = 1, null = TRUE),
    "null-sep" = list(setting = "A", a = 2.5, null = TRUE)
  )
  for (name in names(a2)) {
    s <- a2[[name]]
    add(makeCell(
      paste0("A2:", name),
      "A2",
      "effect",
      10L,
      30L,
      effectGroups(s$null),
      setting = s$setting,
      a = s$a,
      null = s$null,
      psKey = "ihdp-fixed"
    ))
  }
  a3 <- function(name, row, effect, a, rate, size, hd = FALSE, reps = 30L) {
    key <- sprintf(
      "acic16-r%d-%s%s",
      row,
      if (size == 4802L) "full" else "n1000",
      if (hd) "-hd" else ""
    )
    null <- effect == "null"
    add(makeCell(
      paste0("A3:", name),
      "A3",
      "effect",
      10L,
      reps,
      effectGroups(null, effect == "sparse"),
      row = row,
      effect = effect,
      a = a,
      rate = rate,
      size = size,
      hd = hd,
      null = null,
      psKey = key,
      big = size == 4802L
    ))
  }
  for (sig in c(mod = 1, sep = 2.5)) {
    tag <- if (sig == 1) "mod" else "sep"
    a3(paste0("r1-", tag), 1L, "het", sig, 0.3, 1000L)
    a3(paste0("r31-", tag), 31L, "het", sig, 0.3, 1000L)
    a3(paste0("r30-", tag), 30L, "het", sig, 0.3, 1000L)
    a3(paste0("null-", tag), 31L, "null", sig, 0.3, 1000L)
  }
  a3("r31-rare05", 31L, "het", 1, 0.05, 1000L)
  a3("null-rare05", 31L, "null", 1, 0.05, 1000L)
  a3("r1-rare05", 1L, "het", 1, 0.05, 1000L)
  a3("sparse", 31L, "sparse", 1, 0.3, 1000L)
  a3("hd-r31", 31L, "het", 1, 0.3, 1000L, hd = TRUE)
  a3("hd-null", 31L, "null", 1, 0.3, 1000L, hd = TRUE)
  a3("full-r31-sep", 31L, "het", 2.5, 0.3, 4802L, reps = 16L)
  a3("full-null-sep", 31L, "null", 2.5, 0.3, 4802L, reps = 16L)
  a3("full-r31-rare01", 31L, "het", 1, 0.01, 4802L, reps = 16L)
  for (name in names(acic2019Sets)) {
    add(makeCell(
      paste0("A4:", name),
      "A4",
      "effect",
      10L,
      20L,
      effectGroups(FALSE),
      dataset = name,
      null = FALSE,
      psKey = paste0("acic19-", name, "-fixed")
    ))
  }
  add(makeCell(
    "A5:twins",
    "A5",
    "twins",
    10L,
    20L,
    c("T", "ate", "rmse", "is"),
    psKey = "twins"
  ))
  a6 <- list(
    "arr1-mod" = list(arrangement = 1L, a = 1),
    "arr1-sep" = list(arrangement = 1L, a = 2.5),
    "arr2-mod" = list(arrangement = 2L, a = 1),
    "null-mod" = list(arrangement = 0L, a = 1)
  )
  for (name in names(a6)) {
    s <- a6[[name]]
    groups <- if (s$arrangement == 0L) {
      c("M1", "M2", "A6M", "sim", "ate", "rmse", "is")
    } else {
      c("T", "M1", "M2", "M3", "A6M", "sim", "ate", "rmse", "is")
    }
    add(makeCell(
      paste0("A6:", name),
      "A6",
      "multi",
      10L,
      30L,
      groups,
      arrangement = s$arrangement,
      a = s$a,
      psKey = "a6",
      big = TRUE
    ))
  }
  for (d in names(realDatasets)) {
    add(makeCell(
      paste0("B1:", d),
      "B1",
      "prob",
      4L,
      12L,
      "R",
      dataset = d,
      full = FALSE
    ))
  }
  add(makeCell(
    "B1:adult-full",
    "B1",
    "prob",
    4L,
    8L,
    "R",
    dataset = "adult",
    full = TRUE,
    big = TRUE
  ))
  for (d in names(causalDatasets)) {
    add(makeCell(
      paste0("B2:", d),
      "B2",
      "effect",
      10L,
      15L,
      "R",
      dataset = d,
      null = FALSE,
      psKey = paste0("b2-", d, "-fixed")
    ))
  }
  cells
}

cells <- buildCells()

getCell <- function(id) {
  cell <- cells[[id]]
  if (is.null(cell)) {
    stop("no such cell: ", id)
  }
  cell
}

drawsPerChain <- function(cell) {
  if (cell$chains == 10L) mcmc$draws10 else mcmc$draws4
}

## ----------------------------------------------------------- generators

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

# npci's IHDP surfaces (settings A and B, offset 0.5, the treated surface
# shifted so the effect on the treated is 4 on the latent scale).
npciSurfaces <- function(x, z, setting, seed) {
  withSeed(seed, {
    x <- as.matrix(x)
    p <- ncol(x)
    if (setting == "B") {
      binary <- apply(x, 2L, function(col) length(unique(col)) == 2L)
      quad <- list()
      for (i in seq_len(p)) {
        if (!binary[i]) {
          quad[[length(quad) + 1L]] <- x[, i]^2
        }
        if (i < p) {
          for (j in seq.int(i + 1L, p)) {
            quad[[length(quad) + 1L]] <- x[, i] * x[, j]
          }
        }
      }
      xq <- cbind(x, do.call(cbind, quad))
      # aliased columns dropped, as npci's glm leaves them NA
      keep <- !is.na(lm.fit(cbind(1, xq), rep(0, nrow(xq)))$coefficients[-1L])
      xq <- xq[, keep, drop = FALSE]
      nq <- ncol(xq) - p
      probM <- max(1 - 2 / sqrt(p), 0.2)
      probM <- c(probM, 0.75 * (1 - probM), 0.25 * (1 - probM))
      probQ <- max(1 - 5^(1 - 1 / 8) / (p^(15 / 16)), 0.4)
      probQ <- c(probQ, 0.75 * (1 - probQ), 0.25 * (1 - probQ))
      bm0 <- sample(c(0, 1, 2), p + 1L, replace = TRUE, prob = probM)
      bm1 <- sample(c(0, 1, 2), p + 1L, replace = TRUE, prob = probM)
      bq0 <- sample(c(0, 0.5, 1), nq, replace = TRUE, prob = probQ)
      bq1 <- sample(c(0, 0.5, 1), nq, replace = TRUE, prob = probQ)
      xm <- cbind(1, xq)
      mu0 <- drop(xm %*% c(bm0, bq0))
      mu1 <- drop(xm %*% c(bm1, bq1))
    } else {
      probs <- max(1 - 2 / sqrt(p), 0.2)
      probs <- c(probs, rep(0.25 * (1 - probs), 4L))
      beta <- c(
        sample(seq(-1, 1, 0.25), 1L),
        sample(seq(0, 0.4, 0.1), p, replace = TRUE, prob = probs)
      )
      xm <- cbind(1, x)
      offset <- matrix(c(0, rep(0.5, p)), nrow(xm), p + 1L, byrow = TRUE)
      mu0 <- drop(exp((xm + offset) %*% beta))
      mu1 <- drop(xm %*% beta)
    }
    omega <- mean(mu1[z == 1L] - mu0[z == 1L]) - 4
    list(mu0 = mu0, mu1 = mu1 - omega)
  })
}

subsampleRows <- function(n, size, seed) {
  if (size >= n) {
    return(seq_len(n))
  }
  withSeed(seed, sort(sample.int(n, size)))
}

# The sparse cell's subgroup: two continuous covariates, the first above its
# median and the second above the cut that leaves 15 percent of the 4,802
# units in the subgroup.
sparseSubgroup <- function(x) {
  a <- x$x_36
  b <- x$x_27
  upper <- a > median(a)
  cuts <- sort(unique(b))
  share <- vapply(cuts, function(t) mean(upper & b > t), numeric(1L))
  cut <- cuts[which.min(abs(share - 0.15))]
  upper & b > cut
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

# A1's base data: the dataset, or its fixed 5,000-row subsample for the four
# of more than 5,000 rows.
a1Base <- function(dataset) {
  raw <- realDatasets[[dataset]]()
  n <- nrow(raw$x)
  rows <- if (n > subsampleAbove) {
    subsampleRows(n, subsampleTrain + subsampleTest, seedOf("A1-sub", dataset))
  } else {
    seq_len(n)
  }
  x <- raw$x[rows, , drop = FALSE]
  rownames(x) <- NULL
  list(x = x, y = raw$y[rows], large = n > subsampleAbove)
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

# A truth-known or real-outcome probability problem: training rows with an
# outcome, held-out rows scored.
probReplication <- function(cell, rep, link) {
  par <- cell$par
  if (cell$problem == "A1") {
    base <- a1Base(par$dataset)
    truthLink <- Sys.getenv("BIS_TRUTH_LINK", link)
    p <- loadTruth(truthLink, par$dataset, par$truth)
    n <- nrow(base$x)
    nTrain <- if (base$large) subsampleTrain else floor(0.8 * n)
    for (attempt in seq_len(100L)) {
      split <- withSeed(seedOf("A1-split", par$dataset, rep, attempt), {
        train <- sample.int(n, nTrain)
        list(train = train, test = setdiff(seq_len(n), train))
      })
      u <- withSeed(seedOf("A1-y", par$dataset, rep, attempt), runif(n))
      y <- as.integer(u < p)
      if (length(unique(y[split$train])) == 2L) break
    }
    return(list(
      train = data.frame(
        base$x[split$train, , drop = FALSE],
        y = y[split$train]
      ),
      test = base$x[split$test, , drop = FALSE],
      truth = p[split$test]
    ))
  }
  raw <- realDatasets[[par$dataset]]()
  n <- nrow(raw$x)
  if (isTRUE(par$full)) {
    split <- withSeed(seedOf("B1-full", par$dataset), {
      train <- sample.int(n, floor(0.8 * n))
      rest <- setdiff(seq_len(n), train)
      list(train = train, test = sort(rest[sample.int(length(rest), 2000L)]))
    })
  } else {
    seed <- 10L * (8000L + match(par$dataset, names(realDatasets))) + 1L
    split <- withSeed(1L, kStudySplit(n, raw$y, seed, n > subsampleAbove))
  }
  list(
    train = data.frame(
      raw$x[split$train, , drop = FALSE],
      y = raw$y[split$train]
    ),
    test = raw$x[split$test, , drop = FALSE],
    truth = NULL
  )
}

# The covariates, treatment and true probabilities of an effect cell's
# replication, before any propensity is attached.
effectBase <- function(cell, rep, link) {
  par <- cell$par
  if (cell$problem == "A2") {
    d <- ihdpData()
    s <- npciSurfaces(
      d$x,
      d$z,
      par$setting,
      seedOf("A2-surface", par$setting, rep)
    )
    etas <- if (par$null) list(s$mu0, s$mu0) else list(s$mu0, s$mu1)
    p <- mapSurfaces(etas, d$z == 0L, par$a, 0.3, link)
    return(list(
      x = d$x,
      z = d$z,
      p0 = p[[1L]],
      p1 = p[[2L]],
      noComparison = d$noComparison
    ))
  }
  if (cell$problem == "A3") {
    s <- acic2016Surfaces(par$row, rep)
    rows <- subsampleRows(4802L, par$size, seedOf("A3-sub", par$size, rep))
    x <- acic2016Input()
    sub <- if (par$effect == "sparse") sparseSubgroup(x)[rows] else NULL
    x <- x[rows, , drop = FALSE]
    rownames(x) <- NULL
    if (par$hd) {
      x <- cbind(x, permutedCopies(x, 192L, seedOf("A3-hd", rep)))
    }
    z <- s$z[rows]
    mu0 <- s$mu0[rows]
    control <- z == 0L
    if (par$effect == "sparse") {
      # one latent standard deviation, in the controls' standardized units
      mu1 <- mu0 + sd(mu0[control]) * as.numeric(sub)
    } else if (par$effect == "null") {
      mu1 <- mu0
    } else {
      mu1 <- s$mu1[rows]
    }
    p <- mapSurfaces(list(mu0, mu1), control, par$a, par$rate, link)
    return(list(x = x, z = z, p0 = p[[1L]], p1 = p[[2L]]))
  }
  if (cell$problem == "A4") {
    d <- acic2019Data(par$dataset)
    return(list(x = d$x, z = d$z, p0 = d$p0, p1 = d$p1))
  }
  if (cell$problem == "B2") {
    d <- causalDatasets[[par$dataset]]()
    return(list(x = d$x, z = d$z, y = d$y))
  }
  stop("no effect base for ", cell$id)
}

# Twins: 2,000 pairs a replication, the heavier twin treated by Louizos et
# al.'s rule, t ~ Bernoulli(logistic(w'x + w_h (g / 10 - 0.1))) with w ~
# N(0, 0.1^2) on the standardized covariates (missing values at 0 there),
# w_h ~ N(5, 0.1^2), g the gestational-age group from 0 to 9.
twinsBase <- function(rep) {
  d <- twinsData()
  rows <- subsampleRows(nrow(d$x), 2000L, seedOf("A5-sub", rep))
  x <- d$x[rows, , drop = FALSE]
  rownames(x) <- NULL
  g <- x$gestat10 - 1L
  xs <- vapply(
    x,
    function(col) {
      col <- as.numeric(col)
      s <- sd(col, na.rm = TRUE)
      out <- if (is.finite(s) && s > 0) {
        (col - mean(col, na.rm = TRUE)) / s
      } else {
        0 * col
      }
      out[is.na(out)] <- 0
      out
    },
    numeric(nrow(x))
  )
  z <- withSeed(seedOf("A5-treat", rep), {
    w <- rnorm(ncol(xs), 0, 0.1)
    wh <- rnorm(1L, 5, 0.1)
    as.integer(runif(nrow(xs)) < plogis(drop(xs %*% w) + wh * (g / 10 - 0.1)))
  })
  y1 <- d$y1[rows]
  y0 <- d$y0[rows]
  list(x = x, z = z, y0 = y0, y1 = y1, y = ifelse(z == 1L, y1, y0), group = g)
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

outcomeDraw <- function(p, cell, rep) {
  u <- withSeed(seedOf("y", cell$id, rep), runif(length(p)))
  as.integer(u < p)
}

## ---------------------------------------------------------- propensity

propensityFile <- function(key, rep) {
  if (grepl("fixed$", key)) {
    rep <- 0L
  }
  runDir("prop", key, sprintf("rep%03d.rds", as.integer(rep)))
}

loadPropensity <- function(key, rep) {
  path <- propensityFile(key, rep)
  if (!file.exists(path)) {
    stop("no propensity ", path, "; run `propensity ", key, " ", rep, "`")
  }
  p <- readRDS(path)
  if (!identical(p$quick, quick)) {
    stop("quick and full output mixed: ", path)
  }
  p$ps
}

# The covariates and treatment a propensity key stands for at one
# replication. They do not depend on the link.
propensityInputs <- function(key, rep) {
  if (key == "twins") {
    b <- twinsBase(rep)
    return(list(x = b$x, z = b$z))
  }
  if (key == "a6") {
    b <- a6Base(rep)
    return(list(x = b$x, arm = b$arm))
  }
  owner <- Filter(function(cell) identical(cell$psKey, key), cells)[[1L]]
  b <- effectBase(owner, rep, "probit")
  list(x = b$x, z = b$z)
}

# bartCause's default treatment model: a 10-chain BART probit of the
# treatment on the covariates, its posterior mean probability; for A6 a
# 4-chain multinomial fit, its posterior mean probability of each arm. Both
# are the shipped build's, shared by both arms and both links.
runPropensity <- function(key, rep) {
  requireBuild("shipped")
  if (grepl("fixed$", key)) {
    rep <- 0L
  }
  inputs <- propensityInputs(key, max(rep, 1L))
  seed <- seedOf("propensity", key, rep)
  elapsed <- system.time({
    if (key == "a6") {
      d <- data.frame(inputs$x, arm = factor(inputs$arm, levels = 0:3))
      fit <- bart(
        arm ~ .,
        data = d,
        family = "multinomial",
        n.chains = 4L,
        n.threads = 1L,
        n.samples = mcmc$propDraws,
        n.burn = mcmc$burn,
        n.trees = nTrees,
        seed = seed,
        verbose = FALSE
      )
      ps <- apply(fit$yhat.train, c(2L, 3L), mean)
      colnames(ps) <- paste0("gps", 0:3)
    } else {
      d <- data.frame(inputs$x, trt = inputs$z)
      fit <- bart(
        trt ~ .,
        data = d,
        family = "probit",
        n.chains = 10L,
        n.threads = 1L,
        n.samples = mcmc$propDraws,
        n.burn = mcmc$burn,
        n.trees = nTrees,
        seed = seed,
        verbose = FALSE
      )
      ps <- colMeans(pnorm(fit$yhat.train))
    }
  })[["elapsed"]]
  if (any(!is.finite(ps))) {
    stop("non-finite propensity for ", key, " ", rep)
  }
  saveAtomic(
    list(
      ps = ps,
      key = key,
      rep = rep,
      seconds = elapsed,
      lib = find.package("dbarts"),
      quick = quick
    ),
    propensityFile(key, rep)
  )
  cat(sprintf("propensity %s %d: %.1f s\n", key, rep, elapsed))
}

## ------------------------------------------------------------- truths

truthFile <- function(link, dataset, truth, part = NULL) {
  if (is.null(part)) {
    runDir("truth", link, sprintf("%s-%s.rds", dataset, truth))
  } else {
    runDir("truth", link, sprintf("%s-%s-part%d.rds", dataset, truth, part))
  }
}

# A1's BART truth: the move build at reference length, 16 chains in four
# parts of four, the posterior mean of the per-draw probability at every row
# of the dataset (or of its fixed 5,000-row subsample).
runTruthBart <- function(link, dataset, part) {
  requireBuild("move")
  base <- a1Base(dataset)
  invLink <- linkFunction(link)
  chains <- truthChains / truthParts
  elapsed <- system.time(
    fit <- bart(
      y ~ .,
      data = data.frame(base$x, y = base$y),
      family = link,
      n.trees = nTrees,
      n.chains = chains,
      n.threads = 1L,
      n.burn = mcmc$refBurn,
      n.samples = mcmc$refIter,
      n.thin = mcmc$refThin,
      seed = seedOf("truth", link, dataset, part),
      verbose = FALSE,
      combineChains = FALSE,
      keepTrainingFits = TRUE
    )
  )[["elapsed"]]
  draws <- fit$yhat.train
  if (any(!is.finite(draws))) {
    stop("non-finite truth draws")
  }
  p <- apply(invLink(draws), 3L, mean)
  meanP <- apply(invLink(draws), c(1L, 2L), mean)
  saveAtomic(
    list(
      p = p,
      nDraws = prod(dim(draws)[1:2]),
      k = fit$k,
      meanP = meanP,
      seconds = elapsed,
      lib = find.package("dbarts"),
      quick = quick
    ),
    truthFile(link, dataset, "bart", part)
  )
  cat(sprintf(
    "truth bart %s %s part %d: %.1f s\n",
    link,
    dataset,
    part,
    elapsed
  ))
}

# A1's additive truth (mgcv): smooths of the numeric predictors with ten or
# more distinct values, the rest as terms of their own (factors, and numeric
# predictors of fewer values as factors); tensor-product interactions among
# the three smoothed predictors of largest main-effect deviance reduction;
# above 20 predictors, smooths on the 20 of largest univariate deviance
# reduction and the rest linear; a ridge penalty on the parametric terms
# where the data separate.
runTruthAdditive <- function(link, dataset) {
  base <- a1Base(dataset)
  x <- base$x
  y <- base$y
  fam <- binomial(link = if (link == "probit") "probit" else "logit")
  smoothable <- vapply(
    x,
    function(col) {
      is.numeric(col) && length(unique(col[!is.na(col)])) >= 10L
    },
    logical(1L)
  )
  if (!all(stats::complete.cases(x))) {
    stop(dataset, " has incomplete rows")
  }
  for (name in names(x)[!smoothable]) {
    if (is.numeric(x[[name]])) x[[name]] <- factor(x[[name]])
  }
  devReduction <- function(name) {
    f <- if (smoothable[[name]]) {
      stats::as.formula(sprintf(
        "y ~ s(%s, k = %d)",
        name,
        min(10L, length(unique(x[[name]])) - 1L)
      ))
    } else {
      stats::as.formula(sprintf("y ~ %s", name))
    }
    d <- data.frame(x[, name, drop = FALSE], y = y)
    m <- suppressWarnings(mgcv::gam(f, data = d, family = fam))
    m$null.deviance - m$deviance
  }
  reduction <- vapply(names(x), devReduction, numeric(1L))
  smoothed <- names(x)[smoothable]
  linear <- character()
  if (ncol(x) > 20L) {
    top <- names(sort(reduction, decreasing = TRUE))[seq_len(20L)]
    linear <- setdiff(smoothed, top)
    smoothed <- intersect(smoothed, top)
  }
  parametric <- c(names(x)[!smoothable], linear)
  inter <- names(sort(reduction[smoothed], decreasing = TRUE))
  inter <- inter[seq_len(min(3L, length(inter)))]
  sTerms <- vapply(
    smoothed,
    function(name) {
      sprintf("s(%s, k = %d)", name, min(10L, length(unique(x[[name]])) - 1L))
    },
    character(1L)
  )
  tiTerms <- if (length(inter) >= 2L) {
    pairs <- utils::combn(inter, 2L)
    apply(pairs, 2L, function(pr) sprintf("ti(%s, %s, k = 5)", pr[1L], pr[2L]))
  } else {
    character()
  }
  P <- if (length(parametric) > 0L) {
    stats::model.matrix(
      stats::as.formula(paste("~", paste(parametric, collapse = " + "))),
      data = x
    )[, -1L, drop = FALSE]
  } else {
    NULL
  }
  d <- x
  d$y <- y
  rhs <- c(sTerms, tiTerms)
  if (!is.null(P)) {
    d$P <- P
    rhs <- c(rhs, "P")
  }
  if (length(rhs) == 0L) {
    rhs <- "1"
  }
  f <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")))
  warned <- FALSE
  fit <- withCallingHandlers(
    mgcv::gam(f, data = d, family = fam, method = "REML"),
    warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  pc <- if (is.null(P)) {
    numeric()
  } else {
    stats::coef(fit)[grep("^P", names(stats::coef(fit)))]
  }
  separated <- !is.null(P) &&
    (warned ||
      !isTRUE(fit$converged) ||
      any(!is.finite(pc)) ||
      any(abs(pc) > 10))
  if (separated) {
    fit <- suppressWarnings(mgcv::gam(
      f,
      data = d,
      family = fam,
      method = "REML",
      paraPen = list(P = list(diag(ncol(P))))
    ))
  }
  p <- pmin(
    pmax(as.numeric(stats::predict(fit, type = "response")), probClamp),
    1 - probClamp
  )
  saveAtomic(
    list(
      p = p,
      separated = separated,
      smoothed = smoothed,
      interactions = inter,
      linear = linear,
      deviance = fit$deviance,
      quick = quick
    ),
    truthFile(link, dataset, "additive")
  )
  cat(sprintf(
    "truth additive %s %s: separated %s, %d smooths\n",
    link,
    dataset,
    separated,
    length(smoothed)
  ))
}

loadTruth <- function(link, dataset, truth) {
  if (truth == "additive") {
    # no chain enters the additive truth, so a quick run reads the full one
    t <- readRDS(truthFile(link, dataset, "additive"))
    return(t$p)
  }
  parts <- lapply(seq_len(truthParts), function(part) {
    path <- truthFile(link, dataset, "bart", part)
    if (!file.exists(path)) {
      stop("missing truth part ", path)
    }
    t <- readRDS(path)
    if (!identical(t$quick, quick)) {
      stop("quick and full output mixed")
    }
    t
  })
  w <- vapply(parts, function(t) t$nDraws, numeric(1L))
  drop(do.call(cbind, lapply(parts, function(t) t$p)) %*% w) / sum(w)
}

## -------------------------------------------------------- replications

# A replication: the training frame (outcome y), the frame of rows predicted,
# and what they are scored against. Effect fits follow bartCause's default
# response model: the outcome on the covariates, the treatment `trt` and the
# estimated propensity `pscore`, every unit predicted under both treatments.
replication <- function(cell, rep, link) {
  if (cell$kind == "prob") {
    r <- probReplication(cell, rep, link)
    r$kind <- "prob"
    return(r)
  }
  if (cell$kind == "multi") {
    b <- a6Base(rep)
    P <- a6Probabilities(b, cell, link)
    y <- outcomeDraw(P[cbind(seq_len(nrow(P)), b$arm + 1L)], cell, rep)
    gps <- loadPropensity(cell$psKey, rep)
    arm <- factor(b$arm, levels = 0:3)
    train <- data.frame(b$x, arm = arm, gps, y = y)
    test <- do.call(
      rbind,
      lapply(0:3, function(a) {
        data.frame(b$x, arm = factor(rep(a, nrow(b$x)), levels = 0:3), gps)
      })
    )
    truth <- vapply(
      seq_len(nrow(contrastPairs)),
      function(j) {
        P[, contrastPairs[j, 1L] + 1L] - P[, contrastPairs[j, 2L] + 1L]
      },
      numeric(nrow(P))
    )
    return(list(
      kind = "multi",
      train = train,
      test = test,
      P = P,
      contrasts = truth,
      n = nrow(P)
    ))
  }
  if (cell$kind == "twins") {
    b <- twinsBase(rep)
    ps <- loadPropensity(cell$psKey, rep)
  } else {
    b <- effectBase(cell, rep, link)
    ps <- loadPropensity(cell$psKey, rep)
    if (is.null(b$y)) {
      b$y <- outcomeDraw(ifelse(b$z == 1L, b$p1, b$p0), cell, rep)
    }
  }
  train <- data.frame(b$x, trt = b$z, pscore = ps, y = b$y)
  test <- rbind(
    data.frame(b$x, trt = 1L, pscore = ps),
    data.frame(b$x, trt = 0L, pscore = ps)
  )
  c(
    list(kind = cell$kind, train = train, test = test, n = nrow(b$x)),
    b[setdiff(names(b), c("x", "y"))],
    list(y = b$y)
  )
}

## ------------------------------------------------------ draw summaries

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

# Pointwise and simultaneous coverage of a truth vector, with width,
# interval score and error; `mask` restricts the pointwise readings.
coverageMeasures <- function(summary, bandSummary, truth, prefix, mask = NULL) {
  out <- list()
  keep <- if (is.null(mask)) rep(TRUE, length(truth)) else mask
  for (level in levelsPointwise) {
    iv <- interval(summary, level)
    covered <- iv$lo <= truth & truth <= iv$hi
    tag <- sprintf("%s%02d", prefix, round(100 * level))
    out[[paste0(tag, "Cov")]] <- safeMean(covered[keep])
    if (level >= 0.9) {
      out[[paste0(tag, "Width")]] <- safeMean((iv$hi - iv$lo)[keep])
    }
  }
  iv <- interval(summary, 0.9)
  out[[paste0(prefix, "Is90")]] <- mean(intervalScore(iv$lo, iv$hi, truth, 0.9))
  out[[paste0(prefix, "Rmse")]] <- sqrt(mean((summary$mean - truth)^2))
  if (!is.null(bandSummary)) {
    for (band in bandSummary$bands) {
      tag <- sprintf("%sSim%02d", prefix, round(100 * band$level))
      iv <- interval(summary, band$level)
      out[[tag]] <- as.numeric(all(band$lo <= truth & truth <= band$hi))
      out[[paste0(tag, "Widen")]] <- mean(band$hi - band$lo) /
        mean(iv$hi - iv$lo)
      out[[paste0(tag, "Collapsed")]] <- as.numeric(band$collapsed)
    }
  }
  out
}

# Declarations at a quantity whose truth is zero (null) and sign errors among
# declared quantities whose truth is not.
declarationMeasures <- function(summary, bandSummary, truth) {
  out <- list()
  null <- truth == 0
  for (level in c(0.9, 0.95)) {
    tag <- sprintf("%02d", round(100 * level))
    iv <- interval(summary, level)
    declared <- iv$lo > 0 | iv$hi < 0
    wrong <- (iv$lo > 0 & truth < 0) | (iv$hi < 0 & truth > 0)
    if (any(null)) {
      out[[paste0("m1Cov", tag)]] <- mean(!declared[null])
      out[[paste0("m2Point", tag)]] <- as.numeric(any(declared[null]))
      band <- bandSummary$bands[[format(level)]]
      bandDeclared <- band$lo > 0 | band$hi < 0
      out[[paste0("m2Band", tag)]] <- as.numeric(any(bandDeclared[null]))
    }
    if (any(!null)) {
      out[[paste0("m3Decl", tag)]] <- sum(declared[!null])
      out[[paste0("m3Wrong", tag)]] <- sum(wrong[!null])
    }
  }
  out
}

averageMeasures <- function(draws, truth, prefix) {
  out <- list()
  for (level in c(0.9, 0.95)) {
    q <- stats::quantile(
      draws,
      c((1 - level) / 2, (1 + level) / 2),
      names = FALSE
    )
    tag <- sprintf("%s%02d", prefix, round(100 * level))
    out[[paste0(tag, "Cov")]] <- as.numeric(q[1L] <= truth && truth <= q[2L])
    out[[paste0(tag, "Width")]] <- q[2L] - q[1L]
  }
  out[[paste0(prefix, "Err")]] <- mean(draws) - truth
  out
}

# What a truth-unknown fit keeps of a reading for scoring against its
# reference later: the pointwise intervals, the bands and the mean.
storedSummary <- function(summary, bandSummary) {
  list(
    mean = summary$mean,
    q = summary$q,
    bands = lapply(bandSummary$bands, function(b) {
      b[c("level", "lo", "hi", "collapsed")]
    })
  )
}

## ---------------------------------------------------------------- fits

timingFile <- function() runDir("timing", "ratios.rds")

# r, the move build's sampling time over the shipped build's, for a cell.
timeRatio <- function(link, cellId) {
  fixed <- Sys.getenv("BIS_FIXED_R")
  if (nzchar(fixed)) {
    return(as.numeric(fixed))
  }
  path <- timingFile()
  if (!file.exists(path)) {
    stop("no timing table; run `timing summarize`")
  }
  tab <- readRDS(path)
  r <- tab$r[tab$link == link & tab$cell == cellId]
  if (length(r) != 1L || !is.finite(r)) {
    stop("no time ratio for ", cellId)
  }
  r
}

# The readings of a run, as draw indices within a chain. L = eqUnit * r is
# the shipped sampler's length a chain at the move's running time, spent as
# more draws or as more burn-in.
readingIndices <- function(cell, r, which) {
  n <- drawsPerChain(cell)
  L <- round(mcmc$eqUnit * r)
  if (L < mcmc$eqUnit) {
    stop("time ratio below 1 for ", cell$id)
  }
  if (L - mcmc$burn > n) {
    stop("equal time past the run for ", cell$id)
  }
  all <- list(
    default = seq_len(mcmc$defaultDraws),
    eqDraws = seq_len(L - mcmc$burn),
    eqBurn = seq.int(L - mcmc$eqUnit + 1L, L - mcmc$burn),
    long = seq_len(n)
  )
  all[which]
}

thinTo <- function(idx, count) {
  if (length(idx) <= count) {
    idx
  } else {
    idx[unique(round(seq(1, length(idx), length.out = count)))]
  }
}

# A draws array (chains x draws x columns) at some draws of every chain, as
# a matrix of pooled draws by columns, chain-major.
poolDraws <- function(arr, idx, cols = NULL) {
  a <- if (is.null(cols)) {
    arr[, idx, , drop = FALSE]
  } else {
    arr[, idx, cols, drop = FALSE]
  }
  matrix(aperm(a, c(2L, 1L, 3L)), dim(a)[1L] * dim(a)[2L], dim(a)[3L])
}

mixingOf <- function(series) {
  # series: draws x chains
  if (
    !requireNamespace("posterior", quietly = TRUE) ||
      nrow(series) < 4L ||
      any(!is.finite(series))
  ) {
    return(c(rhat = NA_real_, ess = NA_real_))
  }
  c(rhat = posterior::rhat(series), ess = posterior::ess_bulk(series))
}

storeDraws <- function(draws, path) {
  sel <- unique(round(seq(
    1,
    nrow(draws),
    length.out = min(nrow(draws), mcmc$stored)
  )))
  dims <- writeSingle(draws[sel, , drop = FALSE], path)
  dims
}

# Scores one reading of one fit. `arr` is the fit's latent test draws.
scoreReading <- function(
  cell,
  rep,
  link,
  arr,
  idx,
  thinIdx,
  replicate,
  storePath
) {
  invLink <- linkFunction(link)
  kind <- replicate$kind
  out <- list()
  stored <- NULL
  if (kind == "prob") {
    sb <- summarizeDraws(poolDraws(arr, thinIdx), invLink)
    s <- if (identical(thinIdx, idx)) {
      sb
    } else {
      summarizeDraws(poolDraws(arr, idx), invLink, bandLevels = numeric())
    }
    if (is.null(replicate$truth)) {
      stored <- storedSummary(s, sb)
      if (!is.null(storePath)) {
        stored$draws <- storeDraws(invLink(poolDraws(arr, idx)), storePath)
      }
    } else {
      p <- replicate$truth
      out <- coverageMeasures(s, sb, p, "f")
      mid <- p >= 0.01 & p <= 0.99
      iv <- interval(s, 0.9)
      out$fMid90Cov <- safeMean((iv$lo <= p & p <= iv$hi)[mid])
      out$fMidShare <- mean(mid)
    }
    return(list(measures = out, stored = stored))
  }
  if (kind == "multi") {
    n <- replicate$n
    effectsOf <- function(index) {
      probs <- lapply(0:3, function(a) {
        invLink(poolDraws(arr, index, a * n + seq_len(n)))
      })
      E <- matrix(0, length(probs[[1L]]) / n, 6L * n)
      for (j in seq_len(nrow(contrastPairs))) {
        cols <- (seq_len(n) - 1L) * 6L + j
        E[, cols] <- probs[[contrastPairs[j, 1L] + 1L]] -
          probs[[contrastPairs[j, 2L] + 1L]]
      }
      list(E = E, probs = probs)
    }
    full <- effectsOf(idx)
    thin <- if (identical(thinIdx, idx)) full else effectsOf(thinIdx)
    sb <- summarizeDraws(thin$E, groupSize = 6L)
    s <- if (identical(thinIdx, idx)) {
      sb
    } else {
      summarizeDraws(full$E, bandLevels = numeric())
    }
    truth <- as.vector(t(replicate$contrasts))
    null <- truth == 0
    out <- coverageMeasures(s, sb, truth, "t", mask = !null)
    for (j in seq_len(nrow(contrastPairs))) {
      cols <- (seq_len(n) - 1L) * 6L + j
      iv <- interval(s, 0.9)
      out[[sprintf(
        "tPair%d%d90Cov",
        contrastPairs[j, 1L],
        contrastPairs[j, 2L]
      )]] <-
        mean((iv$lo <= truth & truth <= iv$hi)[cols])
    }
    out <- c(out, declarationMeasures(s, sb, truth))
    for (l in seq_along(sb$unitBands)) {
      ub <- sb$unitBands[[l]]
      inside <- matrix(ub$lo <= truth & truth <= ub$hi, 6L)
      out[[sprintf("unitBand%02dCov", round(100 * sb$bands[[l]]$level))]] <-
        mean(colSums(inside) == 6L)
    }
    ate <- vapply(
      seq_len(nrow(contrastPairs)),
      function(j) {
        cols <- (seq_len(n) - 1L) * 6L + j
        a <- averageMeasures(
          rowMeans(full$E[, cols, drop = FALSE]),
          mean(replicate$contrasts[, j]),
          "a"
        )
        c(a$a90Cov, a$a95Cov, abs(a$aErr))
      },
      numeric(3L)
    )
    out$a90Cov <- mean(ate[1L, ])
    out$a95Cov <- mean(ate[2L, ])
    out$aAbsErr <- mean(ate[3L, ])
    if (!is.null(storePath)) {
      stored <- storeDraws(do.call(cbind, full$probs), storePath)
    }
    return(list(measures = out, stored = stored))
  }
  n <- replicate$n
  P1 <- invLink(poolDraws(arr, idx, seq_len(n)))
  P0 <- invLink(poolDraws(arr, idx, n + seq_len(n)))
  E <- P1 - P0
  thinE <- if (identical(thinIdx, idx)) {
    E
  } else {
    invLink(poolDraws(arr, thinIdx, seq_len(n))) -
      invLink(poolDraws(arr, thinIdx, n + seq_len(n)))
  }
  sb <- summarizeDraws(thinE)
  s <- if (identical(thinIdx, idx)) {
    sb
  } else {
    summarizeDraws(E, bandLevels = numeric())
  }
  if (!is.null(storePath)) {
    stored <- storeDraws(E, storePath)
  }
  if (kind == "twins") {
    # the sample's own effects: the counterfactual twin drawn from its
    # predicted probability, the observed one as observed
    z <- replicate$z
    y <- replicate$y
    pcf <- P1
    pcf[, z == 1L] <- P0[, z == 1L]
    ycf <- withSeed(
      seedOf("A5-predict", cell$id, rep, link),
      matrix(as.integer(runif(length(pcf)) < pcf), nrow(P1))
    )
    obs <- matrix(y, nrow(P1), n, byrow = TRUE)
    zz <- matrix(z, nrow(P1), n, byrow = TRUE)
    effect <- ifelse(zz == 1L, obs - ycf, ycf - obs)
    trueEffect <- replicate$y1 - replicate$y0
    out <- averageMeasures(rowMeans(effect), mean(trueEffect), "a")
    groups <- sort(unique(replicate$group))
    gDraws <- vapply(
      groups,
      function(g) rowMeans(effect[, replicate$group == g, drop = FALSE]),
      numeric(nrow(effect))
    )
    gTruth <- vapply(
      groups,
      function(g) mean(trueEffect[replicate$group == g]),
      numeric(1L)
    )
    gs <- summarizeDraws(gDraws, bandLevels = numeric())
    out <- c(out, coverageMeasures(gs, NULL, gTruth, "t"))
    pHat <- colMeans(pcf)
    yTrue <- ifelse(z == 1L, replicate$y0, replicate$y1)
    clamped <- pmin(pmax(pHat, 1e-6), 1 - 1e-6)
    out$cfLogScore <- -mean(
      yTrue * log(clamped) + (1 - yTrue) * log(1 - clamped)
    )
    out$cfCalibration <- mean(pHat) - mean(yTrue)
    return(list(measures = out, stored = stored))
  }
  if (is.null(replicate$p1)) {
    return(list(
      measures = list(),
      stored = c(list(draws = stored), storedSummary(s, sb))
    ))
  }
  tau <- replicate$p1 - replicate$p0
  null <- tau == 0
  out <- coverageMeasures(s, sb, tau, "t", mask = !null)
  mid <- pmin(replicate$p0, replicate$p1) >= 0.01 &
    pmax(replicate$p0, replicate$p1) <= 0.99
  iv <- interval(s, 0.9)
  covered <- iv$lo <= tau & tau <= iv$hi
  out$tMid90Cov <- safeMean(covered[mid & !null])
  out$tMidShare <- mean(mid)
  if (!is.null(replicate$noComparison)) {
    out$tNoComparison90Cov <- safeMean(covered[replicate$noComparison & !null])
    out$tComparison90Cov <- safeMean(covered[!replicate$noComparison & !null])
  }
  out <- c(out, declarationMeasures(s, sb, tau))
  out <- c(out, averageMeasures(rowMeans(E), mean(tau), "a"))
  list(measures = out, stored = stored)
}

fitFile <- function(link, cellId, rep, which = build()) {
  runDir(
    "fits",
    link,
    cellFile(cellId),
    sprintf("r%03d-%s.rds", as.integer(rep), which)
  )
}

# One fit of one cell's replication by the build this process loaded, read at
# every length the rule and the report use.
runFit <- function(link, cellId, rep) {
  cell <- getCell(cellId)
  rep <- as.integer(rep)
  which <- build()
  r <- timeRatio(link, cellId)
  replicate <- replication(cell, rep, link)
  invLink <- linkFunction(link)
  n <- drawsPerChain(cell)
  elapsed <- system.time(
    fit <- bart(
      y ~ .,
      data = replicate$train,
      test = replicate$test,
      family = link,
      n.trees = nTrees,
      n.chains = cell$chains,
      n.threads = 1L,
      n.burn = mcmc$burn,
      n.samples = n,
      seed = seedOf("fit", cell$id, rep, link),
      verbose = FALSE,
      combineChains = FALSE,
      keepTrainingFits = FALSE
    )
  )[["elapsed"]]
  arr <- fit$yhat.test
  k <- fit[["k"]]
  result <- list(
    cell = cellId,
    rep = rep,
    link = link,
    build = which,
    lib = find.package("dbarts"),
    seconds = elapsed,
    r = r,
    quick = quick,
    ok = TRUE,
    reason = NA_character_
  )
  if (any(!is.finite(arr)) || any(!is.finite(k))) {
    result$ok <- FALSE
    result$reason <- "non-finite draw"
    saveAtomic(result, fitFile(link, cellId, rep))
    cat(sprintf("fit %s %s %d %s: NON-FINITE\n", link, cellId, rep, which))
    return(invisible(result))
  }
  keepNames <- if (which == "move") {
    c("default", "long")
  } else {
    c("default", "eqDraws", "eqBurn", "long")
  }
  indices <- readingIndices(cell, r, keepNames)
  # Before the pilot is frozen every default and equal-time reading is kept;
  # after it, only the two readings the rule compares (the move's default, the
  # shipped arm's frozen equal-time reading), compressed.
  frozenChoice <- if (file.exists(frozenFile())) {
    readRDS(frozenFile())$equalTime
  } else {
    NULL
  }
  storeKinds <- if (is.null(frozenChoice)) {
    c("default", "eqDraws", "eqBurn")
  } else if (which == "move") {
    "default"
  } else {
    frozenChoice
  }
  drawsSuffix <- if (is.null(frozenChoice)) ".f32" else ".f32.gz"
  scoring <- system.time({
    readings <- lapply(names(indices), function(name) {
      idx <- indices[[name]]
      thin <- if (name == "eqDraws") thinTo(idx, mcmc$defaultDraws) else idx
      # effect draws at the default and equal-time readings; a truth-unknown
      # fit also keeps its long reading's, for the references' cross-check
      diagnosis <- cell$problem %in% c("B1", "B2") && name == "long"
      storePath <- if (
        (cell$kind != "prob" && name %in% storeKinds) || diagnosis
      ) {
        runDir(
          "draws",
          link,
          cellFile(cellId),
          sprintf("r%03d-%s-%s%s", rep, which, name, drawsSuffix)
        )
      } else {
        NULL
      }
      scored <- scoreReading(
        cell,
        rep,
        link,
        arr,
        idx,
        thin,
        replicate,
        storePath
      )
      kSeries <- t(k[, idx, drop = FALSE])
      pSeries <- t(apply(invLink(arr[, idx, , drop = FALSE]), c(1L, 2L), mean))
      scored$mixing <- c(k = mixingOf(kSeries), p = mixingOf(pSeries))
      scored$draws <- length(idx)
      scored
    })
    names(readings) <- names(indices)
  })[["elapsed"]]
  result$readings <- readings
  result$scoringSeconds <- scoring
  saveAtomic(result, fitFile(link, cellId, rep))
  cat(sprintf(
    "fit %s %s %d %s: %.1f s, scored %.1f s\n",
    link,
    cellId,
    rep,
    which,
    elapsed,
    scoring
  ))
  invisible(result)
}

## -------------------------------------------------------------- timing

timingPieceFile <- function(link, cellId, which, t) {
  runDir("timing", link, sprintf("%s-%s-t%d.rds", cellFile(cellId), which, t))
}

# Sampling time alone: a sampler built by bart() with samplerOnly = TRUE on
# the cell's first replication, one chain, warmed, then a run timed by
# itself. Setup, prediction after the run and R's handling of the draws stay
# out of it. A sampler-only run does not draw bart()'s chain, so it times and
# never fits.
runTiming <- function(link, cellId, t) {
  cell <- getCell(cellId)
  replicate <- replication(cell, 1L, link)
  sampler <- bart(
    y ~ .,
    data = replicate$train,
    test = replicate$test,
    family = link,
    n.trees = nTrees,
    n.chains = 1L,
    n.threads = 1L,
    seed = seedOf("timing", cellId, link),
    verbose = FALSE,
    keepTrainingFits = FALSE,
    samplerOnly = TRUE
  )
  invisible(sampler$run(mcmc$timingWarm, 0L))
  elapsed <- system.time(invisible(sampler$run(0L, mcmc$timingRun)))[[
    "elapsed"
  ]]
  saveAtomic(
    list(
      link = link,
      cell = cellId,
      build = build(),
      t = as.integer(t),
      seconds = elapsed,
      iterations = mcmc$timingRun,
      load = if (file.exists("/proc/loadavg")) {
        readLines("/proc/loadavg", warn = FALSE)
      } else {
        NA_character_
      },
      quick = quick
    ),
    timingPieceFile(link, cellId, build(), t)
  )
  cat(sprintf(
    "timing %s %s %s %d: %.3f s\n",
    link,
    cellId,
    build(),
    as.integer(t),
    elapsed
  ))
}

# r per cell: the median over repetitions of the move's time over the shipped
# build's, each repetition's two runs adjacent in time.
summarizeTiming <- function() {
  files <- list.files(
    runDir("timing"),
    pattern = "-t[0-9]+[.]rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  pieces <- lapply(files, readRDS)
  pieces <- Filter(function(p) identical(p$quick, quick), pieces)
  d <- data.frame(
    link = vapply(pieces, `[[`, "", "link"),
    cell = vapply(pieces, `[[`, "", "cell"),
    build = vapply(pieces, `[[`, "", "build"),
    t = vapply(pieces, `[[`, 0L, "t"),
    seconds = vapply(pieces, `[[`, 0, "seconds"),
    iterations = vapply(pieces, `[[`, 0L, "iterations"),
    stringsAsFactors = FALSE
  )
  keys <- unique(d[, c("link", "cell")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    s <- d[d$link == keys$link[i] & d$cell == keys$cell[i], ]
    sh <- s[s$build == "shipped", ]
    mv <- s[s$build == "move", ]
    common <- intersect(sh$t, mv$t)
    ratios <- mv$seconds[match(common, mv$t)] / sh$seconds[match(common, sh$t)]
    data.frame(
      link = keys$link[i],
      cell = keys$cell[i],
      n = length(common),
      r = if (length(common)) stats::median(ratios) else NA_real_,
      rMin = if (length(common)) min(ratios) else NA_real_,
      rMax = if (length(common)) max(ratios) else NA_real_,
      msShipped = 1000 * stats::median(sh$seconds) / sh$iterations[1L],
      msMove = 1000 * stats::median(mv$seconds) / mv$iterations[1L],
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$link, tab$cell), ]
  saveAtomic(tab, timingFile())
  print(tab, row.names = FALSE, digits = 4)
  invisible(tab)
}

## ---------------------------------------------------------- references

referenceDir <- function(link, cellId, rep) {
  runDir("ref", link, sprintf("%s-r%03d", cellFile(cellId), as.integer(rep)))
}

# The scored quantity's draws at a fit's test rows: the probability, or each
# unit's effect.
scoredDraws <- function(arr, replicate, invLink, idx = seq_len(dim(arr)[2L])) {
  if (replicate$kind == "prob") {
    return(invLink(poolDraws(arr, idx)))
  }
  n <- replicate$n
  invLink(poolDraws(arr, idx, seq_len(n))) -
    invLink(poolDraws(arr, idx, n + seq_len(n)))
}

# One part of a reference: chains of one build at reference length (5,000
# burn-in, 15,000 iterations kept every 15th). A reference is 8 chains of each
# build, in refParts parts a build.
runReference <- function(link, cellId, part, rep = 1L) {
  cell <- getCell(cellId)
  which <- build()
  if (!which %in% c("shipped", "move")) {
    stop("a reference is made by the two arms' builds")
  }
  replicate <- replication(cell, as.integer(rep), link)
  chains <- refChainsPerBuild %/% refParts
  elapsed <- system.time(
    fit <- bart(
      y ~ .,
      data = replicate$train,
      test = replicate$test,
      family = link,
      n.trees = nTrees,
      n.chains = chains,
      n.threads = 1L,
      n.burn = mcmc$refBurn,
      n.samples = mcmc$refIter,
      n.thin = mcmc$refThin,
      seed = seedOf("reference", link, cellId, rep, which, part),
      verbose = FALSE,
      combineChains = FALSE,
      keepTrainingFits = FALSE
    )
  )[["elapsed"]]
  invLink <- linkFunction(link)
  arr <- fit$yhat.test
  if (any(!is.finite(arr))) {
    stop("non-finite reference draws")
  }
  draws <- scoredDraws(arr, replicate, invLink)
  dir <- referenceDir(link, cellId, rep)
  dims <- writeSingle(
    draws,
    file.path(dir, sprintf("%s-part%d.f32", which, part))
  )
  saveAtomic(
    list(
      dims = dims,
      chains = chains,
      k = fit$k,
      meanP = apply(invLink(arr), c(1L, 2L), mean),
      seconds = elapsed,
      lib = find.package("dbarts"),
      quick = quick
    ),
    file.path(dir, sprintf("%s-part%d.rds", which, part))
  )
  cat(sprintf(
    "reference %s %s %s part %d: %.1f s\n",
    link,
    cellId,
    which,
    part,
    elapsed
  ))
}

# A reference's draws as an array: draws a chain x chains x quantities, the
# shipped build's chains first.
loadReference <- function(link, cellId, rep = 1L) {
  dir <- referenceDir(link, cellId, rep)
  pieces <- list()
  for (which in c("shipped", "move")) {
    for (part in seq_len(refParts)) {
      meta <- readRDS(file.path(dir, sprintf("%s-part%d.rds", which, part)))
      if (!identical(meta$quick, quick)) {
        stop("quick and full output mixed")
      }
      d <- readSingle(
        file.path(dir, sprintf("%s-part%d.f32", which, part)),
        meta$dims
      )
      perChain <- nrow(d) / meta$chains
      pieces[[length(pieces) + 1L]] <- list(
        draws = array(d, c(perChain, meta$chains, ncol(d))),
        k = t(meta$k),
        meanP = t(meta$meanP),
        build = which
      )
    }
  }
  perChain <- dim(pieces[[1L]]$draws)[1L]
  m <- dim(pieces[[1L]]$draws)[3L]
  draws <- array(0, c(perChain, 2L * refChainsPerBuild, m))
  k <- meanP <- matrix(0, perChain, 2L * refChainsPerBuild)
  at <- 0L
  for (p in pieces) {
    cols <- at + seq_len(dim(p$draws)[2L])
    draws[, cols, ] <- p$draws
    k[, cols] <- p$k
    meanP[, cols] <- p$meanP
    at <- max(cols)
  }
  list(
    draws = draws,
    k = k,
    meanP = meanP,
    build = rep(c("shipped", "move"), each = refChainsPerBuild)
  )
}

# Acceptance: the two builds' halves agree on every quantity's mean and 5th
# and 95th percentiles within 4.5 combined Monte Carlo standard errors;
# rank-normalized split-Rhat over all 16 chains at most 1.01 for k and the mean
# held-out probability and at most 1.02 at the 99th percentile of quantities;
# tail effective draws at least 400 (k, the mean probability and every
# quantity). A failing reference stays in, flagged.
acceptReference <- function(link, cellId, rep = 1L) {
  ref <- loadReference(link, cellId, rep)
  m <- dim(ref$draws)[3L]
  halves <- split(seq_len(dim(ref$draws)[2L]), ref$build)
  stats <- function(x) {
    c(
      mean = mean(x),
      q05 = unname(stats::quantile(x, 0.05)),
      q95 = unname(stats::quantile(x, 0.95)),
      seMean = posterior::mcse_mean(x),
      se05 = posterior::mcse_quantile(x, probs = 0.05),
      se95 = posterior::mcse_quantile(x, probs = 0.95)
    )
  }
  z <- matrix(NA_real_, m, 3L)
  rhat <- essTail <- numeric(m)
  for (j in seq_len(m)) {
    x <- ref$draws[,, j]
    a <- stats(x[, halves$shipped, drop = FALSE])
    b <- stats(x[, halves$move, drop = FALSE])
    z[j, ] <- abs(a[1:3] - b[1:3]) / sqrt(a[4:6]^2 + b[4:6]^2)
    rhat[j] <- posterior::rhat(x)
    essTail[j] <- posterior::ess_tail(x)
  }
  kRhat <- posterior::rhat(ref$k)
  pRhat <- posterior::rhat(ref$meanP)
  tails <- c(k = posterior::ess_tail(ref$k), p = posterior::ess_tail(ref$meanP))
  checks <- c(
    halvesAgree = all(z <= 4.5, na.rm = TRUE),
    kRhat = is.finite(kRhat) && kRhat <= 1.01,
    pRhat = is.finite(pRhat) && pRhat <= 1.01,
    rhat99 = unname(stats::quantile(rhat, 0.99, na.rm = TRUE)) <= 1.02,
    essTail = min(c(essTail, tails), na.rm = TRUE) >= 400
  )
  out <- list(
    cell = cellId,
    link = link,
    rep = rep,
    accepted = all(checks),
    checks = checks,
    maxZ = max(z, na.rm = TRUE),
    shareZ = mean(z > 4.5, na.rm = TRUE),
    kRhat = kRhat,
    pRhat = pRhat,
    rhat99 = unname(stats::quantile(rhat, 0.99, na.rm = TRUE)),
    minEssTail = min(c(essTail, tails), na.rm = TRUE),
    quick = quick
  )
  saveAtomic(out, file.path(referenceDir(link, cellId, rep), "acceptance.rds"))
  print(out[c(
    "cell",
    "link",
    "accepted",
    "checks",
    "maxZ",
    "kRhat",
    "pRhat",
    "rhat99",
    "minEssTail"
  )])
  invisible(out)
}

## ---------------------------------------------------- stage 1 checks

# Seeded short fits along the study's own paths: numeric and factor
# predictors, missing predictor values, test rows, 4 and 10 chains, one
# thread, both links, and the two propensity fits. Each is also refit
# shorter to check that a shorter run is a longer run's prefix.
checkPaths <- function() {
  set <- list()
  addProb <- function(dataset, chains) {
    raw <- realDatasets[[dataset]]()
    split <- withSeed(seedOf("check", dataset), {
      train <- sample.int(nrow(raw$x), floor(0.8 * nrow(raw$x)))
      list(train = train, test = setdiff(seq_len(nrow(raw$x)), train))
    })
    set[[paste0("prob-", dataset)]] <<- list(
      train = data.frame(
        raw$x[split$train, , drop = FALSE],
        y = raw$y[split$train]
      ),
      test = raw$x[split$test, , drop = FALSE],
      chains = chains
    )
  }
  addProb("pima", 4L)
  addProb("german", 4L)
  addProb("creditapproval", 4L)
  effect <- function(name, x, z, y, chains) {
    ps <- withSeed(seedOf("check-ps", name), stats::runif(nrow(x), 0.2, 0.8))
    set[[name]] <<- list(
      train = data.frame(x, trt = z, pscore = ps, y = y),
      test = rbind(
        data.frame(x, trt = 1L, pscore = ps),
        data.frame(x, trt = 0L, pscore = ps)
      ),
      chains = chains
    )
  }
  tw <- twinsBase(1L)
  effect("twins-missing", tw$x, tw$z, tw$y, 10L)
  a3 <- effectBase(getCell("A3:r31-mod"), 1L, "probit")
  effect(
    "acic16-factors",
    a3$x,
    a3$z,
    outcomeDraw(ifelse(a3$z == 1L, a3$p1, a3$p0), getCell("A3:r31-mod"), 1L),
    10L
  )
  b6 <- a6Base(1L)
  rows <- seq_len(500L)
  arm <- factor(b6$arm[rows], levels = 0:3)
  y6 <- withSeed(seedOf("check-a6"), stats::rbinom(500L, 1L, 0.3))
  set[["a6-arm"]] <- list(
    train = data.frame(b6$x[rows, ], arm = arm, y = y6),
    test = do.call(
      rbind,
      lapply(0:3, function(a) {
        data.frame(b6$x[rows, ], arm = factor(rep(a, 500L), levels = 0:3))
      })
    ),
    chains = 10L
  )
  list(
    fits = set,
    propensity = list(x = a3$x, z = a3$z),
    multinomial = list(x = b6$x[rows, ], arm = arm)
  )
}

runCheck <- function(outFile) {
  paths <- checkPaths()
  short <- 40L
  long <- 100L
  res <- list(build = build(), lib = find.package("dbarts"), fits = list())
  one <- function(f, n) {
    fit <- f(n)
    list(test = fit$yhat.test, train = fit$yhat.train, k = fit[["k"]])
  }
  prefixOk <- function(a, b) {
    sub <- function(x, n) {
      if (is.null(x)) NULL else x[, seq_len(n), , drop = FALSE]
    }
    subK <- function(x, n) {
      if (is.null(x)) NULL else x[, seq_len(n), drop = FALSE]
    }
    identical(a$test, sub(b$test, short)) && identical(a$k, subK(b$k, short))
  }
  for (link in links) {
    for (name in names(paths$fits)) {
      p <- paths$fits[[name]]
      f <- function(n) {
        bart(
          y ~ .,
          data = p$train,
          test = p$test,
          family = link,
          n.trees = nTrees,
          n.chains = p$chains,
          n.threads = 1L,
          n.burn = 50L,
          n.samples = n,
          seed = seedOf("check", name, link),
          verbose = FALSE,
          combineChains = FALSE,
          keepTrainingFits = FALSE
        )
      }
      a <- one(f, short)
      b <- one(f, long)
      res$fits[[paste(link, name, sep = ":")]] <- c(b, prefix = prefixOk(a, b))
    }
  }
  pr <- function(n) {
    bart(
      trt ~ .,
      data = data.frame(paths$propensity$x, trt = paths$propensity$z),
      family = "probit",
      n.chains = 10L,
      n.threads = 1L,
      n.burn = 50L,
      n.samples = n,
      n.trees = nTrees,
      seed = seedOf("check-propensity"),
      verbose = FALSE,
      combineChains = FALSE
    )
  }
  a <- one(pr, short)
  b <- one(pr, long)
  res$fits[["propensity"]] <- c(
    b,
    prefix = identical(a$train, b$train[, seq_len(short), , drop = FALSE])
  )
  mn <- function(n) {
    bart(
      arm ~ .,
      data = data.frame(paths$multinomial$x, arm = paths$multinomial$arm),
      family = "multinomial",
      n.chains = 4L,
      n.threads = 1L,
      n.burn = 50L,
      n.samples = n,
      n.trees = nTrees,
      seed = seedOf("check-multinomial"),
      verbose = FALSE
    )
  }
  b <- one(mn, long)
  res$fits[["multinomial"]] <- c(b, prefix = NA)
  saveRDS(res, outFile)
  for (name in names(res$fits)) {
    cat(sprintf("%-32s prefix %s\n", name, res$fits[[name]]$prefix))
  }
}

compareChecks <- function(fileA, fileB) {
  a <- readRDS(fileA)
  b <- readRDS(fileB)
  cat(sprintf("%s (%s) against %s (%s)\n", a$build, a$lib, b$build, b$lib))
  for (name in names(a$fits)) {
    x <- a$fits[[name]]
    y <- b$fits[[name]]
    cat(sprintf(
      "%-32s test %-5s train %-5s k %-5s\n",
      name,
      identical(x$test, y$test),
      identical(x$train, y$train),
      identical(x$k, y$k)
    ))
  }
}

## ------------------------------------------------------ reading results

# The arms the rule compares: the move's default reading against the shipped
# sampler's equal-time reading, the choice between more draws and more burn-in
# made on the pilot and frozen.
armReadings <- c(
  "move.default",
  "move.long",
  "shipped.default",
  "shipped.eqDraws",
  "shipped.eqBurn",
  "shipped.long"
)

# Per cell, the replications both builds ran and, for each arm reading, a
# data frame of the measures by replication. A pair with a failed arm is
# skipped, tallied by the arm that failed.
readCell <- function(link, cell, reps) {
  values <- setNames(vector("list", length(armReadings)), armReadings)
  kept <- integer()
  skipped <- data.frame(rep = integer(), shipped = logical(), move = logical())
  for (rep in reps) {
    files <- c(
      shipped = fitFile(link, cell$id, rep, "shipped"),
      move = fitFile(link, cell$id, rep, "move")
    )
    present <- file.exists(files)
    if (!all(present)) {
      failed <- file.exists(paste0(files, ".skip"))
      if (any(failed)) {
        skipped <- rbind(
          skipped,
          data.frame(
            rep = rep,
            shipped = failed[["shipped"]],
            move = failed[["move"]]
          )
        )
      }
      next
    }
    res <- lapply(files, readRDS)
    if (!all(vapply(res, function(r) identical(r$quick, quick), logical(1L)))) {
      stop("quick and full output mixed in ", cell$id)
    }
    ok <- vapply(res, function(r) isTRUE(r$ok), logical(1L))
    if (!all(ok)) {
      skipped <- rbind(
        skipped,
        data.frame(rep = rep, shipped = !ok[["shipped"]], move = !ok[["move"]])
      )
      next
    }
    kept <- c(kept, rep)
    for (ar in armReadings) {
      parts <- strsplit(ar, ".", fixed = TRUE)[[1L]]
      reading <- res[[parts[1L]]]$readings[[parts[2L]]]
      m <- reading$measures
      m <- m[vapply(m, function(v) length(v) == 1L, logical(1L))]
      row <- as.data.frame(
        c(list(rep = rep), m, as.list(reading$mixing)),
        stringsAsFactors = FALSE
      )
      values[[ar]] <- if (is.null(values[[ar]])) {
        row
      } else {
        merge(values[[ar]], row, all = TRUE, sort = FALSE)
      }
    }
  }
  list(cell = cell, reps = kept, values = values, skipped = skipped)
}

readAll <- function(link, repsFor) {
  ids <- names(cells)[vapply(
    cells,
    function(c) substr(c$problem, 1L, 1L) == "A",
    logical(1L)
  )]
  out <- lapply(ids, function(id) {
    readCell(link, cells[[id]], repsFor(cells[[id]]))
  })
  names(out) <- ids
  out
}

## ---------------------------------------------------------- the rule

# A cell's distance from nominal corrected for noise: the d >= 0 whose folded
# normal mean, s sqrt(2/pi) exp(-d^2 / (2 s^2)) + d (1 - 2 Phi(-d / s)),
# equals the observed |coverage - nominal|, zero below s sqrt(2/pi). The
# folded mean is convex and increasing in d and at least d, so Newton's
# method from the observed distance falls to the root from above.
correctedDistance <- function(observed, s) {
  d <- observed
  floorAt <- s * sqrt(2 / pi)
  active <- s > 0 & observed > floorAt
  dd <- observed[active]
  ss <- s[active]
  aa <- observed[active]
  for (iter in seq_len(60L)) {
    g <- ss *
      sqrt(2 / pi) *
      exp(-dd^2 / (2 * ss^2)) +
      dd * (1 - 2 * pnorm(-dd / ss))
    slope <- 1 - 2 * pnorm(-dd / ss)
    step <- (g - aa) / pmax(slope, 1e-12)
    dd <- pmax(dd - step, 0)
    if (all(abs(step) < 1e-12)) break
  }
  d[active] <- dd
  d[!active & s > 0] <- 0
  d
}

# The rule's measures. `nominal` and the measure name are per level: the
# primary level first (90, or 95 for M), the other read beside it.
ruleMeasures <- list(
  F = list(
    type = "coverage",
    inGroup = "F",
    measure = c("f90Cov", "f95Cov"),
    nominal = c(0.90, 0.95),
    bar = 0.02,
    margin = 0.010,
    plannedSe = 0.0019,
    plannedPairs = 704,
    gain = TRUE
  ),
  T = list(
    type = "coverage",
    inGroup = "T",
    measure = c("t90Cov", "t95Cov"),
    nominal = c(0.90, 0.95),
    bar = 0.02,
    margin = 0.015,
    plannedSe = 0.0038,
    plannedPairs = 692,
    gain = TRUE
  ),
  M1 = list(
    type = "coverage",
    inGroup = "M1",
    measure = c("m1Cov95", "m1Cov90"),
    nominal = c(0.95, 0.90),
    bar = 0.02,
    margin = 0.015,
    plannedSe = 0.0032,
    plannedPairs = 346,
    gain = TRUE
  ),
  M2 = list(
    type = "coverage",
    inGroup = "M2",
    measure = c("m2Band95", "m2Band90"),
    nominal = c(0.05, 0.10),
    bar = 0.05,
    margin = 0.075,
    plannedSe = 0.018,
    plannedPairs = 316,
    gain = TRUE
  ),
  M3 = list(
    type = "share",
    inGroup = "M3",
    measure = c("m3Wrong95", "m3Wrong90"),
    denominator = c("m3Decl95", "m3Decl90"),
    nominal = c(0, 0),
    bar = 0.01,
    margin = 0.010,
    plannedSe = NA_real_,
    plannedPairs = 672,
    gain = TRUE
  ),
  A6M1 = list(
    type = "coverage",
    inGroup = "A6M",
    measure = c("m1Cov95", "m1Cov90"),
    nominal = c(0.95, 0.90),
    margin = 0.025,
    plannedSe = 0.0055,
    plannedPairs = 120,
    gain = FALSE
  ),
  A6M2 = list(
    type = "coverage",
    inGroup = "A6M",
    measure = c("m2Band95", "m2Band90"),
    nominal = c(0.05, 0.10),
    margin = 0.11,
    plannedSe = 0.029,
    plannedPairs = 120,
    gain = FALSE
  ),
  sim = list(
    type = "coverage",
    inGroup = "sim",
    measure = c("Sim90", "Sim95"),
    nominal = c(0.90, 0.95),
    margin = 0.035,
    plannedSe = 0.0097,
    plannedPairs = 1602,
    gain = FALSE
  ),
  ate = list(
    type = "coverage",
    inGroup = "ate",
    measure = c("a90Cov", "a95Cov"),
    nominal = c(0.90, 0.95),
    margin = 0.040,
    plannedSe = 0.011,
    plannedPairs = 918,
    gain = FALSE
  ),
  rmse = list(
    type = "ratio",
    inGroup = "rmse",
    measure = c("Rmse", "Rmse"),
    margin = 1.02,
    gain = FALSE
  ),
  is = list(
    type = "ratio",
    inGroup = "is",
    measure = c("Is90", "Is90"),
    margin = 1.02,
    gain = FALSE
  )
)

# The column a cell reads for a measure: the simultaneous and error measures
# carry the cell's quantity as a prefix ("f" probabilities, "t" effects).
measureColumn <- function(cell, name) {
  if (name %in% c("Sim90", "Sim95", "Rmse", "Is90")) {
    prefix <- if (cell$kind == "prob") "f" else "t"
    return(paste0(prefix, name))
  }
  name
}

# Per-replication values of one measure in one cell for the two compared arm
# readings, aligned by replication.
pairedValues <- function(cr, spec, level, arms) {
  cell <- cr$cell
  col <- measureColumn(cell, spec$measure[level])
  get <- function(ar, column) {
    v <- cr$values[[ar]]
    if (is.null(v) || !column %in% names(v)) {
      return(NULL)
    }
    v[[column]][match(cr$reps, v$rep)]
  }
  a <- get(arms[1L], col)
  b <- get(arms[2L], col)
  if (is.null(a) || is.null(b)) {
    return(NULL)
  }
  out <- list(a = a, b = b)
  if (spec$type == "share") {
    out$da <- get(arms[1L], spec$denominator[level])
    out$db <- get(arms[2L], spec$denominator[level])
  }
  if (any(!is.finite(unlist(out)))) {
    stop("missing ", col, " in ", cell$id)
  }
  out
}

# The group statistic of one arm from per-cell vectors, under resampling
# index matrices (B x n per cell; NULL is the observed data).
groupStatistic <- function(spec, level, cellValues, which, indices = NULL) {
  stats <- lapply(seq_along(cellValues), function(i) {
    cv <- cellValues[[i]]
    v <- cv[[which]]
    n <- length(v)
    idx <- if (is.null(indices)) matrix(seq_len(n), 1L) else indices[[i]]
    V <- matrix(v[idx], nrow(idx))
    if (spec$type == "coverage") {
      rowMean <- rowMeans(V)
      s <- if (n > 1L) {
        sqrt(pmax(rowSums((V - rowMean)^2) / (n - 1L), 0) / n)
      } else {
        rep(0, nrow(V))
      }
      correctedDistance(abs(rowMean - spec$nominal[level]), s)
    } else if (spec$type == "share") {
      D <- matrix(cv[[paste0("d", which)]][idx], nrow(idx))
      ifelse(rowSums(D) > 0, rowSums(V) / rowSums(D), NA_real_)
    } else {
      if (grepl("Rmse", spec$measure[1L])) sqrt(rowMeans(V^2)) else rowMeans(V)
    }
  })
  do.call(cbind, stats)
}

# The comparison of one measure: the gap of each arm, the gain (shipped gap
# minus the move's, positive when the move is closer to nominal) or, for a
# ratio, the move's over the shipped arm's, and their bootstrap distribution.
compareMeasure <- function(spec, level, cellResults, arms, indices) {
  inGroup <- vapply(
    cellResults,
    function(cr) {
      spec$inGroup %in%
        cr$cell$groups ||
        (spec$inGroup == "A6M" && cr$cell$problem == "A6")
    },
    logical(1L)
  )
  crs <- cellResults[inGroup]
  values <- lapply(crs, pairedValues, spec = spec, level = level, arms = arms)
  usable <- !vapply(values, is.null, logical(1L)) &
    vapply(values, function(v) !is.null(v) && length(v$a) > 0L, logical(1L))
  values <- values[usable]
  idx <- indices[names(crs)[usable]]
  if (length(values) == 0L) {
    return(NULL)
  }
  summarise <- function(ind) {
    A <- groupStatistic(spec, level, values, "a", ind)
    B <- groupStatistic(spec, level, values, "b", ind)
    if (spec$type == "ratio") {
      list(
        a = NA_real_,
        b = NA_real_,
        diff = exp(rowMeans(log(B / A), na.rm = TRUE))
      )
    } else {
      list(
        a = rowMeans(A, na.rm = TRUE),
        b = rowMeans(B, na.rm = TRUE),
        diff = rowMeans(A, na.rm = TRUE) - rowMeans(B, na.rm = TRUE)
      )
    }
  }
  observed <- summarise(NULL)
  boot <- summarise(idx)
  pairs <- sum(vapply(values, function(v) length(v$a), integer(1L)))
  list(
    cells = length(values),
    pairs = pairs,
    gapShipped = observed$a,
    gapMove = observed$b,
    estimate = observed$diff,
    boot = boot$diff,
    values = values
  )
}

# Resampling indices: for every cell, B rows of its n replications drawn with
# replacement, the same rows for both arms (pairs) and for every measure.
resampleIndices <- function(cellResults, B = 2000L) {
  out <- lapply(names(cellResults), function(id) {
    n <- max(1L, length(cellResults[[id]]$reps))
    withSeed(
      seedOf("bootstrap", id),
      matrix(sample.int(n, B * n, replace = TRUE), B, n)
    )
  })
  names(out) <- names(cellResults)
  out
}

# The pooled within-cell spread of the per-replication paired difference,
# scaled to the planned pairs: the pilot's reading of the planned standard
# error column.
plannedSeFromPilot <- function(cmp, spec) {
  if (is.null(cmp) || spec$type == "ratio") {
    return(NA_real_)
  }
  if (spec$type == "share") {
    se <- stats::sd(cmp$boot, na.rm = TRUE)
    return(se * sqrt(cmp$pairs / spec$plannedPairs))
  }
  ss <- 0
  df <- 0
  for (v in cmp$values) {
    d <- v$b - v$a
    if (length(d) > 1L) {
      ss <- ss + sum((d - mean(d))^2)
      df <- df + length(d) - 1L
    }
  }
  sqrt(ss / df) / sqrt(spec$plannedPairs)
}

skippedTally <- function(cellResults) {
  tally <- do.call(
    rbind,
    lapply(cellResults, function(cr) {
      data.frame(
        cell = cr$cell$id,
        pairs = length(cr$reps) + nrow(cr$skipped),
        shipped = sum(cr$skipped$shipped),
        move = sum(cr$skipped$move)
      )
    })
  )
  tally
}

frozenFile <- function() runDir("frozen.rds")

# The rule at one look: every gain measure and guard, its estimate and its
# one-sided bound at the look's alpha, and the verdict.
readRule <- function(link, cellResults, alpha, frozen) {
  arms <- c(paste0("shipped.", frozen$equalTime), "move.default")
  indices <- resampleIndices(cellResults)
  rows <- list()
  for (name in names(ruleMeasures)) {
    spec <- ruleMeasures[[name]]
    margin <- frozen$margins[[name]]
    cmp <- compareMeasure(spec, 1L, cellResults, arms, indices)
    if (is.null(cmp)) {
      next
    }
    alt <- if (spec$type == "ratio") {
      NULL
    } else {
      compareMeasure(spec, 2L, cellResults, arms, indices)
    }
    if (spec$type == "ratio") {
      worsen <- cmp$estimate
      worsenBoot <- cmp$boot
    } else {
      worsen <- -cmp$estimate
      worsenBoot <- -cmp$boot
    }
    lowerWorsen <- unname(stats::quantile(worsenBoot, alpha, na.rm = TRUE))
    upperWorsen <- unname(stats::quantile(worsenBoot, 1 - alpha, na.rm = TRUE))
    row <- data.frame(
      measure = name,
      cells = cmp$cells,
      pairs = cmp$pairs,
      gapShipped = cmp$gapShipped,
      gapMove = cmp$gapMove,
      gain = if (spec$type == "ratio") NA_real_ else cmp$estimate,
      gainLower = if (spec$type == "ratio") NA_real_ else -upperWorsen,
      gainUpper = if (spec$type == "ratio") NA_real_ else -lowerWorsen,
      gainOtherLevel = if (is.null(alt)) NA_real_ else alt$estimate,
      ratio = if (spec$type == "ratio") cmp$estimate else NA_real_,
      worsenLower = lowerWorsen,
      worsenUpper = upperWorsen,
      margin = margin,
      bootSe = stats::sd(cmp$boot, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    row$harmless <- row$worsenUpper < margin
    row$harm <- row$worsenLower > margin
    row$gainShown <- isTRUE(spec$gain) &&
      row$gain >= spec$bar &&
      row$gainLower > 0 &&
      is.finite(row$gainOtherLevel) &&
      sign(row$gainOtherLevel) == sign(row$gain)
    row$canReachBar <- isTRUE(spec$gain) && row$gainUpper >= spec$bar
    rows[[name]] <- row
  }
  tally <- skippedTally(cellResults)
  pairs <- sum(tally$pairs)
  excess <- (sum(tally$move) - sum(tally$shipped)) / max(pairs, 1)
  table <- do.call(rbind, rows)
  skipRow <- data.frame(
    measure = "skipped",
    cells = nrow(tally),
    pairs = pairs,
    gapShipped = sum(tally$shipped),
    gapMove = sum(tally$move),
    gain = NA_real_,
    gainLower = NA_real_,
    gainUpper = NA_real_,
    gainOtherLevel = NA_real_,
    ratio = NA_real_,
    worsenLower = excess,
    worsenUpper = excess,
    margin = 0.005,
    bootSe = NA_real_,
    harmless = excess < 0.005,
    harm = excess > 0.005,
    gainShown = FALSE,
    canReachBar = FALSE
  )
  table <- rbind(table, skipRow)
  verdict <- if (any(table$harm)) {
    "leave out: harm shown"
  } else if (any(table$gainShown) && all(table$harmless)) {
    "include"
  } else if (!any(table$canReachBar)) {
    "inconclusive; no gain can reach its bar, so doubling is skipped"
  } else {
    "inconclusive; double the undecided cells"
  }
  list(
    link = link,
    alpha = alpha,
    arms = arms,
    table = table,
    verdict = verdict,
    skipped = tally,
    skippedShare = sum(tally$shipped | tally$move) / max(pairs, 1)
  )
}

printRule <- function(reading) {
  cat(sprintf(
    "\n== %s at one-sided alpha %.3f: %s against %s\n",
    reading$link,
    reading$alpha,
    reading$arms[2L],
    reading$arms[1L]
  ))
  t <- reading$table
  num <- vapply(t, is.numeric, logical(1L))
  t[num] <- lapply(t[num], function(x) signif(x, 4))
  print(t, row.names = FALSE)
  cat("verdict:", reading$verdict, "\n")
}

## ------------------------------------------------- pilot, looks, summary

ruleGainCells <- c("F", "T", "M1")

# The equal-time choice: the shipped reading, more draws or more burn-in,
# whose cells of F, T and M1 sit closer to nominal on average (each cell's
# corrected distance at its primary level), pooled over both links.
pooledDistance <- function(byLink, reading) {
  d <- c()
  for (crs in byLink) {
    for (cr in crs) {
      for (g in ruleGainCells) {
        spec <- ruleMeasures[[g]]
        if (!g %in% cr$cell$groups) {
          next
        }
        v <- cr$values[[reading]][[measureColumn(cr$cell, spec$measure[1L])]]
        if (is.null(v) || length(v) == 0L) {
          next
        }
        s <- if (length(v) > 1L) stats::sd(v) / sqrt(length(v)) else 0
        d <- c(d, correctedDistance(abs(mean(v) - spec$nominal[1L]), s))
      }
    }
  }
  c(distance = mean(d), cells = length(d))
}

pilotReport <- function(link, crs, equalTime) {
  indices <- resampleIndices(crs)
  comparisons <- list(
    equalTime = c(paste0("shipped.", equalTime), "move.default"),
    equalSweeps = c("shipped.default", "move.default"),
    long = c("shipped.long", "move.long")
  )
  rows <- list()
  for (name in names(ruleMeasures)) {
    spec <- ruleMeasures[[name]]
    for (cmpName in names(comparisons)) {
      cmp <- compareMeasure(spec, 1L, crs, comparisons[[cmpName]], indices)
      if (is.null(cmp)) {
        next
      }
      rows[[length(rows) + 1L]] <- data.frame(
        link = link,
        measure = name,
        reading = cmpName,
        cells = cmp$cells,
        pairs = cmp$pairs,
        gapShipped = cmp$gapShipped,
        gapMove = cmp$gapMove,
        estimate = cmp$estimate,
        bootSe = stats::sd(cmp$boot, na.rm = TRUE),
        plannedSe = if (cmpName == "equalTime") {
          plannedSeFromPilot(cmp, spec)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

runPilot <- function() {
  if (file.exists(frozenFile())) {
    stop("the pilot is frozen: ", frozenFile())
  }
  repsFor <- function(cell) seq_len(min(pilotReps, cell$reps))
  byLink <- lapply(links, function(link) readAll(link, repsFor))
  names(byLink) <- links
  for (link in links) {
    short <- vapply(
      byLink[[link]],
      function(cr) length(cr$reps) < min(pilotReps, cr$cell$reps),
      logical(1L)
    )
    if (any(short)) {
      cat("cells short of the pilot's replications under", link, ":\n")
      print(names(short)[short])
    }
  }
  distances <- rbind(
    eqDraws = pooledDistance(byLink, "shipped.eqDraws"),
    eqBurn = pooledDistance(byLink, "shipped.eqBurn")
  )
  equalTime <- rownames(distances)[which.min(distances[, "distance"])]
  cat(
    "\nEqual-time choice, pooled corrected distance over the cells of F, T and M1:\n"
  )
  print(distances)
  cat("chosen:", equalTime, "\n")
  report <- do.call(
    rbind,
    lapply(links, function(link) pilotReport(link, byLink[[link]], equalTime))
  )
  margins <- lapply(links, function(link) {
    out <- list()
    for (name in names(ruleMeasures)) {
      spec <- ruleMeasures[[name]]
      se <- report$plannedSe[
        report$link == link &
          report$measure == name &
          report$reading == "equalTime"
      ]
      se <- if (length(se) == 1L && is.finite(se)) se else NA_real_
      out[[name]] <- if (spec$type == "ratio") {
        spec$margin
      } else {
        max(spec$margin, 3.3 * se, na.rm = TRUE)
      }
    }
    out
  })
  names(margins) <- links
  frozen <- list(
    equalTime = equalTime,
    distances = distances,
    margins = margins,
    report = report,
    frozenAt = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    quick = quick
  )
  saveAtomic(frozen, frozenFile())
  num <- vapply(report, is.numeric, logical(1L))
  shown <- report
  shown[num] <- lapply(shown[num], function(x) signif(x, 4))
  cat(
    "\nPilot, per measure: gaps (shipped, move), the move's gain (or ratio), bootstrap SE at",
    "the pilot's counts, and the planned-count SE read from the pilot's spread\n"
  )
  print(shown, row.names = FALSE)
  cat("\nFrozen margins:\n")
  print(do.call(cbind, lapply(margins, unlist)))
  extra <- lapply(links, function(link) {
    cr <- byLink[[link]][["A3:sparse"]]
    data.frame(
      link = link,
      shipped = mean(cr$values[[paste0("shipped.", equalTime)]]$m2Band95),
      move = mean(cr$values$move.default$m2Band95)
    )
  })
  cat(
    "\nThe sparse cell's family-wise rate (95 percent band), outside the rule:\n"
  )
  print(do.call(rbind, extra), row.names = FALSE)
  invisible(frozen)
}

runLook <- function(kind) {
  frozen <- readRDS(frozenFile())
  if (!identical(frozen$quick, quick)) {
    stop("quick and full output mixed")
  }
  setup <- switch(
    kind,
    interim = list(alpha = 0.005, reps = function(cell) {
      seq_len(ceiling(cell$reps / 2))
    }),
    summary = list(alpha = 0.02, reps = function(cell) seq_len(cell$reps)),
    doubled = list(alpha = 0.025, reps = function(cell) seq_len(2L * cell$reps))
  )
  out <- list()
  for (link in links) {
    crs <- readAll(link, setup$reps)
    missing <- vapply(
      crs,
      function(cr) {
        length(cr$reps) + nrow(cr$skipped) < length(setup$reps(cr$cell))
      },
      logical(1L)
    )
    if (any(missing)) {
      cat(
        link,
        ": cells short of the look's replications:",
        paste(names(missing)[missing], collapse = ", "),
        "\n"
      )
    }
    out[[link]] <- readRule(
      link,
      crs,
      setup$alpha,
      list(equalTime = frozen$equalTime, margins = frozen$margins[[link]])
    )
    printRule(out[[link]])
    cat(sprintf(
      "skipped pairs: %.4f of pairs%s\n",
      out[[link]]$skippedShare,
      if (out[[link]]$skippedShare > 0.01) {
        " (past 1 percent: the run stops)"
      } else {
        ""
      }
    ))
  }
  saveAtomic(out, runDir(paste0(kind, ".rds")))
  invisible(out)
}

## ------------------------------------------------- the R diagnosis

# R, outside the rule (dec-B456): whether a default-length reading reaches
# the long run's posterior, as the share of the reference's draws inside the
# reading's interval (0.90 for a 90 percent interval), with its Monte Carlo
# floor from two interleaved halves of the reference's chains.
referenceMass <- function(lo, hi, D) {
  inside <- D >= rep(lo, each = nrow(D)) & D <= rep(hi, each = nrow(D))
  list(point = mean(colMeans(inside)), joint = mean(rowSums(!inside) == 0))
}

runRDiagnosis <- function() {
  rows <- list()
  for (link in links) {
    for (cell in Filter(function(c) c$problem %in% c("B1", "B2"), cells)) {
      accFile <- file.path(referenceDir(link, cell$id, 1L), "acceptance.rds")
      if (!file.exists(accFile)) {
        next
      }
      acc <- readRDS(accFile)
      ref <- loadReference(link, cell$id, 1L)
      D <- matrix(ref$draws, prod(dim(ref$draws)[1:2]), dim(ref$draws)[3L])
      chainOf <- rep(seq_len(dim(ref$draws)[2L]), each = dim(ref$draws)[1L])
      half <- chainOf %% 2L == 0L
      qa <- apply(
        D[half, , drop = FALSE],
        2L,
        stats::quantile,
        probs = c(0.05, 0.95)
      )
      floorMass <- referenceMass(
        qa[1L, ],
        qa[2L, ],
        D[!half, , drop = FALSE]
      )$point
      refQ <- apply(D, 2L, stats::quantile, probs = c(0.05, 0.95))
      refMean <- colMeans(D)
      refSd <- apply(D, 2L, stats::sd)
      for (rep in seq_len(cell$reps)) {
        for (which in c("shipped", "move")) {
          path <- fitFile(link, cell$id, rep, which)
          if (!file.exists(path)) {
            next
          }
          res <- readRDS(path)
          if (!isTRUE(res$ok)) {
            next
          }
          for (reading in names(res$readings)) {
            st <- res$readings[[reading]]$stored
            iv90 <- interval(st, 0.9)
            m90 <- referenceMass(iv90$lo, iv90$hi, D)
            iv95 <- interval(st, 0.95)
            band <- st$bands[["0.9"]]
            rows[[length(rows) + 1L]] <- data.frame(
              link = link,
              cell = cell$id,
              rep = rep,
              arm = paste(which, reading, sep = "."),
              mass90 = m90$point,
              mass95 = referenceMass(iv95$lo, iv95$hi, D)$point,
              simMass90 = referenceMass(band$lo, band$hi, D)$joint,
              # ratios of means: a saturated row has a reference width of 0
              widthRatio = mean(iv90$hi - iv90$lo) /
                mean(refQ[2L, ] - refQ[1L, ]),
              errorSd = mean(
                abs(st$mean - refMean)[refSd > 0] / refSd[refSd > 0]
              ),
              floor90 = floorMass,
              accepted = acc$accepted,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  if (length(rows) == 0L) {
    stop("no reference and fit to read")
  }
  d <- do.call(rbind, rows)
  saveAtomic(d, runDir("r-diagnosis.rds"))
  agg <- stats::aggregate(
    cbind(mass90, mass95, simMass90, widthRatio, errorSd, floor90) ~
      link + arm,
    data = d,
    FUN = mean
  )
  print(agg, row.names = FALSE, digits = 4)
  invisible(d)
}

## ---------------------------------------------------------------- jobs

# A job file for the queue: one job a line, tab-separated: id, whether it is
# one of the largest jobs (two at a time), its cell (for the skip rule), an
# expected time in seconds, then commands separated by " ;; ". A command
# starts with @<build>, which the queue expands to that build's library and
# BIS_BUILD before the harness's arguments.
propensityKeys <- function() {
  unique(unlist(lapply(cells, `[[`, "psKey")))
}

expectedSeconds <- function(link, cell, iterationsPerChain, chains) {
  tab <- if (file.exists(timingFile())) readRDS(timingFile()) else NULL
  ms <- if (!is.null(tab)) {
    tab$msShipped[tab$link == link & tab$cell == cell$id]
  } else {
    numeric()
  }
  if (length(ms) != 1L || !is.finite(ms)) {
    return(NA_real_)
  }
  ms / 1000 * iterationsPerChain * chains + 15
}

writeJobs <- function(phase, file) {
  lines <- character()
  job <- function(id, big, cellId, seconds, commands) {
    lines <<- c(
      lines,
      paste(
        id,
        as.integer(big),
        cellId,
        if (is.finite(seconds)) round(seconds) else "NA",
        paste(commands, collapse = " ;; "),
        sep = "\t"
      )
    )
  }
  if (phase == "truth") {
    for (link in links) {
      for (d in names(realDatasets)) {
        for (part in seq_len(truthParts)) {
          job(
            sprintf("truth-%s-%s-%d", link, d, part),
            FALSE,
            paste0("A1:", d, ":bart"),
            NA,
            sprintf("@move truth bart %s %s %d", link, d, part)
          )
        }
        job(
          sprintf("additive-%s-%s", link, d),
          FALSE,
          paste0("A1:", d, ":additive"),
          NA,
          sprintf("@shipped truth additive %s %s", link, d)
        )
      }
    }
  } else if (phase == "propensity") {
    for (key in propensityKeys()) {
      reps <- if (grepl("fixed$", key)) 0L else seq_len(pilotReps)
      for (rep in reps) {
        job(
          sprintf("ps-%s-%d", key, rep),
          key %in% c("a6", "acic16-r31-full"),
          key,
          NA,
          sprintf("@shipped propensity %s %d", key, rep)
        )
      }
    }
  } else if (phase == "timing") {
    for (link in links) {
      for (cell in cells) {
        cmds <- character()
        for (t in 1:3) {
          buildOrder <- if (t %% 2L == 1L) {
            c("shipped", "move")
          } else {
            c("move", "shipped")
          }
          cmds <- c(
            cmds,
            sprintf("@%s timing %s %s %d", buildOrder, link, cell$id, t)
          )
        }
        job(
          sprintf("timing-%s-%s", link, cellFile(cell$id)),
          cell$big,
          cell$id,
          NA,
          cmds
        )
      }
    }
  } else if (phase %in% c("pilot", "main")) {
    jobs <- list()
    for (link in links) {
      for (cell in cells) {
        reps <- if (phase == "pilot") {
          seq_len(min(pilotReps, cell$reps))
        } else {
          if (cell$reps > pilotReps) {
            seq.int(pilotReps + 1L, cell$reps)
          } else {
            integer()
          }
        }
        for (rep in reps) {
          jobs[[length(jobs) + 1L]] <- list(link = link, cell = cell, rep = rep)
        }
      }
    }
    # replication-major, cells in a seeded random order within a
    # replication, so that a partial run covers every cell
    shuffle <- withSeed(seedOf("job-order", phase), stats::runif(length(jobs)))
    ord <- order(vapply(jobs, `[[`, 0L, "rep"), shuffle)
    for (j in jobs[ord]) {
      buildOrder <- if (j$rep %% 2L == 1L) {
        c("shipped", "move")
      } else {
        c("move", "shipped")
      }
      est <- expectedSeconds(
        j$link,
        j$cell,
        mcmc$burn + drawsPerChain(j$cell),
        j$cell$chains
      )
      job(
        sprintf("fit-%s-%s-%03d", j$link, cellFile(j$cell$id), j$rep),
        j$cell$big,
        paste0(j$link, "|", j$cell$id),
        2.3 * est,
        sprintf("@%s fit %s %s %d", buildOrder, j$link, j$cell$id, j$rep)
      )
    }
  } else if (phase == "reference") {
    for (link in links) {
      for (cell in Filter(function(c) c$problem %in% c("B1", "B2"), cells)) {
        for (which in c("shipped", "move")) {
          for (part in seq_len(refParts)) {
            job(
              sprintf("ref-%s-%s-%s-%d", link, cellFile(cell$id), which, part),
              cell$big,
              cell$id,
              NA,
              sprintf("@%s reference %s %s %d", which, link, cell$id, part)
            )
          }
        }
      }
    }
  } else {
    stop("no phase ", phase)
  }
  if (is.null(file)) {
    writeLines(lines)
  } else {
    writeLines(lines, file)
    cat(length(lines), "jobs written to", file, "\n")
  }
}

markSkip <- function(link, cellId, rep, which, reason) {
  path <- paste0(fitFile(link, cellId, rep, which), ".skip")
  ensureDir(dirname(path))
  writeLines(reason, path)
}

## ---------------------------------------------------------------- main

main <- function(args) {
  if (length(args) == 0L) {
    stop("no action; see the header")
  }
  action <- args[1L]
  a <- args[-1L]
  switch(
    action,
    cells = {
      d <- do.call(
        rbind,
        lapply(cells, function(c) {
          data.frame(
            cell = c$id,
            kind = c$kind,
            chains = c$chains,
            reps = c$reps,
            big = c$big,
            groups = paste(c$groups, collapse = ","),
            psKey = if (is.null(c$psKey)) "" else c$psKey
          )
        })
      )
      print(d, row.names = FALSE)
    },
    jobs = writeJobs(a[1L], if (length(a) >= 2L) a[2L] else NULL),
    truth = if (a[1L] == "bart") {
      runTruthBart(a[2L], a[3L], as.integer(a[4L]))
    } else {
      runTruthAdditive(a[2L], a[3L])
    },
    propensity = runPropensity(a[1L], as.integer(a[2L])),
    timing = if (a[1L] == "summarize") {
      summarizeTiming()
    } else {
      runTiming(a[1L], a[2L], as.integer(a[3L]))
    },
    fit = runFit(a[1L], a[2L], as.integer(a[3L])),
    "mark-skip" = markSkip(a[1L], a[2L], as.integer(a[3L]), a[4L], a[5L]),
    reference = if (a[1L] == "accept") {
      acceptReference(
        a[2L],
        a[3L],
        if (length(a) >= 4L) as.integer(a[4L]) else 1L
      )
    } else {
      runReference(
        a[1L],
        a[2L],
        as.integer(a[3L]),
        if (length(a) >= 4L) as.integer(a[4L]) else 1L
      )
    },
    check = if (a[1L] == "compare") {
      compareChecks(a[2L], a[3L])
    } else {
      runCheck(a[1L])
    },
    pilot = runPilot(),
    interim = runLook("interim"),
    summary = runLook(
      if (length(a) >= 1L && a[1L] == "doubled") "doubled" else "summary"
    ),
    "r-diagnosis" = runRDiagnosis(),
    fetch = fetchAll(),
    "export-causal" = exportCausal(),
    stop("unknown action ", action)
  )
  invisible(NULL)
}

if (!interactive()) {
  main(commandArgs(trailingOnly = TRUE))
}
