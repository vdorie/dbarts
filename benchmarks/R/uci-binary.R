#!/usr/bin/env Rscript

# Sixteen binary-outcome datasets from the UCI Machine Learning Repository,
# fetched on demand, cached on disk and verified against a recorded sha256.
#
# WHY. benchmarks/R/binary-hyperprior.R scores the binary end-node hyperprior
# on six small datasets that ship with R and its recommended packages. The
# real-data column is the one that could move the package default, and six
# datasets cannot settle it. These sixteen widen that column to the breadth
# the original k-sensitivity work had: two hundred to fifty thousand rows,
# three to sixty predictors, positive rates from 0.085 to 0.65, seven
# datasets whose predictors are wholly or partly factors, and three with an
# imbalance strong enough to matter.
#
# WHAT EACH ONE IS. Rows and predictors are after the per-dataset cleaning
# described in its function; "rate" is the fraction of rows with y == 1.
#
#  name            id   rows  pred   rate  notes
#  adult            2  48842    14  0.239  census income, 8 factors
#  bank           222  45211    16  0.117  term deposit take-up, 9 factors
#  magic          159  19020    10  0.352  Cherenkov telescope, all numeric
#  mushroom        73   8124    21  0.482  all factors
#  spambase        94   4601    57  0.394  widest all-numeric
#  banknote       267   1372     4  0.445  wavelet features
#  german         144   1000    20  0.300  credit risk, 13 factors
#  tictactoe      101    958     9  0.653  all factors
#  transfusion    176    748     4  0.238  donor return
#  creditapproval  27    653    15  0.453  mixed types, complete cases
#  wdbc            17    569    30  0.373  breast cancer, diagnostic
#  climate        252    540    18  0.085  simulation crashes, imbalanced
#  ionosphere      52    351    33  0.641  radar returns
#  haberman        43    306     3  0.265  surgical survival
#  cleveland       45    297    13  0.461  heart disease, complete cases
#  sonar          151    208    60  0.534  widest p against smallest n
#
# SOURCES AND CITATIONS. Every file comes from archive.ics.uci.edu over https
# and is checked against the sha256 recorded in `uciSources` below; a mismatch
# is an error, not a warning, because a silently re-cut upstream file would
# move every score that depends on it. The repository's own citation policy
# asks that the repository be cited alongside the dataset:
#
#   Kelly, M., Longjohn, R. and Nottingham, K. The UCI Machine Learning
#   Repository. https://archive.ics.uci.edu
#
#  adult           Becker, B. and Kohavi, R. (1996). Adult. UCI ML Repository.
#  bank            Moro, S., Rita, P. and Cortez, P. (2014). Bank Marketing.
#  magic           Bock, R. (2004). MAGIC Gamma Telescope.
#  mushroom        (1981). Mushroom. From Lincoff, Audubon Field Guide.
#  spambase        Hopkins, M. et al. (1999). Spambase.
#  banknote        Lohweg, V. (2012). Banknote Authentication.
#  german          Hofmann, H. (1994). Statlog (German Credit Data).
#  tictactoe       Aha, D. (1991). Tic-Tac-Toe Endgame.
#  transfusion     Yeh, I. (2008). Blood Transfusion Service Center.
#  creditapproval  Quinlan, J. R. (1987). Credit Approval.
#  wdbc            Wolberg, W., Mangasarian, O., Street, N. and Street, W.
#                  (1993). Breast Cancer Wisconsin (Diagnostic).
#  climate         Lucas, D. et al. (2013). Climate Model Simulation Crashes.
#  ionosphere      Sigillito, V. et al. (1989). Ionosphere.
#  haberman        Haberman, S. (1976). Haberman's Survival.
#  cleveland       Janosi, A., Steinbrunn, W., Pfisterer, M. and Detrano, R.
#                  (1989). Heart Disease; the processed Cleveland file.
#  sonar           Sejnowski, T. and Gorman, R. (1988). Connectionist Bench
#                  (Sonar, Mines vs. Rocks).
#
# CACHE. Files land in $DBARTS_BENCH_DATA, or in tools::R_user_dir("dbarts",
# "cache") when that is unset. The directory is created on first use and is
# not part of the package or the repository; deleting it costs one download.
# Nothing here is run at package build, check or install time.
#
# Usage:
#   source("benchmarks/R/uci-binary.R")
#   d <- uciBinary$sonar()          # data.frame, 0/1 column y, typed x
#   uciBinaryFetchAll()             # warm the cache, report the table

## ------------------------------------------------------------------ cache

uciBase <- "https://archive.ics.uci.edu/ml/machine-learning-databases"

# One row per file, not per dataset: adult is two files, bank is a zip.
uciSources <- list(
  "adult.data" = list(
    url = file.path(uciBase, "adult/adult.data"),
    sha256 = "5b00264637dbfec36bdeaab5676b0b309ff9eb788d63554ca0a249491c86603d"
  ),
  "adult.test" = list(
    url = file.path(uciBase, "adult/adult.test"),
    sha256 = "a2a9044bc167a35b2361efbabec64e89d69ce82d9790d2980119aac5fd7e9c05"
  ),
  "bank.zip" = list(
    url = file.path(uciBase, "00222/bank.zip"),
    sha256 = "99d7e8eb12401ed278b793984423915411ea8df099e1795f9fefe254f513fe5e"
  ),
  "magic04.data" = list(
    url = file.path(uciBase, "magic/magic04.data"),
    sha256 = "e9314b7ebd4b4b59a3b3d65f7316663963777b16a46786877651dbbaa640b36a"
  ),
  "agaricus-lepiota.data" = list(
    url = file.path(uciBase, "mushroom/agaricus-lepiota.data"),
    sha256 = "e65d082030501a3ebcbcd7c9f7c71aa9d28fdfff463bf4cf4716a3fe13ac360e"
  ),
  "spambase.data" = list(
    url = file.path(uciBase, "spambase/spambase.data"),
    sha256 = "b1ef93de71f97714d3d7d4f58fc9f718da7bbc8ac8a150eff2778616a8097b12"
  ),
  "data_banknote_authentication.txt" = list(
    url = file.path(uciBase, "00267/data_banknote_authentication.txt"),
    sha256 = "d0539aaed2139ba7a587b3e34fb345ce503ff7d5d33dbf9912d8e195ce425cb9"
  ),
  "german.data" = list(
    url = file.path(uciBase, "statlog/german/german.data"),
    sha256 = "b21f3d81db8071257d5ff1deaeba1fd4303b62712e6fcc9715c7a86202cb5871"
  ),
  "tic-tac-toe.data" = list(
    url = file.path(uciBase, "tic-tac-toe/tic-tac-toe.data"),
    sha256 = "a996a756484d0cf8cbad1680d823aadbfcff8ae2d78c450e6999362d580301c5"
  ),
  "transfusion.data" = list(
    url = file.path(uciBase, "blood-transfusion/transfusion.data"),
    sha256 = "96c8e1091b9c037bcaf25a19b24b49d07771cd88689ffd273e056e9e8845ffe7"
  ),
  "crx.data" = list(
    url = file.path(uciBase, "credit-screening/crx.data"),
    sha256 = "fff49bc186cbddb3ace7371d40d9fbbb3af4f126019c13ff3f562249b1454f4d"
  ),
  "wdbc.data" = list(
    url = file.path(uciBase, "breast-cancer-wisconsin/wdbc.data"),
    sha256 = "d606af411f3e5be8a317a5a8b652b425aaf0ff38ca683d5327ffff94c3695f4a"
  ),
  "pop_failures.dat" = list(
    url = file.path(uciBase, "00252/pop_failures.dat"),
    sha256 = "8a892cbb056072e84749aed7409d6f17c7b4a1c899fc98ce52dc2af8b727ff9c"
  ),
  "ionosphere.data" = list(
    url = file.path(uciBase, "ionosphere/ionosphere.data"),
    sha256 = "46d52186b84e20be52918adb93e8fb9926b34795ff7504c24350ae0616a04bbd"
  ),
  "haberman.data" = list(
    url = file.path(uciBase, "haberman/haberman.data"),
    sha256 = "b4b7a32586a5668f9f4d6dc8be9d1bc8cd4822523affb1f6b5bfc350681ef3e2"
  ),
  "processed.cleveland.data" = list(
    url = file.path(uciBase, "heart-disease/processed.cleveland.data"),
    sha256 = "a74b7efa387bc9d108d7d0115d831fe9b414b29ae7124f331b622b4efa0427c8"
  ),
  "sonar.all-data" = list(
    url = file.path(
      uciBase,
      "undocumented/connectionist-bench/sonar/sonar.all-data"
    ),
    sha256 = "e90434cdbf00fcf93ffa911fe447ae25606979658e60f1d32e155c3b5240234d"
  )
)

uciCacheDir <- function() {
  dir <- Sys.getenv("DBARTS_BENCH_DATA")
  if (!nzchar(dir)) {
    dir <- tools::R_user_dir("dbarts", "cache")
  }
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  dir
}

# tools::sha256sum is recent; shasum is everywhere this harness runs.
uciDigest <- function(path) {
  if (exists("sha256sum", where = asNamespace("tools"))) {
    return(unname(get("sha256sum", asNamespace("tools"))(path)))
  }
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", out[1L])
}

# Downloads once, verifies always, and never leaves a partial or wrong file
# in the cache: a failed check deletes the file so the next call refetches.
uciFile <- function(key) {
  source <- uciSources[[key]]
  if (is.null(source)) {
    stop("no such UCI file: ", key)
  }
  path <- file.path(uciCacheDir(), key)
  if (!file.exists(path)) {
    temp <- paste0(path, ".part")
    status <- tryCatch(
      utils::download.file(
        source$url,
        temp,
        method = "curl",
        extra = "-sSfL --retry 3 --max-time 600",
        mode = "wb",
        quiet = TRUE
      ),
      error = function(e) {
        unlink(temp)
        stop("download failed for ", source$url, ": ", conditionMessage(e))
      }
    )
    if (!identical(as.integer(status), 0L) || !file.exists(temp)) {
      unlink(temp)
      stop("download failed for ", source$url)
    }
    file.rename(temp, path)
  }
  got <- uciDigest(path)
  if (!identical(got, source$sha256)) {
    unlink(path)
    stop(sprintf(
      paste0(
        "sha256 mismatch for %s\n  expected %s\n  got      %s\n",
        "the cached copy has been removed; the upstream file may have changed"
      ),
      key,
      source$sha256,
      got
    ))
  }
  path
}

# bank ships its csv inside a zip; extracted once, next to the zip.
uciZipMember <- function(key, member) {
  archive <- uciFile(key)
  path <- file.path(uciCacheDir(), member)
  if (!file.exists(path)) {
    utils::unzip(archive, files = member, exdir = uciCacheDir())
  }
  path
}

## -------------------------------------------------------------- helpers

# A 0/1 outcome named y in the last column, predictors ahead of it, and no
# stray row names, which is the shape the hyperprior harness splits.
uciFrame <- function(x, y) {
  y <- as.integer(y)
  if (!all(y %in% c(0L, 1L))) {
    stop("outcome is not 0/1")
  }
  out <- data.frame(x, y = y, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

# Columns with one level carry no rule and are dropped; which ones these are
# is dataset-specific and noted where it happens.
dropConstant <- function(x) {
  keep <- vapply(
    x,
    function(column) length(unique(column[!is.na(column)])) > 1L,
    logical(1L)
  )
  x[, keep, drop = FALSE]
}

asFactors <- function(x, columns) {
  for (column in columns) {
    x[[column]] <- factor(x[[column]])
  }
  x
}

## ------------------------------------------------------------- datasets

uciBinary <- list(
  # Census income. The training and test files are one sample split by the
  # donor; both are used. The test file carries a banner line and a trailing
  # period on every label. "?" is kept as a factor level rather than dropping
  # the 3,620 rows that carry one, which would change the positive rate.
  adult = function() {
    columns <- c(
      "age",
      "workclass",
      "fnlwgt",
      "education",
      "education.num",
      "marital.status",
      "occupation",
      "relationship",
      "race",
      "sex",
      "capital.gain",
      "capital.loss",
      "hours.per.week",
      "native.country",
      "income"
    )
    read <- function(path, skip) {
      utils::read.table(
        path,
        sep = ",",
        header = FALSE,
        skip = skip,
        col.names = columns,
        strip.white = TRUE,
        colClasses = "character",
        na.strings = character(0L)
      )
    }
    d <- rbind(read(uciFile("adult.data"), 0L), read(uciFile("adult.test"), 1L))
    d <- d[nzchar(d$income), ]
    y <- as.integer(sub("\\.$", "", d$income) == ">50K")
    numeric <- c(
      "age",
      "fnlwgt",
      "education.num",
      "capital.gain",
      "capital.loss",
      "hours.per.week"
    )
    x <- d[, setdiff(columns, "income"), drop = FALSE]
    for (column in numeric) {
      x[[column]] <- as.numeric(x[[column]])
    }
    x <- asFactors(x, setdiff(names(x), numeric))
    uciFrame(x, y)
  },

  # Portuguese bank telemarketing, the full 45,211-row file. The outcome
  # column is itself called y upstream; it becomes the 0/1 y here.
  bank = function() {
    d <- utils::read.table(
      uciZipMember("bank.zip", "bank-full.csv"),
      sep = ";",
      header = TRUE,
      quote = "\"",
      stringsAsFactors = TRUE
    )
    y <- as.integer(d$y == "yes")
    uciFrame(d[, setdiff(names(d), "y"), drop = FALSE], y)
  },

  # Cherenkov telescope showers; hadron is the positive class, which is the
  # minority and, per the donor, the costly error.
  magic = function() {
    columns <- c(
      "fLength",
      "fWidth",
      "fSize",
      "fConc",
      "fConc1",
      "fAsym",
      "fM3Long",
      "fM3Trans",
      "fAlpha",
      "fDist",
      "class"
    )
    d <- utils::read.table(
      uciFile("magic04.data"),
      sep = ",",
      header = FALSE,
      col.names = columns
    )
    uciFrame(
      d[, setdiff(columns, "class"), drop = FALSE],
      as.integer(d$class == "h")
    )
  },

  # Poisonous is the positive class. "?" in stalk.root is a level; veil.type
  # is constant and is dropped, which is why this is 21 predictors and not 22.
  mushroom = function() {
    columns <- c(
      "class",
      "cap.shape",
      "cap.surface",
      "cap.color",
      "bruises",
      "odor",
      "gill.attachment",
      "gill.spacing",
      "gill.size",
      "gill.color",
      "stalk.shape",
      "stalk.root",
      "stalk.surface.above.ring",
      "stalk.surface.below.ring",
      "stalk.color.above.ring",
      "stalk.color.below.ring",
      "veil.type",
      "veil.color",
      "ring.number",
      "ring.type",
      "spore.print.color",
      "population",
      "habitat"
    )
    d <- utils::read.table(
      uciFile("agaricus-lepiota.data"),
      sep = ",",
      header = FALSE,
      col.names = columns,
      colClasses = "character",
      na.strings = character(0L)
    )
    x <- dropConstant(d[, setdiff(columns, "class"), drop = FALSE])
    uciFrame(asFactors(x, names(x)), as.integer(d$class == "p"))
  },

  # Word and character frequencies plus three run-length summaries; spam is
  # the positive class.
  spambase = function() {
    d <- utils::read.table(uciFile("spambase.data"), sep = ",", header = FALSE)
    names(d) <- c(sprintf("f%02d", seq_len(ncol(d) - 1L)), "spam")
    uciFrame(d[, seq_len(ncol(d) - 1L), drop = FALSE], d$spam)
  },

  # Wavelet moments of banknote images; forgery is class 1 upstream.
  banknote = function() {
    d <- utils::read.table(
      uciFile("data_banknote_authentication.txt"),
      sep = ",",
      header = FALSE,
      col.names = c("variance", "skewness", "curtosis", "entropy", "class")
    )
    uciFrame(d[, 1:4, drop = FALSE], d$class)
  },

  # Statlog German credit, the A-coded file. Bad risk (upstream 2) is the
  # positive class, so the rate is the 0.30 the cost matrix is built around.
  german = function() {
    columns <- c(
      "status",
      "duration",
      "history",
      "purpose",
      "amount",
      "savings",
      "employment",
      "installment.rate",
      "personal",
      "debtors",
      "residence",
      "property",
      "age",
      "other.plans",
      "housing",
      "existing.credits",
      "job",
      "dependents",
      "telephone",
      "foreign.worker",
      "risk"
    )
    d <- utils::read.table(
      uciFile("german.data"),
      header = FALSE,
      col.names = columns,
      stringsAsFactors = TRUE
    )
    uciFrame(
      d[, setdiff(columns, "risk"), drop = FALSE],
      as.integer(d$risk == 2L)
    )
  },

  # Nine board squares, three levels each; a win for x is positive.
  tictactoe = function() {
    d <- utils::read.table(
      uciFile("tic-tac-toe.data"),
      sep = ",",
      header = FALSE,
      col.names = c(
        paste0("s", 1:9),
        "class"
      ),
      stringsAsFactors = TRUE
    )
    uciFrame(
      d[, paste0("s", 1:9), drop = FALSE],
      as.integer(d$class == "positive")
    )
  },

  # Recency, frequency, monetary and time; the outcome is a March 2007
  # donation. Monetary is 250 times frequency by construction and is kept,
  # since a duplicated column is a fair thing for a tree prior to face.
  transfusion = function() {
    d <- utils::read.csv(uciFile("transfusion.data"), header = TRUE)
    names(d) <- c("recency", "frequency", "monetary", "time", "donated")
    uciFrame(d[, 1:4, drop = FALSE], d$donated)
  },

  # Anonymized credit applications, mixed types with "?" scattered over
  # seven columns; 37 incomplete rows are dropped rather than coded, because
  # here the missingness is thin enough that a level would be mostly noise.
  creditapproval = function() {
    d <- utils::read.table(
      uciFile("crx.data"),
      sep = ",",
      header = FALSE,
      col.names = c(sprintf("a%02d", 1:15), "class"),
      na.strings = "?",
      colClasses = c(
        "factor",
        "numeric",
        "numeric",
        "factor",
        "factor",
        "factor",
        "factor",
        "numeric",
        "factor",
        "factor",
        "numeric",
        "factor",
        "factor",
        "numeric",
        "numeric",
        "character"
      )
    )
    d <- d[stats::complete.cases(d), ]
    x <- d[, sprintf("a%02d", 1:15), drop = FALSE]
    x <- asFactors(x, names(x)[vapply(x, is.factor, logical(1L))])
    uciFrame(x, as.integer(d$class == "+"))
  },

  # Thirty cell-nucleus features; malignant is positive. Column 1 is an
  # identifier and is dropped.
  wdbc = function() {
    d <- utils::read.table(uciFile("wdbc.data"), sep = ",", header = FALSE)
    names(d) <- c("id", "diagnosis", sprintf("f%02d", 1:30))
    uciFrame(
      d[, sprintf("f%02d", 1:30), drop = FALSE],
      as.integer(d$diagnosis == "M")
    )
  },

  # Parallel ocean model runs over eighteen scaled parameters; the positive
  # class is a crash, which is 46 of 540 rows. Study and Run index the Latin
  # hypercube and are not predictors.
  climate = function() {
    d <- utils::read.table(uciFile("pop_failures.dat"), header = TRUE)
    x <- d[, setdiff(names(d), c("Study", "Run", "outcome")), drop = FALSE]
    uciFrame(x, as.integer(d$outcome == 0L))
  },

  # Radar returns; "good" (a structure detected) is positive. The second
  # column is zero for every row and is dropped, leaving 33 predictors.
  ionosphere = function() {
    d <- utils::read.table(
      uciFile("ionosphere.data"),
      sep = ",",
      header = FALSE
    )
    names(d) <- c(sprintf("v%02d", 1:34), "class")
    x <- dropConstant(d[, sprintf("v%02d", 1:34), drop = FALSE])
    uciFrame(x, as.integer(d$class == "g"))
  },

  # Age, year of operation and positive axillary nodes; death within five
  # years is positive.
  haberman = function() {
    d <- utils::read.table(
      uciFile("haberman.data"),
      sep = ",",
      header = FALSE,
      col.names = c("age", "year", "nodes", "status")
    )
    uciFrame(d[, 1:3, drop = FALSE], as.integer(d$status == 2L))
  },

  # The processed Cleveland heart-disease file. The upstream outcome is a
  # severity 0-4; any disease is positive. Six rows carry "?" in ca or thal
  # and are dropped. cp, restecg, slope and thal are unordered codes.
  cleveland = function() {
    columns <- c(
      "age",
      "sex",
      "cp",
      "trestbps",
      "chol",
      "fbs",
      "restecg",
      "thalach",
      "exang",
      "oldpeak",
      "slope",
      "ca",
      "thal",
      "num"
    )
    d <- utils::read.table(
      uciFile("processed.cleveland.data"),
      sep = ",",
      header = FALSE,
      col.names = columns,
      na.strings = "?"
    )
    d <- d[stats::complete.cases(d), ]
    x <- d[, setdiff(columns, "num"), drop = FALSE]
    x <- asFactors(x, c("cp", "restecg", "slope", "thal"))
    uciFrame(x, as.integer(d$num > 0L))
  },

  # Sixty energy bands from a sonar return; a mine is positive. The widest
  # predictor set against the smallest sample in this set.
  sonar = function() {
    d <- utils::read.table(uciFile("sonar.all-data"), sep = ",", header = FALSE)
    names(d) <- c(sprintf("b%02d", 1:60), "class")
    uciFrame(
      d[, sprintf("b%02d", 1:60), drop = FALSE],
      as.integer(d$class == "M")
    )
  }
)

## ---------------------------------------------------------------- driver

# Warms the cache and prints what the harness will see, which is also how the
# table in this header was produced.
uciBinaryFetchAll <- function(names. = names(uciBinary)) {
  rows <- lapply(names., function(name) {
    d <- uciBinary[[name]]()
    x <- d[, setdiff(names(d), "y"), drop = FALSE]
    data.frame(
      name = name,
      rows = nrow(d),
      predictors = ncol(x),
      factors = sum(vapply(x, is.factor, logical(1L))),
      rate = round(mean(d$y), 3),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  print(out, row.names = FALSE)
  invisible(out)
}

if (identical(environment(), globalenv()) && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0L && args[1L] == "fetch") {
    uciBinaryFetchAll(if (length(args) > 1L) args[-1L] else names(uciBinary))
  }
}
