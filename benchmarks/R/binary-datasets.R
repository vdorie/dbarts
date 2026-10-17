# The real binary datasets two harnesses share: the six from R and its
# recommended packages, then the sixteen UCI ones of benchmarks/R/uci-binary.R,
# in a fixed order, and the k study's subsampling constants. A dataset's seed
# in benchmarks/R/binary-hyperprior.R and benchmarks/R/binary-interval-study.R
# is a function of its position here, so the order is part of both studies'
# record: append, never reorder. The sourcing script defines harnessDir(), the
# directory uci-binary.R is sourced from on first use.

## ------------------------------------------------------------- real data

# Six binary-outcome datasets from R and its recommended packages: small to
# moderate n, base rates from 0.21 to 0.40, numeric and factor predictors,
# one dichotomized survival outcome. The sixteen UCI datasets are appended
# below, after these, so the seed of an existing dataset does not move.
realDatasets <- list(
  pima = function() {
    data("Pima.tr", package = "MASS", envir = environment())
    data("Pima.te", package = "MASS", envir = environment())
    d <- rbind(Pima.tr, Pima.te)
    list(x = d[, setdiff(names(d), "type")], y = as.integer(d$type == "Yes"))
  },
  biopsy = function() {
    data("biopsy", package = "MASS", envir = environment())
    d <- biopsy[complete.cases(biopsy), ]
    list(
      x = d[, paste0("V", 1:9)],
      y = as.integer(d$class == "malignant")
    )
  },
  infert = function() {
    d <- datasets::infert
    list(
      x = d[, c(
        "education",
        "age",
        "parity",
        "induced",
        "spontaneous"
      )],
      y = as.integer(d$case)
    )
  },
  kyphosis = function() {
    d <- rpart::kyphosis
    list(
      x = d[, c("Age", "Number", "Start")],
      y = as.integer(d$Kyphosis == "present")
    )
  },
  pbc = function() {
    d <- survival::pbc
    keep <- c(
      "age",
      "sex",
      "ascites",
      "hepato",
      "spiders",
      "edema",
      "bili",
      "chol",
      "albumin",
      "copper",
      "alk.phos",
      "ast",
      "trig",
      "platelet",
      "protime",
      "stage"
    )
    d <- d[!is.na(d$trt), c("status", keep)]
    d <- d[complete.cases(d), ]
    list(x = d[, keep], y = as.integer(d$status == 2L))
  },
  birthwt = function() {
    data("birthwt", package = "MASS", envir = environment())
    d <- birthwt
    d$race <- factor(d$race)
    d$ftv <- factor(pmin(d$ftv, 2L))
    list(
      x = d[, c("age", "lwt", "race", "smoke", "ptl", "ht", "ui", "ftv")],
      y = as.integer(d$low)
    )
  }
)

# The sixteen UCI datasets, in descending size. benchmarks/R/uci-binary.R
# holds the urls, the sha256 of every file and the per-dataset cleaning; it is
# sourced the first time a UCI block asks for data, so a simulated block or a
# block on one of the six above needs no cache and no network. The names are
# listed here rather than read from that file because `blocks` and the seed of
# a real block are both functions of this order, and neither should depend on
# whether anything has been downloaded yet.
uciDatasetNames <- c(
  "adult",
  "bank",
  "magic",
  "mushroom",
  "spambase",
  "banknote",
  "german",
  "tictactoe",
  "transfusion",
  "creditapproval",
  "wdbc",
  "climate",
  "ionosphere",
  "haberman",
  "cleveland",
  "sonar"
)

uciLoaded <- FALSE

ensureUci <- function() {
  if (!uciLoaded) {
    source(file.path(harnessDir(), "uci-binary.R"))
    uciLoaded <<- TRUE
  }
}

uciLoader <- function(name) {
  force(name)
  function() {
    ensureUci()
    d <- get("uciBinary")[[name]]()
    list(x = d[, setdiff(names(d), "y"), drop = FALSE], y = d$y)
  }
}

for (uciName in uciDatasetNames) {
  realDatasets[[uciName]] <- uciLoader(uciName)
}
rm(uciName)

# Above this many rows a split takes a fixed 4,000 training rows and 1,000
# held-out rows instead of 80/20 of everything, which is what keeps adult and
# bank inside the same per-fit budget as the rest.
subsampleAbove <- 5000L
subsampleTrain <- 4000L
subsampleTest <- 1000L
