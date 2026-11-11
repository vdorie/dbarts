# The real-data side of pricing (docs/plans/k-move-pricing.md, The real-data
# guard and Step 3; the sealed rule and list, k-move-pricing-datasets.md):
# the committed loader and its covariate-only pre-fit check, the sealed
# separation scores, the R2 and R3 assignment by the rule's own code, the
# guard's references (the k study's acceptance, each chain from its own draw
# from the prior), its scored fits, internal pilot and reading, the held-out
# score, and step 3's trigger selection and test. Sourced by
# k-move-pricing.R after the k study's harness.
#
# What is sealed (runDir("sealed")) is never printed: the outcomes, the
# prediction rows' outcomes and the separation scores. No file of these data,
# or derived from them, is committed or shipped; the archive (KMP_ARCHIVE)
# is private and outside the tree.

## ----------------------------------------------------------------- the list

# The 24 passers of the sealed list, in frame order, each with its file's
# source, the outcome named by the suite's default target or the donor, the
# columns dropped as not predictors, the columns the donor lists as
# categorical (set as factors; every other column keeps the file's type), and
# the listing's rows and predictors. The OpenML file ids and md5 sums are the
# suites' records, read from OpenML's metadata; SECOM and Vertebral Column
# enter only if step 3 runs (their second screens first).
passers <- data.frame(
  name = c(
    "APSFailure",
    "coil2000_insurance_policies",
    "credit_card_clients_default",
    "NATICUSdroid",
    "online_shoppers_intention",
    "polish_companies_bankruptcy",
    "qsar-biodeg",
    "seismic-bumps",
    "jm1",
    "ilpd",
    "ozone-level-8hr",
    "wilt",
    "Internet-Advertisements",
    "irish",
    "skin-segmentation",
    "spect",
    "breast_cancer",
    "hepatitis",
    "Parkinsons",
    "Diabetic_Retinopathy_Debrecen",
    "vote",
    "Toxicity",
    "Bone marrow transplant: children",
    "Mammographic Mass"
  ),
  frame = c(
    2L,
    8L,
    10L,
    24L,
    25L,
    26L,
    27L,
    28L,
    30L,
    51L,
    54L,
    62L,
    64L,
    70L,
    77L,
    97L,
    101L,
    126L,
    134L,
    136L,
    150L,
    152L,
    193L,
    200L
  ),
  source = c(
    rep("openml", 15L),
    rep("pmlb", 3L),
    "uci",
    "uci",
    "openml",
    "uci",
    "uci",
    "uci"
  ),
  id = c(
    "46908",
    "46916",
    "46919",
    "46969",
    "46947",
    "46950",
    "46952",
    "46956",
    "46979",
    "1480",
    "1487",
    "40983",
    "40978",
    "451",
    "1502",
    "spect",
    "breast_cancer",
    "hepatitis",
    "174",
    "329",
    "56",
    "728",
    "565",
    "161"
  ),
  fileId = c(
    "22125219",
    "22125227",
    "22125230",
    "22125386",
    "22125258",
    "22125261",
    "22125263",
    "22125267",
    "22125396",
    "1590565",
    "1592279",
    "18151926",
    "18140371",
    "52563",
    "1592294",
    "",
    "",
    "",
    "",
    "",
    "56",
    "",
    "",
    ""
  ),
  md5 = c(
    "3896ad0c77f3c96cc9eafad55b268863",
    "e4fd3362bf8126f724475ea2f318c123",
    "14af3e5dd1fa9af94c8d0f34c0068711",
    "7d290894bf8f7454baf61d6998aa2d56",
    "793c609105b8eacb5ed5acdcce64c00d",
    "7bb8841d3c0426943845f6a7e1997e83",
    "48f4696362b73acc569fc13a1510b8c8",
    "4eb73b1e0e9b5233d1b0b3f457af662b",
    "879ec0e136c85f164775651f2eec099f",
    "155295549a459387fb71312caf1b8360",
    "6e6b7aeb382e30ae79149e4e5845cf2d",
    "aa7c0aaf3fd671f06bbf6a68be4754f7",
    "4f8829d20c886c86de4bc06c4eeb5cc6",
    "b7f36a00bc8f9e2ab73d530b6f3ccb3d",
    "5dd4d1f27c2b1c6361c18dfda6fb61d6",
    "",
    "",
    "",
    "",
    "",
    "3c16059c5b92f6551f720f97d0eccc09",
    "",
    "",
    ""
  ),
  target = c(
    "AirPressureSystemFailure",
    "MobileHomePolicy",
    "DefaultOnPaymentNextMonth",
    "Malware",
    "Revenue",
    "company_bankrupt",
    "Biodegradable",
    "HighEnergySeismicBump",
    "defects",
    "Class",
    "Class",
    "class",
    "class",
    "Leaving_Certificate",
    "Class",
    "target",
    "target",
    "target",
    "status",
    "Class",
    "Class",
    "Class",
    "survival_status",
    "Severity"
  ),
  drop = c(
    rep("", 10L),
    "Date",
    rep("", 7L),
    "name",
    rep("", 4L),
    "BI-RADS"
  ),
  factors = c(
    rep("", 4L),
    "Month;VisitorType",
    rep("", 11L),
    "age;menopause;tumor-size;inv-nodes;breast-quad",
    rep("", 5L),
    "DonorABO;RecipientABO;CMVstatus;Disease;HLAmatch;Antigen;Allele;HLAgrI;Recipientageint",
    "Shape;Margin"
  ),
  rows = c(
    76000L,
    9822L,
    30000L,
    7491L,
    12330L,
    5910L,
    1054L,
    2584L,
    10885L,
    583L,
    2534L,
    4839L,
    3279L,
    500L,
    245057L,
    267L,
    286L,
    155L,
    195L,
    1151L,
    435L,
    171L,
    187L,
    961L
  ),
  predictors = c(
    170L,
    85L,
    23L,
    86L,
    17L,
    64L,
    41L,
    15L,
    21L,
    10L,
    72L,
    5L,
    1558L,
    5L,
    3L,
    22L,
    9L,
    19L,
    22L,
    19L,
    16L,
    1203L,
    38L,
    4L
  ),
  stringsAsFactors = FALSE
)
r2Size <- 12L

# The smoke (KMP_SMOKE=1) replaces the list with three stand-in datasets
# simulated here, so that no passer's file is read before its step.
smokeGuard <- nzchar(Sys.getenv("KMP_SMOKE"))
if (smokeGuard) {
  passers <- data.frame(
    name = c("standin a", "standin b", "standin c"),
    frame = 1:3,
    source = "uci",
    id = c("sa", "sb", "sc"),
    fileId = "",
    md5 = "",
    target = "outcome",
    drop = c("", "rowid", ""),
    factors = c("", "", "g"),
    rows = c(150L, 200L, 120L),
    predictors = c(4L, 4L, 5L),
    stringsAsFactors = FALSE
  )
  r2Size <- 2L
}

# The stand-ins' files, written to KMP_LOADER_SOURCE: covariates, an id
# column the loader drops, a factor, missing values and a binary outcome.
writeStandins <- function() {
  if (!smokeGuard) {
    stop("stand-ins belong to the smoke")
  }
  dir <- ensureDir(Sys.getenv("KMP_LOADER_SOURCE"))
  for (i in seq_len(nrow(passers))) {
    p <- passers[i, ]
    d <- withSeed(seedOf("standin", p$name), {
      n <- p$rows
      x <- data.frame(
        x1 = runif(n),
        x2 = rnorm(n),
        x3 = runif(n),
        x4 = rnorm(n)
      )
      if (p$factors == "g") {
        x$g <- sample(c("u", "v", "w"), n, replace = TRUE)
      }
      x$x2[sample.int(n, 3L)] <- NA
      eta <- 3 * (x$x1 - 0.5) + x$x4 + c(-1, 0, 1.5)[i]
      out <- data.frame(
        x,
        outcome = ifelse(rbinom(n, 1L, plogis(eta)) == 1L, "yes", "no")
      )
      if (p$drop == "rowid") {
        out$rowid <- seq_len(n)
      }
      out
    })
    utils::write.csv(
      d,
      file.path(dir, paste0(fileTag(p$name), ".csv")),
      row.names = FALSE
    )
  }
  cat("stand-ins written:", nrow(passers), "\n")
}
fitRowCap <- 4000L
guardFitsStart <- 40L
step3Fits <- 20L # step 3's default fits a library a dataset, no pilot raise

# The fits a library a dataset reads: the pilot's count on R2; step 3's 20
# on an R3 dataset (and a quick run's 4 either way before any pilot).
guardFitCount <- function(name) {
  pilotPath <- runDir("guard", fileTag(name), "pilot.rds")
  if (file.exists(pilotPath)) {
    return(readChecked(pilotPath)$fits)
  }
  a <- tryCatch(assignment(), error = function(e) NULL)
  if (quick) {
    4L
  } else if (!is.null(a) && name %in% a$r3) {
    step3Fits
  } else {
    guardFitsStart
  }
}
guardFitsRaised <- 80L
guardSpreadRaise <- 0.95
guardMargin <- 1
guardResamples <- function() if (quick) 200L else 2000L
guardIndependent <- 100L
quantileGrid <- (seq_len(50L) - 0.5) / 50
guardEstimands <- c("mean", "sd", "q05", "q95", "prob")
guardGroups <- c("all", "interior", "nearEnd")

regimes$G <- regime(
  "G",
  "real",
  "prob",
  kmpLink,
  4L,
  c(1L, 1L),
  c(1L, 1L),
  NA_integer_,
  refLooks = 1L
)
guardCell <- function(name) paste0("G:", gsub("[^A-Za-z0-9_.-]", "_", name))

sealedPath <- function(...) runDir("sealed", ...)
passerRow <- function(name) {
  i <- match(name, passers$name)
  if (is.na(i)) {
    stop("not a passer of the sealed list: ", name)
  }
  passers[i, ]
}
fileTag <- function(name) gsub("[^A-Za-z0-9_.-]", "_", name)

## ------------------------------------------------- the loader and pre-fit

archiveDir <- function() {
  a <- Sys.getenv("KMP_ARCHIVE")
  if (!nzchar(a)) {
    stop("KMP_ARCHIVE must name the private archive, outside the tree")
  }
  ensureDir(a)
}

sourceUrl <- function(p) {
  switch(
    p$source,
    openml = sprintf("https://www.openml.org/data/v1/download/%s", p$fileId),
    pmlb = sprintf(
      "https://github.com/EpistasisLab/pmlb/raw/master/datasets/%s/%s.tsv.gz",
      p$id,
      p$id
    ),
    uci = sprintf("https://archive.ics.uci.edu/static/public/%s/data.csv", p$id)
  )
}

# The file, from the archive when it is there, else fetched once and
# archived; its sha256 recorded at the first download and checked at every
# later load. KMP_LOADER_SOURCE names a directory holding the files under
# their archive names instead of the network (the smoke's stand-ins).
fetchPasser <- function(p) {
  ext <- switch(p$source, openml = "arff", pmlb = "tsv.gz", uci = "csv")
  local <- file.path(archiveDir(), paste0(fileTag(p$name), ".", ext))
  shaPath <- sealedPath("sha256", paste0(fileTag(p$name), ".txt"))
  if (!file.exists(local)) {
    stand <- Sys.getenv("KMP_LOADER_SOURCE")
    if (nzchar(stand)) {
      file.copy(file.path(stand, basename(local)), local)
    } else {
      utils::download.file(sourceUrl(p), local, mode = "wb", quiet = TRUE)
    }
  }
  sha <- sha256Of(local)
  if (file.exists(shaPath)) {
    if (!identical(readLines(shaPath, warn = FALSE), sha)) {
      stop("the archived file of ", p$name, " does not match its sha256")
    }
  } else {
    ensureDir(dirname(shaPath))
    writeLines(sha, shaPath)
  }
  md5 <- if (nzchar(p$md5)) unname(tools::md5sum(local)) else NA_character_
  list(
    path = local,
    sha256 = sha,
    md5Matches = if (nzchar(p$md5)) md5 == p$md5 else NA
  )
}

# foreign::read.arff refuses ARFF's integer type; the header is rewritten to
# numeric (same values) in a temporary copy, the archived file untouched.
readArff <- function(path) {
  lines <- readLines(path, warn = FALSE)
  d <- grep("^[[:space:]]*@data", lines, ignore.case = TRUE)[1L]
  if (is.na(d)) {
    stop("no @data line in ", path)
  }
  hdr <- seq_len(d - 1L)
  lines[hdr] <- sub(
    "^([[:space:]]*@attribute[[:space:]]+('[^']*'|\"[^\"]*\"|[^[:space:]]+)[[:space:]]+)integer\\b",
    "\\1numeric",
    lines[hdr],
    ignore.case = TRUE,
    perl = TRUE
  )
  tmp <- tempfile(fileext = ".arff")
  on.exit(unlink(tmp))
  writeLines(lines, tmp)
  foreign::read.arff(tmp)
}

parsePasser <- function(p, path) {
  switch(
    p$source,
    openml = readArff(path),
    pmlb = utils::read.delim(gzfile(path), check.names = FALSE),
    uci = utils::read.csv(
      path,
      check.names = FALSE,
      na.strings = c("", "NA", "?")
    )
  )
}

# The pre-fit check (the rule's section 4): read, record the sha256, remove
# the outcome immediately after parsing (kept sealed, never printed), keep
# the predictors as the source lists them, set the donor's categorical
# columns as factors, drop a column with no observed value; then the report
# from covariates only: rows and predictors against the listing, missing
# values per column, columns after factor expansion as bart() expands them.
runLoad <- function(name) {
  p <- passerRow(name)
  f <- fetchPasser(p)
  raw <- parsePasser(p, f$path)
  if (!p$target %in% names(raw)) {
    stop("no outcome column ", p$target, " in ", name)
  }
  y <- raw[[p$target]]
  raw[[p$target]] <- NULL
  saveAtomic(
    list(y = y, quick = quick),
    sealedPath("outcomes", paste0(fileTag(name), ".rds"))
  )
  rm(y)
  dropCols <- strsplit(p$drop, ";", fixed = TRUE)[[1L]]
  droppedPresent <- intersect(dropCols, names(raw))
  raw <- raw[, setdiff(names(raw), dropCols), drop = FALSE]
  facCols <- strsplit(p$factors, ";", fixed = TRUE)[[1L]]
  facMatched <- intersect(facCols, names(raw))
  for (cl in facMatched) {
    raw[[cl]] <- factor(raw[[cl]])
  }
  for (cl in names(raw)) {
    if (is.character(raw[[cl]])) {
      raw[[cl]] <- factor(raw[[cl]])
    }
  }
  empty <- names(raw)[vapply(raw, function(v) all(is.na(v)), NA)]
  raw <- raw[, setdiff(names(raw), empty), drop = FALSE]
  expanded <- sum(vapply(
    raw,
    function(v) if (is.factor(v) && nlevels(v) > 2L) nlevels(v) else 1L,
    0L
  ))
  missing <- vapply(raw, function(v) sum(is.na(v)), 0L)
  report <- list(
    name = name,
    frame = p$frame,
    sha256 = f$sha256,
    md5Matches = f$md5Matches,
    rows = nrow(raw),
    rowsListed = p$rows,
    predictors = ncol(raw) + length(empty),
    predictorsListed = p$predictors,
    droppedAsNotPredictors = droppedPresent,
    droppedEmpty = empty,
    donorFactorsUnmatched = setdiff(facCols, facMatched),
    expandedColumns = expanded,
    missingByColumn = missing[missing > 0L],
    rowsWithMissing = sum(stats::complete.cases(raw) == FALSE),
    pass = nrow(raw) == p$rows &&
      (is.na(p$predictors) || ncol(raw) + length(empty) == p$predictors) &&
      expanded <= 2000L,
    quick = quick
  )
  saveAtomic(raw, sealedPath("covariates", paste0(fileTag(name), ".rds")))
  saveAtomic(report, runDir("prefit", paste0(fileTag(name), ".rds")))
  cat(sprintf(
    "prefit %s: %d rows (listed %d), %d predictors (listed %d), %d columns expanded: %s\n",
    name,
    report$rows,
    p$rows,
    report$predictors,
    p$predictors,
    expanded,
    if (report$pass) "pass" else "FAIL"
  ))
  invisible(report)
}

## -------------------------------------- the sealed scores and assignment

sealedRng <- function(seed) {
  set.seed(
    seed,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
}

# A passer's rows (the rule's section 6, item 4): the outcome attached, a
# row with a missing outcome dropped, then m = min(n, 5,000) rows drawn, four
# fifths fitted and the rest prediction rows. The outcome is coded 1 for the
# minority class.
passerRows <- function(name) {
  p <- passerRow(name)
  x <- readRDS(sealedPath("covariates", paste0(fileTag(name), ".rds")))
  y <- readRDS(sealedPath("outcomes", paste0(fileTag(name), ".rds")))$y
  y <- if (is.factor(y)) as.character(y) else y
  keep <- !is.na(y)
  x <- x[keep, , drop = FALSE]
  y <- y[keep]
  lev <- names(sort(table(y)))
  if (length(lev) != 2L) {
    stop(name, " has not two outcome classes")
  }
  y01 <- as.integer(y == lev[1L])
  n <- nrow(x)
  sealedRng(20261109L + p$frame)
  m <- min(n, 5000L)
  rows <- sort(sample.int(n, m))
  fit <- sort(rows[sample.int(m, (4L * m) %/% 5L)])
  pred <- setdiff(rows, fit)
  sealedRng(20262109L + p$frame)
  folds <- sample(rep_len(1:5, length(fit)))
  list(
    x = x,
    y = y01,
    fit = fit,
    pred = pred,
    folds = folds,
    droppedMissingOutcome = sum(!keep)
  )
}

# The k study's separation score with the rule's folds in place of its keyed
# seed: 5-fold cross-validated rpart AUC on the fitted rows.
sealedScore <- function(x, y, folds) {
  pred <- numeric(length(y))
  d <- data.frame(x, .y = factor(y, levels = 0:1))
  for (f in 1:5) {
    test <- folds == f
    train <- d[!test, , drop = FALSE]
    if (length(unique(train$.y)) < 2L) {
      pred[test] <- mean(as.integer(as.character(train$.y)))
      next
    }
    fit <- rpart::rpart(
      .y ~ .,
      data = train,
      method = "class",
      control = rpart::rpart.control(xval = 0L)
    )
    pred[test] <- stats::predict(fit, d[test, , drop = FALSE], type = "prob")[,
      2L
    ]
  }
  auc(pred, y)
}

# Every passer's score, written sealed; the pre-fit check must have passed on
# every passer first (the draw needs every score).
runSealedScores <- function() {
  reports <- lapply(passers$name, function(nm) {
    readChecked(runDir("prefit", paste0(fileTag(nm), ".rds")))
  })
  failed <- passers$name[!vapply(reports, `[[`, NA, "pass")]
  if (length(failed)) {
    stop(
      "the pre-fit check failed on ",
      paste(failed, collapse = ", "),
      ": the census reads on before any score"
    )
  }
  scores <- vapply(
    passers$name,
    function(nm) {
      r <- passerRows(nm)
      sealedScore(r$x[r$fit, , drop = FALSE], r$y[r$fit], r$folds)
    },
    0
  )
  saveAtomic(
    list(score = scores, frame = passers$frame, quick = quick),
    sealedPath("scores.rds")
  )
  cat("sealed scores written for", length(scores), "passers\n")
}

# The allocation by the rule's own code: R2 one from each of 12 blocks by
# score, R3 every other passer. Only names are written and printed.
runAssign <- function() {
  s <- readChecked(sealedPath("scores.rds"))
  ord <- names(s$score)[order(s$score, s$frame)]
  drawBlocks <- function(o, k) {
    b <- ceiling(k * seq_along(o) / length(o))
    o[vapply(
      seq_len(k),
      function(j) {
        i <- which(b == j)
        i[sample.int(length(i), 1L)]
      },
      1L
    )]
  }
  sealedRng(20261108L)
  r2 <- drawBlocks(ord, r2Size)
  r3 <- ord[!ord %in% r2]
  out <- list(r2 = r2, r3 = sort(r3), quick = quick)
  saveAtomic(out, runDir("assignment.rds"))
  cat("R2:", paste(r2, collapse = "; "), "\n")
  cat("R3:", paste(sort(r3), collapse = "; "), "\n")
  invisible(out)
}

assignment <- function() readChecked(runDir("assignment.rds"))

# A dataset is fitted only by the step that reads it: R2 by the guard, R3 by
# step 3 once its whole form is committed (step3 select).
guardAllowed <- function(name) {
  a <- assignment()
  if (name %in% a$r2) {
    return(TRUE)
  }
  if (name %in% a$r3 && file.exists(runDir("step3", "commit.rds"))) {
    return(TRUE)
  }
  stop(name, " is not read by any step that has started")
}

# The guard's design: fitted rows (at most 4,000) and prediction rows, whose
# outcomes are written sealed apart.
runGuardDesign <- function(name) {
  guardAllowed(name)
  r <- passerRows(name)
  fit <- r$fit
  if (length(fit) > fitRowCap) {
    stop("more fitted rows than the cap")
  }
  x <- r$x
  rownames(x) <- NULL
  d <- list(
    regime = "G",
    cell = guardCell(name),
    link = kmpLink,
    rep = 1L,
    task = "prob",
    dataset = name,
    x = x[fit, , drop = FALSE],
    xh = x[r$pred, , drop = FALSE],
    y = r$y[fit],
    droppedMissingOutcome = r$droppedMissingOutcome,
    quick = quick
  )
  rownames(d$x) <- NULL
  rownames(d$xh) <- NULL
  saveAtomic(
    list(y = r$y[r$pred], quick = quick),
    sealedPath("heldout", paste0(fileTag(name), ".rds"))
  )
  saveAtomic(d, designFile(d$cell, kmpLink, 1L))
  cat(sprintf(
    "guard design %s: %d fitted, %d prediction rows\n",
    name,
    nrow(d$x),
    nrow(d$xh)
  ))
}

## ---------------------------------------------------- the guard's references

# The guard's reference parts: one chain a part, started from its own draw
# from the model's prior (trees, leaves and k), then the k study's form: refBurn
# sweeps of burn-in, refFirst recorded every refThin, refTail unthinned. A
# guard cell's parts take this in place of the k study's referenceSampler,
# which the reference jobs call; every other cell keeps the k study's.
kstudyReferenceSampler <- referenceSampler

# A one-chain sampler returns its fits as rows x draws; the k study's binder,
# which a continuation round also calls, takes rows x draws x chains.
# The smoke forces its first stand-in up the reference ladder: every round
# before the first past refBlockMax's next opens the next (extend while it
# stays within refBlockMax, then one escalation to the shipped half alone),
# so the one-chain parts run extendReferencePart on both kinds of round; the
# round after decides as written.
if (kmpSmoke) {
  kstudyAcceptReference <- acceptReference
  kstudyGuardDecision <- guardDecision
  forcedLadder <- FALSE
  acceptReference <- function(cell, link, rep, round) {
    first <- guardCell(passers$name[1L])
    nextL <- refLengths()[round + 2L]
    forcedLadder <<- identical(cell, first) &&
      !is.na(nextL) &&
      refLengths()[round + 1L] <= mcmc$refBlockMax
    kstudyAcceptReference(cell, link, rep, round)
  }
  guardDecision <- function(cur, prev, nextL, projected, escalating) {
    if (forcedLadder) {
      return(if (nextL <= mcmc$refBlockMax) "extend" else "escalate")
    }
    kstudyGuardDecision(cur, prev, nextL, projected, escalating)
  }
}

kstudyAbind3 <- abind3
abind3 <- function(a, b) {
  as3d <- function(m) {
    if (length(dim(m)) == 2L) array(m, c(dim(m), 1L)) else m
  }
  kstudyAbind3(as3d(a), as3d(b))
}
referenceSampler <- function(design, link, task, seed, part, which, dir) {
  if (!startsWith(design$cell, "G:")) {
    return(kstudyReferenceSampler(design, link, task, seed, part, which, dir))
  }
  s <- bart(
    y ~ .,
    data = data.frame(design$x, y = design$y),
    test = design$xh,
    family = link,
    n.trees = nTrees,
    n.chains = refPartChains,
    n.threads = 1L,
    n.burn = mcmc$refBurn,
    n.samples = mcmc$refFirst,
    n.thin = mcmc$refThin,
    seed = seed,
    verbose = FALSE,
    samplerOnly = TRUE,
    keepTrainingFits = TRUE
  )
  k0 <- withSeed(seedOf("reference-start", seed), priorStart(s))
  r1 <- s$run()
  saveSampler(s, samplerFile(dir, which, part))
  ctrl <- s$control
  ctrl@n.thin <- 1L
  s$setControl(ctrl)
  r2 <- s$run(0L, mcmc$refTail)
  pick <- function(r) if (is.null(r$test)) r$train else abind3(r$train, r$test)
  list(
    record = aperm(pick(r1), c(3L, 2L, 1L)),
    tail = aperm(pick(r2), c(3L, 2L, 1L)),
    k = r1$k,
    kTail = r2$k,
    chains = refPartChains,
    kStart = k0
  )
}

# The k study's acceptance with posteriordb's draw count read at the
# accepting end (the call with `project`), and agreement on the probability's
# posterior mean at every row beside the latent fit's mean and percentiles.
kstudyReferenceChecks <- referenceChecks
referenceChecks <- function(ref, design, link, labels, project = NULL) {
  out <- kstudyReferenceChecks(ref, design, link, labels, project)
  A <- poolDraws(ref$shipped$record, seq_len(dim(ref$shipped$record)[2L]))
  B <- poolDraws(ref$move$record, seq_len(dim(ref$move$record)[2L]))
  pz <- vapply(
    seq_len(ncol(A)),
    function(j) {
      a <- matrix(plogis(A[, j]), ncol = dim(ref$shipped$record)[1L])
      b <- matrix(plogis(B[, j]), ncol = dim(ref$move$record)[1L])
      abs(mean(a) - mean(b)) /
        sqrt(posterior::mcse_mean(a)^2 + posterior::mcse_mean(b)^2)
    },
    0
  )
  out$checks <- c(out$checks, probAgree = all(pz <= 4.5, na.rm = TRUE))
  if (!is.null(project)) {
    draws <- nrow(B) + nrow(evenRows(A, nrow(B)))
    out$checks <- c(out$checks, drawCount = draws >= tenth(10000))
  } else {
    out$checks <- c(out$checks, drawCount = TRUE)
  }
  out$probZMax <- max(pz, na.rm = TRUE)
  out
}

guardReference <- function(name) {
  cell <- guardCell(name)
  acc <- finalAcceptance(cell, kmpLink, 1L)
  ref <- loadReference(cell, kmpLink, 1L)
  S <- poolDraws(ref$shipped$record, seq_len(dim(ref$shipped$record)[2L]))
  M <- poolDraws(ref$move$record, seq_len(dim(ref$move$record)[2L]))
  S <- evenRows(S, nrow(M))
  list(
    acc = acc,
    pooled = rbind(S, M),
    shipped = S,
    move = M,
    k = c(
      as.vector(evenRows(ref$shipped$k, nrow(as.matrix(ref$move$k)))),
      as.vector(ref$move$k)
    ),
    chainsS = dim(ref$shipped$record)[1L],
    chainsM = dim(ref$move$record)[1L]
  )
}

## ------------------------------------------------------ the scored fits

guardFitPath <- function(name, which, i, tag = "") {
  runDir(
    "guard",
    fileTag(name),
    sprintf(
      "%s%s-fit%03d.rds",
      which,
      if (nzchar(tag)) paste0("-", tag) else "",
      as.integer(i)
    )
  )
}

# The estimands at every fitted and prediction row: the latent fit's mean,
# standard deviation, 5th and 95th percentiles, and the probability's mean;
# beside them k's mean and percentiles, the average probability over the
# prediction rows, the 90 percent band and the quantile grid (for the
# Wasserstein distance).
fitEstimands <- function(L, k, nFit) {
  q <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    t(matrixStats::colQuantiles(L, probs = c(0.05, 0.95), type = 7L))
  } else {
    apply(L, 2L, stats::quantile, c(0.05, 0.95), names = FALSE)
  }
  P <- plogis(L)
  est <- rbind(
    mean = colMeans(L),
    sd = apply(L, 2L, stats::sd),
    q05 = q[1L, ],
    q95 = q[2L, ],
    prob = colMeans(P)
  )
  band <- summarizeDraws(L, bandLevels = 0.9)$bands[["0.9"]]
  grid <- apply(L, 2L, stats::quantile, quantileGrid, names = FALSE)
  list(
    est = est,
    k = c(
      mean = mean(k),
      q05 = stats::quantile(k, 0.05, names = FALSE),
      q95 = stats::quantile(k, 0.95, names = FALSE)
    ),
    avgProbPred = if (ncol(L) > nFit) {
      mean(rowMeans(P[, -seq_len(nFit), drop = FALSE]))
    } else {
      NA_real_
    },
    bandLo = band$lo,
    bandHi = band$hi,
    grid = grid
  )
}

runGuardFits <- function(name, from, to, tag = "") {
  guardAllowed(name)
  cell <- guardCell(name)
  design <- blind(loadDesign(cell, kmpLink, 1L))
  chains <- 4L
  n <- if (kmpBuild == "shipped") equalTimeDraws(cell) else mcmc$defaultDraws
  for (i in seq.int(from, to)) {
    path <- guardFitPath(name, kmpBuild, i, tag)
    if (file.exists(path)) {
      next
    }
    seed <- taggedSeedOf("guard-fit", name, kmpBuild, tag, i)
    res <- fitBart(design, kmpLink, chains, c(mcmc$burn, n), seed)
    L <- poolDraws(res$arr, seq_len(n))
    e <- fitEstimands(L, as.vector(res$k), nrow(design$x))
    # each row's effective draws (draws x chains), the median over rows
    essRow <- vapply(
      seq_len(dim(res$arr)[3L]),
      function(j) {
        m <- t(res$arr[, seq_len(n), j])
        c(posterior::ess_bulk(m), posterior::ess_tail(m))
      },
      numeric(2L)
    )
    e$essPerSecond <- c(
      bulk = stats::median(essRow[1L, ]),
      tail = stats::median(essRow[2L, ])
    ) /
      res$seconds
    gridPath <- sub("[.]rds$", "-grid.f32", path)
    writeSingle(e$grid, gridPath)
    e$grid <- NULL
    e$gridDims <- c(length(quantileGrid), ncol(L))
    saveAtomic(
      c(
        e,
        list(
          name = name,
          build = kmpBuild,
          i = as.integer(i),
          tag = tag,
          drawsPerChain = n,
          seconds = res$seconds,
          quick = quick
        )
      ),
      path
    )
    cat(sprintf(
      "guard fit %s %s %d: %d draws a chain, %.1f s\n",
      name,
      kmpBuild,
      as.integer(i),
      n,
      res$seconds
    ))
  }
}

## ----------------------------------------------------- errors and reading

# The reference's value of each estimand at each row, and its unit: the
# expected squared error of the same estimate from 100 independent draws of
# the pooled reference, by 2,000 resamples of 100 draws.
referenceValues <- function(M) {
  q <- apply(M, 2L, stats::quantile, c(0.05, 0.95), names = FALSE)
  rbind(
    mean = colMeans(M),
    sd = apply(M, 2L, stats::sd),
    q05 = q[1L, ],
    q95 = q[2L, ],
    prob = colMeans(plogis(M))
  )
}

referenceUnits <- function(pooled, values, name) {
  B <- guardResamples()
  acc <- matrix(0, nrow(values), ncol(values))
  withSeed(seedOf("guard-unit", name), {
    for (b in seq_len(B)) {
      d <- pooled[sample.int(nrow(pooled), guardIndependent), , drop = FALSE]
      acc <- acc + (referenceValues(d) - values)^2
    }
  })
  acc / B
}

guardLabels <- function(values) {
  values["prob", ] >= interiorRange[1L] & values["prob", ] <= interiorRange[2L]
}

loadGuardFits <- function(name, which, n, tag = "") {
  lapply(seq_len(n), function(i) {
    p <- guardFitPath(name, which, i, tag)
    if (file.exists(p)) readChecked(p) else NULL
  })
}

# Each fit's error on each row estimand and grouping against a reference
# value, in units of the error of 100 independent draws.
fitErrors <- function(fits, values, units, labels, nFit = NA_integer_) {
  groups <- list(
    all = rep(TRUE, length(labels)),
    interior = labels,
    nearEnd = !labels
  )
  if (!is.na(nFit)) {
    groups <- list(prediction = seq_along(labels) > nFit)
  }
  t(vapply(
    fits,
    function(f) {
      e <- (f$est - values)^2 / units
      unlist(lapply(groups, function(g) {
        if (!any(g)) {
          return(stats::setNames(rep(NA_real_, nrow(e)), rownames(e)))
        }
        rowMeans(e[, g, drop = FALSE])
      }))
    },
    numeric(length(guardEstimands) * length(groups))
  ))
}

skewness <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 3L || stats::sd(v) == 0) {
    return(NA_real_)
  }
  mean((v - mean(v))^3) / stats::sd(v)^3
}

# The internal pilot: after 40 fits a library, each library's spread and
# skewness of the error over its own fits on every row estimand and grouping;
# no difference is computed. A spread past 0.95 raises the dataset to 80.
runGuardPilot <- function(name) {
  guardAllowed(name)
  if (!file.exists(runDir("step3", "commit.rds"))) {
    stop("the guard's pilot is read after step 3 is committed")
  }
  r <- guardReference(name)
  if (!isTRUE(r$acc$accepted)) {
    cat("guard pilot", name, ": reference not accepted, no pilot\n")
    return(invisible(NULL))
  }
  values <- referenceValues(r$pooled)
  units <- referenceUnits(r$pooled, values, name)
  labels <- guardLabels(values)
  n0 <- if (quick) 4L else guardFitsStart
  per <- lapply(c("shipped", "move"), function(which) {
    fits <- Filter(Negate(is.null), loadGuardFits(name, which, n0))
    if (length(fits) < n0) {
      stop("the pilot needs ", n0, " fits of ", which)
    }
    E <- fitErrors(fits, values, units, labels)
    rbind(s = apply(E, 2L, stats::sd), skew = apply(E, 2L, skewness))
  })
  names(per) <- c("shipped", "move")
  raise <- any(vapply(
    per,
    function(m) any(m["s", ] > guardSpreadRaise, na.rm = TRUE),
    NA
  ))
  out <- list(
    name = name,
    spreads = per,
    fits = if (raise) {
      if (quick) 8L else guardFitsRaised
    } else {
      n0
    },
    quick = quick
  )
  saveAtomic(out, runDir("guard", fileTag(name), "pilot.rds"))
  cat("guard pilot", name, ":", out$fits, "fits a library\n")
  invisible(out)
}

# The pilot's raised count, the fits past the first 40 (none when the pilot
# kept 40).
runGuardMore <- function(name) {
  path <- runDir("guard", fileTag(name), "pilot.rds")
  if (!file.exists(path)) {
    cat("guard", name, ": no pilot (reference not accepted), no more fits\n")
    return(invisible(NULL))
  }
  p <- readChecked(path)
  n0 <- if (quick) 4L else guardFitsStart
  if (p$fits > n0) {
    runGuardFits(name, n0 + 1L, p$fits)
  } else {
    cat("guard", name, ": no fits past", n0, "\n")
  }
}

welch <- function(a, b, level = 0.975) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  va <- stats::var(a) / length(a)
  vb <- stats::var(b) / length(b)
  se <- sqrt(va + vb)
  df <- (va + vb)^2 / (va^2 / (length(a) - 1) + vb^2 / (length(b) - 1))
  D <- mean(a) - mean(b)
  tq <- stats::qt(level, df)
  c(D = D, se = se, df = df, lo = D - tq * se, hi = D + tq * se)
}

# One dataset's reading: per row estimand and grouping, D (the move's mean
# error minus the shipped fits') on Welch's t, harmless when its upper
# one-sided 0.025 bound is below 1, harm shown when D is at least 1 with its
# lower bound above 0; D again against each library's half alone; the error
# ratio with its bootstrap interval; each library's skew; and beside, framing
# nothing, k, the average probability, the band's reference mass, the
# reference mass inside each fit's interval, the Wasserstein distance.
guardDatasetReading <- function(name, tag = "") {
  r <- guardReference(name)
  acc <- r$acc
  base <- list(
    name = name,
    accepted = isTRUE(acc$accepted),
    acceptance = acc[c(
      "decision",
      "checks",
      "checksBlockBefore",
      "history",
      "kRhat",
      "kEss",
      "maxZ",
      "rhat99"
    )]
  )
  if (!base$accepted) {
    return(base)
  }
  nFit <- nrow(loadDesign(guardCell(name), kmpLink, 1L)$x)
  n <- guardFitCount(name)
  values <- referenceValues(r$pooled)
  units <- referenceUnits(r$pooled, values, name)
  labels <- guardLabels(values)
  vS <- referenceValues(r$shipped)
  vM <- referenceValues(r$move)
  fits <- lapply(c("shipped", "move"), function(w) {
    Filter(Negate(is.null), loadGuardFits(name, w, n, tag))
  })
  names(fits) <- c("shipped", "move")
  if (any(lengths(fits) < 2L)) {
    stop("guard ", name, ": fewer than two fits of a library to read")
  }
  errs <- lapply(
    fits,
    fitErrors,
    values = values,
    units = units,
    labels = labels
  )
  errS <- lapply(fits, fitErrors, values = vS, units = units, labels = labels)
  errM <- lapply(fits, fitErrors, values = vM, units = units, labels = labels)
  tests <- colnames(errs$move)
  tab <- do.call(
    rbind,
    lapply(tests, function(tt) {
      w <- welch(errs$move[, tt], errs$shipped[, tt])
      wS <- welch(errS$move[, tt], errS$shipped[, tt])
      wM <- welch(errM$move[, tt], errM$shipped[, tt])
      boot <- withSeed(seedOf("guard-ratio", name, tt), {
        replicate(bootstrapResamples, {
          a <- sample(errs$move[, tt], replace = TRUE)
          b <- sample(errs$shipped[, tt], replace = TRUE)
          mean(a, na.rm = TRUE) / mean(b, na.rm = TRUE)
        })
      })
      data.frame(
        test = tt,
        D = w[["D"]],
        se = w[["se"]],
        lo = w[["lo"]],
        hi = w[["hi"]],
        harmless = is.finite(w[["hi"]]) && w[["hi"]] < guardMargin,
        harm = is.finite(w[["lo"]]) && w[["D"]] >= guardMargin && w[["lo"]] > 0,
        harmlessVsShippedHalf = is.finite(wS[["hi"]]) &&
          wS[["hi"]] < guardMargin,
        harmVsMoveHalf = is.finite(wM[["lo"]]) &&
          wM[["D"]] >= guardMargin &&
          wM[["lo"]] > 0,
        errMove = mean(errs$move[, tt], na.rm = TRUE),
        errShipped = mean(errs$shipped[, tt], na.rm = TRUE),
        ratio = mean(errs$move[, tt], na.rm = TRUE) /
          mean(errs$shipped[, tt], na.rm = TRUE),
        ratioLo = stats::quantile(boot, 0.025, na.rm = TRUE, names = FALSE),
        ratioHi = stats::quantile(boot, 0.975, na.rm = TRUE, names = FALSE),
        skewMove = skewness(errs$move[, tt]),
        skewShipped = skewness(errs$shipped[, tt])
      )
    })
  )
  # a grouping with no rows on this dataset reads nothing
  tab$harmless[!is.finite(tab$D)] <- NA
  tab$ratioFlag <- tab$ratioLo > costBar
  tab$harmStands <- tab$harm & tab$harmVsMoveHalf
  tab$flagged <- tab$harmless & !tab$harmlessVsShippedHalf
  # beside, framing nothing
  P <- r$pooled
  inside <- function(f) {
    lo <- f$est["q05", ]
    hi <- f$est["q95", ]
    m <- colMeans(sweep(P, 2L, lo, ">=") & sweep(P, 2L, hi, "<="))
    c(all = mean(m), interior = mean(m[labels]), nearEnd = mean(m[!labels]))
  }
  bandMass <- function(f) {
    mean(apply(
      sweep(P, 2L, f$bandLo, ">=") & sweep(P, 2L, f$bandHi, "<="),
      1L,
      all
    ))
  }
  refGrid <- apply(P, 2L, stats::quantile, quantileGrid, names = FALSE)
  refSd <- values["sd", ]
  w1 <- function(f, which) {
    g <- readSingle(
      sub("[.]rds$", "-grid.f32", guardFitPath(name, which, f$i, tag)),
      f$gridDims
    )
    mean(colMeans(abs(g - refGrid)) / refSd)
  }
  beside <- lapply(names(fits), function(w) {
    list(
      mass = t(vapply(fits[[w]], inside, numeric(3L))),
      bandMass = vapply(fits[[w]], bandMass, 0),
      avgProbPred = vapply(fits[[w]], `[[`, 0, "avgProbPred"),
      k = t(vapply(fits[[w]], `[[`, numeric(3L), "k")),
      wasserstein = vapply(fits[[w]], w1, 0, which = w),
      essPerSecond = t(vapply(fits[[w]], `[[`, numeric(2L), "essPerSecond"))
    )
  })
  names(beside) <- names(fits)
  if (ncol(P) > nFit) {
    predOnly <- lapply(
      fits,
      fitErrors,
      values = values,
      units = units,
      labels = labels,
      nFit = nFit
    )
    beside$predictionRows <- data.frame(
      test = colnames(predOnly$move),
      errMove = colMeans(predOnly$move, na.rm = TRUE),
      errShipped = colMeans(predOnly$shipped, na.rm = TRUE)
    )
  }
  avgPred <- NULL
  if (ncol(P) > nFit) {
    # the average probability over the prediction rows, each half's chains
    # (draws x chains) for its R-hat and effective draws, reported beside
    avgOf <- function(rec) {
      a <- poolDraws(rec, seq_len(dim(rec)[2L]))
      matrix(
        rowMeans(plogis(a[, -seq_len(nFit), drop = FALSE])),
        ncol = dim(rec)[1L]
      )
    }
    ref0 <- loadReference(guardCell(name), kmpLink, 1L)
    mS <- avgOf(ref0$shipped$record)
    mM <- avgOf(ref0$move$record)
    both <- cbind(evenRows(mS, nrow(mM)), mM)
    avgPred <- c(
      mean = mean(both),
      rhat = posterior::rhat(both),
      essTail = posterior::ess_tail(both),
      essBulk = posterior::ess_bulk(both)
    )
    rm(ref0)
  }
  lag1 <- vapply(
    seq_len(ncol(r$pooled)),
    function(j) stats::acf(r$move[, j], lag.max = 1L, plot = FALSE)$acf[2L],
    0
  )
  c(
    base,
    list(
      fits = vapply(fits, length, 0L),
      table = tab,
      harmless = all(tab$harmless, na.rm = TRUE),
      harmShown = any(tab$harm),
      harmStands = any(tab$harmStands),
      beside = beside,
      reference = list(
        avgProbPred = avgPred,
        k = stats::quantile(r$k, c(0.05, 0.5, 0.95), names = FALSE),
        lag1Mean = mean(lag1, na.rm = TRUE)
      ),
      rowsInterior = sum(labels),
      rowsNearEnd = sum(!labels)
    )
  )
}

# The guard over datasets (Reading over datasets): on every row estimand and
# grouping, harmless only when every accepted dataset reads harmless on its
# own bound and the one-sided 0.025 bound of a one-sample t over the accepted
# datasets' D is below 1; a confirmed harm on any dataset is a harm. Beside:
# the exact bound on the share harmed, the Wilcoxon signed-rank test of D
# below 1, the sign count. Run only after step 3's commit (step3 select).
commitMd5 <- function() {
  unname(tools::md5sum(runDir("step3", "commit.rds")))
}

runGuardRead <- function(excludeDeterministic = FALSE) {
  if (!file.exists(runDir("step3", "commit.rds"))) {
    stop("no R2 reading is opened before step 3 is committed")
  }
  a <- assignment()
  names2 <- a$r2
  if (excludeDeterministic) {
    s <- readChecked(sealedPath("scores.rds"))$score
    names2 <- names2[s[names2] < deterministicScore]
  }
  per <- lapply(names2, guardDatasetReading)
  names(per) <- names2
  acc <- Filter(function(p) isTRUE(p$accepted), per)
  J <- length(acc)
  tests <- if (J) acc[[1L]]$table$test else character()
  over <- do.call(
    rbind,
    lapply(tests, function(tt) {
      D <- vapply(acc, function(p) p$table$D[p$table$test == tt], 0)
      each <- vapply(acc, function(p) p$table$harmless[p$table$test == tt], NA)
      harm <- vapply(
        acc,
        function(p) p$table$harmStands[p$table$test == tt],
        NA
      )
      # a dataset with no rows in this grouping has nothing to read on it
      has <- is.finite(D)
      D <- D[has]
      each <- each[has]
      harm <- harm[has]
      J <- length(D)
      tUpper <- if (J >= 2L) {
        mean(D) + stats::qt(0.975, J - 1L) * stats::sd(D) / sqrt(J)
      } else {
        NA_real_
      }
      wil <- if (J >= 2L) {
        suppressWarnings(
          stats::wilcox.test(D, mu = guardMargin, alternative = "less")$p.value
        )
      } else {
        NA_real_
      }
      data.frame(
        test = tt,
        accepted = J,
        everyHarmless = all(each),
        tUpper = tUpper,
        harmless = if (J == 0L) {
          NA
        } else {
          J >= 2L && all(each) && tUpper < guardMargin
        },
        harmShownOn = sum(harm),
        shareHarmedUpper = if (all(each)) 1 - 0.025^(1 / J) else NA_real_,
        wilcoxonP = wil,
        signsPositive = sum(D > 0)
      )
    })
  )
  J <- length(acc)
  verdict <- if (J == 0L) {
    "unresolved: no reference accepted"
  } else if (any(over$harmShownOn > 0)) {
    "harm shown: confirm on fresh fits before it frames anything"
  } else if (J == 1L) {
    "unresolved: one accepted dataset, no claim over datasets"
  } else if (all(over$harmless, na.rm = TRUE)) {
    "harmless"
  } else {
    "unresolved"
  }
  out <- list(
    datasets = per,
    over = over,
    verdict = verdict,
    excludeDeterministic = excludeDeterministic,
    step3CommitMd5 = commitMd5(),
    quick = quick
  )
  saveAtomic(
    out,
    runDir(
      "readings",
      if (excludeDeterministic) "guard-no995.rds" else "guard.rds"
    )
  )
  print(over, row.names = FALSE, digits = 3)
  cat("guard:", verdict, "\n")
  invisible(out)
}

# The held-out score, framing nothing, once the guard's reading is written:
# the prediction rows' outcomes unsealed, each fit's log score and Brier
# score, per dataset the mean difference move minus shipped with its Welch
# interval, on every R2 dataset.
runHeldout <- function() {
  if (!file.exists(runDir("readings", "guard.rds"))) {
    stop("the prediction rows stay sealed until the guard's reading is written")
  }
  a <- assignment()
  rows <- lapply(a$r2, function(name) {
    y <- readRDS(sealedPath("heldout", paste0(fileTag(name), ".rds")))$y
    nFit <- nrow(loadDesign(guardCell(name), kmpLink, 1L)$x)
    sc <- lapply(c("shipped", "move"), function(w) {
      fits <- Filter(Negate(is.null), loadGuardFits(name, w, guardFitsRaised))
      t(vapply(
        fits,
        function(f) {
          p <- pmin(pmax(f$est["prob", -seq_len(nFit)], 1e-12), 1 - 1e-12)
          c(
            log = mean(y * log(p) + (1 - y) * log(1 - p)),
            brier = mean((p - y)^2)
          )
        },
        numeric(2L)
      ))
    })
    wl <- welch(sc[[2L]][, "log"], sc[[1L]][, "log"])
    wb <- welch(sc[[2L]][, "brier"], sc[[1L]][, "brier"])
    data.frame(
      name = name,
      logDiff = wl[["D"]],
      logLo = wl[["lo"]],
      logHi = wl[["hi"]],
      brierDiff = wb[["D"]],
      brierLo = wb[["lo"]],
      brierHi = wb[["hi"]]
    )
  })
  tab <- do.call(rbind, rows)
  saveAtomic(
    list(table = tab, quick = quick),
    runDir("readings", "heldout.rds")
  )
  print(tab, row.names = FALSE, digits = 3)
}

## ------------------------------------------------------------------ step 3

# The trigger's features of every step 2 replication, from its data alone:
# the separation score on the fit's training rows and the smaller class's
# count.
triggerFeatures <- function(cell, rep, design = NULL) {
  if (is.null(design)) {
    design <- blind(loadDesign(cell, kmpLink, rep))
  }
  fd <- fitData(design)
  train <- fd$train
  y <- train$y
  train$y <- NULL
  if (design$task == "loop") {
    seen <- setdiff(seq_along(y), design$mask)
    train <- train[seen, , drop = FALSE]
    y <- y[seen]
  }
  c(
    score = separationScore(train, y, paste("trigger", cell, rep, seedTag)),
    minority = min(sum(y == 1L), sum(y == 0L))
  )
}

runStep3Features <- function() {
  rows <- list()
  for (id in step2Regimes) {
    n <- lookCount(id, 2L)
    for (rep in seq_len(n)) {
      if (!file.exists(step2ScorePath(id, rep))) {
        next
      }
      f <- triggerFeatures(id, rep)
      rows[[length(rows) + 1L]] <- data.frame(
        regime = id,
        rep = rep,
        score = f[["score"]],
        minority = f[["minority"]]
      )
    }
  }
  tab <- do.call(rbind, rows)
  saveAtomic(list(table = tab, quick = quick), runDir("step3", "features.rds"))
  cat("step 3 features:", nrow(tab), "replications\n")
}

triggerCandidates <- function() {
  ts <- seq(0.70, 0.95, by = 0.01)
  es <- c(10L, 20L, 40L, 80L)
  single <- c(
    lapply(ts, function(t) list(kind = "score", t = t, e = NA)),
    lapply(es, function(e) list(kind = "count", t = NA, e = e))
  )
  union <- unlist(
    lapply(ts, function(t) {
      lapply(es, function(e) list(kind = "union", t = t, e = e))
    }),
    recursive = FALSE
  )
  c(single, union)
}

fires <- function(cand, f) {
  s <- !is.na(cand$t) & f$score >= cand$t
  e <- !is.na(cand$e) & f$minority <= cand$e
  switch(cand$kind, score = s, count = e, union = s | e)
}

# Pooled gain over a set of replications: d(shipped at equal time) minus
# d(move) on a row class, each from the coverage pooled over the set; its
# one-sided 0.025 upper bound by a bootstrap over the set's replications.
pooledGainUpper <- function(sc, cl, key) {
  if (length(sc) == 0L) {
    return(-Inf)
  }
  get1 <- function(arm, field) {
    vapply(
      sc,
      function(s) {
        v <- s$coverage[[arm]][[cl]]
        if (is.null(v)) NA_real_ else v[[field]]
      },
      0
    )
  }
  n <- get1("move.default", "n")
  cm <- get1("move.default", "cov90")
  cs <- get1("shipped.eq", "cov90")
  ok <- is.finite(n) & n > 0
  if (!any(ok)) {
    return(-Inf)
  }
  g <- function(i) {
    abs(sum(cs[i]) / sum(n[i]) - 0.9) - abs(sum(cm[i]) / sum(n[i]) - 0.9)
  }
  w <- which(ok)
  b <- withSeed(seedOf("trigger-boot", key), {
    replicate(bootstrapResamples, g(w[sample.int(length(w), replace = TRUE)]))
  })
  stats::quantile(b, 0.975, names = FALSE)
}

# The mechanical selection: among candidates whose spared set's pooled-gain
# bound is below 0.02 in every A cell, the one firing on the fewest of B's
# replications; ties by the single feature over the union, then firing on
# more of A, then a score candidate over a count one, then the smaller t,
# then the larger e. Written with the rest of step 3 as its commit.
runStep3Select <- function(r3Size = NA_integer_) {
  feats <- readChecked(runDir("step3", "features.rds"))$table
  scores <- lapply(seq_len(nrow(feats)), function(i) {
    readChecked(step2ScorePath(feats$regime[i], feats$rep[i]))
  })
  aCells <- list(
    S1 = c("fInt", "fTail"),
    S2 = c("fInt", "fTail"),
    S3 = c("fInt", "fTail"),
    S4 = c("fInt", "fTail")
  )
  rows <- list()
  cands <- triggerCandidates()
  for (ci in seq_along(cands)) {
    cand <- cands[[ci]]
    fire <- fires(cand, feats)
    ok <- TRUE
    for (id in names(aCells)) {
      spared <- which(feats$regime == id & !fire)
      for (cl in aCells[[id]]) {
        u <- pooledGainUpper(scores[spared], cl, paste(ci, id, cl))
        if (u >= step2Bar) {
          ok <- FALSE
        }
      }
    }
    rows[[ci]] <- data.frame(
      ci = ci,
      kind = cand$kind,
      t = cand$t,
      e = cand$e,
      ok = ok,
      firesB = sum(fire & feats$regime %in% setB),
      firesA = sum(fire & feats$regime %in% setA)
    )
  }
  tab <- do.call(rbind, rows)
  okTab <- tab[tab$ok, , drop = FALSE]
  if (nrow(okTab) == 0L) {
    trigger <- NULL
  } else {
    o <- order(
      okTab$firesB,
      okTab$kind == "union",
      -okTab$firesA,
      okTab$kind != "score",
      ifelse(is.na(okTab$t), Inf, okTab$t),
      -ifelse(is.na(okTab$e), -Inf, okTab$e)
    )
    trigger <- cands[[okTab$ci[o[1L]]]]
  }
  commit <- list(
    trigger = trigger,
    candidates = tab,
    needRule = "accepted reference; 20 shipped fits at equal time hold a reference mass on a row class at least 0.02 from 0.90 (t bound excluding 0) while 20 move fits are within 0.02; not needing when the shipped fits' bound is within 0.02 on every row class",
    refutationRule = "a confirmed miss, or the exact one-sided 0.025 upper bound on specificity below one half",
    r3Size = r3Size,
    quick = quick
  )
  saveAtomic(commit, runDir("step3", "commit.rds"))
  str(trigger)
  invisible(commit)
}

# A real dataset's need of the move (step 3): from the reference mass inside
# each of the first 20 fits' 90 percent intervals on the row classes.
needsMove <- function(reading) {
  if (!isTRUE(reading$accepted)) {
    return("unresolved")
  }
  n20 <- function(m) m[seq_len(min(20L, nrow(m))), , drop = FALSE]
  ms <- n20(reading$beside$shipped$mass)
  mm <- n20(reading$beside$move$mass)
  verdict <- vapply(
    c("interior", "nearEnd"),
    function(cl) {
      a <- ms[, cl] - 0.9
      b <- mm[, cl] - 0.9
      a <- a[is.finite(a)]
      b <- b[is.finite(b)]
      if (length(a) < 2L || length(b) < 2L) {
        return("none")
      }
      ta <- stats::t.test(a)$conf.int
      tb <- stats::t.test(b)$conf.int
      shippedOff <- abs(mean(a)) >= step2Bar && (ta[1L] > 0 || ta[2L] < 0)
      moveWithin <- tb[1L] > -step2Bar && tb[2L] < step2Bar
      shippedWithin <- ta[1L] > -step2Bar && ta[2L] < step2Bar
      if (shippedOff && moveWithin) {
        "needs"
      } else if (shippedWithin) {
        "not"
      } else {
        "neither"
      }
    },
    ""
  )
  verdict <- verdict[verdict != "none"]
  if (any(verdict == "needs")) {
    "needs"
  } else if (length(verdict) && all(verdict == "not")) {
    "not"
  } else {
    "unresolved"
  }
}

clopper <- function(x, n) {
  if (n == 0L) {
    return(c(NA_real_, NA_real_))
  }
  c(
    if (x == 0L) 0 else stats::qbeta(0.025, x, n - x + 1),
    if (x == n) 1 else stats::qbeta(0.975, x + 1, n - x)
  )
}

# Step 3's reading: the trigger on each dataset's training rows, the need of
# each, the refutation rule, and sensitivity and specificity with exact
# intervals at every count; R2's datasets join when the commit predates any R2
# reading (the guard reading saw this commit file, by its md5). The trigger's
# separation score is each passer's sealed one.
runStep3Read <- function() {
  commit <- readChecked(runDir("step3", "commit.rds"))
  a <- assignment()
  # R2 joins when the guard was read under this commit (the md5 of the
  # commit file it saw), so the commit predates every R2 reading
  gf <- runDir("readings", "guard.rds")
  guardOpened <- file.exists(gf) &&
    !identical(readChecked(gf)$step3CommitMd5, commitMd5())
  sealed <- readChecked(sealedPath("scores.rds"))$score
  names3 <- c(a$r3, if (!guardOpened) a$r2)
  rows <- lapply(names3, function(name) {
    cell <- guardCell(name)
    if (!file.exists(acceptanceFile(cell, kmpLink, 1L))) {
      return(NULL)
    }
    d <- blind(loadDesign(cell, kmpLink, 1L))
    # the sealed separation score, the trigger's own feature
    f <- c(
      score = unname(sealed[name]),
      minority = min(sum(d$y), sum(1L - d$y))
    )
    rd <- guardDatasetReading(name)
    data.frame(
      name = name,
      set = if (name %in% a$r2) "R2" else "R3",
      fires = if (is.null(commit$trigger)) {
        NA
      } else {
        fires(commit$trigger, as.list(f))
      },
      need = needsMove(rd)
    )
  })
  tab <- do.call(rbind, Filter(Negate(is.null), rows))
  needing <- tab[tab$need == "needs", , drop = FALSE]
  notNeeding <- tab[tab$need == "not", , drop = FALSE]
  sens <- c(x = sum(needing$fires), n = nrow(needing))
  spec <- c(x = sum(!notNeeding$fires), n = nrow(notNeeding))
  specCi <- clopper(spec[["x"]], spec[["n"]])
  missed <- needing$name[!needing$fires]
  refuted <- length(missed) > 0L || (!is.na(specCi[2L]) && specCi[2L] < 0.5)
  out <- list(
    table = tab,
    sensitivity = sens,
    sensitivityCi = clopper(sens[["x"]], sens[["n"]]),
    specificity = spec,
    specificityCi = specCi,
    missesToConfirm = missed,
    refuted = refuted,
    r2Joined = !guardOpened,
    quick = quick
  )
  saveAtomic(out, runDir("readings", "step3.rds"))
  print(tab, row.names = FALSE)
  cat("step 3: misses", length(missed), "refuted", refuted, "\n")
  invisible(out)
}

# The trigger's sensitivity on fresh replications of A's cells and its
# specificity on B's, drawn after the commit (a run directory of their own,
# its step 2 scores and features made there; KMP_STEP3_COMMIT names the
# committed trigger's file).
runStep3Fresh <- function() {
  commit <- readRDS(Sys.getenv("KMP_STEP3_COMMIT"))
  feats <- readChecked(runDir("step3", "features.rds"))$table
  fire <- fires(commit$trigger, feats)
  inA <- feats$regime %in% setA
  out <- list(
    sensitivity = c(x = sum(fire[inA]), n = sum(inA)),
    sensitivityCi = clopper(sum(fire[inA]), sum(inA)),
    specificity = c(x = sum(!fire[!inA]), n = sum(!inA)),
    specificityCi = clopper(sum(!fire[!inA]), sum(!inA)),
    quick = quick
  )
  saveAtomic(out, runDir("readings", "step3-fresh.rds"))
  str(out)
}
