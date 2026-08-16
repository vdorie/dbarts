#!/usr/bin/env Rscript

# Permanent exact-posterior gate for the change move (docs/design/
# change-move-balance.md). The move redraws a node's split variable from the
# prior and a new rule over the descendant-valid good set; its acceptance now
# carries the proposal-density ratio, composed per side (an ordinal side
# contributes its counted good-set/interval ratio; a categorical side proposes
# from the node prior, whose density cancels outright), so it satisfies
# detailed balance. The omitted ratio was a CGM-lineage defect: for a
# cross-variable change at the root of a single-split "stump" it is exactly
# a_{v'}/a_v (a_v = number of available cuts of the split variable), which left
# the stationary distribution carrying the rule prior SQUARED, ~ (p_var / a_v)^2
# instead of p_var / a_v, biasing the root toward low-cardinality variables when
# cut counts are unequal (it cancels when they are equal).
#
# Three scenarios: MAIN (unequal ordinal cut counts, the original defect),
# CONTROL (equal ordinal cut counts, must stay clean), and MIXED (ordinal vs
# categorical: cross-type changes exercise BOTH sides of the correction, so a
# kernel whose correction branches on only one variable's type fails it).
#
# Three arms, all sharing one integrated-likelihood / tree-prior machinery:
#   engine       : a long single-tree dbarts run, root-split frequencies;
#   exact        : the true single-tree posterior over trees (region DP that
#                  sums the CGM prior x integrated likelihood over every
#                  subtree, empty-leaf splits excluded exactly as the engine's
#                  -1e7 branch veto excludes them);
#   wrong-target : the same DP with the rule-selection term (p_var / a_v)
#                  SQUARED at each internal node - the predicted defective
#                  stationary distribution. This is a leading-order (stump-
#                  exact) approximation; birth/death still target the correct
#                  pi, so the realized engine bias should sit BETWEEN exact and
#                  wrong-target.
# CONTROL: x2 given the same cut count as x1 (equal cardinality) - the bias
# cancels and the engine must match the exact posterior, validating the DP.
#
# Calibration matched to the engine exactly (else the comparison is void):
#   * useQuantiles cuts: x1 with 20 distinct values -> 19 cuts, x2 with 3 ->
#     2 cuts (verified: inducedNumCuts = numUnique - 1, data.hpp).
#   * range scaling: yScaled = (y - min y)/(max y - min y) - 0.5, offset NULL.
#   * fixed(1) residual variance -> internal sigma = 1 / range, known exactly.
#   * constant-leaf conjugate marginal with priorPrecision = (k/scale)^2,
#     scale = node.scale / sqrt(ntree) = 0.5, k = 2 (model.hpp
#     ConstantGaussianLeaf::logIntegratedLikelihood, reproduced verbatim).
#
# Usage: Rscript change-balance.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

# three engine arms x two scenarios share this budget; 300k keeps the variant
# arms' |z| < 4 a tight test while the current arm's failure stays z >> 10
nKept <- if (quick) 100000L else 300000L
batchSize <- if (quick) 25000L else 50000L
nThin <- 20L
nBurn <- 2000L
engineSeed <- 20260708L
maxDepthMain <- 6L
maxDepthCtrl <- 6L

base <- 0.95
power <- 2.0
kLeaf <- 2.0
nodeScale <- 0.5 # gaussian node.scale default; scale = nodeScale / sqrt(ntree=1)

logSumExp <- function(v) {
  m <- max(v)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(v - m)))
}

# ---- integrated likelihood, verbatim from ConstantGaussianLeaf ----

makeLogIL <- function(residVar) {
  priorPrecision <- (kLeaf / nodeScale)^2
  function(n, S, SS) {
    if (n == 0) {
      return(0)
    }
    posteriorPrecision <- n / residVar
    mean <- S / n
    centeredSumOfSquares <- SS - S * mean
    0.5 *
      log(priorPrecision / (priorPrecision + posteriorPrecision)) -
      0.5 * centeredSumOfSquares / residVar -
      0.5 *
        ((priorPrecision * mean) * (posteriorPrecision * mean)) /
        (priorPrecision + posteriorPrecision)
  }
}

# ---- exact / wrong-target posterior over single trees (region DP) ----
#
# A node's reachable set is a rectangle in (x1 sorted-value index, x2 index).
# splitInterval is ancestor-only, so region [a1,b1]x[a2,b2] has a_x1 = b1 - a1
# available cuts (b1 - a1 >= 1 to be available), likewise x2. Splitting x1 at
# the cut between values c and c+1 (a1 <= c <= b1 - 1) sends [a1,c] left and
# [c+1,b1] right; empty children are skipped (engine's -1e7 leaf veto). M(.)
# is the log of the summed prior x likelihood over every subtree rooted in the
# region at that depth; ruleMult = 2 squares the rule-selection term for the
# wrong-target arm. The root layer is kept split-by-split for the marginals.

buildPosterior <- function(
  vals1,
  vals2,
  cellN,
  cellS,
  cellSS,
  residVar,
  maxDepth,
  computeWrong = TRUE
) {
  K1 <- length(vals1)
  K2 <- length(vals2)
  cuts1 <- if (K1 >= 2L) (vals1[-K1] + vals1[-1L]) / 2 else numeric(0)
  cuts2 <- if (K2 >= 2L) (vals2[-K2] + vals2[-1L]) / 2 else numeric(0)
  logIL <- makeLogIL(residVar)

  # 2-D prefix sums for O(1) rectangle sufficient statistics
  pfx <- function(m) {
    p <- matrix(0, K1 + 1L, K2 + 1L)
    p[-1L, -1L] <- m
    p <- t(apply(p, 1L, cumsum))
    apply(p, 2L, cumsum)
  }
  PN <- pfx(cellN)
  PS <- pfx(cellS)
  PSS <- pfx(cellSS)
  query <- function(a1, b1, a2, b2) {
    idx <- function(P) {
      P[b1 + 1L, b2 + 1L] - P[a1, b2 + 1L] - P[b1 + 1L, a2] + P[a1, a2]
    }
    list(n = idx(PN), S = idx(PS), SS = idx(PSS))
  }

  runArm <- function(ruleMult) {
    memo <- new.env(parent = emptyenv())
    M <- function(a1, b1, a2, b2, depth) {
      key <- paste(a1, b1, a2, b2, depth, sep = ",")
      hit <- memo[[key]]
      if (!is.null(hit)) {
        return(hit)
      }
      st <- query(a1, b1, a2, b2)
      lLeaf <- logIL(st$n, st$S, st$SS)
      avail1 <- b1 > a1
      avail2 <- b2 > a2
      numAvail <- avail1 + avail2
      growth <- if (numAvail > 0L) base / (1 + depth)^power else 0
      terms <- log(1 - growth) + lLeaf # stop-here (leaf) contribution
      if (numAvail > 0L && depth < maxDepth) {
        logGrowth <- log(growth)
        if (avail1) {
          aV <- b1 - a1
          ruleLP <- ruleMult * (-log(numAvail) - log(aV))
          for (c in a1:(b1 - 1L)) {
            nL <- query(a1, c, a2, b2)$n
            nR <- query(c + 1L, b1, a2, b2)$n
            if (nL == 0 || nR == 0) {
              next
            }
            terms <- c(
              terms,
              logGrowth +
                ruleLP +
                M(a1, c, a2, b2, depth + 1L) +
                M(c + 1L, b1, a2, b2, depth + 1L)
            )
          }
        }
        if (avail2) {
          aV <- b2 - a2
          ruleLP <- ruleMult * (-log(numAvail) - log(aV))
          for (c in a2:(b2 - 1L)) {
            nL <- query(a1, b1, a2, c)$n
            nR <- query(a1, b1, c + 1L, b2)$n
            if (nL == 0 || nR == 0) {
              next
            }
            terms <- c(
              terms,
              logGrowth +
                ruleLP +
                M(a1, b1, a2, c, depth + 1L) +
                M(a1, b1, c + 1L, b2, depth + 1L)
            )
          }
        }
      }
      res <- logSumExp(terms)
      memo[[key]] <- res
      res
    }

    # root layer, kept split-by-split
    numAvail0 <- (K1 > 1L) + (K2 > 1L)
    growth0 <- base / (1 + 0)^power
    stFull <- query(1L, K1, 1L, K2)
    labels <- "root-only"
    variable <- 0L
    cutIndex <- NA_integer_
    logNum <- log(1 - growth0) + logIL(stFull$n, stFull$S, stFull$SS)
    addSplit <- function(v, ci, ln) {
      labels[[length(labels) + 1L]] <<- paste0("x", v, ".c", ci)
      variable[[length(variable) + 1L]] <<- v
      cutIndex[[length(cutIndex) + 1L]] <<- ci
      logNum[[length(logNum) + 1L]] <<- ln
    }
    if (K1 > 1L) {
      ruleLP <- ruleMult * (-log(numAvail0) - log(K1 - 1L))
      for (c in 1:(K1 - 1L)) {
        addSplit(
          1L,
          c,
          log(growth0) +
            ruleLP +
            M(1L, c, 1L, K2, 1L) +
            M(c + 1L, K1, 1L, K2, 1L)
        )
      }
    }
    if (K2 > 1L) {
      ruleLP <- ruleMult * (-log(numAvail0) - log(K2 - 1L))
      for (c in 1:(K2 - 1L)) {
        addSplit(
          2L,
          c,
          log(growth0) +
            ruleLP +
            M(1L, K1, 1L, c, 1L) +
            M(1L, K1, c + 1L, K2, 1L)
        )
      }
    }
    logZ <- logSumExp(logNum)
    data.frame(
      label = labels,
      variable = variable,
      cutIndex = cutIndex,
      prob = exp(logNum - logZ),
      stringsAsFactors = FALSE
    )
  }

  out <- list(exact = runArm(1), cuts1 = cuts1, cuts2 = cuts2)
  if (computeWrong) {
    out$wrong <- runArm(2)
  }
  out
}

# marginals P(root-only), P(root=x1), P(root=x2) from a root table
rootMarginals <- function(tab) {
  c(
    rootOnly = sum(tab$prob[tab$variable == 0L]),
    x1 = sum(tab$prob[tab$variable == 1L]),
    x2 = sum(tab$prob[tab$variable == 2L])
  )
}

# ---- exact posterior, mixed ordinal x categorical (region DP) ----
#
# x1 ordinal with K1 values (K1 - 1 cuts), x2 categorical with L levels. A
# node's reachable set is an x1 value interval [a1,b1] times a level subset S
# (ancestor-filtered and occupancy-blind, as the engine's reachableCategories).
# Ordinal splits as in buildPosterior. A categorical split assigns directions:
# each of the 2^|S| - 2 nonempty proper subsets D of S is an equally likely
# rule under the node prior, sending D right and S \ D left (ordered - the
# mirrored assignment is a distinct rule, as in the engine's gauge). Empty
# children are skipped (the engine's -1e7 branch veto). Subsets are bitmasks
# over level indices; no wrong-target arm (the defective correction here mixes
# omitted terms with out-of-domain ones, so no clean closed form exists).

buildMixedPosterior <- function(
  K1,
  L,
  cellN,
  cellS,
  cellSS,
  residVar,
  maxDepth
) {
  logIL <- makeLogIL(residVar)
  fullMask <- 2L^L - 1L
  levelBits <- 2L^(seq_len(L) - 1L)
  maskCols <- lapply(0:fullMask, function(S) which(bitwAnd(S, levelBits) != 0L))
  properSubsets <- lapply(0:fullMask, function(S) {
    if (S < 3L) {
      return(integer(0))
    } # fewer than two levels: no split
    cand <- seq_len(S - 1L)
    cand[bitwAnd(cand, S) == cand]
  })

  query <- function(a1, b1, S) {
    cols <- maskCols[[S + 1L]]
    list(
      n = sum(cellN[a1:b1, cols]),
      S = sum(cellS[a1:b1, cols]),
      SS = sum(cellSS[a1:b1, cols])
    )
  }

  memo <- new.env(parent = emptyenv())
  M <- function(a1, b1, S, depth) {
    key <- paste(a1, b1, S, depth, sep = ",")
    hit <- memo[[key]]
    if (!is.null(hit)) {
      return(hit)
    }
    st <- query(a1, b1, S)
    sizeS <- length(maskCols[[S + 1L]])
    avail1 <- b1 > a1
    avail2 <- sizeS >= 2L
    numAvail <- avail1 + avail2
    growth <- if (numAvail > 0L) base / (1 + depth)^power else 0
    terms <- log(1 - growth) + logIL(st$n, st$S, st$SS)
    if (numAvail > 0L && depth < maxDepth) {
      logGrowth <- log(growth)
      if (avail1) {
        ruleLP <- -log(numAvail) - log(b1 - a1)
        for (c in a1:(b1 - 1L)) {
          if (query(a1, c, S)$n == 0 || query(c + 1L, b1, S)$n == 0) {
            next
          }
          terms <- c(
            terms,
            logGrowth +
              ruleLP +
              M(a1, c, S, depth + 1L) +
              M(c + 1L, b1, S, depth + 1L)
          )
        }
      }
      if (avail2) {
        ruleLP <- -log(numAvail) - log(2^sizeS - 2)
        for (D in properSubsets[[S + 1L]]) {
          left <- S - D # D is a subset of S, so mask difference is set difference
          if (query(a1, b1, left)$n == 0 || query(a1, b1, D)$n == 0) {
            next
          }
          terms <- c(
            terms,
            logGrowth +
              ruleLP +
              M(a1, b1, left, depth + 1L) +
              M(a1, b1, D, depth + 1L)
          )
        }
      }
    }
    res <- logSumExp(terms)
    memo[[key]] <- res
    res
  }

  # root layer, kept split-by-split
  numAvail0 <- (K1 > 1L) + (L > 1L)
  growth0 <- base / (1 + 0)^power
  stFull <- query(1L, K1, fullMask)
  labels <- "root-only"
  variable <- 0L
  logNum <- log(1 - growth0) + logIL(stFull$n, stFull$S, stFull$SS)
  addSplit <- function(v, lab, ln) {
    labels[[length(labels) + 1L]] <<- lab
    variable[[length(variable) + 1L]] <<- v
    logNum[[length(logNum) + 1L]] <<- ln
  }
  if (K1 > 1L) {
    ruleLP <- -log(numAvail0) - log(K1 - 1L)
    for (c in 1:(K1 - 1L)) {
      addSplit(
        1L,
        paste0("x1.c", c),
        log(growth0) +
          ruleLP +
          M(1L, c, fullMask, 1L) +
          M(c + 1L, K1, fullMask, 1L)
      )
    }
  }
  if (L > 1L) {
    ruleLP <- -log(numAvail0) - log(2^L - 2)
    for (D in properSubsets[[fullMask + 1L]]) {
      addSplit(
        2L,
        paste0("x2.d", D),
        log(growth0) + ruleLP + M(1L, K1, fullMask - D, 1L) + M(1L, K1, D, 1L)
      )
    }
  }
  logZ <- logSumExp(logNum)
  data.frame(
    label = labels,
    variable = variable,
    prob = exp(logNum - logZ),
    stringsAsFactors = FALSE
  )
}

# ---- engine arm ----

runEngine <- function(x, y, nCutsVec) {
  ctl <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 1L,
    useQuantiles = TRUE,
    keepTrees = TRUE,
    n.samples = batchSize,
    n.burn = nBurn,
    n.thin = nThin,
    updateState = TRUE,
    seed = engineSeed,
    n.cuts = nCutsVec
  )
  s <- dbarts(
    x,
    y,
    control = ctl,
    tree.prior = cgm(power, base),
    node.prior = normal(kLeaf),
    resid.prior = fixed(1)
  )
  stopifnot(is.null(s$data@offset))
  nBatch <- as.integer(ceiling(nKept / batchSize))
  var <- integer(0)
  val <- numeric(0)
  first <- TRUE
  for (b in seq_len(nBatch)) {
    if (first) {
      invisible(s$run(nBurn, batchSize))
      first <- FALSE
    } else {
      invisible(s$run(0L, batchSize))
    }
    tr <- s$getTrees()
    roots <- tr[!duplicated(tr$sample), ]
    var <- c(var, roots$var)
    val <- c(val, roots$value)
  }
  list(var = var, val = val)
}

# ---- batch-means MC error and z-scores ----

batchMeanSE <- function(indicator, nBatch = 400L) {
  N <- length(indicator)
  bl <- N %/% nBatch
  x <- indicator[seq_len(bl * nBatch)]
  bm <- colMeans(matrix(x, nrow = bl))
  sd(bm) / sqrt(nBatch)
}

# ---- one scenario end to end ----

runScenario <- function(name, x1, x2, y, maxDepth, nCutsVec) {
  vals1 <- sort(unique(x1))
  vals2 <- sort(unique(x2))
  K1 <- length(vals1)
  K2 <- length(vals2)
  yRange <- max(y) - min(y)
  yScaled <- (y - min(y)) / yRange - 0.5
  residVar <- (1 / yRange)^2

  i1 <- match(x1, vals1)
  i2 <- match(x2, vals2)
  cellN <- matrix(0, K1, K2)
  cellS <- matrix(0, K1, K2)
  cellSS <- matrix(0, K1, K2)
  for (o in seq_along(y)) {
    cellN[i1[o], i2[o]] <- cellN[i1[o], i2[o]] + 1
    cellS[i1[o], i2[o]] <- cellS[i1[o], i2[o]] + yScaled[o]
    cellSS[i1[o], i2[o]] <- cellSS[i1[o], i2[o]] + yScaled[o]^2
  }

  post <- buildPosterior(vals1, vals2, cellN, cellS, cellSS, residVar, maxDepth)
  eng <- runEngine(cbind(x1 = x1, x2 = x2), y, nCutsVec)

  exM <- rootMarginals(post$exact)
  wrM <- rootMarginals(post$wrong)
  N <- length(eng$var)
  enM <- c(
    rootOnly = mean(eng$var == -1L),
    x1 = mean(eng$var == 1L),
    x2 = mean(eng$var == 2L)
  )
  se <- c(
    rootOnly = batchMeanSE(eng$var == -1L),
    x1 = batchMeanSE(eng$var == 1L),
    x2 = batchMeanSE(eng$var == 2L)
  )

  cat(sprintf(
    "\n=== %s (n=%d, x1 %d cuts, x2 %d cuts, engine N=%d) ===\n",
    name,
    length(y),
    K1 - 1L,
    K2 - 1L,
    N
  ))
  cat(sprintf(
    "%-10s %10s %10s %10s %10s\n",
    "root",
    "engine",
    "exact",
    "wrong",
    "MCse"
  ))
  for (q in c("rootOnly", "x1", "x2")) {
    cat(sprintf(
      "%-10s %10.4f %10.4f %10.4f %10.5f\n",
      q,
      enM[q],
      exM[q],
      wrM[q],
      se[q]
    ))
  }
  zEx <- (enM - exM) / se
  zWr <- (enM - wrM) / se
  cat(sprintf(
    "z(uncond) vs exact : rootOnly %+.1f  x1 %+.1f  x2 %+.1f\n",
    zEx["rootOnly"],
    zEx["x1"],
    zEx["x2"]
  ))
  cat(sprintf(
    "z(uncond) vs wrong : rootOnly %+.1f  x1 %+.1f  x2 %+.1f\n",
    zWr["rootOnly"],
    zWr["x1"],
    zWr["x2"]
  ))

  # Conditional on the root being split: the change move governs only the
  # among-variables/among-cuts distribution, and for a stump the squared-rule
  # wrong-target prediction is exact - so here the three arms bracket cleanly.
  # (Unconditionally the wrong arm over-inflates root-only because squaring the
  # rule prior also squashes the split/no-split balance, which the real defect
  # leaves to the correct birth/death moves.)
  condX2 <- function(m) m["x2"] / (m["x1"] + m["x2"])
  enC <- condX2(enM)
  exC <- condX2(exM)
  wrC <- condX2(wrM)
  split <- eng$var != -1L
  seC <- batchMeanSE(eng$var[split] == 2L)
  zExC <- (enC - exC) / seC
  zWrC <- (enC - wrC) / seC
  frac <- (enC - exC) / (wrC - exC)
  cat(sprintf(
    paste0(
      "P(root=x2 | split): engine %.4f  exact %.4f  wrong %.4f",
      "  (MCse %.5f)\n"
    ),
    enC,
    exC,
    wrC,
    seC
  ))
  cat(sprintf(
    paste0(
      "  z vs exact %+.1f, z vs wrong %+.1f; engine sits %.0f%%",
      " of the way exact->wrong\n"
    ),
    zExC,
    zWrC,
    100 * frac
  ))

  # cut distribution within each root variable (engine vs exact vs wrong)
  cutReport <- function(v, cuts) {
    idx <- match(eng$val[eng$var == v], cuts)
    engCut <- tabulate(idx, nbins = length(cuts)) / sum(eng$var == v)
    exCut <- post$exact$prob[post$exact$variable == v]
    exCut <- exCut / sum(exCut)
    wrCut <- post$wrong$prob[post$wrong$variable == v]
    wrCut <- wrCut / sum(wrCut)
    cat(sprintf(
      "cut dist x%d (engine/exact/wrong) max|eng-exact|=%.4f%s\n",
      v,
      max(abs(engCut - exCut)),
      if (length(cuts) <= 3L) "" else " (per-cut arrays suppressed)"
    ))
    if (length(cuts) <= 3L) {
      cat("  engine:", sprintf("%.3f", engCut), "\n")
      cat("  exact :", sprintf("%.3f", exCut), "\n")
      cat("  wrong :", sprintf("%.3f", wrCut), "\n")
    }
  }
  cutReport(1L, post$cuts1)
  cutReport(2L, post$cuts2)

  list(
    engine = enM,
    exact = exM,
    wrong = wrM,
    se = se,
    zEx = zEx,
    zWr = zWr,
    condEngine = enC,
    condExact = exC,
    condWrong = wrC,
    condSe = seC,
    condZex = zExC,
    condZwr = zWrC
  )
}

# ---- truncation check ----

truncationCheck <- function(x1, x2, y, depths) {
  vals1 <- sort(unique(x1))
  vals2 <- sort(unique(x2))
  K1 <- length(vals1)
  K2 <- length(vals2)
  yRange <- max(y) - min(y)
  yScaled <- (y - min(y)) / yRange - 0.5
  residVar <- (1 / yRange)^2
  i1 <- match(x1, vals1)
  i2 <- match(x2, vals2)
  cellN <- matrix(0, K1, K2)
  cellS <- matrix(0, K1, K2)
  cellSS <- matrix(0, K1, K2)
  for (o in seq_along(y)) {
    cellN[i1[o], i2[o]] <- cellN[i1[o], i2[o]] + 1
    cellS[i1[o], i2[o]] <- cellS[i1[o], i2[o]] + yScaled[o]
    cellSS[i1[o], i2[o]] <- cellSS[i1[o], i2[o]] + yScaled[o]^2
  }
  prev <- NULL
  for (d in depths) {
    p <- buildPosterior(
      vals1,
      vals2,
      cellN,
      cellS,
      cellSS,
      residVar,
      d,
      computeWrong = FALSE
    )
    m <- rootMarginals(p$exact)
    if (!is.null(prev)) {
      cat(sprintf(
        "  depth %d -> %d: max root-marginal shift %.2e\n",
        prev$d,
        d,
        max(abs(m - prev$m))
      ))
    }
    prev <- list(d = d, m = m)
  }
}

# ---- the mixed scenario end to end ----

mixedCellStats <- function(x1, x2f, y) {
  vals1 <- sort(unique(x1))
  K1 <- length(vals1)
  L <- nlevels(x2f)
  yRange <- max(y) - min(y)
  yScaled <- (y - min(y)) / yRange - 0.5
  i1 <- match(x1, vals1)
  i2 <- as.integer(x2f)
  cellN <- matrix(0, K1, L)
  cellS <- matrix(0, K1, L)
  cellSS <- matrix(0, K1, L)
  for (o in seq_along(y)) {
    cellN[i1[o], i2[o]] <- cellN[i1[o], i2[o]] + 1
    cellS[i1[o], i2[o]] <- cellS[i1[o], i2[o]] + yScaled[o]
    cellSS[i1[o], i2[o]] <- cellSS[i1[o], i2[o]] + yScaled[o]^2
  }
  list(
    K1 = K1,
    L = L,
    cellN = cellN,
    cellS = cellS,
    cellSS = cellSS,
    residVar = (1 / yRange)^2
  )
}

runMixedScenario <- function(name, x1, x2f, y, maxDepth) {
  cs <- mixedCellStats(x1, x2f, y)
  post <- buildMixedPosterior(
    cs$K1,
    cs$L,
    cs$cellN,
    cs$cellS,
    cs$cellSS,
    cs$residVar,
    maxDepth
  )
  eng <- runEngine(data.frame(x1 = x1, x2 = x2f), y, c(100L, 100L))

  exM <- rootMarginals(post)
  N <- length(eng$var)
  enM <- c(
    rootOnly = mean(eng$var == -1L),
    x1 = mean(eng$var == 1L),
    x2 = mean(eng$var == 2L)
  )
  se <- c(
    rootOnly = batchMeanSE(eng$var == -1L),
    x1 = batchMeanSE(eng$var == 1L),
    x2 = batchMeanSE(eng$var == 2L)
  )

  cat(sprintf(
    "\n=== %s (n=%d, x1 %d cuts, x2 %d levels, engine N=%d) ===\n",
    name,
    length(y),
    cs$K1 - 1L,
    cs$L,
    N
  ))
  cat(sprintf("%-10s %10s %10s %10s\n", "root", "engine", "exact", "MCse"))
  for (q in c("rootOnly", "x1", "x2")) {
    cat(sprintf("%-10s %10.4f %10.4f %10.5f\n", q, enM[q], exM[q], se[q]))
  }
  zEx <- (enM - exM) / se
  cat(sprintf(
    "z(uncond) vs exact : rootOnly %+.1f  x1 %+.1f  x2 %+.1f\n",
    zEx["rootOnly"],
    zEx["x1"],
    zEx["x2"]
  ))

  condX2 <- function(m) m["x2"] / (m["x1"] + m["x2"])
  enC <- condX2(enM)
  exC <- condX2(exM)
  split <- eng$var != -1L
  seC <- batchMeanSE(eng$var[split] == 2L)
  zExC <- (enC - exC) / seC
  cat(sprintf(
    "P(root=x2 | split): engine %.4f  exact %.4f  (MCse %.5f)  z %+.1f\n",
    enC,
    exC,
    seC,
    zExC
  ))

  list(
    engine = enM,
    exact = exM,
    se = se,
    zEx = zEx,
    condEngine = enC,
    condExact = exC,
    condSe = seC,
    condZex = zExC
  )
}

mixedTruncationCheck <- function(x1, x2f, y, depths) {
  cs <- mixedCellStats(x1, x2f, y)
  prev <- NULL
  for (d in depths) {
    p <- buildMixedPosterior(
      cs$K1,
      cs$L,
      cs$cellN,
      cs$cellS,
      cs$cellSS,
      cs$residVar,
      d
    )
    m <- rootMarginals(p)
    if (!is.null(prev)) {
      cat(sprintf(
        "  depth %d -> %d: max root-marginal shift %.2e\n",
        prev$d,
        d,
        max(abs(m - prev$m))
      ))
    }
    prev <- list(d = d, m = m)
  }
}

# ================= MAIN scenario =================

set.seed(101L)
n <- 100L
x1 <- as.double(rep(1:20, each = 5L)) # 20 distinct -> 19 cuts
x2 <- as.double(rep(1:3, length.out = n)) # 3 distinct -> 2 cuts
yMain <- 0.2 * as.vector(scale(x1)) + rnorm(n) # weak x1 signal, x2 noise

cat("truncation check (exact root marginals vs depth), MAIN:\n")
truncationCheck(x1, x2, yMain, c(4L, 5L, 6L))

resMain <- runScenario(
  "MAIN: unequal cut counts (19 vs 2)",
  x1,
  x2,
  yMain,
  maxDepthMain,
  c(100L, 100L)
)

# ================= CONTROL scenario (equal cut counts) =================

set.seed(202L)
x1c <- as.double(rep(1:20, each = 5L)) # 20 distinct -> 19 cuts
x2c <- as.double(rep(1:20, times = 5L)) # 20 distinct -> 19 cuts, independent
yCtrl <- 0.2 * as.vector(scale(x1c)) + rnorm(n) # signal on x1c only

cat("\ntruncation check (exact root marginals vs depth), CONTROL:\n")
truncationCheck(x1c, x2c, yCtrl, c(5L, 6L))

resCtrl <- runScenario(
  "CONTROL: equal cut counts (19 vs 19)",
  x1c,
  x2c,
  yCtrl,
  maxDepthCtrl,
  c(100L, 100L)
)

# ================= MIXED scenario (ordinal vs categorical) =================
# x2 is a 4-level factor exactly aligned with x1's step blocks, so both
# variables encode the same 4-level response and every mixed tree shape has a
# likelihood-neutral cross-type partner: change proposals swapping x1 cuts for
# x2 level subsets (and back) at STACKED nodes flow freely, and the acceptance
# there is pure proposal-density correction. A one-sided kernel mis-weights
# exactly those transitions - ord->cat omits the constrained old ordinal
# side's counted ratio (|Valid| < |SI| whenever a same-variable descendant
# split sits below), cat->ord feeds a categorical column through the ordinal
# counters - so the x1-vs-x2 root balance drifts from the exact posterior.

set.seed(303L)
nMix <- 120L
x1m <- as.double(rep(1:10, each = 12L)) # 10 distinct -> 9 cuts
blockMix <- c(0, 0, 1, 1, 1, 2, 2, 2, 3, 3) # steps at x1 cuts 2, 5, 8
x2m <- factor(letters[blockMix[x1m] + 1L]) # a..d, one level per block
yMix <- blockMix[x1m] + rnorm(nMix)

cat("\ntruncation check (exact root marginals vs depth), MIXED:\n")
mixedTruncationCheck(x1m, x2m, yMix, c(4L, 5L, 6L))

resMix <- runMixedScenario(
  "MIXED: ordinal (9 cuts) vs categorical (4 levels)",
  x1m,
  x2m,
  yMix,
  6L
)

# ================= verdict =================
# Gate: the engine must MATCH the exact posterior (|z| < 4) on MAIN and MIXED
# and sit far from the wrong (squared-rule-prior) target on MAIN; CONTROL
# (equal cut counts) must match the exact posterior too (|z| < 4) - the depth-6
# DP-truncation residual is ~1e-7 (see truncationCheck), far below what a z = 4
# gate demands, so the equal-cut arm is held to the same bound as MAIN/MIXED.

cat("\n================ VERDICT ================\n")
cat(sprintf(
  "%-26s %10s %10s %10s %10s %8s\n",
  "scenario",
  "P(x2|spl)",
  "exact",
  "wrong",
  "z-exact",
  "gate"
))
mainZ <- resMain$condZex
mainOk <- abs(mainZ) < 4
cat(sprintf(
  "%-26s %10.4f %10.4f %10.4f %+10.1f %8s\n",
  "MAIN (unequal 19 vs 2)",
  resMain$condEngine,
  resMain$condExact,
  resMain$condWrong,
  mainZ,
  if (mainOk) "PASS" else "FAIL"
))
cat(sprintf("  (engine z vs wrong target %+.1f)\n", resMain$condZwr))
ctrlZ <- resCtrl$condZex
ctrlOk <- abs(ctrlZ) < 4
cat(sprintf(
  "%-26s %10.4f %10.4f %10.4f %+10.1f %8s\n",
  "CONTROL (equal 19 vs 19)",
  resCtrl$condEngine,
  resCtrl$condExact,
  resCtrl$condWrong,
  ctrlZ,
  if (ctrlOk) "match" else "MISMATCH"
))
mixZ <- resMix$condZex
mixOk <- abs(mixZ) < 4
cat(sprintf(
  "%-26s %10.4f %10.4f %10s %+10.1f %8s\n",
  "MIXED (ordinal vs categ.)",
  resMix$condEngine,
  resMix$condExact,
  "-",
  mixZ,
  if (mixOk) "PASS" else "FAIL"
))

cat(sprintf(
  "\nBALANCE GATE: %s (MAIN matches exact: %s; CONTROL clean: %s; MIXED matches exact: %s)\n",
  if (mainOk && ctrlOk && mixOk) "PASS" else "FAIL",
  if (mainOk) "yes" else "NO",
  if (ctrlOk) "yes" else "NO",
  if (mixOk) "yes" else "NO"
))

# Exit nonzero on failure so a runner/CI can catch it (matches bd-balance.R).
if (!(mainOk && ctrlOk && mixOk)) {
  quit(status = 1L, save = "no")
}
