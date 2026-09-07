#!/usr/bin/env Rscript

# Exact-posterior gate for the bartcore BCF two-forest sampler under the
# LATENT binary families (docs/design/bcf.md). The gaussian twin of this
# gate is bcf-exact.R; four things are different under probit and logistic
# and they are what this one covers. Sigma is pinned at exactly 1 with no
# draw, so the sigma quadrature disappears. The leaf-scale map anchors on
# the link's latent scale (1 under probit, pi/sqrt(3) under logistic) at a
# half-Cauchy scale of 1, not gaussian's 2. The latent refresh runs against
# the COMBINED location a mu + b_z tau rather than either forest's fits.
# And the glue is drawn under both links, so all three glue modes apply.
#
# The design is one ordinal predictor with two cells, n.cuts = 1 and z
# balanced within each cell, so each single-tree forest realizes exactly two
# trees (after the one split each child holds a single cell with an empty
# cut interval, where the CGM growth probability is zero) and the joint
# space is four configurations. Conditional on a configuration the leaf
# parameters couple only through the cells, so the leaf integral factorizes
# over the bipartite mu-leaf/tau-leaf graph into blocks of dimension at most
# three; each block is integrated by Newton to the mode plus an adaptive
# Gauss-Hermite product rule. The glue axes take a tangent substitution,
# under which the Cauchy measure on a is uniform in t; a truncated tensor
# grid is not accurate enough here, and its error is a systematic offset the
# gate could not see.
#
# Modes: (1) glue fixed a = 1, b0 = 0, b1 = 1 - match E[mu], E[tau];
# (2a) a free (Cauchy(0, sd.control)), b fixed - match E[a mu], E[tau];
# (2b) b0/b1 free (N(0, bPriorVariance)), a fixed - match E[mu],
# E[(b1-b0)tau]. Every mode also matches E[F(eta)] at the four (cell, z)
# groups, the reported probability surface, which is ridge-invariant.
# Two arms extend the coverage: an AGGREGATED logistic arm (two weighted
# rows per group in place of 400 unit rows) gates the frequency-weight
# channel probit cannot reach, and a K = 3, mode-1 arm at each link
# restores the depth-decay and cut-selection coverage a two-cell design is
# blind to (the only interior node is the root, where base / (1 + depth)^
# power equals base for every power).
#
# Gating statistic: a per-channel batch-means z at |z| <= 4, not an
# absolute tolerance - the binary leaf posteriors are 0.13 to 0.31 wide by
# channel against the gaussian gate's 0.0034, so an absolute bound carries
# no fixed meaning across the two. Standard errors are formed per seed and
# pooled as sqrt(sum se_s^2) / S, never batched across a chain seam, and the
# batch LENGTH is raised until the batch means decorrelate: several channels
# here mix far too slowly for a fixed batch count to report their own error
# honestly. About 80 tests at |z| <= 4 is a family-wise false-failure rate
# near 5e-3, accepted in advance. Failure means the backfit, the glue draw,
# the latent refresh or the calibration map is wrong.
#
# Usage: Rscript bcf-latent-exact.R [quick]

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

zBound <- 4
nBurn <- if (quick) 5000L else 20000L
nKept <- if (quick) 12500L else 100000L
nThin <- if (quick) 5L else 10L
nSeeds <- 3L
# mode 2b's E[mu] couples to the slowly mixing (b0 + b1) / (b1 - b0) ratio,
# the bscale ridge the gaussian gate met with thin 200 and eight seeds; this
# gate does the same, at the kept-draw count that holds the per-seed SWEEP
# budget equal to the other modes' (100000 kept draws at thin 200 would be
# 2e7 sweeps a seed and eight seeds a link).
nKept2b <- if (quick) 750L else 5000L
nThin2b <- 200L
nSeeds2b <- if (quick) 4L else 8L
# glue quadrature: an open trapezoid in t under a = sd.control tan(t) and
# b = sqrt(bPriorVariance) tan(t), plus the refinement the oracle refuses to
# run without
aN <- 101L
aNFine <- 401L
bN <- 81L
bNFine <- 121L
refineTolerance <- 1e-6
# the design keeps more than one configuration alive; the guard is stated
# from below on the SECOND largest weight at the realized counts, because a
# cap on the largest would refuse a third of realizations of this design
runnerUpFloor <- 0.02
ghNodes <- 12L

# ---- design ----

dataSeed <- 20260907L
nPerCell <- 200L # even, so z is balanced within each cell
muBase <- 0.95
muPower <- 2
tauBase <- 0.25
tauPower <- 3
sdControl <- 1 # the latent families' defaultAmplitudePriorScale
sdModerate <- 1
bVar <- 0.5
nodeScaleDivisor <- 0.674 # the half-normal median of b1 - b0 ~ N(0, 1)

# leaf values, one per cell, held small enough that every group probability
# is well inside the unit interval under both links (which is also what
# keeps 0 < s_g < n_g at every group for the aggregated arm)
muValue2 <- c(0.0, 0.3)
tauValue2 <- c(0.5, 0.8)
muValue3 <- c(-0.3, 0.0, 0.3)
tauValue3 <- c(0.4, 0.6, 0.8)

# ---- links ----

# Each link supplies the latent-scale anchor the calibration map reads, the
# log link probabilities the binomial integrand sums, and the first two
# derivatives of that integrand in the index (used by Newton alone; the
# quadrature needs only the log probabilities).
probitLink <- list(
  name = "probit",
  anchor = 1,
  logs = function(eta) {
    list(lp = pnorm(eta, log.p = TRUE), lq = pnorm(-eta, log.p = TRUE))
  },
  derivs = function(eta, s, m) {
    logDensity <- dnorm(eta, log = TRUE)
    r1 <- exp(logDensity - pnorm(eta, log.p = TRUE))
    r0 <- exp(logDensity - pnorm(-eta, log.p = TRUE))
    list(
      gradient = s * r1 - (m - s) * r0,
      curvature = -s * (eta * r1 + r1^2) + (m - s) * (eta * r0 - r0^2)
    )
  },
  linkinv = pnorm
)
logisticLink <- list(
  name = "logistic",
  anchor = pi / sqrt(3), # the logistic law's sd, latentScaleAnchor's value
  logs = function(eta) {
    list(lp = plogis(eta, log.p = TRUE), lq = plogis(-eta, log.p = TRUE))
  },
  derivs = function(eta, s, m) {
    p <- plogis(eta)
    list(gradient = s - m * p, curvature = -m * p * (1 - p))
  },
  linkinv = plogis
)

# ---- Gauss-Hermite rules ----

# Golub-Welsch on the physicists' Hermite recurrence: weight exp(-y^2).
ghRule <- function(n) {
  i <- seq_len(n - 1L)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(i, i + 1L)] <- sqrt(i / 2)
  jacobi[cbind(i + 1L, i)] <- sqrt(i / 2)
  eg <- eigen(jacobi, symmetric = TRUE)
  ord <- order(eg$values)
  list(node = eg$values[ord], weight = sqrt(pi) * eg$vectors[1L, ord]^2)
}

# The d-fold product rule, with the exp(+|y|^2) that undoes the Gaussian
# weight folded into the log weights once.
ghProduct <- function(rule, d) {
  nodes <- as.matrix(expand.grid(rep(list(rule$node), d)))
  logWeight <- rowSums(
    as.matrix(expand.grid(rep(list(log(rule$weight)), d)))
  )
  list(node = nodes, logWeight = logWeight + rowSums(nodes^2))
}

ghCache <- local({
  rule <- ghRule(ghNodes)
  lapply(1:4, function(d) ghProduct(rule, d))
})

# ---- tree enumeration and the bipartite block decomposition ----

# CGM prior over trees on cut interval [loCut, hiCut]: a leaf with
# probability 1 - base / (1 + depth)^power, forced once the interval empties
# (the growth probability there is exactly zero, so the forced leaf costs no
# prior mass), else a uniformly drawn cut splits the cell range. Leaves are
# contiguous cell ranges.
enumerate <- function(loCell, hiCell, loCut, hiCut, depth, base, power) {
  growth <- if (hiCut >= loCut) base / (1 + depth)^power else 0
  result <- list(list(
    leaves = list(c(loCell, hiCell)),
    logPrior = log(1 - growth)
  ))
  if (hiCut < loCut) {
    return(result)
  }
  for (j in loCut:hiCut) {
    lefts <- enumerate(loCell, j, loCut, j - 1L, depth + 1L, base, power)
    rights <- enumerate(j + 1L, hiCell, j + 1L, hiCut, depth + 1L, base, power)
    for (left in lefts) {
      for (right in rights) {
        result[[length(result) + 1L]] <- list(
          leaves = c(left$leaves, right$leaves),
          logPrior = log(growth) -
            log(hiCut - loCut + 1) +
            left$logPrior +
            right$logPrior
        )
      }
    }
  }
  result
}

leafAssignment <- function(tree, K) {
  leafOf <- integer(K)
  for (li in seq_along(tree$leaves)) {
    rng <- tree$leaves[[li]]
    leafOf[rng[1L]:rng[2L]] <- li
  }
  leafOf
}

# Connected components of the bipartite graph whose vertices are the mu- and
# tau-leaves and whose edges are the cells: the likelihood factorizes over
# them, so the leaf integral is a product of block integrals.
blockPartition <- function(muLeaf, tauLeaf, K) {
  pMu <- max(muLeaf)
  parent <- seq_len(pMu + max(tauLeaf))
  root <- function(i) {
    while (parent[i] != i) {
      i <- parent[i]
    }
    i
  }
  for (c in seq_len(K)) {
    rootMu <- root(muLeaf[c])
    rootTau <- root(pMu + tauLeaf[c])
    if (rootMu != rootTau) {
      parent[rootMu] <- rootTau
    }
  }
  label <- vapply(seq_along(parent), root, integer(1L))
  lapply(unique(label), function(r) {
    list(params = which(label == r), cells = which(label[muLeaf] == r))
  })
}

# Everything about a block that the glue does not move: which parameters it
# holds, which (cell, z) groups enter its likelihood, their sufficient
# statistics, and the two design pieces the glue scales (the index is
# a * mu_l(c) + b_z * tau_m(c), linear in the block's parameters).
blockDesign <- function(block, muLeaf, tauLeaf, pMu, sCz, nCz) {
  params <- block$params
  d <- length(params)
  nGroup <- 2L * length(block$cells)
  designMu <- matrix(0, nGroup, d)
  designTau0 <- matrix(0, nGroup, d)
  designTau1 <- matrix(0, nGroup, d)
  cellOf <- integer(nGroup)
  zOf <- integer(nGroup)
  gi <- 0L
  for (c in block$cells) {
    for (zg in 1:2) {
      gi <- gi + 1L
      designMu[gi, match(muLeaf[c], params)] <- 1
      column <- match(pMu + tauLeaf[c], params)
      if (zg == 1L) {
        designTau0[gi, column] <- 1
      } else {
        designTau1[gi, column] <- 1
      }
      cellOf[gi] <- c
      zOf[gi] <- zg
    }
  }
  index <- cbind(cellOf, zOf)
  list(
    params = params,
    d = d,
    isMu = params <= pMu,
    designMu = designMu,
    designTau0 = designTau0,
    designTau1 = designTau1,
    cellOf = cellOf,
    zOf = zOf,
    successes = sCz[index],
    trials = nCz[index]
  )
}

# The joint tree space: every (mu tree, tau tree) pair, its CGM prior mass
# and its block decomposition.
makeConfigurations <- function(K, sCz, nCz) {
  nCuts <- K - 1L
  muTrees <- enumerate(1L, K, 1L, nCuts, 0L, muBase, muPower)
  tauTrees <- enumerate(1L, K, 1L, nCuts, 0L, tauBase, tauPower)
  configurations <- list()
  for (mt in seq_along(muTrees)) {
    for (tt in seq_along(tauTrees)) {
      muLeaf <- leafAssignment(muTrees[[mt]], K)
      tauLeaf <- leafAssignment(tauTrees[[tt]], K)
      pMu <- max(muLeaf)
      blocks <- blockPartition(muLeaf, tauLeaf, K)
      configurations[[length(configurations) + 1L]] <- list(
        muLeaf = muLeaf,
        tauLeaf = tauLeaf,
        pMu = pMu,
        logPrior = muTrees[[mt]]$logPrior + tauTrees[[tt]]$logPrior,
        design = lapply(
          blocks,
          blockDesign,
          muLeaf = muLeaf,
          tauLeaf = tauLeaf,
          pMu = pMu,
          sCz = sCz,
          nCz = nCz
        )
      )
    }
  }
  configurations
}

# ---- one block's marginal and posterior means ----

# The log integrand is strictly concave (log F and log(1 - F) are concave in
# a linear index, the leaf prior strictly so), so Newton from the origin
# reaches the mode; the quadrature is then a product Gauss-Hermite rule in
# the mode's own curvature metric. At 12 nodes per dimension the dim-3 block
# moves 2.6e-13 against a 200-node reference, so the leaf rule is not what
# the refinement check below has to certify.

# The undamped Newton iteration is what runs: on a strictly concave objective
# it converges here in a handful of steps, and it is the inner loop of the
# whole gate. `damp` is the fallback the caller reaches for when it does not,
# and costs one objective evaluation per step.
newtonMode <- function(
  design,
  precision,
  des,
  link,
  logIntegrand,
  damp,
  start
) {
  d <- ncol(design)
  theta <- start
  value <- if (damp) logIntegrand(theta) else NA_real_
  for (it in seq_len(if (damp) 200L else 50L)) {
    dv <- link$derivs(as.vector(design %*% theta), des$successes, des$trials)
    gradient <- as.vector(crossprod(design, dv$gradient)) - precision * theta
    hessian <- crossprod(design * dv$curvature, design) - diag(precision, d)
    step <- tryCatch(solve(hessian, gradient), error = function(e) NULL)
    if (is.null(step) || !all(is.finite(step))) {
      return(NULL)
    }
    # the Newton step, not the gradient, is the scale-free convergence test:
    # a glue value of 100 scales the gradient by 100 and no absolute
    # gradient bound then means the same thing across the axis
    if (max(abs(step)) <= 1e-10 * (1 + max(abs(theta)))) {
      return(theta)
    }
    if (!damp) {
      theta <- theta - step
      next
    }
    damping <- 1
    repeat {
      candidate <- theta - damping * step
      candidateValue <- logIntegrand(candidate)
      if (is.finite(candidateValue) && candidateValue >= value) {
        break
      }
      damping <- damping / 2
      if (damping < 1e-10) {
        return(theta)
      }
    }
    # a damped ascent that has stopped climbing is at the mode as closely as
    # this conditioning allows
    if (candidateValue - value <= 1e-12 * (1 + abs(value))) {
      return(candidate)
    }
    theta <- candidate
    value <- candidateValue
  }
  # the damped iteration only ever climbs, so its last point is the best
  # available; the undamped one has no such guarantee and reports failure
  if (damp) theta else NULL
}

blockQuadrature <- function(des, a, bz, sdVec, link, start) {
  d <- des$d
  design <- a * des$designMu + bz[1L] * des$designTau0 + bz[2L] * des$designTau1
  precision <- 1 / sdVec^2
  logNormalizer <- -0.5 * d * log(2 * pi) - sum(log(sdVec))
  successes <- des$successes
  failures <- des$trials - successes
  logIntegrand <- function(theta) {
    lg <- link$logs(as.vector(design %*% theta))
    sum(successes * lg$lp) +
      sum(failures * lg$lq) -
      0.5 * sum(precision * theta^2) +
      logNormalizer
  }
  theta <- newtonMode(design, precision, des, link, logIntegrand, FALSE, start)
  if (is.null(theta)) {
    theta <- newtonMode(
      design,
      precision,
      des,
      link,
      logIntegrand,
      TRUE,
      numeric(d)
    )
  }
  if (is.null(theta)) {
    stop("a block's Newton iteration did not reach its mode")
  }
  dv <- link$derivs(as.vector(design %*% theta), des$successes, des$trials)
  hessian <- crossprod(design * dv$curvature, design) - diag(precision, d)
  cholSigma <- t(chol(solve(-hessian)))
  rule <- ghCache[[d]]
  draws <- rep(theta, each = nrow(rule$node)) +
    sqrt(2) * rule$node %*% t(cholSigma)
  eta <- draws %*% t(design)
  lg <- link$logs(eta)
  logLik <- as.vector(
    matrix(lg$lp, nrow(draws)) %*%
      successes +
      matrix(lg$lq, nrow(draws)) %*% failures
  )
  logWeight <- rule$logWeight +
    logLik -
    0.5 * as.vector((draws^2) %*% precision) +
    logNormalizer
  maxLogWeight <- max(logWeight)
  weight <- exp(logWeight - maxLogWeight)
  total <- sum(weight)
  list(
    mode = theta,
    logMarginal = maxLogWeight +
      log(total) +
      0.5 * d * log(2) +
      sum(log(diag(cholSigma))),
    mean = as.vector(crossprod(draws, weight)) / total,
    probability = as.vector(
      crossprod(matrix(exp(lg$lp), nrow(draws)), weight)
    ) /
      total
  )
}

# ---- the exact posterior over the enumerated joint space ----

# glueGrid columns: a, b0, b1, logPrior. Every reported quantity is a
# posterior mean over the configuration mixture and the glue grid.
exactLatentBCF <- function(glueGrid, configurations, link, K, scales) {
  maxLogWeight <- -Inf
  total <- 0
  acc <- list(
    mu = numeric(K),
    tau = numeric(K),
    aMu = numeric(K),
    bTau = numeric(K),
    probability = matrix(0, K, 2L)
  )
  configurationWeight <- numeric(length(configurations))
  # the glue grid moves one coordinate at a time, so the previous point's mode
  # is a warm start that takes the Newton iteration from a handful of steps to
  # two or three - this is the gate's inner loop
  warm <- lapply(configurations, function(cfg) {
    lapply(cfg$design, function(des) numeric(des$d))
  })
  for (gi in seq_len(nrow(glueGrid))) {
    a <- glueGrid[gi, 1L]
    bz <- glueGrid[gi, 2:3]
    gluePrior <- glueGrid[gi, 4L]
    for (ci in seq_along(configurations)) {
      cfg <- configurations[[ci]]
      logMarginal <- 0
      mu <- numeric(K)
      tau <- numeric(K)
      probability <- matrix(0, K, 2L)
      for (bi in seq_along(cfg$design)) {
        des <- cfg$design[[bi]]
        block <- blockQuadrature(
          des,
          a,
          bz,
          ifelse(des$isMu, scales[1L], scales[2L]),
          link,
          warm[[ci]][[bi]]
        )
        warm[[ci]][[bi]] <- block$mode
        logMarginal <- logMarginal + block$logMarginal
        for (j in seq_len(des$d)) {
          p <- des$params[j]
          if (p <= cfg$pMu) {
            mu[cfg$muLeaf == p] <- block$mean[j]
          } else {
            tau[cfg$tauLeaf == (p - cfg$pMu)] <- block$mean[j]
          }
        }
        probability[cbind(des$cellOf, des$zOf)] <- block$probability
      }
      logWeight <- cfg$logPrior + gluePrior + logMarginal
      if (logWeight > maxLogWeight) {
        rescale <- exp(maxLogWeight - logWeight)
        total <- total * rescale
        acc <- lapply(acc, function(v) v * rescale)
        configurationWeight <- configurationWeight * rescale
        maxLogWeight <- logWeight
      }
      weight <- exp(logWeight - maxLogWeight)
      total <- total + weight
      configurationWeight[ci] <- configurationWeight[ci] + weight
      acc$mu <- acc$mu + weight * mu
      acc$tau <- acc$tau + weight * tau
      acc$aMu <- acc$aMu + weight * a * mu
      acc$bTau <- acc$bTau + weight * (bz[2L] - bz[1L]) * tau
      acc$probability <- acc$probability + weight * probability
    }
  }
  result <- lapply(acc, function(v) v / total)
  result$weights <- configurationWeight / total
  result
}

# ---- the glue axes ----

# Both axes take a tangent substitution on an open trapezoid in t. Under
# a = sd.control tan(t) the Cauchy measure is uniform in t, so every point
# carries the same log prior; the normal b axis keeps its Jacobian. The
# truncated tensor grid the gaussian gate uses does not meet the bar here:
# at 45 x 45 on [-3.5, 3.5] the b integral is off by 1.1e-3.
openTrapezoid <- function(n) {
  h <- pi / (n + 1)
  list(t = -pi / 2 + seq_len(n) * h, h = h)
}

fixedGlue <- matrix(c(1, 0, 1, 0), nrow = 1L)

aGlueGrid <- function(n) {
  rule <- openTrapezoid(n)
  cbind(sdControl * tan(rule$t), 0, 1, log(rule$h / pi))
}

bGlueGrid <- function(n) {
  rule <- openTrapezoid(n)
  b <- sqrt(bVar) * tan(rule$t)
  logWeight <- dnorm(tan(rule$t), log = TRUE) +
    log1p(tan(rule$t)^2) +
    log(rule$h)
  pair <- expand.grid(i0 = seq_len(n), i1 = seq_len(n))
  cbind(1, b[pair$i0], b[pair$i1], logWeight[pair$i0] + logWeight[pair$i1])
}

# The refinement self-check runs on the GLUE axis, where the error is:
# refining the leaf rule would certify nothing at 1e-13 by 12 nodes.
refineOrQuit <- function(name, coarse, fine) {
  gap <- abs(unlist(coarse) - unlist(fine))
  moved <- which.max(gap)
  cat(sprintf(
    "  refinement %-11s moves %.2e at %s\n",
    name,
    gap[moved],
    names(gap)[moved]
  ))
  if (gap[moved] > refineTolerance) {
    cat("FAIL: the glue quadrature is not converged; refusing to gate\n")
    quit(status = 1L)
  }
  invisible(NULL)
}

# ---- sampler collection ----

# One chain. The glue and each forest's fits are current-state reads, so the
# draws are taken one at a time with `thin` sweeps between them.
samplerFit <- function(seed, arm, updateA, updateB, ndpost, thin) {
  set.seed(seed)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 1L,
    n.cuts = arm$K - 1L,
    updateState = FALSE
  )
  # cgm() is internal, reached the way every sibling gate reaches it: as an
  # unevaluated argument, resolved inside dbarts()'s own frame
  host <- if (is.null(arm$weights)) {
    dbarts(
      arm$x,
      arm$y,
      control = control,
      family = arm$link$name,
      tree.prior = cgm(muPower, muBase)
    )
  } else {
    dbarts(
      arm$x,
      arm$y,
      weights = arm$weights,
      control = control,
      family = arm$link$name,
      tree.prior = cgm(muPower, muBase)
    )
  }
  bc <- dbarts:::bartcoreBCFSampler(
    host,
    arm$z,
    family = arm$link$name,
    n.trees.treatment = 1L,
    treatment.base = tauBase,
    treatment.power = tauPower,
    sd.control = sdControl,
    sd.moderate = sdModerate,
    b.prior.variance = bVar,
    update.a = updateA,
    update.b = updateB
  )
  bartcoreRun(bc, nBurn, 1L)
  K <- arm$K
  muDraws <- matrix(0, ndpost, K)
  tauDraws <- matrix(0, ndpost, K)
  aDraws <- numeric(ndpost)
  b0Draws <- numeric(ndpost)
  b1Draws <- numeric(ndpost)
  for (d in seq_len(ndpost)) {
    bartcoreRun(bc, 0L, thin)
    muDraws[d, ] <- bartcoreForestFits(bc, 0L)[arm$repObs, 1L]
    tauDraws[d, ] <- bartcoreForestFits(bc, 1L)[arm$repObs, 1L]
    glue <- bartcoreForestAmplitudes(bc)[, 1L]
    aDraws[d] <- glue[1L]
    b0Draws[d] <- glue[2L]
    b1Draws[d] <- glue[3L]
  }
  channels <- list()
  for (c in seq_len(K)) {
    channels[[sprintf("mu[%d]", c)]] <- muDraws[, c]
    channels[[sprintf("tau[%d]", c)]] <- tauDraws[, c]
    channels[[sprintf("a mu[%d]", c)]] <- aDraws * muDraws[, c]
    channels[[sprintf("(b1-b0) tau[%d]", c)]] <- (b1Draws - b0Draws) *
      tauDraws[, c]
    for (zg in 1:2) {
      bDraws <- if (zg == 1L) b0Draws else b1Draws
      channels[[sprintf("F(eta[%d,%d])", c, zg - 1L)]] <- arm$link$linkinv(
        aDraws * muDraws[, c] + bDraws * tauDraws[, c]
      )
    }
  }
  channels
}

# A batch-means standard error whose batch LENGTH is raised until the batch
# means decorrelate, with any residual lag-1 correlation charged to the
# estimate as an AR(1) inflation. A fixed 400 batches is not honest here:
# mode 2a's chain is metastable (its (a, mu) state sits in one place for
# 1e5 kept draws at a time) and the K = 3 arms switch tree partitions
# slowly, so at 400 batches those channels understate their own error by up
# to 30x - measured against the spread of independent seeds - and the gate
# would fail a correct sampler. The lag-1 autocorrelation of the batch means
# is reported beside every z, so a mixing-limited channel says so; the
# correlation is capped at 0.95 rather than allowed to make the gate vacuous.
batchStats <- function(v) {
  best <- NULL
  for (nBatches in c(400L, 200L, 100L, 50L, 25L)) {
    if (length(v) %/% nBatches < 10L) {
      next
    }
    len <- (length(v) %/% nBatches) * nBatches
    bm <- colMeans(matrix(v[seq_len(len)], ncol = nBatches))
    acf1 <- if (sd(bm) > 0) cor(bm[-nBatches], bm[-1L]) else 0
    best <- list(bm = bm, nBatches = nBatches, acf1 = acf1)
    if (abs(acf1) < 0.1) {
      break
    }
  }
  if (is.null(best)) {
    stop("too few kept draws for a batch-means standard error")
  }
  r <- min(max(best$acf1, 0), 0.95)
  list(
    mean = mean(v),
    se = sd(best$bm) / sqrt(best$nBatches) * sqrt((1 + r) / (1 - r)),
    acf1 = best$acf1
  )
}

# Seeds are pooled as the mean of per-seed means with sqrt(sum se_s^2) / S,
# never batched across a chain seam, and floored by the spread of the seed
# means themselves - the one estimator that sees a state no single chain
# leaves. The floor is purely protective: it can only widen the interval.
samplerChannels <- function(arm, updateA, updateB, ndpost, thin, seeds) {
  fits <- lapply(
    seq_len(seeds),
    samplerFit,
    arm = arm,
    updateA = updateA,
    updateB = updateB,
    ndpost = ndpost,
    thin = thin
  )
  out <- list()
  for (nm in names(fits[[1L]])) {
    stats <- lapply(fits, function(f) batchStats(f[[nm]]))
    means <- vapply(stats, function(s) s$mean, numeric(1L))
    within <- sqrt(sum(vapply(stats, function(s) s$se^2, numeric(1L)))) / seeds
    between <- if (seeds >= 3L) sd(means) / sqrt(seeds) else 0
    out[[nm]] <- list(
      mean = mean(means),
      se = max(within, between),
      acf1 = max(abs(vapply(stats, function(s) s$acf1, numeric(1L))))
    )
  }
  out
}

# ---- the arms ----

# The design's data, one realization at a pinned seed shared by the links: the
# same uniforms are thresholded at each link's own group probabilities.
makeArm <- function(link, K, muValue, tauValue, aggregate = FALSE) {
  cell <- rep(seq_len(K), each = nPerCell)
  z <- rep(rep(c(0, 1), nPerCell / 2L), K)
  index <- muValue[cell] + z * tauValue[cell]
  set.seed(dataSeed)
  y <- as.double(runif(length(index)) < link$linkinv(index))
  group <- cbind(cell, ifelse(z != 0, 2L, 1L))
  sCz <- matrix(0, K, 2L)
  nCz <- matrix(0, K, 2L)
  for (i in seq_along(y)) {
    sCz[group[i, 1L], group[i, 2L]] <- sCz[group[i, 1L], group[i, 2L]] + y[i]
    nCz[group[i, 1L], group[i, 2L]] <- nCz[group[i, 1L], group[i, 2L]] + 1
  }
  weights <- NULL
  if (aggregate) {
    # Two rows per group, y = 1 with count s and y = 0 with count n - s. The
    # oracle reads the group sufficient statistics alone, so this is the same
    # target bit for bit as the unit-row arm, and both fits must meet it.
    cell <- rep(seq_len(K), each = 4L)
    z <- rep(c(0, 0, 1, 1), K)
    y <- rep(c(1, 0, 1, 0), K)
    group <- cbind(cell, rep(c(1L, 1L, 2L, 2L), K))
    weights <- ifelse(y == 1, sCz[group], nCz[group] - sCz[group])
  }
  x <- matrix(as.double(cell), ncol = 1L)
  # uniform cuts separate the cells exactly (dbarts n.cuts = K - 1)
  cuts <- min(x) + seq_len(K - 1L) * (max(x) - min(x)) / K
  stopifnot(identical(findInterval(x[, 1L], cuts) + 1L, cell))
  list(
    link = link,
    K = K,
    x = x,
    y = as.double(y),
    z = as.double(z),
    weights = weights,
    repObs = vapply(seq_len(K), function(c) which(cell == c)[1L], integer(1L)),
    sCz = sCz,
    nCz = nCz,
    scales = c(link$anchor, sdModerate * link$anchor / nodeScaleDivisor),
    configurations = makeConfigurations(K, sCz, nCz)
  )
}

# ---- reporting ----

anyFailure <- FALSE
worstZ <- list()

reportChannel <- function(mode, label, fitted, exact) {
  z <- (fitted$mean - exact) / fitted$se
  failed <- !is.finite(z) || abs(z) > zBound
  cat(sprintf(
    "  %-16s %+9.5f %+9.5f %9.2e %+7.2f %+6.2f%s\n",
    label,
    fitted$mean,
    exact,
    fitted$se,
    z,
    fitted$acf1,
    if (failed) " <- FAIL" else ""
  ))
  worstZ[[mode]] <<- max(worstZ[[mode]], abs(z))
  failed
}

# `matched` names the mode's two leaf channels, each a vector over cells; the
# probability surface is matched in every mode.
reportMode <- function(name, arm, fitted, exact, matched, elapsed) {
  cat(sprintf("%s (%.1f s)\n", name, elapsed))
  cat(sprintf(
    "  %-16s %9s %9s %9s %7s %6s\n",
    "quantity",
    "sampler",
    "exact",
    "se",
    "z",
    "acf1"
  ))
  worstZ[[name]] <<- 0
  failed <- FALSE
  for (label in names(matched)) {
    for (c in seq_len(arm$K)) {
      key <- sprintf("%s[%d]", label, c)
      failed <- reportChannel(name, key, fitted[[key]], matched[[label]][c]) ||
        failed
    }
  }
  for (c in seq_len(arm$K)) {
    for (zg in 1:2) {
      key <- sprintf("F(eta[%d,%d])", c, zg - 1L)
      failed <- reportChannel(
        name,
        key,
        fitted[[key]],
        exact$probability[c, zg]
      ) ||
        failed
    }
  }
  failed
}

runMode <- function(
  name,
  arm,
  exact,
  matched,
  updateA,
  updateB,
  ndpost,
  thin,
  seeds
) {
  fitted <- NULL
  elapsed <- system.time(
    fitted <- samplerChannels(arm, updateA, updateB, ndpost, thin, seeds)
  )[["elapsed"]]
  reportMode(name, arm, fitted, exact, matched, elapsed)
}

# ---- the design, its realized counts and the runner-up guard ----

cat(sprintf(
  "BCF latent exact gate (%s): %d kept draws thin %d over %d seeds\n",
  if (quick) "quick" else "full",
  nKept,
  nThin,
  nSeeds
))
cat(sprintf(
  "mode 2b: %d kept draws thin %d over %d seeds\n",
  nKept2b,
  nThin2b,
  nSeeds2b
))

arms <- list(
  probit = makeArm(probitLink, 2L, muValue2, tauValue2),
  logistic = makeArm(logisticLink, 2L, muValue2, tauValue2),
  aggregated = makeArm(logisticLink, 2L, muValue2, tauValue2, aggregate = TRUE),
  probit3 = makeArm(probitLink, 3L, muValue3, tauValue3),
  logistic3 = makeArm(logisticLink, 3L, muValue3, tauValue3)
)

cat("\nconfiguration weights at the realized counts\n")
fixedStart <- proc.time()[["elapsed"]]
exactFixed <- list()
for (nm in names(arms)) {
  arm <- arms[[nm]]
  exactFixed[[nm]] <- exactLatentBCF(
    fixedGlue,
    arm$configurations,
    arm$link,
    arm$K,
    arm$scales
  )
  weights <- sort(exactFixed[[nm]]$weights, decreasing = TRUE)
  cat(sprintf(
    "  %-11s s %-14s %2d configurations, top %s\n",
    nm,
    paste(as.vector(t(arm$sCz)), collapse = " "),
    length(weights),
    paste(sprintf("%.4f", weights[1:3]), collapse = " ")
  ))
  if (weights[2L] < runnerUpFloor) {
    cat(sprintf(
      "FAIL: %s runner-up weight %.4f under %.2f; no tree channel\n",
      nm,
      weights[2L],
      runnerUpFloor
    ))
    quit(status = 1L)
  }
}

cat(sprintf(
  "  fixed-glue oracles: %.1f s\n",
  proc.time()[["elapsed"]] - fixedStart
))

# ---- the free-glue oracles and their refinement ----

cat("\nglue quadrature\n")
oracleStart <- proc.time()[["elapsed"]]
exactFreeA <- list()
exactFreeB <- list()
for (nm in c("probit", "logistic")) {
  arm <- arms[[nm]]
  oracle <- function(grid) {
    exactLatentBCF(grid, arm$configurations, arm$link, arm$K, arm$scales)
  }
  exactFreeA[[nm]] <- oracle(aGlueGrid(aN))
  refineOrQuit(sprintf("a/%s", nm), exactFreeA[[nm]], oracle(aGlueGrid(aNFine)))
  exactFreeB[[nm]] <- oracle(bGlueGrid(bN))
  refineOrQuit(sprintf("b/%s", nm), exactFreeB[[nm]], oracle(bGlueGrid(bNFine)))
}
cat(sprintf(
  "  free-glue oracles and refinement: %.1f s\n",
  proc.time()[["elapsed"]] - oracleStart
))

# ---- run the modes ----

cat("\n")
for (nm in c("probit", "logistic")) {
  arm <- arms[[nm]]
  anyFailure <- runMode(
    sprintf("mode 1 (fixed glue), %s", nm),
    arm,
    exactFixed[[nm]],
    list(mu = exactFixed[[nm]]$mu, tau = exactFixed[[nm]]$tau),
    FALSE,
    FALSE,
    nKept,
    nThin,
    nSeeds
  ) ||
    anyFailure
  anyFailure <- runMode(
    sprintf("mode 2a (a free), %s", nm),
    arm,
    exactFreeA[[nm]],
    list(`a mu` = exactFreeA[[nm]]$aMu, tau = exactFreeA[[nm]]$tau),
    TRUE,
    FALSE,
    nKept,
    nThin,
    nSeeds
  ) ||
    anyFailure
  anyFailure <- runMode(
    sprintf("mode 2b (b free), %s", nm),
    arm,
    exactFreeB[[nm]],
    list(
      mu = exactFreeB[[nm]]$mu,
      `(b1-b0) tau` = exactFreeB[[nm]]$bTau
    ),
    FALSE,
    TRUE,
    nKept2b,
    nThin2b,
    nSeeds2b
  ) ||
    anyFailure
}

# the trial-count channel: the same oracle, two weighted rows per group
anyFailure <- runMode(
  "mode 1, logistic aggregated counts",
  arms$aggregated,
  exactFixed$aggregated,
  list(mu = exactFixed$aggregated$mu, tau = exactFixed$aggregated$tau),
  FALSE,
  FALSE,
  nKept,
  nThin,
  nSeeds
) ||
  anyFailure

# depth decay and cut selection, which a two-cell design cannot see: at two
# cells the only interior node is the root, where base / (1 + depth)^power
# equals base for every power
for (nm in c("probit3", "logistic3")) {
  arm <- arms[[nm]]
  anyFailure <- runMode(
    sprintf("mode 1, K = 3, %s", arm$link$name),
    arm,
    exactFixed[[nm]],
    list(mu = exactFixed[[nm]]$mu, tau = exactFixed[[nm]]$tau),
    FALSE,
    FALSE,
    nKept,
    nThin,
    nSeeds
  ) ||
    anyFailure
}

cat("\nworst |z| by configuration\n")
for (nm in names(worstZ)) {
  cat(sprintf("  %-38s %5.2f\n", nm, worstZ[[nm]]))
}

if (anyFailure) {
  cat("\nFAIL: latent BCF sampler deviates from the exact posterior\n")
  quit(status = 1L)
}
cat("\nOK: latent BCF sampler matches the exact posterior\n")
