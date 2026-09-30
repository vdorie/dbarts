#!/usr/bin/env Rscript

# Exact enumeration of the monotone single-tree posterior over multi-split
# trees, with an R prototype of the birth/death + leaf Gibbs chain
# (docs/plans/monotone-exact-birth-death.md), under both monotone priors.
# "leaf": the CGM tree prior and, given a tree T, iid (c-inflated) normal
# leaves restricted to the monotone cone C(T), normalized per tree by
# Z_T = P(unconstrained leaves in C(T)) = e(P_T) / L!, e the number of linear
# extensions of the leaf order the constraints impose. A tree's posterior
# weight is p(T) exp(sum of unconstrained leaf log-marginals) P_post(C(T)) / Z_T.
# "joint": the CGM tree prior and iid leaves conditioned on the cone together,
# unnormalized, so the weight drops the 1 / Z_T.
#
# Each design is a grid of cells on one to three predictors with 10 rows per
# cell and sigma fixed; its tree structures are enumerated exactly. Draws are
# compared per root rule, the split variable and cut at the root (the only way
# to change it is through the root-only tree, which some designs visit rarely,
# so a chain can hold one root rule for its whole run), each group against its
# conditional law, by a Hotelling T^2 test of the batch-mean cell frequencies
# (cells under 0.2% pooled into the omitted reference cell). A group fails at
# p < 1e-4.
#
# Usage: Rscript monotone-exact-enumeration.R [quick] [mode] [prior] [design ...]
#   prior  leaf, joint (default: both, each run naming its prior); engine and
#          prototype modes only, the others take "leaf"
#   mode   engine (default): the installed dbarts sampler; exits 1 on failure
#          prototype: the R prototype of the corrected move; exits 1 on failure
#          prototype-old: the prototype with the conditional-d ratio the
#            engine shipped with (reproduces the engine before the fix)
#          hybrid: the prototype with every move whose finer tree's pair
#            components hold 3 or more leaves switched to the count-free
#            Barker move (docs/design/monotone-barker-hybrid.md): two-coin
#            acceptance, the coin an exact linear extension by Huber's
#            coupling from the past; exits 1 on failure
#          hybrid-nocoin: hybrid with the coin always heads (drops theta);
#            must fail
#          zcheck: Z_T = e(P_T) / L! against a brute-force monotonicity check
#          unnormalized: the exact posterior against the one without 1 / Z_T
#            (the BART prior conditioned on every tree being monotone); no sampling
#   design c1 c2 c3 cN cM cF (default: all)

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
args <- setdiff(args, "quick")
modes <- c(
  "engine",
  "prototype",
  "prototype-old",
  "hybrid",
  "hybrid-nocoin",
  "zcheck",
  "unnormalized"
)
mode <- intersect(args, modes)
mode <- if (length(mode)) mode[1L] else "engine"
priors <- intersect(c("leaf", "joint"), args)
if (!length(priors)) {
  priors <- c("leaf", "joint")
}
if (!mode %in% c("engine", "prototype")) {
  priors <- "leaf"
}
designs <- list(
  # one constrained predictor, four cells: every tree is a chain
  c1 = list(nc = 4L, dirs = 1L, mu = c(0, 0.2, 0.4, 0.6), sigma = 0.6),
  # two constrained predictors
  c2 = list(
    nc = c(3L, 2L),
    dirs = c(1L, 1L),
    mu = c(0, 0.2, 0.4, 0.3, 0.5, 0.7),
    sigma = 1.5
  ),
  # one constrained and one free predictor
  c3 = list(
    nc = c(3L, 2L),
    dirs = c(1L, 0L),
    mu = c(0, 0.2, 0.4, 0.3, 0.5, 0.7),
    sigma = 1.0
  ),
  # the two halves of an x1 split change level at different x2 cuts, so the
  # posterior puts about a quarter of its mass on N-shaped leaf orders
  cN = list(
    nc = c(2L, 3L),
    dirs = c(1L, 0L),
    mu = c(0, 3, 0.6, 3, 0.6, 3.6),
    sigma = 0.6,
    # a tree rooted on x2 holds its subtree roots for whole runs and carries
    # no N-shaped order, so only the x1-rooted group is tested
    groups = "x1 cut 1"
  ),
  # cN's cells reflected along x1, with x1 decreasing: the same orders reached
  # through a decreasing axis, the move's and the redraws' other direction
  cM = list(
    nc = c(2L, 3L),
    dirs = c(-1L, 0L),
    mu = c(3, 0, 3, 0.6, 3.6, 0.6),
    sigma = 0.6,
    groups = "x1 cut 1"
  ),
  # x1 constrained and x2 a free two-level factor, x1 rising at one level and
  # falling at the other: the engine splits x2 by a level subset, which the
  # enumeration reads as its one cut, so the order relates leaves through a
  # level set
  cF = list(
    nc = c(3L, 2L),
    dirs = c(1L, 0L),
    mu = c(0, 0.3, 0.6, 0.6, 0.3, 0),
    sigma = 1.0,
    factor = 2L
  )
)
chosen <- intersect(args, names(designs))
if (!length(chosen)) {
  chosen <- names(designs)
}
nDraws <- if (quick) 300000L else 900000L
nChains <- 4L
nBurn <- 2000L
pFail <- 1e-4

base <- 0.95
power <- 2
kLeaf <- 2
nodeScale <- 0.5
cScale <- sqrt(pi / (pi - 1))

# ---- trees, leaf orders and their counts -------------------------------------

regionKey <- function(lo, hi) paste(lo, hi, sep = "-", collapse = ",")
ruleKey <- function(lo, hi, v, cut) {
  paste0("[", regionKey(lo, hi), "|", v, ":", cut, "]")
}
treeKey <- function(rules) {
  if (length(rules)) paste(sort(rules), collapse = "") else "ROOT"
}

# every tree on the region [lo, hi] (cell indices, inclusive) at depth d, with
# its CGM log prior: grow with base / (1 + d)^power when a split is available,
# the variable uniform over those with a cut, the cut uniform over its cuts
enumerateTrees <- function(lo, hi, d = 0) {
  avail <- which(hi > lo)
  pSplit <- if (length(avail)) base / (1 + d)^power else 0
  out <- list(list(
    rules = character(0),
    logPrior = log1p(-pSplit),
    leaves = list(list(lo = lo, hi = hi))
  ))
  for (v in avail) {
    for (cut in lo[v]:(hi[v] - 1L)) {
      hiLeft <- hi
      hiLeft[v] <- cut
      loRight <- lo
      loRight[v] <- cut + 1L
      lefts <- enumerateTrees(lo, hiLeft, d + 1)
      rights <- enumerateTrees(loRight, hi, d + 1)
      for (l in lefts) {
        for (r in rights) {
          out[[length(out) + 1L]] <- list(
            rules = c(ruleKey(lo, hi, v, cut), l$rules, r$rules),
            logPrior = log(pSplit) -
              log(length(avail)) -
              log(hi[v] - lo[v]) +
              l$logPrior +
              r$logPrior,
            leaves = c(l$leaves, r$leaves)
          )
        }
      }
    }
  }
  out
}

# R[j, k] when mu_j <= mu_k is imposed directly: along a constrained axis j's
# box ends one cell below k's start, and the boxes share a cell on every other
# axis (the pairs the engine's monotoneNeighborBounds constrains)
leafRelation <- function(leaves, dirs) {
  n <- length(leaves)
  rel <- matrix(FALSE, n, n)
  for (j in seq_len(n)) {
    for (k in seq_len(n)) {
      if (j == k) {
        next
      }
      a <- leaves[[j]]
      b <- leaves[[k]]
      for (i in which(dirs != 0L)) {
        if (a$hi[i] + 1L != b$lo[i]) {
          next
        }
        other <- setdiff(seq_along(dirs), i)
        shared <- all(
          pmax(a$lo[other], b$lo[other]) <= pmin(a$hi[other], b$hi[other])
        )
        if (shared) {
          if (dirs[i] > 0L) rel[j, k] <- TRUE else rel[k, j] <- TRUE
        }
      }
    }
  }
  rel
}

# the down-sets of the order, in layers by size, as logical membership rows
downSetLayers <- function(rel) {
  n <- nrow(rel)
  layers <- list(matrix(FALSE, 1L, n))
  for (size in seq_len(n)) {
    prev <- layers[[size]]
    nxt <- list()
    for (r in seq_len(nrow(prev))) {
      set <- prev[r, ]
      for (x in which(!set)) {
        if (all(set[rel[, x]])) {
          grown <- set
          grown[x] <- TRUE
          nxt[[paste(which(grown), collapse = ".")]] <- grown
        }
      }
    }
    layers[[size + 1L]] <- do.call(rbind, unname(nxt))
  }
  layers
}

# a quantity propagated over the down-set lattice: value(I + x) accumulates
# step(x, value(I)) over every x minimal outside I; returns value(all leaves)
propagate <- function(rel, start, step) {
  layers <- downSetLayers(rel)
  vals <- list(start)
  for (size in seq_len(nrow(rel))) {
    prev <- layers[[size]]
    cur <- layers[[size + 1L]]
    curKeys <- apply(cur, 1L, function(s) paste(which(s), collapse = "."))
    nextVals <- vector("list", nrow(cur))
    for (r in seq_len(nrow(prev))) {
      set <- prev[r, ]
      for (x in which(!set)) {
        if (!all(set[rel[, x]])) {
          next
        }
        grown <- set
        grown[x] <- TRUE
        idx <- match(paste(which(grown), collapse = "."), curKeys)
        add <- step(x, vals[[r]])
        nextVals[[idx]] <- if (is.null(nextVals[[idx]])) {
          add
        } else {
          nextVals[[idx]] + add
        }
      }
    }
    vals <- nextVals
  }
  vals[[1L]]
}

countExtensions <- function(rel) propagate(rel, 1, function(x, v) v)

# P(independent N(m_k, s_k^2) leaves respect the order): the same lattice
# recursion carrying G_I(x) = P(leaves of I fall below x in an order the
# constraints admit), the new minimum's density integrated against it by the
# trapezoid rule at steps min(s) / 100 and / 200, Richardson-extrapolated
orderProbability <- function(rel, m, s) {
  atStep <- function(h) {
    grid <- seq(min(m - 10 * s), max(m + 10 * s), by = h)
    step <- function(x, g) {
      f <- dnorm(grid, m[x], s[x]) * g
      c(0, cumsum((f[-1L] + f[-length(f)]) * h / 2))
    }
    g <- propagate(rel, rep(1, length(grid)), step)
    g[length(g)]
  }
  # trapezoid error is O(h^2): one Richardson step removes it
  coarse <- atStep(min(s) / 100)
  fine <- atStep(min(s) / 200)
  (4 * fine - coarse) / 3
}

# ---- designs: data, enumeration, exact law -----------------------------------

buildDesign <- function(spec, nPer = 10L, seed = 11L) {
  set.seed(seed)
  nc <- spec$nc
  cells <- as.matrix(expand.grid(lapply(nc, seq_len)))
  cellOf <- rep(seq_len(nrow(cells)), each = nPer)
  x <- matrix(
    as.double(cells[cellOf, ]),
    ncol = length(nc),
    dimnames = list(NULL, paste0("x", seq_along(nc)))
  )
  y <- spec$mu[cellOf] + rnorm(length(cellOf), sd = 0.3)
  yRange <- max(y) - min(y)
  z <- (y - min(y)) / yRange - 0.5
  residVar <- (spec$sigma / yRange)^2
  sumZ <- as.numeric(tapply(z, cellOf, sum))
  count <- as.numeric(table(cellOf))
  cellIndex <- function(leaf) {
    box <- as.matrix(expand.grid(lapply(seq_along(nc), function(i) {
      leaf$lo[i]:leaf$hi[i]
    })))
    drop((box - 1) %*% cumprod(c(1, nc[-length(nc)]))) + 1
  }
  trees <- enumerateTrees(rep(1L, length(nc)), as.integer(nc))
  for (t in seq_along(trees)) {
    tr <- trees[[t]]
    tr$key <- treeKey(tr$rules)
    tr$leafKeys <- vapply(tr$leaves, function(l) regionKey(l$lo, l$hi), "")
    tr$sumZ <- vapply(tr$leaves, function(l) sum(sumZ[cellIndex(l)]), 0)
    tr$n <- vapply(tr$leaves, function(l) sum(count[cellIndex(l)]), 0)
    tr$rel <- leafRelation(tr$leaves, spec$dirs)
    tr$related <- rowSums(tr$rel) + colSums(tr$rel) > 0
    tr$tau <- nodeScale / kLeaf * ifelse(tr$related, cScale, 1)
    prec <- tr$n / residVar + 1 / tr$tau^2
    tr$m <- (tr$sumZ / residVar) / prec
    tr$s <- sqrt(1 / prec)
    tr$base <- -0.5 *
      log(1 + tr$n * tr$tau^2 / residVar) +
      0.5 * (tr$sumZ / residVar)^2 / prec
    tr$e <- countExtensions(tr$rel)
    tr$logZ <- log(tr$e) - lfactorial(length(tr$leaves))
    tr$logPostCone <- log(orderProbability(tr$rel, tr$m, tr$s))
    tr$avail <- lapply(tr$leaves, function(l) which(l$hi > l$lo))
    trees[[t]] <- tr
  }
  names(trees) <- vapply(trees, `[[`, "", "key")
  # "leaf" divides by Z_T, "joint" does not
  lawOf <- function(normalized) {
    logW <- vapply(
      trees,
      function(tr) {
        tr$logPrior +
          sum(tr$base) +
          tr$logPostCone -
          if (normalized) tr$logZ else 0
      },
      0
    )
    law <- exp(logW - max(logW))
    law / sum(law)
  }
  list(
    spec = spec,
    x = x,
    y = y,
    trees = trees,
    law = lawOf(TRUE),
    laws = list(leaf = lawOf(TRUE), joint = lawOf(FALSE)),
    residVar = residVar
  )
}

rootRule <- function(keys) {
  out <- rep("", length(keys))
  for (i in seq_along(keys)) {
    if (keys[i] == "ROOT") {
      next
    }
    rules <- regmatches(keys[i], gregexpr("\\[[^]]*\\]", keys[i]))[[1L]]
    lens <- vapply(
      rules,
      function(r) {
        reg <- sub("^\\[(.*)\\|.*$", "\\1", r)
        b <- matrix(
          as.integer(unlist(strsplit(strsplit(reg, ",")[[1L]], "-"))),
          2L
        )
        sum(b[2L, ] - b[1L, ])
      },
      0
    )
    out[i] <- sub(
      "^.*\\|([0-9]+):([0-9]+)\\]$",
      "x\\1 cut \\2",
      rules[which.max(lens)]
    )
  }
  out
}

# ---- hybrid: the count-free Barker move -------------------------------------

# the transitive closure of an order given as a relation matrix
closeOrder <- function(rel) {
  n <- nrow(rel)
  for (k in seq_len(n)) {
    rel <- rel | outer(rel[, k], rel[k, ], `&`)
  }
  rel
}

# connected components of an order's comparability graph, as labels
orderComponents <- function(rel) {
  n <- nrow(rel)
  adj <- rel | t(rel)
  label <- rep(0L, n)
  for (s in seq_len(n)) {
    if (label[s]) {
      next
    }
    label[s] <- s
    stack <- s
    while (length(stack)) {
      x <- stack[1L]
      stack <- stack[-1L]
      new <- which(adj[x, ] & !label)
      label[new] <- s
      stack <- c(stack, new)
    }
  }
  label
}

# One uniform linear extension of a closed order (closed[j, k] when j < k),
# by Huber's bounding chain for the Karzanov-Khachiyan chain and his
# non-Markovian coupling from the past (Discrete Math. 306, 2006, sections 3
# and 4), doubling the block length. Returns the elements in order.
huberExtension <- function(closed) {
  n <- nrow(closed)
  if (n == 1L) {
    return(1L)
  }
  # relabel so the identity is an extension: fewer predecessors first
  ord <- order(colSums(closed))
  p <- closed[ord, ord, drop = FALSE]
  state <- new.env()
  run <- function(iv, cv, advance) {
    r <- rep(n, n)
    who <- c(rep(0L, n - 1L), 1L)
    inPlay <- 1L
    x <- state$x
    for (t in seq_along(iv)) {
      i <- iv[t]
      cc <- cv[t]
      a <- who[i]
      b <- who[i + 1L]
      if (advance) {
        cx <- if (b > 0L && x[i] == b) 1L - cc else cc
        if (cx == 1L && !p[x[i], x[i + 1L]]) x[c(i, i + 1L)] <- x[c(i + 1L, i)]
      }
      if (cc == 1L) {
        if (a > 0L && b > 0L) {
          if (!p[a, b]) {
            r[a] <- i + 1L
            r[b] <- i
            who[i] <- b
            who[i + 1L] <- a
          }
        } else if (a > 0L) {
          r[a] <- i + 1L
          who[i] <- 0L
          who[i + 1L] <- a
        } else if (b > 0L) {
          r[b] <- i
          who[i + 1L] <- 0L
          who[i] <- b
        }
      }
      if (inPlay < n && who[n] == 0L) {
        inPlay <- inPlay + 1L
        r[inPlay] <- n
        who[n] <- inPlay
      }
    }
    state$x <- x
    state$r <- r
    inPlay == n
  }
  base <- max(16L, n * n)
  blocks <- list()
  level <- 0L
  repeat {
    len <- base * 2L^level
    blocks[[level + 1L]] <- list(
      iv = sample.int(n - 1L, len, replace = TRUE),
      cv = sample.int(2L, len, replace = TRUE) - 1L
    )
    if (run(blocks[[level + 1L]]$iv, blocks[[level + 1L]]$cv, FALSE)) {
      break
    }
    level <- level + 1L
  }
  x <- integer(n)
  x[state$r] <- seq_len(n)
  state$x <- x
  for (l in rev(seq_len(level))) {
    run(blocks[[l]]$iv, blocks[[l]]$cv, TRUE)
  }
  ord[state$x]
}

# the coin: one exact uniform extension of the finer tree's components that
# hold the pair, drawn per component and riffled uniformly; heads when c2
# immediately follows c1
adjacencyCoin <- function(tr, c1, c2) {
  held <- unique(tr$comp[c(c1, c2)])
  parts <- lapply(held, function(h) {
    members <- which(tr$comp == h)
    members[huberExtension(tr$closed[members, members, drop = FALSE])]
  })
  if (length(parts) == 1L) {
    seq <- parts[[1L]]
  } else {
    a <- length(parts[[1L]])
    b <- length(parts[[2L]])
    slots <- sort(sample.int(a + b, a))
    seq <- integer(a + b)
    seq[slots] <- parts[[1L]]
    seq[-slots] <- parts[[2L]]
  }
  match(c2, seq) == match(c1, seq) + 1L
}

# Barker's acceptance by the two-coin algorithm (Goncalves, Latuszynski and
# Roberts 2017, section 2) for a pair whose finer tree's weight carries the
# unknown factor theta: birth accepts with c theta / (1 + c theta), a death
# with 1 / (1 + c theta). Returns the decision and the coins drawn.
twoCoin <- function(logC, coin, birth) {
  draws <- 0L
  repeat {
    if (runif(1L) < stats::plogis(logC)) {
      draws <- draws + 1L
      if (coin()) {
        return(list(accept = birth, draws = draws))
      }
    } else {
      return(list(accept = !birth, draws = draws))
    }
  }
}

# ---- samplers ----------------------------------------------------------------

engineKeys <- function(design, nDraw, seed, prior) {
  suppressPackageStartupMessages(library(dbarts))
  nc <- design$spec$nc
  dirs <- design$spec$dirs
  block <- 5000L
  ctl <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 1L,
    keepTrees = TRUE,
    n.samples = block,
    updateState = TRUE,
    seed = seed,
    n.cuts = as.integer(nc - 1L)
  )
  mono <- setNames(dirs, colnames(design$x))[dirs != 0L]
  # a factor axis's cells are its levels, in order
  factorAxis <- design$spec$factor
  x <- as.data.frame(design$x)
  for (v in factorAxis) {
    x[[v]] <- factor(letters[x[[v]]], levels = letters[seq_len(nc[v])])
  }
  sampler <- dbarts(
    y ~ .,
    data.frame(y = design$y, x),
    control = ctl,
    tree.prior = cgm(power, base),
    leaf.prior = normal(kLeaf),
    family = gaussian(sigma = fixed(design$spec$sigma^2)),
    monotone = monotone(mono, prior = prior)
  )
  cuts <- lapply(nc, function(n) seq(1, n, length.out = n + 1L)[-c(1L, n + 1L)])
  # a two-level factor's split is its one cut; the engine's left child holds
  # the higher level when the lower one goes right
  keyOf <- function(var, value, directions) {
    i <- 0L
    rules <- character(0)
    walk <- function(lo, hi) {
      i <<- i + 1L
      if (var[i] == -1L) {
        return(invisible())
      }
      v <- var[i]
      swap <- v %in% factorAxis
      cut <- if (swap) 1L else match(round(value[i], 8), round(cuts[[v]], 8))
      swap <- swap && startsWith(directions[i], "R")
      rules <<- c(rules, ruleKey(lo, hi, v, cut))
      hiLeft <- hi
      hiLeft[v] <- cut
      loRight <- lo
      loRight[v] <- cut + 1L
      if (swap) {
        walk(loRight, hi)
        walk(lo, hiLeft)
      } else {
        walk(lo, hiLeft)
        walk(loRight, hi)
      }
    }
    walk(rep(1L, length(nc)), as.integer(nc))
    treeKey(rules)
  }
  keys <- character(0)
  for (b in seq_len(ceiling(nDraw / block))) {
    sampler$run(if (b == 1L) nBurn else 0L, block)
    tr <- sampler$getTrees()
    bySample <- split(tr, tr$sample)[as.character(sort(unique(tr$sample)))]
    keys <- c(
      keys,
      vapply(
        bySample,
        function(d) {
          keyOf(d$var, d$value, if (is.null(d$directions)) NA else d$directions)
        },
        ""
      )
    )
  }
  keys[seq_len(nDraw)]
}

drawTruncated <- function(m, s, a, b) {
  za <- (a - m) / s
  zb <- (b - m) / s
  u <- runif(1L)
  if (za > 0) {
    qa <- pnorm(za, lower.tail = FALSE)
    qb <- pnorm(zb, lower.tail = FALSE)
    m + s * qnorm(qa - u * (qa - qb), lower.tail = FALSE)
  } else {
    pa <- pnorm(za)
    pb <- pnorm(zb)
    m + s * qnorm(pa + u * (pb - pa))
  }
}

# The prototype chain on (T, M): birth/death with the engine's proposal (leaf
# uniform over birthable leaves, rule from the prior, birth probability 1 at
# the root, 0.5 otherwise, 0 with nothing birthable; death uniform over nodes
# with two leaf children), the touched leaves integrated out given the rest
# and redrawn exactly on acceptance, then a Gibbs sweep over the leaves.
# exact = FALSE divides each touched marginal by its conditional prior cone
# mass d and drops the Z_T ratio, as the engine did before the fix; the "joint"
# prior drops the Z_T ratio alone.
prototypeKeys <- function(
  design,
  nDraw,
  seed,
  exact,
  switchAt = Inf,
  coin = TRUE,
  prior = "leaf"
) {
  set.seed(seed)
  trees <- design$trees
  dirs <- design$spec$dirs
  splitRule <- function(rule) {
    p <- regmatches(rule, regexec("^\\[(.*)\\|([0-9]+):([0-9]+)\\]$", rule))[[
      1L
    ]]
    b <- matrix(
      as.integer(unlist(strsplit(strsplit(p[2L], ",")[[1L]], "-"))),
      2L
    )
    list(
      lo = b[1L, ],
      hi = b[2L, ],
      v = as.integer(p[3L]),
      cut = as.integer(p[4L])
    )
  }
  childKeys <- function(r) {
    hiLeft <- r$hi
    hiLeft[r$v] <- r$cut
    loRight <- r$lo
    loRight[r$v] <- r$cut + 1L
    c(regionKey(r$lo, hiLeft), regionKey(loRight, r$hi))
  }
  for (t in seq_along(trees)) {
    tr <- trees[[t]]
    tr$birthable <- which(lengths(tr$avail) > 0L)
    tr$closed <- closeOrder(tr$rel)
    tr$comp <- orderComponents(tr$rel)
    tr$nog <- Filter(
      function(rule) all(childKeys(splitRule(rule)) %in% tr$leafKeys),
      tr$rules
    )
    trees[[t]] <- tr
  }
  pBirth <- function(tr) {
    if (!length(tr$birthable)) {
      0
    } else if (!length(tr$rules)) {
      1
    } else {
      0.5
    }
  }
  bounds <- function(tr, k, mu, skip = integer(0)) {
    below <- setdiff(which(tr$rel[, k]), skip)
    above <- setdiff(which(tr$rel[k, ]), skip)
    c(
      if (length(below)) max(mu[tr$leafKeys[below]]) else -Inf,
      if (length(above)) min(mu[tr$leafKeys[above]]) else Inf
    )
  }
  logMass <- function(ab, m, s) {
    if (!(ab[1L] < ab[2L])) {
      return(-Inf)
    }
    if ((ab[1L] - m) / s > 0) {
      log(
        pnorm((ab[1L] - m) / s, lower.tail = FALSE) -
          pnorm((ab[2L] - m) / s, lower.tail = FALSE)
      )
    } else {
      log(pnorm((ab[2L] - m) / s) - pnorm((ab[1L] - m) / s))
    }
  }
  oneLeaf <- function(tr, k, mu, skip) {
    ab <- bounds(tr, k, mu, skip)
    out <- tr$base[k]
    if (tr$related[k]) {
      out <- out + logMass(ab, tr$m[k], tr$s[k])
      if (!exact) out <- out - logMass(ab, 0, tr$tau[k])
    }
    out
  }
  pairMass <- function(abL, abU, mL, sL, mU, sU) {
    lo <- max(abU[1L], abL[1L])
    if (lo >= abU[2L]) {
      return(0)
    }
    f <- function(u) {
      dnorm(u, mU, sU) *
        pmax(pnorm(pmin(abL[2L], u), mL, sL) - pnorm(abL[1L], mL, sL), 0)
    }
    lower <- max(lo, mU - 12 * sU)
    upper <- min(abU[2L], mU + 12 * sU)
    if (lower >= upper) {
      return(0)
    }
    integrate(f, lower, upper, rel.tol = 1e-10)$value
  }
  twoLeaves <- function(tr, kL, kU, mu) {
    abL <- bounds(tr, kL, mu, kU)
    abU <- bounds(tr, kU, mu, kL)
    out <- tr$base[kL] +
      tr$base[kU] +
      log(pairMass(abL, abU, tr$m[kL], tr$s[kL], tr$m[kU], tr$s[kU]))
    if (!exact) {
      out <- out - log(pairMass(abL, abU, 0, tr$tau[kL], 0, tr$tau[kU]))
    }
    out
  }
  branch <- function(tr, r, mu) {
    k <- match(childKeys(r), tr$leafKeys)
    if (dirs[r$v] > 0L) {
      return(twoLeaves(tr, k[1L], k[2L], mu))
    }
    if (dirs[r$v] < 0L) {
      return(twoLeaves(tr, k[2L], k[1L], mu))
    }
    oneLeaf(tr, k[1L], mu, k[2L]) + oneLeaf(tr, k[2L], mu, k[1L])
  }
  logZRatio <- function(from, to) {
    if (exact && prior == "leaf") from$logZ - to$logZ else 0
  }
  # the pair (c1 lower, c2 upper) in the finer tree fine, and m, the size of
  # fine's components that hold it
  pairOf <- function(fine, r) {
    k <- match(childKeys(r), fine$leafKeys)
    if (dirs[r$v] < 0L) {
      k <- rev(k)
    }
    list(c1 = k[1L], c2 = k[2L], m = sum(fine$comp %in% fine$comp[k]))
  }
  # switched moves decide by Barker's two-coin algorithm; the rest by MH
  decide <- function(logR, fine, from, to, r, birth) {
    pr <- pairOf(fine, r)
    if (pr$m < switchAt) {
      return(isTRUE(log(runif(1L)) < logR))
    }
    if (!is.finite(logR) && logR < 0) {
      return(FALSE)
    }
    logRest <- logR - (from$logZ - to$logZ)
    # a birth's ratio is c theta, a death's 1 / (c theta)
    logC <- if (birth) logRest + log(pr$m) else log(pr$m) - logRest
    flip <- if (coin) {
      function() adjacencyCoin(fine, pr$c1, pr$c2)
    } else {
      function() TRUE
    }
    out <- twoCoin(logC, flip, birth)
    stats$switched <- stats$switched + 1L
    stats$draws <- stats$draws + out$draws
    stats$accB <- stats$accB + stats::plogis(logR)
    stats$accMH <- stats$accMH + min(1, exp(logR))
    out$accept
  }
  stats <- new.env()
  stats$switched <- stats$draws <- 0L
  stats$accB <- stats$accMH <- 0
  cur <- trees[["ROOT"]]
  mu <- setNames(0, cur$leafKeys)
  keys <- character(nDraw)
  for (it in seq_len(nBurn + nDraw)) {
    from <- cur
    pb <- pBirth(from)
    if (runif(1L) < pb) {
      l <- from$birthable[sample.int(length(from$birthable), 1L)]
      leaf <- from$leaves[[l]]
      av <- from$avail[[l]]
      v <- av[sample.int(length(av), 1L)]
      cut <- leaf$lo[v] + sample.int(leaf$hi[v] - leaf$lo[v], 1L) - 1L
      r <- list(lo = leaf$lo, hi = leaf$hi, v = v, cut = cut)
      to <- trees[[treeKey(c(from$rules, ruleKey(leaf$lo, leaf$hi, v, cut)))]]
      logQ <- log(1 - pBirth(to)) -
        log(length(to$nog)) -
        log(pb) +
        log(length(from$birthable)) +
        log(length(av)) +
        log(leaf$hi[v] - leaf$lo[v])
      logR <- to$logPrior -
        from$logPrior +
        logQ +
        branch(to, r, mu) -
        oneLeaf(from, l, mu, integer(0)) +
        logZRatio(from, to)
      if (decide(logR, to, from, to, r, TRUE)) {
        k <- match(childKeys(r), to$leafKeys)
        mu[to$leafKeys[k]] <- 0
        abA <- bounds(to, k[1L], mu, k[2L])
        abB <- bounds(to, k[2L], mu, k[1L])
        repeat {
          a <- drawTruncated(to$m[k[1L]], to$s[k[1L]], abA[1L], abA[2L])
          b <- drawTruncated(to$m[k[2L]], to$s[k[2L]], abB[1L], abB[2L])
          if (dirs[v] == 0L || dirs[v] * (b - a) >= 0) break
        }
        mu[to$leafKeys[k]] <- c(a, b)
        mu <- mu[to$leafKeys]
        cur <- to
      }
    } else {
      rule <- from$nog[sample.int(length(from$nog), 1L)]
      r <- splitRule(rule)
      to <- trees[[treeKey(setdiff(from$rules, rule))]]
      l <- match(regionKey(r$lo, r$hi), to$leafKeys)
      logQ <- log(pBirth(to)) -
        log(length(to$birthable)) -
        log(length(to$avail[[l]])) -
        log(r$hi[r$v] - r$lo[r$v]) -
        log(1 - pb) +
        log(length(from$nog))
      muTo <- mu
      muTo[to$leafKeys[l]] <- 0
      logR <- to$logPrior -
        from$logPrior +
        logQ +
        oneLeaf(to, l, muTo, integer(0)) -
        branch(from, r, mu) +
        logZRatio(from, to)
      if (decide(logR, from, from, to, r, FALSE)) {
        ab <- bounds(to, l, muTo)
        muTo[to$leafKeys[l]] <- drawTruncated(to$m[l], to$s[l], ab[1L], ab[2L])
        mu <- muTo[to$leafKeys]
        cur <- to
      }
    }
    for (k in seq_along(cur$leafKeys)) {
      ab <- bounds(cur, k, mu)
      mu[k] <- drawTruncated(cur$m[k], cur$s[k], ab[1L], ab[2L])
    }
    if (it > nBurn) keys[it - nBurn] <- cur$key
  }
  attr(keys, "switched") <- as.list(stats)
  keys
}

# ---- the test ----------------------------------------------------------------

# Per root-rule group: batch means of the cell indicators over draws in that
# group (50 batches per chain), cells below 0.2% of the group's conditional
# law pooled into the omitted reference cell, Hotelling T^2 on the rest.
testGroups <- function(design, chains, label) {
  law <- design$law
  keys <- names(law)
  group <- rootRule(keys)
  failed <- FALSE
  chainGroups <- lapply(chains, rootRule)
  tested <- setdiff(sort(unique(group)), "")
  if (!is.null(design$spec$groups)) {
    tested <- intersect(tested, design$spec$groups)
  }
  for (g in tested) {
    inGroup <- Map(function(k, r) k[r == g], chains, chainGroups)
    inGroup <- inGroup[lengths(inGroup) >= 5000L]
    if (!length(inGroup)) {
      next
    }
    lawG <- ifelse(group == g, law, 0)
    lawG <- lawG / sum(lawG)
    kept <- which(lawG >= 0.002)
    if (length(kept) == sum(lawG > 0)) {
      kept <- kept[-which.min(lawG[kept])]
    }
    nb <- 50L
    means <- do.call(
      rbind,
      lapply(inGroup, function(k) {
        len <- length(k) %/% nb * nb
        batch <- split(k[seq_len(len)], rep(seq_len(nb), each = len / nb))
        t(vapply(
          batch,
          function(b) {
            as.numeric(table(factor(b, levels = keys[kept]))) / length(b)
          },
          numeric(length(kept))
        ))
      })
    )
    diff <- colMeans(means) - lawG[kept]
    # floored at the iid multinomial variance, so a kept cell no chain
    # visited counts against the sampler rather than making T2 singular
    nGroup <- sum(lengths(inGroup))
    covMean <- cov(means) /
      nrow(means) +
      diag(lawG[kept] * (1 - lawG[kept]) / nGroup, length(kept))
    t2 <- drop(crossprod(diff, solve(covMean, diff)))
    p <- length(kept)
    nbt <- nrow(means)
    pValue <- pf(
      (nbt - p) / ((nbt - 1) * p) * t2,
      p,
      nbt - p,
      lower.tail = FALSE
    )
    fail <- pValue < pFail
    failed <- failed || fail
    cat(sprintf(
      "  %-4s root %-8s: %d chain(s), %7d draws, %2d cells, T2 %8.1f, p %.2g%s\n",
      label,
      g,
      length(inGroup),
      sum(lengths(inGroup)),
      p,
      t2,
      pValue,
      if (fail) "  <- FAIL" else ""
    ))
  }
  failed
}

# ---- zcheck: Z_T = e(P_T) / L! against a brute-force monotonicity check ------

zcheck <- function() {
  set.seed(5L)
  worst <- 0
  for (cfg in list(
    list(nc = c(3L, 3L), dirs = c(1L, 0L)),
    list(nc = c(3L, 2L), dirs = c(1L, 1L)),
    list(nc = c(3L, 3L), dirs = c(1L, -1L)),
    list(nc = c(2L, 2L, 2L), dirs = c(1L, 1L, 1L))
  )) {
    nc <- cfg$nc
    cells <- as.matrix(expand.grid(lapply(nc, seq_len)))
    trees <- enumerateTrees(rep(1L, length(nc)), nc)
    nLeaves <- vapply(trees, function(tr) length(tr$leaves), 0L)
    for (t in sample(which(nLeaves >= 4L), 8L)) {
      leaves <- trees[[t]]$leaves
      leafOf <- apply(cells, 1L, function(x) {
        which(vapply(leaves, function(l) all(x >= l$lo & x <= l$hi), TRUE))
      })
      nSim <- 200000L
      draws <- matrix(rnorm(nSim * length(leaves)), nSim)
      ok <- rep(TRUE, nSim)
      for (i in which(cfg$dirs != 0L)) {
        for (r in seq_len(nrow(cells))) {
          up <- cells[r, ]
          up[i] <- up[i] + 1L
          if (up[i] > nc[i]) {
            next
          }
          r2 <- which(apply(cells, 1L, function(z) all(z == up)))
          ok <- ok &
            cfg$dirs[i] * (draws[, leafOf[r2]] - draws[, leafOf[r]]) >= 0
        }
      }
      z <- countExtensions(leafRelation(leaves, cfg$dirs)) /
        factorial(length(leaves))
      zStat <- (mean(ok) - z) / sqrt(z * (1 - z) / nSim)
      worst <- max(worst, abs(zStat))
      cat(sprintf(
        "  %s dirs %s, %d leaves: e/L! %.5f, simulated %.5f, z %5.2f\n",
        paste(nc, collapse = "x"),
        paste(cfg$dirs, collapse = ","),
        length(leaves),
        z,
        mean(ok),
        zStat
      ))
    }
  }
  worst
}

# ---- unnormalized: the posterior without the per-tree 1 / Z_T --------------

# the unnormalized prior p_CGM(T) prod phi 1{M in C(T)} weights each tree by
# Z_T more than the documented one, so its law is the exact law times Z_T
compareUnnormalized <- function(name) {
  design <- buildDesign(designs[[name]])
  trees <- design$trees
  constrained <- which(design$spec$dirs != 0L)
  leaves <- vapply(trees, function(tr) length(tr$leaves), 0)
  splits <- vapply(
    trees,
    function(tr) {
      v <- as.integer(sub(".*\\|([0-9]+):.*", "\\1", tr$rules))
      sum(v %in% constrained)
    },
    0
  )
  logZ <- vapply(trees, `[[`, 0, "logZ")
  unnorm <- design$law * exp(logZ - max(logZ))
  unnorm <- unnorm / sum(unnorm)
  summ <- function(w) {
    sprintf(
      "leaves %.3f, constrained splits %.3f, root-only %.3f",
      sum(w * leaves),
      sum(w * splits),
      w[["ROOT"]]
    )
  }
  cat(sprintf(
    "%s: %d structures
  documented   %s
  unnormalized %s
  TV %.3f
",
    name,
    length(trees),
    summ(design$law),
    summ(unnorm),
    sum(abs(design$law - unnorm)) / 2
  ))
}

# ---- main --------------------------------------------------------------------

# the lattice quadrature against the closed form for two leaves
local({
  worst <- 0
  for (s in c(0.2, 0.05, 0.005)) {
    for (gap in c(0, -2, 3)) {
      m <- c(0.1, 0.1 + gap * s)
      rel <- matrix(c(FALSE, FALSE, TRUE, FALSE), 2L)
      exact <- pnorm((m[2L] - m[1L]) / (sqrt(2) * s))
      worst <- max(worst, abs(orderProbability(rel, m, c(s, s)) / exact - 1))
    }
  }
  cat(sprintf("quadrature self-check: worst relative error %.1e\n", worst))
  if (worst > 1e-6) stop("order-probability quadrature is inaccurate")
})

if (mode == "zcheck") {
  worst <- zcheck()
  cat(sprintf("zcheck: worst |z| %.2f\n", worst))
  if (worst > 4) {
    quit(status = 1L)
  }
  quit(status = 0L)
}

if (mode == "unnormalized") {
  for (name in chosen) {
    compareUnnormalized(name)
  }
  quit(status = 0L)
}

anyFailure <- FALSE
for (prior in priors) {
  for (name in chosen) {
    design <- buildDesign(designs[[name]])
    design$law <- design$laws[[prior]]
    cat(sprintf(
      "%s: %d structures, %s draws over %d chains (%s, prior \"%s\")\n",
      name,
      length(design$trees),
      format(nDraws, big.mark = ","),
      nChains,
      mode,
      prior
    ))
    perChain <- nDraws %/% nChains
    chains <- lapply(seq_len(nChains), function(s) {
      switch(
        mode,
        engine = engineKeys(design, perChain, s, prior),
        prototype = prototypeKeys(
          design,
          perChain,
          s,
          exact = TRUE,
          prior = prior
        ),
        "prototype-old" = prototypeKeys(design, perChain, s, exact = FALSE),
        hybrid = prototypeKeys(
          design,
          perChain,
          s,
          exact = TRUE,
          switchAt = 3L
        ),
        "hybrid-nocoin" = prototypeKeys(
          design,
          perChain,
          s,
          exact = TRUE,
          switchAt = 3L,
          coin = FALSE
        )
      )
    })
    if (startsWith(mode, "hybrid")) {
      st <- Reduce(
        function(a, b) Map(`+`, a, b),
        lapply(chains, attr, "switched")
      )
      cat(sprintf(
        "  %d switched moves, %.2f coins each, acceptance %.3f of MH's\n",
        st$switched,
        st$draws / st$switched,
        st$accB / st$accMH
      ))
    }
    anyFailure <- testGroups(design, chains, name) || anyFailure
  }
}

if (anyFailure) {
  cat("\nFAIL: the sampler does not target the documented monotone posterior\n")
  if (!mode %in% c("prototype-old", "hybrid-nocoin")) quit(status = 1L)
} else {
  cat("\nOK: every design matches the enumerated posterior\n")
}
