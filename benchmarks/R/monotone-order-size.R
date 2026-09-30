#!/usr/bin/env Rscript

# Size of the monotone leaf order in fitted trees, the input to the count
# limit in docs/plans/monotone-exact-birth-death.md (Decision). The exact
# normalizer Z_T = e(P_T) / L! is counted by a DP over the down-sets of the
# leaf order P_T, one component at a time; its work is about (down-sets) x
# (leaves) per component. A fit's trees are read with getTrees(); each tree's
# order is built with the engine's adjacency rule (j < k when j's box ends
# where k's starts on a constrained axis and the boxes overlap on every other
# axis), closed transitively, split into components, and its down-sets are
# counted by memoized recursion.
#
# Per fit it prints the sampler's time per sweep, the most leaves any kept tree
# has, and, as median and maximum over kept sweeps: the largest whole-tree
# down-set count (the product over components) in the sweep, with the share of
# sweeps where it exceeds 2^22; the largest component count; and the sweep's
# count work, summed over its trees, of a move whose touched leaf is uniform:
# sum over components of (|C| / L) D(C) |C|, in down-set x leaf units.
#
# closure mode checks which measures deaths never raise, over every death of
# random guillotine trees, and prints the two-by-two tree whose every death
# raises the largest component's count (a per-component budget of 4 traps it
# away from the root).
#
# Usage: Rscript monotone-order-size.R nTrees n p nConstrained [seed]
#        Rscript monotone-order-size.R closure
#   x is uniform on [0, 1]^p; y increases in the first nConstrained axes and
#   has 2 sin(2 pi x_free) (1 + x1) in the first free one; 1000 burn-in and
#   200 kept draws of the installed engine.

args <- commandArgs(trailingOnly = TRUE)
budget <- 2^22
# a component whose recursion outgrows this many states is reported as Inf
memoCap <- 2^23

# ---- trees, orders and counts ------------------------------------------------

# a tree is list() for a leaf or list(var, cut, left, right); left is x <= cut
fromPreorder <- function(var, value) {
  i <- 0L
  build <- function() {
    i <<- i + 1L
    if (var[i] == -1L) {
      return(list())
    }
    node <- list(var[i], value[i])
    node[[3L]] <- build()
    node[[4L]] <- build()
    node
  }
  build()
}

leafBoxes <- function(tree, p) {
  out <- list()
  walk <- function(node, lo, hi) {
    if (!length(node)) {
      out[[length(out) + 1L]] <<- list(lo = lo, hi = hi)
      return(invisible())
    }
    v <- node[[1L]]
    hiLeft <- hi
    hiLeft[v] <- node[[2L]]
    walk(node[[3L]], lo, hiLeft)
    loRight <- lo
    loRight[v] <- node[[2L]]
    walk(node[[4L]], loRight, hi)
  }
  walk(tree, rep(-Inf, p), rep(Inf, p))
  out
}

# R[j, k] when leaf j lies below leaf k, transitively closed; every
# constrained axis increasing
leafOrder <- function(boxes, nc) {
  n <- length(boxes)
  p <- length(boxes[[1L]]$lo)
  rel <- matrix(FALSE, n, n)
  for (j in seq_len(n)) {
    for (k in seq_len(n)) {
      if (j == k) {
        next
      }
      a <- boxes[[j]]
      b <- boxes[[k]]
      for (i in seq_len(nc)) {
        if (a$hi[i] != b$lo[i]) {
          next
        }
        other <- setdiff(seq_len(p), i)
        if (
          all(pmax(a$lo[other], b$lo[other]) < pmin(a$hi[other], b$hi[other]))
        ) {
          rel[j, k] <- TRUE
        }
      }
    }
  }
  repeat {
    closed <- rel | (rel %*% rel > 0)
    if (all(closed == rel)) {
      return(rel)
    }
    rel <- closed
  }
}

components <- function(rel) {
  adj <- rel | t(rel)
  n <- nrow(adj)
  label <- seq_len(n)
  repeat {
    nextLabel <- vapply(
      seq_len(n),
      function(j) min(label[adj[j, ] | seq_len(n) == j]),
      integer(1L)
    )
    if (all(nextLabel == label)) {
      return(unname(split(seq_len(n), label)))
    }
    label <- nextLabel
  }
}

# down-sets of S: those without x (nothing above x) plus those with x (all of
# its down-closure), each a down-set of what is left
countDownSets <- function(rel) {
  memo <- new.env(hash = TRUE)
  states <- 0L
  count <- function(set) {
    if (!length(set)) {
      return(1)
    }
    key <- paste(set, collapse = ".")
    hit <- memo[[key]]
    if (!is.null(hit)) {
      return(hit)
    }
    states <<- states + 1L
    if (states > memoCap) {
      stop("memo cap")
    }
    x <- set[1L]
    above <- set[rel[x, set]]
    below <- set[rel[set, x]]
    result <- count(setdiff(set, c(x, above))) +
      count(setdiff(set, c(x, below)))
    assign(key, result, envir = memo)
    result
  }
  tryCatch(count(seq_len(nrow(rel))), error = function(e) Inf)
}

orderSize <- function(tree, p, nc) {
  rel <- leafOrder(leafBoxes(tree, p), nc)
  parts <- components(rel)
  counts <- vapply(
    parts,
    function(ix) countDownSets(rel[ix, ix, drop = FALSE]),
    numeric(1L)
  )
  sizes <- lengths(parts)
  c(
    leaves = nrow(rel),
    whole = prod(counts),
    largest = max(counts),
    work = sum(sizes / nrow(rel) * counts * sizes)
  )
}

# ---- fits --------------------------------------------------------------------

fitSizes <- function(nTrees, n, p, nc, seed) {
  suppressPackageStartupMessages(library(dbarts))
  set.seed(seed)
  x <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
  f <- 2 * rowSums(x[, seq_len(nc), drop = FALSE])
  if (p > nc) {
    f <- f + 2 * sin(2 * pi * x[, nc + 1L]) * (1 + x[, 1L])
  }
  y <- f + rnorm(n, 0, 0.3)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = nTrees,
    n.samples = 200L,
    keepTrees = TRUE,
    updateState = FALSE,
    seed = seed + 1L
  )
  monotone <- setNames(rep(1L, nc), colnames(x)[seq_len(nc)])
  sampler <- dbarts(x, y, control = control, monotone = monotone)
  seconds <- system.time(invisible(sampler$run(1000L, 200L)))[["elapsed"]]
  trees <- sampler$getTrees()
  keys <- unique(trees[c("sample", "tree")])
  sizes <- t(vapply(
    seq_len(nrow(keys)),
    function(r) {
      d <- trees[
        trees$sample == keys$sample[r] & trees$tree == keys$tree[r],
      ]
      orderSize(fromPreorder(d$var, d$value), p, nc)
    },
    numeric(4L)
  ))
  perSweep <- function(column, f) tapply(sizes[, column], keys$sample, f)
  whole <- perSweep("whole", max)
  largest <- perSweep("largest", max)
  q <- function(v) {
    sprintf("%.2g / %.2g", stats::median(v), max(v))
  }
  cat(sprintf(
    paste0(
      "trees %d n %d p %d nc %d seed %d: %.2f ms/sweep, leaves max %d; ",
      "per sweep med/max: whole-tree %s (over 2^22 in %.0f%%), largest ",
      "component %s, count work %s\n"
    ),
    nTrees,
    n,
    p,
    nc,
    seed,
    1000 * seconds / 1200,
    max(sizes[, "leaves"]),
    q(whole),
    100 * mean(whole > budget),
    q(largest),
    q(perSweep("work", sum))
  ))
}

# ---- closure under deaths ----------------------------------------------------

randomTree <- function(nLeaves, p) {
  tree <- list()
  paths <- list(integer(0L))
  boxes <- list(list(lo = rep(0, p), hi = rep(1, p)))
  while (length(paths) < nLeaves) {
    j <- sample.int(length(paths), 1L)
    box <- boxes[[j]]
    v <- sample.int(p, 1L)
    cuts <- (1:15 / 16)[1:15 / 16 > box$lo[v] & 1:15 / 16 < box$hi[v]]
    if (!length(cuts)) {
      next
    }
    cut <- cuts[sample.int(length(cuts), 1L)]
    node <- list(v, cut, list(), list())
    if (length(paths[[j]])) {
      tree[[paths[[j]]]] <- node
    } else {
      tree <- node
    }
    left <- box
    left$hi[v] <- cut
    right <- box
    right$lo[v] <- cut
    paths <- c(paths[-j], list(c(paths[[j]], 3L), c(paths[[j]], 4L)))
    boxes <- c(boxes[-j], list(left, right))
  }
  tree
}

# paths to the internal nodes whose children are both leaves
nogPaths <- function(node, path = integer(0L)) {
  if (!length(node)) {
    return(list())
  }
  if (!length(node[[3L]]) && !length(node[[4L]])) {
    return(list(path))
  }
  c(nogPaths(node[[3L]], c(path, 3L)), nogPaths(node[[4L]], c(path, 4L)))
}

closureCheck <- function(nTrees = 300L, seed = 5L) {
  set.seed(seed)
  raised <- c(whole = 0L, largest = 0L, work = 0L)
  deaths <- 0L
  for (r in seq_len(nTrees)) {
    p <- sample(2:3, 1L)
    nc <- sample.int(p, 1L)
    tree <- randomTree(sample(6:16, 1L), p)
    before <- orderSize(tree, p, nc)
    for (path in nogPaths(tree)) {
      dead <- tree
      if (length(path)) {
        dead[[path]] <- list()
      } else {
        dead <- list()
      }
      after <- orderSize(dead, p, nc)
      deaths <- deaths + 1L
      for (m in names(raised)) {
        raised[[m]] <- raised[[m]] + (after[[m]] > before[[m]] * (1 + 1e-12))
      }
    }
  }
  cat(sprintf(
    "%d deaths of %d random trees raise: whole-tree %d, largest component %d, move work %d\n",
    deaths,
    nTrees,
    raised[["whole"]],
    raised[["largest"]],
    raised[["work"]]
  ))

  # x1 constrained, x2 free: x1 split, then x2 on both sides; A < R1, B < R2
  trap <- list(
    1L,
    0.5,
    list(2L, 0.5, list(), list()),
    list(2L, 0.5, list(), list())
  )
  show <- function(label, tree) {
    s <- orderSize(tree, 2L, 1L)
    cat(sprintf(
      "%s: whole-tree %g, largest component %g\n",
      label,
      s[["whole"]],
      s[["largest"]]
    ))
  }
  show("two-by-two tree", trap)
  lower <- trap
  lower[[3L]] <- list()
  show("  after the lower death", lower)
  upper <- trap
  upper[[4L]] <- list()
  show("  after the upper death", upper)
}

if (identical(args[1L], "closure")) {
  closureCheck()
} else {
  a <- as.integer(args)
  fitSizes(a[1L], a[2L], a[3L], a[4L], if (length(a) > 4L) a[5L] else 1L)
}
